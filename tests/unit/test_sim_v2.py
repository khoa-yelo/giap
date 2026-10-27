"""gippasim v2: 100-kb genome with signals and precomputed layers (issue M2-05; GIPPA_SPEC.md §11.3).

Examples 1-10 and the edge cases of docs/specs/M2-05-simulator-v2.md, the round trip through
M2-03's ``parse_bakta_json``, M2-02's ``load_precomputed`` and M1-05's ``read_gff3``,
backward compatibility with the M0 and M3-18 outputs, determinism and Hypothesis properties.

Conventions these tests pin down where the spec is silent (see the PR description):

* the plasmid record of ``SIM100K`` is named ``p1`` in the spec's example 1 but M3-18 (and its
  merged tests) name plasmids ``plasmid_<n>``; either id is accepted, the record must be the
  second one, 8,000 bp, circular and ``type: plasmid``;
* example 2's "gene count" is the number of ``cds``/``pseudogene`` genes (M3-18 keeps sORFs out of
  ``truth/genes.tsv``; test_bakta_json_sim.py pins that);
* ``truth/tus.tsv`` has at least ``tu_id, contig, start, end, strand, genes`` (§9.3 ``tus.tsv``
  names, ``;``-joined locus tags), ``promoter`` (``true``/``false``) and ``tss`` (1-based TSS on
  the TU strand; ``NA`` allowed when ``promoter`` is false); coordinates 1-based closed;
* ``truth/{promoters,terminators,rbs,regions}.tsv`` hold the columns of the matching M2-02 layer
  plus ``true_element`` (``true``/``false``) and exactly the rows of ``precomputed/<layer>.tsv``
  (the layer = true elements + false positives);
* ``layer_false_positive_rate`` r: the expected number of false rows of a layer is r times its
  number of true rows (pooled calibration in tests/calibration/test_sim_v2_calibration.py);
* signal strength s is the per-element probability of planting the consensus (``0`` = none,
  ``1`` = every TU/gene); the ``-10`` box ``TATAAT`` ends 1-12 nt upstream of the TSS, after
  ``TTGACA`` and a 17-nt spacer; the RBS spacer is the number of nt between the last ``G`` of
  ``AGGAGG`` and the first base of the start codon (5-10);
* ``truth/features.tsv`` gains a trailing ``anticodon`` column (DNA alphabet) and rows of type
  ``crispr`` and ``oric``; the M3-18 columns stay first;
* the bundled sequences live in ``gippasim/data/`` next to the package.
"""

from __future__ import annotations

import csv
import dataclasses
import hashlib
import itertools
import json
import re
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import gippasim
from gippa.annotate.bakta import BaktaAnnotation, parse_bakta_json
from gippa.annotate.precomputed import PRECOMPUTED_COLUMNS, load_precomputed
from gippa.config import resolve_config
from gippa.io.gff3 import read_gff3
from gippa.model.coords import Interval
from gippasim import SimConfig, simulate

ANTI_SD_TAIL = "GATCACCTCCTTA"
STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
LAYERS = ("operon_pairs", "promoters", "tx_rates", "terminators", "rbs", "regions")
TRUTH = ("genes", "tus", "promoters", "terminators", "rbs", "regions", "features")
TRUTH_LAYERS = ("promoters", "terminators", "rbs", "regions")
PROMOTER = re.compile(r"(?=(TTGACA[ACGT]{17}TATAAT))")
COMP = str.maketrans("ACGTN", "TGCAN")
REGION_TYPES = {"prophage", "plasmid", "bgc", "crispr", "amr"}
FEATURE_KEYS = {
    "type", "sequence", "start", "stop", "strand", "locus", "gene", "product", "start_type",
    "nt", "aa", "db_xrefs", "pseudo", "edge", "anti_codon", "amino_acid", "id", "hypothetical",
}  # fmt: skip
SEQUENCE_KEYS = {"id", "type", "topology", "length", "nt"}

#: Output of the merged M3-18 simulator (main at bd06a0c) for ``M318_CONFIG`` and seed 3. M2-05's
#: fields default to off, so the genome, gene truth and parsed gene model must not change.
M318_CONFIG = {
    "genome_length": 60_000,
    "n_plasmids": 1,
    "plasmid_length": (6_000, 9_000),
    "origin_crossing_gene": True,
    "n_rrna_operons": 1,
    "n_trna": 4,
    "n_pseudogenes": 2,
    "n_sorfs": 3,
}
M318_SHA256 = {
    "genome.fasta": "82985e5898ad793f89410da740e6c23fc3acf6256ebf1b845bbbb22b2789f62a",
    "truth/genes.tsv": "8b35e4e29b19d65ba3941630c7b9407a92ef5b34f8a96ef2c13c3f8bea55c62e",
    # main at 42b9412 (review decision 2: no ``anticodon`` column without v2 output)
    "truth/features.tsv": "8149a88ee76437fbaf79f779ab19f2f6131f8138e00c3974a970ed77703b4c5b",
}
M318_MODEL_SHA256 = "00267c904e4528bf7a477963782fd87dfe01788808e1dc4ccae340869130a8ee"


# ---------------------------------------------------------------- helpers


def sim100k() -> tuple[Any, int]:
    for mod in (gippasim, gippasim.simulate):
        if hasattr(mod, "SIM100K"):
            return mod.SIM100K, mod.SIM100K_SEED
    raise AssertionError("gippasim exports no SIM100K / SIM100K_SEED")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def header(path: Path) -> list[str]:
    with path.open() as fh:
        return fh.readline().rstrip("\n").split("\t")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """(header line without '>', sequence) per record."""
    out: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            out.append((line[1:], []))
        elif line.strip():
            out[-1][1].append(line.strip().upper())
    return [(h, "".join(s)) for h, s in out]


def is_true(v: str) -> bool:
    assert v in ("true", "false"), v
    return v == "true"


def internal_stop(nt: str) -> bool:
    return any(nt[i : i + 3] in STOPS for i in range(0, len(nt) - 3, 3))


def gc(s: str) -> float:
    return (s.count("G") + s.count("C")) / len(s)


def tree(root: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()
    }


def model_digest(ann: BaktaAnnotation) -> str:
    m = ann.model
    rows: list[list[Any]] = [
        ["rep", r.id, r.length, r.topology, r.kind] for r in m.replicons.values()
    ]
    rows += [
        ["gene", g.locus_tag, g.feature_type, g.interval.contig, g.interval.start,
         g.interval.end, g.interval.strand, g.aa, g.nt, g.start_codon]
        for g in m.genes
    ]  # fmt: skip
    rows += [
        ["nc", f.id, f.type, f.interval.contig, f.interval.start, f.interval.end,
         f.interval.strand, f.anticodon, f.amino_acid]
        for f in m.features
    ]  # fmt: skip
    return hashlib.sha256(json.dumps(rows).encode()).hexdigest()


class Run:
    """One simulator output directory, parsed."""

    def __init__(self, res: Any, cfg: Any) -> None:
        self.res, self.cfg = res, cfg
        self.out = Path(res.outdir)
        self.data: dict[str, Any] = json.loads((self.out / "bakta_like.json").read_text())
        self.ann = parse_bakta_json(self.out / "bakta_like.json")
        self.reps = self.ann.model.replicons
        self.seqs = dict(self.ann.sequences)

    def truth(self, name: str) -> list[dict[str, str]]:
        return read_rows(self.out / "truth" / f"{name}.tsv")

    def layer(self, name: str) -> list[dict[str, str]]:
        return read_rows(self.out / "precomputed" / f"{name}.tsv")

    def coding(self) -> list[Any]:
        return [g for g in self.ann.model.genes]

    def upstream(self, contig: str, start0: int, end0: int, strand: str, n: int) -> str:
        """``n`` nt upstream of [start0, end0) in transcription orientation (circular wrap)."""
        s = self.seqs[contig]
        doubled = s + s + s
        L = len(s)
        if strand == "+":
            return doubled[L + start0 - n : L + start0]
        e = end0 % L if end0 > L else end0
        return revcomp(doubled[L + e : L + e + n])

    def slice(self, iv: Interval) -> str:
        s = self.seqs[iv.contig]
        x = (s + s)[iv.start : iv.end]
        return revcomp(x) if iv.strand == "-" else x


def run_sim(cfg: Any, seed: int, outdir: Path) -> Run:
    return Run(simulate(cfg, seed, outdir), cfg)


@pytest.fixture(scope="module")
def r100(tmp_path_factory: pytest.TempPathFactory) -> Run:
    cfg, seed = sim100k()
    return run_sim(cfg, seed, tmp_path_factory.mktemp("sim100k"))


@pytest.fixture(scope="module")
def rbs1(tmp_path_factory: pytest.TempPathFactory) -> Run:
    cfg, seed = sim100k()
    cfg = dataclasses.replace(cfg, rbs_strength=1.0)
    return run_sim(cfg, seed, tmp_path_factory.mktemp("rbs1"))


def bundled() -> Path:
    return Path(gippasim.__file__).parent / "data"


# ---------------------------------------------------------------- config surface


def test_new_fields_default_off() -> None:
    cfg = SimConfig()
    assert cfg.real_rrna_trna is False
    assert (cfg.crispr, cfg.oric) == (False, False)
    assert (cfg.prophage_bp, cfg.bgc_bp) == (0, 0)
    assert cfg.prophage_gc_shift == pytest.approx(-0.08)
    assert (cfg.rbs_strength, cfg.promoter_strength, cfg.terminator_strength) == (0.0, 0.0, 0.0)
    assert cfg.layer_false_positive_rate == pytest.approx(0.05)
    assert cfg.write_precomputed is False


def test_sim100k_constants() -> None:
    cfg, seed = sim100k()
    assert seed == 20261014
    assert cfg == SimConfig(
        genome_length=100_000,
        n_plasmids=1,
        plasmid_length=(8_000, 8_000),
        n_rrna_operons=1,
        n_trna=10,
        real_rrna_trna=True,
        crispr=True,
        oric=True,
        prophage_bp=12_000,
        bgc_bp=10_000,
        n_pseudogenes=2,
        n_sorfs=3,
        rbs_strength=0.8,
        promoter_strength=0.8,
        terminator_strength=0.8,
        write_precomputed=True,
    )
    hash(cfg)  # stays hashable (frozen dataclass)


def test_simresult_paths(r100: Run) -> None:
    res, out = r100.res, r100.out
    assert Path(res.gff3) == out / "annotation.gff3" and Path(res.gff3).is_file()
    assert Path(res.config_yaml) == out / "config.yaml" and Path(res.config_yaml).is_file()
    assert set(res.precomputed) == set(LAYERS)
    for k, p in res.precomputed.items():
        assert Path(p) == out / "precomputed" / f"{k}.tsv" and Path(p).is_file(), k
    assert set(res.truth) == set(TRUTH)
    for k, p in res.truth.items():
        assert Path(p) == out / "truth" / f"{k}.tsv" and Path(p).is_file(), k


# ---------------------------------------------------------------- example 1: records, speed


def test_example1_fasta_records(r100: Run) -> None:
    recs = read_fasta(r100.out / "genome.fasta")
    assert len(recs) == 2
    (h1, s1), (h2, s2) = recs
    assert h1.split()[0] == "chr1" and "[topology=circular]" in h1
    assert len(s1) == 100_000
    assert h2.split()[0] in ("p1", "plasmid_1") and "[topology=circular]" in h2
    assert len(s2) == 8_000
    assert set(s1 + s2) <= set("ACGT")


@pytest.mark.slow  # wall clock: flaky on loaded shared nodes, so not in the every-commit suite
def test_example1_under_1s(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    best = float("inf")
    for i in range(3):
        t0 = time.perf_counter()
        simulate(cfg, seed, tmp_path / str(i))
        best = min(best, time.perf_counter() - t0)
    assert best < 1.0, best


# ---------------------------------------------------------------- example 2: Bakta-like JSON


def test_example2_parse_matches_gene_truth(r100: Run) -> None:
    truth = {r["locus_tag"]: r for r in r100.truth("genes")}
    genes = {
        g.locus_tag: g for g in r100.ann.model.genes if g.feature_type in ("cds", "pseudogene")
    }
    assert len(genes) == len(truth)
    assert set(genes) == set(truth)
    for locus, row in truth.items():
        iv = genes[locus].interval
        assert (iv.start + 1, iv.end, iv.strand) == (
            int(row["start"]),
            int(row["end"]),
            row["strand"],
        ), locus
        assert genes[locus].start_codon == row["start_codon"], locus


def test_example2_bakta_field_names(r100: Run) -> None:
    d = r100.data
    assert set(d) <= {"sequences", "features", "version"}
    assert d["version"] == {"bakta": "sim", "db": {"version": "sim", "type": "sim"}}
    ids = []
    for s in d["sequences"]:
        assert SEQUENCE_KEYS <= set(s) <= SEQUENCE_KEYS | {"complete"}, sorted(s)
        assert len(s["nt"]) == s["length"]
        ids.append(s["id"])
    for f in d["features"]:
        assert set(f) <= FEATURE_KEYS, (f.get("locus"), sorted(set(f) - FEATURE_KEYS))
        assert f["sequence"] in ids
        assert isinstance(f["start"], int) and isinstance(f["stop"], int)
    for f in d["features"]:
        if f["type"] in ("cds", "sorf"):
            for k in ("locus", "gene", "product", "start_type", "nt", "aa", "db_xrefs"):
                assert k in f, (f["locus"], k)
            assert isinstance(f["db_xrefs"], list)


def test_example2_parser_reads_everything(r100: Run) -> None:
    ann = r100.ann
    assert dict(ann.skipped_types) == {}
    assert (ann.bakta_version, ann.db_version, ann.db_type) == ("sim", "sim", "unknown")
    fasta = {h.split()[0]: s for h, s in read_fasta(r100.out / "genome.fasta")}
    assert dict(ann.sequences) == fasta
    kinds = [(r.id, r.length, r.topology, r.kind) for r in ann.model.replicons.values()]
    assert kinds[0] == ("chr1", 100_000, "circular", "chromosome")
    assert kinds[1][1:] == (8_000, "circular", "plasmid")
    assert len(kinds) == 2
    # every coding feature's nt is the genome slice; aa is its translation's protein
    for g in ann.model.genes:
        assert g.nt == r100.slice(g.interval), g.locus_tag


def test_plasmid_has_1_to_3_genes(r100: Run) -> None:
    pid = list(r100.reps)[1]
    n = sum(g.interval.contig == pid for g in r100.ann.model.genes)
    assert 1 <= n <= 3, n
    seq = next(s for s in r100.data["sequences"] if s["id"] == pid)
    assert seq["type"] == "plasmid"


def test_gff3_matches_json(r100: Run) -> None:
    feats = read_gff3(r100.out / "annotation.gff3", r100.reps)
    ids = {f.attributes.get("locus_tag") or f.attributes.get("ID") for f in feats}
    assert {g.locus_tag for g in r100.ann.model.genes} <= ids


# ---------------------------------------------------------------- example 3: 16S and anti-SD


def rrna16s_nts(run: Run) -> list[str]:
    return [
        f["nt"]
        for f in run.data["features"]
        if f["type"] == "rRNA" and ("16S" in f.get("product", "") or f.get("gene") == "rrs")
    ]


def test_example3_16s_tail(r100: Run) -> None:
    nts = rrna16s_nts(r100)
    assert len(nts) == 1
    assert nts[0].endswith(ANTI_SD_TAIL)
    nc = [f for f in r100.ann.model.features if f.type == "rrna" and "16S" in f.product]
    assert len(nc) == 1 and nc[0].nt == nts[0]
    assert r100.slice(nc[0].interval) == nts[0]


def test_example3_anti_sd(r100: Run) -> None:
    layers = pytest.importorskip("gippa.annotate.bakta_layers", reason="M2-04 not merged yet")
    assert layers.anti_sd(rrna16s_nts(r100)) == ("ACCTCCTTA", ())


def test_bundled_real_sequences() -> None:
    d = bundled()
    rrn, trn = d / "ecoli_rrnB.fa", d / "ecoli_trna.fa"
    assert rrn.is_file() and trn.is_file() and (d / "PROVENANCE.md").is_file()
    prov = (d / "PROVENANCE.md").read_text()
    assert "MG1655" in prov and "NC_000913" in prov
    r, t = read_fasta(rrn), read_fasta(trn)
    assert sum(len(s) for _, s in r + t) < 20_000
    assert len(t) >= 10
    assert any(s.endswith(ANTI_SD_TAIL) for _, s in r)
    for _, s in r + t:
        assert set(s) <= set("ACGT")


def test_real_rrna_trna_inserted_verbatim(r100: Run) -> None:
    pool = {s for _, s in read_fasta(bundled() / "ecoli_rrnB.fa")}
    pool |= {s for _, s in read_fasta(bundled() / "ecoli_trna.fa")}
    nc = [f for f in r100.ann.model.features if f.type in ("rrna", "trna")]
    assert sum(f.type == "rrna" for f in nc) == 3
    assert sum(f.type == "trna" for f in nc) >= 10
    for f in nc:
        assert r100.slice(f.interval) in pool, (f.id, f.type)
    assert any(f.nt in pool for f in nc if f.type == "rrna" and f.nt)


def test_no_gene_overlaps_rna_crispr_oric(r100: Run) -> None:
    m = r100.ann.model
    for f in m.features:
        hits = [x for x in m.overlapping(f.interval) if x is not f]
        genes = [getattr(x, "locus_tag", None) for x in hits if hasattr(x, "locus_tag")]
        assert not genes, (f.id, f.type, genes)


# ---------------------------------------------------------------- example 4: features truth


def test_example4_features_truth(r100: Run) -> None:
    rows = r100.truth("features")
    cols = header(r100.out / "truth" / "features.tsv")
    assert cols[:6] == ["feature_id", "type", "contig", "start", "end", "strand"]
    assert "anticodon" in cols
    trna = [r for r in rows if r["type"] == "tRNA"]
    assert len(trna) >= 10
    for r in trna:
        assert re.fullmatch("[ACGT]{3}", r["anticodon"]), r
    types = Counter(r["type"] for r in rows)
    assert types["crispr"] == 1
    assert types["oric"] == 1
    # truth tRNA anticodons are the JSON's
    json_ac = {f["locus"] if "locus" in f else f["id"]: f.get("anti_codon") for f in
               r100.data["features"] if f["type"] == "tRNA"}  # fmt: skip
    for r in trna:
        assert json_ac[r["feature_id"]].upper().replace("U", "T") == r["anticodon"]


def test_example4_crispr_and_oric_in_json(r100: Run) -> None:
    rows = {r["type"]: r for r in r100.truth("features") if r["type"] in ("crispr", "oric")}
    nc = {f.type: f for f in r100.ann.model.features if f.type in ("crispr", "oric")}
    assert set(nc) == {"crispr", "oric"}
    for t in ("crispr", "oric"):
        iv, r = nc[t].interval, rows[t]
        assert (iv.contig, iv.start + 1, iv.end) == (r["contig"], int(r["start"]), int(r["end"]))


def test_example4_oric_next_to_mnmG(r100: Run) -> None:
    oric = next(f for f in r100.ann.model.features if f.type == "oric")
    assert oric.interval.contig == "chr1"
    L = r100.reps["chr1"].length

    def dist(g: Any) -> int:
        a, b = g.interval, oric.interval
        d = max(a.start - b.end, b.start - a.end, 0)
        return min(d, L - max(a.end, b.end) + min(a.start, b.start))

    genes = [g for g in r100.ann.model.genes if g.interval.contig == "chr1"]
    nearest = min(genes, key=dist)
    assert nearest.gene == "mnmG", (nearest.locus_tag, nearest.gene)
    named = [g for g in r100.ann.model.genes if g.gene == "mnmG"]
    assert len(named) == 1


# ---------------------------------------------------------------- example 5: promoters


def tus(run: Run) -> list[dict[str, str]]:
    rows = run.truth("tus")
    cols = header(run.out / "truth" / "tus.tsv")
    for c in ("tu_id", "contig", "start", "end", "strand", "genes", "promoter", "tss"):
        assert c in cols, c
    return rows


def check_promoter(run: Run, tu: dict[str, str]) -> None:
    genes = {g.locus_tag: g for g in run.ann.model.genes}
    members = [genes[x] for x in tu["genes"].split(";")]
    contig, strand, tss = tu["contig"], tu["strand"], int(tu["tss"])
    if strand == "+":
        window = run.upstream(contig, tss - 1, tss, "+", 50)
        first = min(g.interval.start for g in members)
        assert tss - 1 <= first, tu["tu_id"]
    else:
        window = run.upstream(contig, tss - 1, tss, "-", 50)
        assert tss >= max(g.interval.end for g in members), tu["tu_id"]
    # the -10 box ends 1-12 nt upstream of the TSS (window[-1] is position -1)
    ok = [m for m in PROMOTER.finditer(window) if 0 <= len(window) - (m.start() + 29) <= 11]
    assert ok, (tu["tu_id"], strand, window)


def test_example5_promoter_motifs(r100: Run) -> None:
    rows = tus(r100)
    promoted = [r for r in rows if is_true(r["promoter"])]
    assert promoted
    for r in promoted:
        check_promoter(r100, r)


def test_tus_cover_gene_truth(r100: Run) -> None:
    rows = tus(r100)
    by_tu = {r["tu_id"]: set(r["genes"].split(";")) for r in rows}
    genes = r100.truth("genes")
    loci = {g.locus_tag for g in r100.ann.model.genes}
    for members in by_tu.values():
        assert members <= loci
    for g in genes:
        if g["tu_id"]:
            assert g["locus_tag"] in by_tu[g["tu_id"]], g


# ---------------------------------------------------------------- example 6: RBS


def rbs_spacer(run: Run, g: Any) -> int | None:
    iv = g.interval
    up = run.upstream(iv.contig, iv.start, iv.end, iv.strand, 20)
    for sp in range(5, 11):
        if up[len(up) - sp - 6 : len(up) - sp] == "AGGAGG":
            return sp
    return None


def test_example6_every_gene_has_rbs(rbs1: Run) -> None:
    missing = [g.locus_tag for g in rbs1.ann.model.genes if rbs_spacer(rbs1, g) is None]
    assert not missing, missing


def test_true_rbs_rows_match_genome(r100: Run) -> None:
    genes = {g.locus_tag: g for g in r100.ann.model.genes}
    true_rows = [r for r in r100.truth("rbs") if is_true(r["true_element"])]
    assert true_rows
    for r in true_rows:
        sp = int(r["spacer_bp"])
        assert 5 <= sp <= 10, r
        g = genes[r["locus_tag"]]
        assert rbs_spacer(r100, g) == sp, r
        assert r["start_codon"] == g.start_codon


# ---------------------------------------------------------------- example 7: prophage GC


def true_regions(run: Run, kind: str) -> list[dict[str, str]]:
    return [
        r for r in run.truth("regions") if r["region_type"] == kind and is_true(r["true_element"])
    ]


def test_example7_prophage_gc_shift(r100: Run) -> None:
    (pro,) = true_regions(r100, "prophage")
    chr1 = r100.seqs["chr1"]
    s, e = int(pro["start"]) - 1, int(pro["end"])
    assert pro["contig"] == "chr1" and 0 <= s < e <= len(chr1)
    diff = gc(chr1[:s] + chr1[e:]) - gc(chr1[s:e])
    assert 0.06 <= diff <= 0.10, diff


@pytest.mark.parametrize(("kind", "attr"), [("prophage", "prophage_bp"), ("bgc", "bgc_bp")])
def test_blocks_sized_and_made_of_genes(r100: Run, kind: str, attr: str) -> None:
    (reg,) = true_regions(r100, kind)
    want = getattr(r100.cfg, attr)
    s, e = int(reg["start"]) - 1, int(reg["end"])
    assert 0.9 * want <= e - s <= 1.1 * want, (kind, e - s)
    genes = {g.locus_tag: g for g in r100.ann.model.genes}
    listed = reg["genes"].split(";")
    assert listed and set(listed) <= set(genes)
    for x in listed:
        iv = genes[x].interval
        assert iv.contig == reg["contig"] and s <= iv.start and iv.end <= e, (kind, x)


# ---------------------------------------------------------------- example 8: pseudogenes, sORFs


def test_example8_pseudogenes_and_sorfs(r100: Run) -> None:
    pseudo = [f for f in r100.data["features"] if f["type"] == "cds" and f.get("pseudo")]
    assert len(pseudo) == 2
    for f in pseudo:
        nt = f["nt"]
        assert len(nt) % 3 != 0 or internal_stop(nt), f["locus"]
    sorfs = [f for f in r100.data["features"] if f["type"] == "sorf"]
    assert len(sorfs) == 3
    for f in sorfs:
        assert len(f["aa"].rstrip("*")) < 30, f["locus"]
    for g in r100.ann.model.genes:  # everything else is a clean CDS
        if g.feature_type == "pseudogene":
            continue
        nt = g.nt or ""
        assert len(nt) % 3 == 0 and nt[:3] in STARTS and nt[-3:] in STOPS, g.locus_tag
        assert not internal_stop(nt), g.locus_tag


# ---------------------------------------------------------------- example 9: layers, config


@pytest.mark.parametrize("key", LAYERS)
def test_example9_layers_load(r100: Run, key: str) -> None:
    df = load_precomputed(key, r100.out / "precomputed" / f"{key}.tsv", r100.reps)
    assert len(df) > 0, key
    assert header(r100.out / "precomputed" / f"{key}.tsv")[: len(PRECOMPUTED_COLUMNS[key])] == list(
        PRECOMPUTED_COLUMNS[key]
    )


def test_example9_config_yaml(r100: Run) -> None:
    path = r100.out / "config.yaml"
    cfg = resolve_config(path)  # validates against GIPPA's config model
    want = {"bakta_json": r100.out / "bakta_like.json"}
    want |= {k: r100.out / "precomputed" / f"{k}.tsv" for k in LAYERS}
    pre = cfg.annotation.precomputed
    for key, target in want.items():
        got = getattr(pre, key)
        assert got is not None, key
        got = Path(got) if Path(got).is_absolute() else path.parent / got
        assert got.resolve() == target.resolve(), key


@pytest.mark.parametrize("key", TRUTH_LAYERS)
def test_truth_lists_layer_rows(r100: Run, key: str) -> None:
    cols = list(PRECOMPUTED_COLUMNS[key])
    assert set(cols) | {"true_element"} <= set(header(r100.out / "truth" / f"{key}.tsv"))
    truth = r100.truth(key)
    layer = r100.layer(key)
    assert Counter(tuple(r[c] for c in cols) for r in truth) == Counter(
        tuple(r[c] for c in cols) for r in layer
    )
    assert all(r["true_element"] in ("true", "false") for r in truth)


def test_true_promoter_rows_are_promoter_tus(r100: Run) -> None:
    promoted = [r for r in tus(r100) if is_true(r["promoter"])]
    true_rows = [r for r in r100.truth("promoters") if is_true(r["true_element"])]
    assert len(true_rows) == len(promoted)
    sites: dict[tuple[str, str], list[int]] = {}
    for r in promoted:
        sites.setdefault((r["contig"], r["strand"]), []).append(int(r["tss"]))
    for r in true_rows:
        near = sites.get((r["contig"], r["strand"]), [])
        assert any(abs(int(r["end3p"]) - t) <= 50 for t in near), r


def test_tx_rates_cover_promoter_tss(r100: Run) -> None:
    have = {(r["contig"], r["strand"], r["tss"]) for r in r100.layer("tx_rates")}
    for r in tus(r100):
        if is_true(r["promoter"]):
            assert (r["contig"], r["strand"], r["tss"]) in have, r["tu_id"]


def test_true_terminators_follow_tu_ends(r100: Run) -> None:
    ends: dict[tuple[str, str], list[int]] = {}
    genes = {g.locus_tag: g for g in r100.ann.model.genes}
    for tu in tus(r100):
        ivs = [genes[x].interval for x in tu["genes"].split(";")]
        end = max(i.end for i in ivs) if tu["strand"] == "+" else min(i.start for i in ivs) + 1
        ends.setdefault((tu["contig"], tu["strand"]), []).append(end)
    true_rows = [r for r in r100.truth("terminators") if is_true(r["true_element"])]
    assert true_rows
    for r in true_rows:
        s, e = int(r["start"]), int(r["end"])
        cands = ends.get((r["contig"], r["strand"]), [])
        if r["strand"] == "+":
            assert any(0 <= s - x <= 150 for x in cands), r
        else:
            assert any(0 <= x - e <= 150 for x in cands), r
        seq = r100.slice(Interval(r["contig"], s - 1, e, r["strand"]))  # type: ignore[arg-type]
        assert "TTTT" in seq, (r, seq)


def test_operon_pairs(r100: Run) -> None:
    rows = r100.layer("operon_pairs")
    genes = {g.locus_tag: g for g in r100.ann.model.genes}
    for r in rows:
        a, b = genes[r["gene_a"]].interval, genes[r["gene_b"]].interval
        assert (a.contig, a.strand) == (b.contig, b.strand), r
        assert 0.0 <= float(r["uniop_p"]) <= 1.0
        assert r["source"] in ("uniop", "distance_model")
    p = {frozenset((r["gene_a"], r["gene_b"])): float(r["uniop_p"]) for r in rows}
    within = []
    for tu in tus(r100):
        m = tu["genes"].split(";")
        for x, y in itertools.pairwise(m):
            assert frozenset((x, y)) in p, (tu["tu_id"], x, y)
            within.append(frozenset((x, y)))
    assert within
    between = [v for k, v in p.items() if k not in set(within)]
    assert between
    mean_w = sum(p[k] for k in within) / len(within)
    assert mean_w > sum(between) / len(between) + 0.3


def test_regions_vocabulary(r100: Run) -> None:
    for r in r100.truth("regions"):
        assert r["region_type"] in REGION_TYPES, r


# ---------------------------------------------------------------- strength 0 / signals off


def test_zero_strength_plants_nothing(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    cfg = dataclasses.replace(cfg, rbs_strength=0.0, promoter_strength=0.0, terminator_strength=0.0)
    run = run_sim(cfg, seed, tmp_path)
    assert not any(is_true(r["promoter"]) for r in tus(run))
    for key in ("promoters", "terminators", "rbs"):
        assert not any(is_true(r["true_element"]) for r in run.truth(key)), key


def test_zero_false_positive_rate(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    run = run_sim(dataclasses.replace(cfg, layer_false_positive_rate=0.0), seed, tmp_path)
    for key in TRUTH_LAYERS:
        assert all(is_true(r["true_element"]) for r in run.truth(key)), key


def test_false_positives_present_at_high_rate(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    run = run_sim(dataclasses.replace(cfg, layer_false_positive_rate=0.5), seed, tmp_path)
    for key in ("promoters", "terminators", "rbs"):
        assert any(not is_true(r["true_element"]) for r in run.truth(key)), key


# ---------------------------------------------------------------- example 10: determinism


def test_example10_same_seed_byte_identical(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    simulate(cfg, seed, tmp_path / "a")
    simulate(cfg, seed, tmp_path / "b")
    a, b = tree(tmp_path / "a"), tree(tmp_path / "b")
    want = {"genome.fasta", "bakta_like.json", "annotation.gff3", "config.yaml"}
    want |= {f"precomputed/{k}.tsv" for k in LAYERS} | {f"truth/{k}.tsv" for k in TRUTH}
    assert want <= set(a)
    a_cfg, b_cfg = a.pop("config.yaml"), b.pop("config.yaml")
    assert a == b
    # config.yaml may hold absolute paths to its own directory; otherwise identical
    assert a_cfg.replace(str(tmp_path / "a").encode(), b"@") == b_cfg.replace(
        str(tmp_path / "b").encode(), b"@"
    )


def test_different_seed_differs(tmp_path: Path) -> None:
    cfg, seed = sim100k()
    simulate(cfg, seed, tmp_path / "a")
    simulate(cfg, seed + 1, tmp_path / "b")
    for name in ("genome.fasta", "truth/tus.tsv", "precomputed/promoters.tsv", "truth/rbs.tsv"):
        assert (tmp_path / "a" / name).read_bytes() != (tmp_path / "b" / name).read_bytes(), name


# ---------------------------------------------------------------- backward compatibility


def test_m0_defaults_write_no_layers(tmp_path: Path) -> None:
    res = simulate(SimConfig(), 1, tmp_path)
    assert not (tmp_path / "precomputed").exists()
    assert not dict(getattr(res, "precomputed", {}) or {})


def test_m318_output_unchanged(tmp_path: Path) -> None:
    simulate(SimConfig(**M318_CONFIG), 3, tmp_path)
    for name, digest in M318_SHA256.items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest, name
    assert model_digest(parse_bakta_json(tmp_path / "bakta_like.json")) == M318_MODEL_SHA256


def test_explicit_defaults_equal_implicit(tmp_path: Path) -> None:
    new = {
        "real_rrna_trna": False, "crispr": False, "oric": False, "prophage_bp": 0, "bgc_bp": 0,
        "rbs_strength": 0.0, "promoter_strength": 0.0, "terminator_strength": 0.0,
        "write_precomputed": False,
    }  # fmt: skip
    simulate(SimConfig(**M318_CONFIG), 5, tmp_path / "a")
    simulate(SimConfig(**M318_CONFIG, **new), 5, tmp_path / "b")
    assert tree(tmp_path / "a") == tree(tmp_path / "b")


# ---------------------------------------------------------------- properties


def small_cfg(**kw: Any) -> Any:
    base = {
        "genome_length": 60_000,
        "n_plasmids": 1,
        "plasmid_length": (6_000, 8_000),
        "n_rrna_operons": 1,
        "n_trna": 10,
        "real_rrna_trna": True,
        "crispr": True,
        "oric": True,
        "prophage_bp": 6_000,
        "bgc_bp": 5_000,
        "n_pseudogenes": 2,
        "n_sorfs": 2,
        "rbs_strength": 1.0,
        "promoter_strength": 1.0,
        "terminator_strength": 1.0,
        "write_precomputed": True,
    }
    return SimConfig(**(base | kw))


@settings(max_examples=5, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(seed=st.integers(min_value=0, max_value=2**32 - 1))
def test_property_full_strength_invariants(seed: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        run = run_sim(small_cfg(), seed, Path(d))
        assert dict(run.ann.skipped_types) == {}
        for g in run.ann.model.genes:
            assert rbs_spacer(run, g) is not None, g.locus_tag
            if g.feature_type != "pseudogene":
                nt = g.nt or ""
                assert nt[:3] in STARTS and nt[-3:] in STOPS and not internal_stop(nt)
        rows = tus(run)
        assert rows and all(is_true(r["promoter"]) for r in rows)
        for r in rows:
            check_promoter(run, r)
        for key in LAYERS:
            load_precomputed(key, Path(d) / "precomputed" / f"{key}.tsv", run.reps)
        m = run.ann.model
        for f in m.features:
            assert not [x for x in m.overlapping(f.interval) if hasattr(x, "locus_tag")], f.id


@settings(max_examples=4, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    seed=st.integers(min_value=0, max_value=2**32 - 1),
    strength=st.sampled_from([0.0, 0.3, 0.8]),
    fp=st.sampled_from([0.0, 0.05, 0.3]),
)
def test_property_deterministic_and_truth_consistent(seed: int, strength: float, fp: float) -> None:
    cfg = small_cfg(
        rbs_strength=strength,
        promoter_strength=strength,
        terminator_strength=strength,
        layer_false_positive_rate=fp,
    )
    with tempfile.TemporaryDirectory() as d:
        simulate(cfg, seed, Path(d) / "a")
        simulate(cfg, seed, Path(d) / "b")
        a, b = tree(Path(d) / "a"), tree(Path(d) / "b")
        a.pop("config.yaml", None), b.pop("config.yaml", None)
        assert a == b
        out = Path(d) / "a"
        for key in TRUTH_LAYERS:
            cols = list(PRECOMPUTED_COLUMNS[key])
            truth = read_rows(out / "truth" / f"{key}.tsv")
            layer = read_rows(out / "precomputed" / f"{key}.tsv")
            assert Counter(tuple(r[c] for c in cols) for r in truth) == Counter(
                tuple(r[c] for c in cols) for r in layer
            )
            if fp == 0.0:
                assert all(r["true_element"] == "true" for r in truth), key
        if strength == 0.0:
            assert not any(r["promoter"] == "true" for r in read_rows(out / "truth" / "tus.tsv"))
