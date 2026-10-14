"""giapsim genome features for mapping (issue M3-18; GIAP_SPEC.md §11.3).

Examples 1-10 of docs/specs/M3-18-sim-genome-features.md, the round trip through M3-01's
``gene_model_from_bakta_json`` and M1-05's ``read_gff3``, truth tables, determinism and speed.

Conventions these tests pin down where the spec is silent (see the PR description):

* a JSON feature's ID is its ``locus`` (falling back to ``id``); ``truth/features.tsv`` rows use it;
* rRNA subunits are recognised by ``product``/``gene``/``name`` (``16S``/``rrs``, ``23S``/``rrl``,
  ``5S``/``rrf``); an operon is 16S, 23S, 5S on one strand in that transcription order and holds
  no tRNA, so ``n_trna`` is the exact tRNA count;
* ``true_start`` is the 1-based genome coordinate of the first base of the true start codon in
  transcription direction (``< start`` on ``+``, ``> end`` on ``-``); empty cells are ``""``,
  ``NA`` or ``.``;
* ``truth/strain_differences.tsv`` follows VCF orientation: ``ref`` is the allele in
  ``reference.fasta`` at ``pos_reference``, ``alt`` the allele in ``genome.fasta`` (the sample) at
  ``pos_sample``; ``ins`` means ``alt`` is longer than ``ref``. Indels are left-normalized and
  anchored VCF-style (``pos`` is the base before the event, both alleles include it; review
  decision 3, pinned in test_sim_features_review.py); ``-`` or ``.`` (an empty allele) are still
  tolerated by the helpers here.
"""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.io.gff3 import read_gff3
from giap.model.coords import Interval, Replicon, split_at_origin
from giap.model.genes import Gene, GeneModel, NcFeature, gene_model_from_bakta_json

BASES = "TCAG"
AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {
    a + b + c: AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}
STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
ANTI_SD_TAIL = "GATCACCTCCTTA"
EMPTY = {"", "NA", "."}
COMP = str.maketrans("ACGT", "TGCA")

#: SHA-256 of the M0 simulator's output for ``SimConfig()`` and seed 1 (computed on main at
#: d048f27 before M3-18). New fields default to off, so these must not change (example 1).
M0_SHA256 = {
    "genome.fasta": "698373a45c787caa64c01d7050dd1d2a7b31e03771e4e389056ca510c6df1971",
    "bakta_like.json": "d312420eb3d404cb1abbb6d07844cccbc83c0ce5e9c912d3a02aa17da6c49064",
    "truth/genes.tsv": "9c2c98152eed218d3855ccddf478bf49ea42bc69727fce466111f704e0387b0f",
}

FEATURES_COLUMNS = [
    "feature_id",
    "type",
    "contig",
    "start",
    "end",
    "strand",
    "true_start",
    "pseudo_reason",
]
NOVEL_COLUMNS = ["orf_id", "class", "contig", "start", "end", "strand", "aa"]
STRAIN_COLUMNS = [
    "contig",
    "pos_sample",
    "pos_reference",
    "ref",
    "alt",
    "type",
    "gene",
    "aa_change",
]


# ---------------------------------------------------------------- helpers


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def translate(nt: str) -> str:
    return "".join(CODE[nt[i : i + 3]] for i in range(0, len(nt) - len(nt) % 3, 3))


def protein(nt: str) -> str:
    """Translation of a CDS without its stop, initiator forced to M (as in the M0 JSON)."""
    aa = translate(nt)
    return "M" + aa[1:-1] if aa.endswith("*") else "M" + aa[1:]


def internal_stops(nt: str) -> list[int]:
    return [i for i in range(3, len(nt) - 3, 3) if nt[i : i + 3] in STOPS]


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """``{id: (header, sequence)}`` in file order."""
    out: dict[str, tuple[str, str]] = {}
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                out[header.split()[0]] = (header, "".join(chunks))
            header, chunks = line[1:], []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        out[header.split()[0]] = (header, "".join(chunks))
    return out


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def ftype(f: dict[str, Any]) -> str:
    return str(f["type"]).lower()


def fid(f: dict[str, Any]) -> str:
    return str(f.get("locus") or f.get("id") or "")


def fcontig(f: dict[str, Any]) -> str:
    return str(f.get("sequence", f.get("contig")))


def by_type(data: dict[str, Any], t: str) -> list[dict[str, Any]]:
    return [f for f in data["features"] if ftype(f) == t]


def is_pseudo(f: dict[str, Any]) -> bool:
    return bool(f.get("pseudo")) or "pseudogene" in f


def rrna_kind(f: dict[str, Any]) -> str | None:
    text = " ".join(str(f.get(k) or "") for k in ("product", "gene", "name", "subtype"))
    for kind, keys in (("16S", ("16S", "rrs")), ("23S", ("23S", "rrl")), ("5S", ("5S", "rrf"))):
        if any(k in text for k in keys):
            return kind
    return None


def seq_of(f: dict[str, Any], genomes: dict[str, str]) -> str:
    """Feature sequence on its own strand; Bakta-style ``start > stop`` wraps the origin."""
    g = genomes[fcontig(f)]
    s, e = int(f["start"]), int(f["stop"])
    seg = g[s - 1 :] + g[:e] if s > e else g[s - 1 : e]
    return revcomp(seg) if f["strand"] == "-" else seg


def replicons_of(data: dict[str, Any]) -> dict[str, Replicon]:
    out = {}
    for s in data["sequences"]:
        kind = s.get("type") if s.get("type") in ("chromosome", "plasmid") else "unknown"
        out[s["id"]] = Replicon(s["id"], int(s["length"]), s["topology"], kind)
    return out


def overlap_bp(a: Interval, b: Interval, reps: dict[str, Replicon]) -> int:
    if a.contig != b.contig:
        return 0
    rep = reps[a.contig]
    return sum(
        max(0, min(x.end, y.end) - max(x.start, y.start))
        for x in split_at_origin(a, rep)
        for y in split_at_origin(b, rep)
    )


class Run:
    """One simulation: parsed JSON, genome sequences, replicons and the GIAP gene model."""

    def __init__(self, out: Path) -> None:
        self.out = out
        self.data: dict[str, Any] = json.loads((out / "bakta_like.json").read_text())
        self.fasta = read_fasta(out / "genome.fasta")
        self.genomes = {k: v[1] for k, v in self.fasta.items()}
        self.reps = replicons_of(self.data)
        self.model: GeneModel = gene_model_from_bakta_json(out / "bakta_like.json", self.reps)

    def truth(self, name: str) -> tuple[list[str], list[dict[str, str]]]:
        return read_tsv(self.out / "truth" / name)

    def coding(self) -> list[dict[str, Any]]:
        return [f for f in self.data["features"] if ftype(f) in ("cds", "sorf")]

    def feature(self, ident: str) -> dict[str, Any]:
        hits = [f for f in self.data["features"] if fid(f) == ident]
        assert len(hits) == 1, (ident, len(hits))
        return hits[0]

    def interval(self, f: dict[str, Any]) -> Interval:
        rep = self.reps[fcontig(f)]
        s, e = int(f["start"]), int(f["stop"])
        end = e + rep.length if s > e else e
        return Interval(fcontig(f), s - 1, end, f["strand"])


def strain_config(**kw: Any) -> Any:
    import importlib

    import giapsim

    cls = getattr(giapsim, "StrainConfig", None)
    if cls is None:
        cls = importlib.import_module("giapsim.strain").StrainConfig
    return cls(**kw)


def simulate(out: Path, seed: int = 1, **kw: Any) -> Run:
    import giapsim

    giapsim.simulate(giapsim.SimConfig(**kw), seed=seed, outdir=out)
    return Run(out)


def all_options(genome_length: int) -> dict[str, Any]:
    return {
        "genome_length": genome_length,
        "n_plasmids": 2,
        "plasmid_length": (6_000, 12_000),
        "linear_contig_length": 20_000,
        "origin_crossing_gene": True,
        "n_rrna_operons": 2,
        "n_trna": 8,
        "n_pseudogenes": 3,
        "n_sorfs": 4,
        "n_nterm_extensions": 3,
        "n_novel_orfs": (("intergenic", 2), ("cds_outframe", 2), ("antisense", 2)),
        "strain": strain_config(snv_rate_per_kb=1.0, indel_rate_per_kb=0.2),
    }


@pytest.fixture(scope="module")
def mk(tmp_path_factory: pytest.TempPathFactory) -> Any:
    cache: dict[str, Run] = {}

    def make(name: str, seed: int = 1, **kw: Any) -> Run:
        if name not in cache:
            cache[name] = simulate(tmp_path_factory.mktemp(name), seed, **kw)
        return cache[name]

    return make


@pytest.fixture(scope="module")
def full(mk: Any) -> Run:
    run: Run = mk("full", seed=3, **all_options(100_000))
    return run


# ---------------------------------------------------------------- example 1: M0 unchanged


def test_example1_defaults_byte_identical_to_m0(tmp_path: Path) -> None:
    import giapsim

    giapsim.simulate(giapsim.SimConfig(), seed=1, outdir=tmp_path)
    for name, digest in M0_SHA256.items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest, name


def test_new_config_fields_default_off() -> None:
    import giapsim

    cfg = giapsim.SimConfig()
    assert cfg.n_plasmids == 0
    assert tuple(cfg.plasmid_length) == (5_000, 60_000)
    assert cfg.linear_contig_length == 0
    assert cfg.origin_crossing_gene is False
    assert (cfg.n_rrna_operons, cfg.n_trna, cfg.n_pseudogenes, cfg.n_sorfs) == (0, 0, 0, 0)
    assert cfg.n_nterm_extensions == 0
    assert dict(cfg.n_novel_orfs) == {}
    assert cfg.strain is None


# ---------------------------------------------------------------- example 2: sORFs


def test_example2_sorfs(mk: Any) -> None:
    run = mk("sorfs", n_sorfs=3)
    sorfs = by_type(run.data, "sorf")
    assert len(sorfs) == 3
    for f in sorfs:
        nt = seq_of(f, run.genomes)
        assert f["nt"] == nt, fid(f)
        assert len(nt) % 3 == 0 and len(nt) == int(f["stop"]) - int(f["start"]) + 1
        assert nt[:3] in STARTS and nt[-3:] in STOPS, fid(f)
        assert not internal_stops(nt), fid(f)
        assert len(f["aa"].rstrip("*")) < 30, fid(f)
        assert f["aa"].rstrip("*") == protein(nt), fid(f)
    genes = [g for g in run.model.genes if g.feature_type == "sorf"]
    assert len(genes) == 3
    _, rows = run.truth("features.tsv")
    assert sorted(r["feature_id"] for r in rows if r["type"].lower() == "sorf") == sorted(
        fid(f) for f in sorfs
    )


# ---------------------------------------------------------------- example 3: pseudogenes


def test_example3_pseudogenes(mk: Any) -> None:
    run = mk("pseudo", n_pseudogenes=2)
    pseudo = [f for f in by_type(run.data, "cds") if is_pseudo(f)]
    assert len(pseudo) == 2
    assert all(f.get("pseudo") is True for f in pseudo)
    _, rows = run.truth("features.tsv")
    reasons = {r["feature_id"]: r["pseudo_reason"] for r in rows if r["pseudo_reason"] not in EMPTY}
    assert sorted(reasons) == sorted(fid(f) for f in pseudo)
    for f in pseudo:
        nt = seq_of(f, run.genomes)
        assert f["nt"] == nt, fid(f)
        reason = reasons[fid(f)]
        assert reason in {"internal_stop", "frameshift"}
        if reason == "internal_stop":
            assert len(nt) % 3 == 0 and internal_stops(nt), fid(f)
        else:
            assert len(nt) % 3 != 0, fid(f)  # a 1-nt indel leaves the ORF out of frame
    assert sum(g.feature_type == "pseudogene" for g in run.model.genes) == 2


def test_pseudogenes_cover_both_reasons(mk: Any) -> None:
    run = mk("pseudo_many", seed=4, genome_length=200_000, n_pseudogenes=12)
    _, rows = run.truth("features.tsv")
    reasons = Counter(r["pseudo_reason"] for r in rows if r["pseudo_reason"] not in EMPTY)
    assert sum(reasons.values()) == 12
    assert set(reasons) == {"internal_stop", "frameshift"}


# ---------------------------------------------------------------- example 4: rRNA / tRNA


def test_example4_rrna_operon(mk: Any) -> None:
    run = mk("rrna", n_rrna_operons=1)
    rrnas = by_type(run.data, "rrna")
    kinds = Counter(rrna_kind(f) for f in rrnas)
    assert kinds == Counter({"16S": 1, "23S": 1, "5S": 1})
    ssu = next(f for f in rrnas if rrna_kind(f) == "16S")
    assert seq_of(ssu, run.genomes).endswith(ANTI_SD_TAIL)
    assert len({f["strand"] for f in rrnas}) == 1
    order = sorted(rrnas, key=lambda f: int(f["start"]), reverse=rrnas[0]["strand"] == "-")
    assert [rrna_kind(f) for f in order] == ["16S", "23S", "5S"]
    assert sum(x.type == "rrna" for x in run.model.features) == 3
    assert not by_type(run.data, "trna")


def test_rrna_operon_count_and_16s_tail(mk: Any) -> None:
    run = mk("rrna3", seed=2, genome_length=200_000, n_rrna_operons=3)
    rrnas = by_type(run.data, "rrna")
    assert Counter(rrna_kind(f) for f in rrnas) == Counter({"16S": 3, "23S": 3, "5S": 3})
    for f in rrnas:
        if rrna_kind(f) == "16S":
            assert seq_of(f, run.genomes).endswith(ANTI_SD_TAIL), fid(f)


def test_trna_genes_with_anticodons(mk: Any) -> None:
    run = mk("trna", n_trna=6)
    trnas = by_type(run.data, "trna")
    assert len(trnas) == 6
    for f in trnas:
        ac = str(f.get("anti_codon", ""))
        assert len(ac) == 3 and set(ac) <= set("ACGT"), fid(f)
        assert ac in seq_of(f, run.genomes), fid(f)
        assert f["strand"] in {"+", "-"}
    assert sum(x.type == "trna" for x in run.model.features) == 6
    _, rows = run.truth("features.tsv")
    assert sum(r["type"].lower() == "trna" for r in rows) == 6


# ---------------------------------------------------------------- example 5: N-terminal extensions


def test_example5_nterm_extensions(mk: Any) -> None:
    run = mk("nterm", n_nterm_extensions=2)
    _, rows = run.truth("features.tsv")
    ext = [r for r in rows if r["true_start"] not in EMPTY]
    assert len(ext) == 2
    for r in ext:
        f = run.feature(r["feature_id"])
        assert ftype(f) == "cds" and not is_pseudo(f)
        g = run.genomes[fcontig(f)]
        ts, s, e = int(r["true_start"]), int(f["start"]), int(f["stop"])
        if f["strand"] == "+":
            assert ts < s
            full_nt = g[ts - 1 : e]
        else:
            assert ts > e
            full_nt = revcomp(g[s - 1 : ts])
        extra = len(full_nt) - len(f["nt"])
        assert extra % 3 == 0 and 5 <= extra // 3 <= 30, (r["feature_id"], extra)
        assert full_nt.endswith(f["nt"])
        assert full_nt[:3] in STARTS and f["nt"][:3] in STARTS
        assert not internal_stops(full_nt), r["feature_id"]
        assert protein(full_nt).endswith(f["aa"].rstrip("*")[1:])


# ---------------------------------------------------------------- example 6: origin-crossing gene


def test_example6_origin_crossing_gene(mk: Any) -> None:
    run = mk("origin", origin_crossing_gene=True)
    wrapped = [f for f in run.data["features"] if int(f["start"]) > int(f["stop"])]
    assert len(wrapped) == 1
    f = wrapped[0]
    assert ftype(f) == "cds" and not is_pseudo(f)
    rep = run.reps[fcontig(f)]
    assert rep.topology == "circular"
    nt = seq_of(f, run.genomes)
    assert f["nt"] == nt
    assert nt[:3] in STARTS and nt[-3:] in STOPS and not internal_stops(nt)
    gene = run.model.by_locus(fid(f))
    assert gene.interval.end > rep.length
    assert gene.interval.crosses_origin(rep)
    assert gene.interval.start == int(f["start"]) - 1
    assert gene.interval.end == int(f["stop"]) + rep.length


# ---------------------------------------------------------------- example 7: linear contig


def test_example7_linear_contig(mk: Any) -> None:
    run = mk("linear", linear_contig_length=20_000)
    seqs = run.data["sequences"]
    assert len(seqs) == 2
    assert seqs[0]["topology"] == "circular"
    lin = seqs[1]
    assert lin["topology"] == "linear" and int(lin["length"]) == 20_000
    assert len(run.genomes[lin["id"]]) == 20_000
    assert "[topology=circular]" not in run.fasta[lin["id"]][0]
    on_lin = [f for f in run.data["features"] if fcontig(f) == lin["id"]]
    for f in on_lin:
        assert 1 <= int(f["start"]) <= int(f["stop"]) <= 20_000, fid(f)


# ---------------------------------------------------------------- plasmids


def test_plasmids(mk: Any) -> None:
    run = mk("plasmids", n_plasmids=2, plasmid_length=(5_000, 9_000))
    plasmids = [s for s in run.data["sequences"] if s.get("type") == "plasmid"]
    assert len(plasmids) == 2
    assert len(run.data["sequences"]) == 3 and len(run.fasta) == 3
    assert len({s["id"] for s in run.data["sequences"]}) == 3
    for s in plasmids:
        assert s["topology"] == "circular"
        assert 5_000 <= int(s["length"]) <= 9_000
        assert len(run.genomes[s["id"]]) == int(s["length"])
        assert "[topology=circular]" in run.fasta[s["id"]][0]
        assert run.reps[s["id"]].kind == "plasmid"
    assert run.data["sequences"][0]["type"] == "chromosome"


# ---------------------------------------------------------------- example 8: novel ORFs


def _novel(run: Run, cls: str) -> list[dict[str, str]]:
    header, rows = run.truth("novel_orfs.tsv")
    assert header == NOVEL_COLUMNS
    return [r for r in rows if r["class"] == cls]


def _check_orf(run: Run, r: dict[str, str]) -> tuple[Interval, str]:
    g = run.genomes[r["contig"]]
    s, e = int(r["start"]), int(r["end"])
    assert 1 <= s < e <= len(g)
    nt = g[s - 1 : e] if r["strand"] == "+" else revcomp(g[s - 1 : e])
    assert len(nt) % 3 == 0 and nt[:3] in STARTS and nt[-3:] in STOPS, r["orf_id"]
    assert not internal_stops(nt), r["orf_id"]
    assert r["aa"].rstrip("*")[1:] == protein(nt)[1:], r["orf_id"]
    annotated = {(fcontig(f), int(f["start"]), int(f["stop"])) for f in run.data["features"]}
    assert (r["contig"], s, e) not in annotated, r["orf_id"]
    assert all(r["aa"].rstrip("*") != f.get("aa", "").rstrip("*") for f in run.coding())
    return Interval(r["contig"], s - 1, e, r["strand"]), nt  # type: ignore[arg-type]


def test_example8_antisense_orfs(mk: Any) -> None:
    run = mk("antisense", n_novel_orfs=(("antisense", 2),))
    rows = _novel(run, "antisense")
    assert len(rows) == 2 and len(run.truth("novel_orfs.tsv")[1]) == 2
    for r in rows:
        iv, _ = _check_orf(run, r)
        hits = [
            f
            for f in by_type(run.data, "cds")
            if f["strand"] != r["strand"] and overlap_bp(iv, run.interval(f), run.reps)
        ]
        assert hits, r["orf_id"]


def test_novel_orf_classes(mk: Any) -> None:
    run = mk(
        "novel", seed=5, n_novel_orfs=(("intergenic", 3), ("cds_outframe", 3), ("antisense", 2))
    )
    _, all_rows = run.truth("novel_orfs.tsv")
    assert Counter(r["class"] for r in all_rows) == {
        "intergenic": 3,
        "cds_outframe": 3,
        "antisense": 2,
    }
    assert len({r["orf_id"] for r in all_rows}) == len(all_rows)
    for r in _novel(run, "intergenic"):
        iv, _ = _check_orf(run, r)
        touching = [f for f in run.data["features"] if overlap_bp(iv, run.interval(f), run.reps)]
        assert not touching, r["orf_id"]
    for r in _novel(run, "cds_outframe"):
        iv, _ = _check_orf(run, r)
        ok = False
        for f in by_type(run.data, "cds"):
            civ = run.interval(f)
            if f["strand"] != r["strand"] or not overlap_bp(iv, civ, run.reps):
                continue
            shift = iv.start - civ.start if r["strand"] == "+" else civ.end - iv.end
            ok |= shift % 3 != 0
        assert ok, r["orf_id"]


# ---------------------------------------------------------------- example 9: strain differences


def _strain_rows(run: Run) -> list[dict[str, str]]:
    header, rows = run.truth("strain_differences.tsv")
    assert header == STRAIN_COLUMNS
    return rows


def _allele(a: str) -> str:
    return "" if a in {"-", "."} else a


def test_example9_strain_snvs(mk: Any) -> None:
    run = mk(
        "snv",
        genome_length=100_000,
        strain=strain_config(snv_rate_per_kb=1.0, indel_rate_per_kb=0.0),
    )
    rows = _strain_rows(run)
    assert 60 <= len(rows) <= 140
    assert {r["type"] for r in rows} == {"snv"}
    ref = {k: v[1] for k, v in read_fasta(run.out / "reference.fasta").items()}
    assert set(ref) == set(run.genomes)
    expected: dict[tuple[str, int], tuple[str, str]] = {}
    for r in rows:
        assert r["pos_sample"] == r["pos_reference"]
        assert len(r["ref"]) == len(r["alt"]) == 1 and r["ref"] != r["alt"]
        expected[(r["contig"], int(r["pos_sample"]))] = (r["ref"], r["alt"])
    assert len(expected) == len(rows)
    diffs = {}
    for c, s in run.genomes.items():
        assert len(ref[c]) == len(s)
        for i, (a, b) in enumerate(zip(ref[c], s, strict=True)):
            if a != b:
                diffs[(c, i + 1)] = (a, b)
    assert diffs == expected
    for r in rows:  # ``gene`` names the CDS hit by a SNV and is empty outside genes
        pos0 = int(r["pos_sample"]) - 1
        point = Interval(r["contig"], pos0, pos0 + 1, "+")
        hits = [
            g
            for g in run.model.overlapping(point)
            if isinstance(g, Gene) and g.feature_type == "cds"
        ]
        if not hits:
            assert r["gene"] in EMPTY
        elif len(hits) == 1:
            assert r["gene"] == hits[0].locus_tag


def test_strain_reference_annotation_is_lifted(mk: Any) -> None:
    run = mk(
        "indel",
        seed=6,
        genome_length=100_000,
        strain=strain_config(snv_rate_per_kb=0.0, indel_rate_per_kb=0.5),
    )
    rows = _strain_rows(run)
    assert rows and {r["type"] for r in rows} <= {"ins", "del"}
    ref = {k: v[1] for k, v in read_fasta(run.out / "reference.fasta").items()}
    shift = Counter[str]()
    for r in rows:
        a, b = _allele(r["ref"]), _allele(r["alt"])
        assert (r["type"] == "ins") == (len(b) > len(a))
        assert (r["type"] == "del") == (len(b) < len(a))
        p_ref, p_s = int(r["pos_reference"]) - 1, int(r["pos_sample"]) - 1
        assert ref[r["contig"]][p_ref : p_ref + len(a)] == a
        assert run.genomes[r["contig"]][p_s : p_s + len(b)] == b
        shift[r["contig"]] += len(a) - len(b)
    for c, s in run.genomes.items():
        assert len(ref[c]) - len(s) == shift[c], c
    rdata = json.loads((run.out / "reference_bakta_like.json").read_text())
    rreps = replicons_of(rdata)
    assert {k: r.length for k, r in rreps.items()} == {k: len(v) for k, v in ref.items()}
    gene_model_from_bakta_json(run.out / "reference_bakta_like.json", rreps)
    for f in rdata["features"]:
        if ftype(f) in ("cds", "sorf"):
            assert f["nt"] == seq_of(f, ref), fid(f)


# ---------------------------------------------------------------- round trip with GIAP's readers


def test_full_run_read_by_gene_model(full: Run) -> None:
    coding = full.coding()
    assert len(full.model.genes) == len(coding)
    assert sum(g.feature_type == "sorf" for g in full.model.genes) == 4
    assert sum(g.feature_type == "pseudogene" for g in full.model.genes) == 3
    assert sum(x.type == "rrna" for x in full.model.features) == 6
    assert sum(x.type == "trna" for x in full.model.features) == 8
    for f in coding:
        g = full.model.by_locus(fid(f))
        assert g.interval == full.interval(f)
        assert g.aa == f["aa"].rstrip("*")
    kinds = Counter(r.kind for r in full.reps.values())
    assert kinds["chromosome"] == 1 and kinds["plasmid"] == 2
    assert sum(r.topology == "linear" for r in full.reps.values()) == 1


def test_full_run_json_sequences_match_fasta(full: Run) -> None:
    assert [s["id"] for s in full.data["sequences"]] == list(full.fasta)
    for s in full.data["sequences"]:
        assert int(s["length"]) == len(full.genomes[s["id"]])


def test_full_run_gff3_readable(full: Run) -> None:
    gff = read_gff3(full.out / "annotation.gff3", full.reps)
    have = {(f.interval.contig, f.interval.start, f.interval.end, f.interval.strand) for f in gff}
    for f in full.data["features"]:
        if ftype(f) in ("cds", "sorf", "rrna", "trna"):
            iv = full.interval(f)
            assert (iv.contig, iv.start, iv.end, iv.strand) in have, fid(f)


def test_full_run_features_truth(full: Run) -> None:
    header, rows = full.truth("features.tsv")
    assert header[: len(FEATURES_COLUMNS)] == FEATURES_COLUMNS  # + ref_pseudo (review decision 2)
    ids = [r["feature_id"] for r in rows]
    assert len(ids) == len(set(ids))
    types = Counter(r["type"].lower() for r in rows)
    assert types["rrna"] == 6 and types["trna"] == 8 and types["sorf"] == 4
    assert sum(r["pseudo_reason"] not in EMPTY for r in rows) == 3
    assert sum(r["true_start"] not in EMPTY for r in rows) == 3
    for r in rows:
        f = full.feature(r["feature_id"])
        assert r["contig"] == fcontig(f) and r["strand"] == f["strand"]
        iv = full.interval(f)  # ``end`` unwrapped (> length across the origin; decision 5)
        assert (int(r["start"]), int(r["end"])) == (iv.start + 1, iv.end)


def test_full_run_novel_orfs_and_strain_tables(full: Run) -> None:
    _, rows = full.truth("novel_orfs.tsv")
    assert Counter(r["class"] for r in rows) == {"intergenic": 2, "cds_outframe": 2, "antisense": 2}
    assert _strain_rows(full)
    assert (full.out / "reference.fasta").is_file()
    assert (full.out / "reference_bakta_like.json").is_file()


# ---------------------------------------------------------------- determinism


def _tree(root: Path) -> dict[str, bytes]:
    return {
        str(p.relative_to(root)): p.read_bytes() for p in sorted(root.rglob("*")) if p.is_file()
    }


def test_same_seed_byte_identical_directories(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(**all_options(50_000))
    giapsim.simulate(cfg, seed=7, outdir=tmp_path / "a")
    giapsim.simulate(cfg, seed=7, outdir=tmp_path / "b")
    a, b = _tree(tmp_path / "a"), _tree(tmp_path / "b")
    assert {"genome.fasta", "bakta_like.json", "annotation.gff3", "reference.fasta"} <= set(a)
    assert {"truth/features.tsv", "truth/novel_orfs.tsv", "truth/strain_differences.tsv"} <= set(a)
    assert a == b


def test_different_seed_differs(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(**all_options(50_000))
    giapsim.simulate(cfg, seed=7, outdir=tmp_path / "a")
    giapsim.simulate(cfg, seed=8, outdir=tmp_path / "b")
    for name in ("genome.fasta", "truth/features.tsv", "truth/strain_differences.tsv"):
        assert (tmp_path / "a" / name).read_bytes() != (tmp_path / "b" / name).read_bytes()


# ---------------------------------------------------------------- properties


def _check_invariants(run: Run) -> None:
    items: list[Gene | NcFeature] = [*run.model.genes, *run.model.features]
    for x in run.model.features:  # rRNA/tRNA never overlap anything else
        others = [y for y in run.model.overlapping(x.interval) if y is not x]
        assert not others, (x.id, [getattr(y, "locus_tag", getattr(y, "id", "")) for y in others])
    genes = list(run.model.genes)
    for i, a in enumerate(genes):  # coding features keep M0's <= 20 bp overlaps
        for b in genes[i + 1 :]:
            assert overlap_bp(a.interval, b.interval, run.reps) <= 20, (a.locus_tag, b.locus_tag)
    for f in run.coding():
        nt = seq_of(f, run.genomes)
        assert f["nt"] == nt, fid(f)
        if is_pseudo(f):
            continue
        assert len(nt) % 3 == 0 and nt[:3] in STARTS and nt[-3:] in STOPS, fid(f)
        assert not internal_stops(nt), fid(f)
        assert f["aa"].rstrip("*") == protein(nt), fid(f)
    assert len(items) == len(run.coding()) + len(run.model.features)


@settings(
    max_examples=6,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
@given(seed=st.integers(min_value=0, max_value=2**32 - 1))
def test_property_no_unrequested_overlaps_and_clean_cds(seed: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        _check_invariants(simulate(Path(d), seed, **all_options(60_000)))


def test_invariants_on_full_run(full: Run) -> None:
    _check_invariants(full)


# ---------------------------------------------------------------- speed (example 10)


def test_100kb_all_options_under_1s(tmp_path: Path) -> None:
    import giapsim

    cfg = giapsim.SimConfig(**all_options(100_000))
    best = float("inf")
    for i in range(3):
        t0 = time.perf_counter()
        giapsim.simulate(cfg, seed=9, outdir=tmp_path / str(i))
        best = min(best, time.perf_counter() - t0)
    # loose 10x guard; the 1 s budget is enforced in benchmarks/test_giapsim.py
    assert best < 10.0, best


@pytest.mark.slow
def test_example10_5mb_all_options_under_30s(tmp_path: Path) -> None:
    import giapsim

    opts = all_options(5_000_000) | {"n_rrna_operons": 7, "n_trna": 86, "n_sorfs": 40}
    cfg = giapsim.SimConfig(**opts)
    t0 = time.perf_counter()
    giapsim.simulate(cfg, seed=10, outdir=tmp_path)
    # loose 10x guard; the 30 s budget is enforced in benchmarks/test_giapsim.py
    assert time.perf_counter() - t0 < 300.0
    run = Run(tmp_path)
    assert len(run.genomes[run.data["sequences"][0]["id"]]) == 5_000_000
    assert sum(x.type == "trna" for x in run.model.features) == 86
