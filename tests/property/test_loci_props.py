"""Properties of locus classification on gippasim genomes (issue M3-14; GIPPA_SPEC.md §7.2.6).

Acceptance properties of docs/specs/M3-14-loci-and-assignment.md:

- every peptide window cut in frame from inside a CDS classifies ``cds_inframe`` with that gene
  (origin-crossing genes included); the same window on the other strand classifies ``antisense``;
- classification is invariant to the order of genes and features in the model;

plus invariants that follow from the normative definitions:

- mirror symmetry: reverse-complementing the replicon and mirroring the model and the hit leaves
  (class, gene) unchanged (the plus and minus rules are mirror images);
- agreement with simulator truth (M3-18 / M2-05): planted sORFs → ``sorf``, pseudogenes →
  ``pseudogene``, rRNA/tRNA → ``rna_overlap``, the region between the true and the annotated start
  of an N-terminal extension → ``nterm_extension`` with that gene, and novel ORFs (interior, stop
  codon excluded) → their truth class (``intergenic``, ``cds_outframe``, ``antisense``). The spec
  says "every", so agreement must be 100 %; windows touching a second annotated feature are skipped
  because priorities legitimately change their class;
- ``annotate_peptide_hits`` is invariant to the row order of its input.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from gippa.mapping.hits import HIT_COLUMNS
from gippa.model.coords import Interval, Replicon, split_at_origin
from gippa.model.genes import Gene, GeneModel, NcFeature, gene_model_from_bakta_json

PROPS = settings(
    max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)  # fmt: skip
COMP = str.maketrans("ACGTN", "TGCAN")
RNA = {"rrna", "trna", "tmrna", "ncrna"}


class Sim:
    def __init__(self, out: Path) -> None:
        data = json.loads((out / "bakta_like.json").read_text())
        self.out = out
        self.reps = {
            s["id"]: Replicon(
                s["id"],
                int(s["length"]),
                s["topology"],
                s.get("type") if s.get("type") in ("chromosome", "plasmid") else "unknown",
            )
            for s in data["sequences"]
        }
        self.genome = _fasta(out / "genome.fasta")
        self.model: GeneModel = gene_model_from_bakta_json(out / "bakta_like.json", self.reps)

    def truth(self, name: str) -> list[dict[str, str]]:
        with (self.out / "truth" / name).open(newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif line:
            out[cur].append(line.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def _simulate(out: Path, seed: int, **kw: Any) -> Sim:
    import gippasim

    gippasim.simulate(gippasim.SimConfig(**kw), seed=seed, outdir=out)
    return Sim(out)


FEATURES_KW: dict[str, Any] = {
    "genome_length": 100_000,
    "n_plasmids": 1,
    "plasmid_length": (8_000, 12_000),
    "linear_contig_length": 8_000,
    "origin_crossing_gene": True,
    "n_rrna_operons": 1,
    "n_trna": 4,
    "n_pseudogenes": 4,
    "n_sorfs": 4,
    "n_nterm_extensions": 4,
    "n_novel_orfs": (("intergenic", 3), ("cds_outframe", 3), ("antisense", 3)),
}


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Sim:
    return _simulate(tmp_path_factory.mktemp("loci_sim"), 14, **FEATURES_KW)


@pytest.fixture(scope="module")
def sim_plain(tmp_path_factory: pytest.TempPathFactory) -> Sim:
    """No origin-crossing feature: used by the mirror property."""
    kw = {k: v for k, v in FEATURES_KW.items() if k != "origin_crossing_gene"}
    return _simulate(tmp_path_factory.mktemp("loci_plain"), 3, **kw)


def classify(i: Interval, gm: GeneModel, genome: dict[str, str]) -> tuple[str, str | None]:
    from gippa.mapping.loci import classify_locus

    cls, gene = classify_locus(i, gm, genome)
    return str(cls), (None if gene is None or pd.isna(gene) else str(gene))


def others(gm: GeneModel, i: Interval, me: str | None) -> list[Gene | NcFeature]:
    """Annotated items overlapping ``i`` other than ``me`` (any strand)."""
    out = []
    for x in gm.overlapping(i):
        tag = x.locus_tag if isinstance(x, Gene) else x.id
        if tag != me:
            out.append(x)
    return out


def window(g: Gene, rep: Replicon, off_codons: int, n_codons: int) -> Interval:
    """A window of ``n_codons`` in frame inside ``g``, ``off_codons`` from its 5' end."""
    gs, ge = g.interval.start, g.interval.end
    s = gs + 3 * off_codons if g.interval.strand == "+" else ge - 3 * (off_codons + n_codons)
    e = s + 3 * n_codons
    if s >= rep.length:
        s, e = s - rep.length, e - rep.length
    return Interval(g.interval.contig, s, e, g.interval.strand)


def flip(i: Interval) -> Interval:
    return Interval(i.contig, i.start, i.end, "-" if i.strand == "+" else "+")


STOPS = {"TAA", "TAG", "TGA"}


def _rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def _nterm_possible(i: Interval, sim: Sim) -> bool:
    """Oracle: could ``i`` be an N-terminal extension of some CDS on its strand (no stop)?"""
    rep = sim.reps[i.contig]
    seq, n = sim.genome[i.contig], rep.length
    circ = rep.topology == "circular"
    for g in sim.model.genes:
        if g.feature_type != "cds" or g.interval.contig != i.contig:
            continue
        if g.interval.strand != i.strand:
            continue
        for d in (0, n, -n) if circ else (0,):
            s, e = i.start + d, i.end + d
            if i.strand == "+":
                gs = g.interval.start
                if not (s < gs and (gs - s) % 3 == 0 and gs - s <= n):
                    continue
                codons = [seq[(x + k) % n] for x in range(e, gs, 3) for k in range(3)]
            else:
                ge = g.interval.end
                if not (e > ge and (e - ge) % 3 == 0 and e - ge <= n):
                    continue
                codons = [seq[(x + k) % n] for x in range(ge, s, 3) for k in range(3)]
                codons = list(_rc("".join(codons)))
            nt = "".join(codons)
            if not any(nt[k : k + 3] in STOPS for k in range(0, len(nt) - 2, 3)):
                return True
    return False


# --- acceptance properties --------------------------------------------------------------------


@PROPS
@given(data=st.data())
def test_inframe_window_is_cds_inframe_and_flip_is_antisense(sim: Sim, data: st.DataObject) -> None:
    cdss = [g for g in sim.model.genes if g.feature_type == "cds" and g.interval.length_bp >= 60]
    g = data.draw(st.sampled_from(cdss), label="gene")
    n_cod = g.interval.length_bp // 3
    k = data.draw(st.integers(7, min(30, n_cod)), label="aa")
    off = data.draw(st.integers(0, n_cod - k), label="offset")
    rep = sim.reps[g.interval.contig]
    w = window(g, rep, off, k)
    near = others(sim.model, w, g.locus_tag)
    assume(not any(isinstance(x, NcFeature) and x.type in RNA for x in near))
    assume(not any(isinstance(x, Gene) and x.feature_type == "pseudogene" for x in near))
    assert classify(w, sim.model, sim.genome) == ("cds_inframe", g.locus_tag)
    rc = flip(w)
    # on the other strand nothing may outrank antisense: no gene there, no N-terminal extension
    assume(not any(isinstance(x, Gene) and x.interval.strand == rc.strand for x in near))
    assume(not _nterm_possible(rc, sim))
    assert classify(rc, sim.model, sim.genome) == ("antisense", None)


def _random_hit(data: st.DataObject, sim: Sim) -> Interval:
    contig = data.draw(st.sampled_from(sorted(sim.reps)), label="contig")
    rep = sim.reps[contig]
    k = data.draw(st.integers(7, 30), label="aa")
    hi = rep.length - 1 if rep.topology == "circular" else rep.length - 3 * k
    s = data.draw(st.integers(0, hi), label="start")
    strand = data.draw(st.sampled_from(["+", "-"]), label="strand")
    return Interval(contig, s, s + 3 * k, strand)  # type: ignore[arg-type]


def _near_gene_hit(data: st.DataObject, sim: Sim) -> Interval:
    """Hits near annotated features (where the rules interact), any frame and strand."""
    items = [*sim.model.genes, *sim.model.features]
    x = data.draw(st.sampled_from(items), label="item")
    rep = sim.reps[x.interval.contig]
    k = data.draw(st.integers(7, 30), label="aa")
    s = x.interval.start + data.draw(st.integers(-120, x.interval.length_bp + 60), label="d")
    s %= rep.length
    if rep.topology == "linear":
        s = min(s, rep.length - 3 * k)
    strand = data.draw(st.sampled_from(["+", "-"]), label="strand")
    return Interval(x.interval.contig, s, s + 3 * k, strand)  # type: ignore[arg-type]


@PROPS
@given(data=st.data())
def test_invariant_to_feature_order(sim: Sim, data: st.DataObject) -> None:
    perm_g = data.draw(st.permutations(sim.model.genes), label="genes")
    perm_f = data.draw(st.permutations(sim.model.features), label="features")
    shuffled = GeneModel(sim.reps, tuple(perm_g), tuple(perm_f))
    for _ in range(5):
        h = _near_gene_hit(data, sim) if data.draw(st.booleans()) else _random_hit(data, sim)
        assert classify(h, shuffled, sim.genome) == classify(h, sim.model, sim.genome), h


def _mirror_iv(i: Interval, n: int) -> Interval:
    strand = {"+": "-", "-": "+"}.get(i.strand, i.strand)
    return Interval(i.contig, n - i.end, n - i.start, strand)  # type: ignore[arg-type]


@pytest.fixture(scope="module")
def mirrored(sim_plain: Sim) -> tuple[GeneModel, dict[str, str]]:
    import dataclasses

    gm = sim_plain.model
    n = {c: r.length for c, r in sim_plain.reps.items()}
    genes = tuple(
        dataclasses.replace(g, interval=_mirror_iv(g.interval, n[g.interval.contig]))
        for g in gm.genes
    )
    feats = tuple(
        dataclasses.replace(f, interval=_mirror_iv(f.interval, n[f.interval.contig]))
        for f in gm.features
    )
    genome = {c: _rc(s) for c, s in sim_plain.genome.items()}
    return GeneModel(sim_plain.reps, genes, feats), genome


@PROPS
@given(data=st.data())
def test_mirror_symmetry(
    sim_plain: Sim, mirrored: tuple[GeneModel, dict[str, str]], data: st.DataObject
) -> None:
    h = (
        _near_gene_hit(data, sim_plain)
        if data.draw(st.booleans())
        else _random_hit(data, sim_plain)
    )
    rep = sim_plain.reps[h.contig]
    assume(h.end <= rep.length)  # mirroring an origin-crossing hit is not a plain interval
    gm2, genome2 = mirrored
    assert classify(_mirror_iv(h, rep.length), gm2, genome2) == classify(
        h, sim_plain.model, sim_plain.genome
    )


# --- agreement with simulator truth -----------------------------------------------------------


def _interior(i: Interval, trim_start: int, trim_stop: int) -> Interval:
    """Drop ``trim_start`` nt at the 5' end and ``trim_stop`` nt at the 3' end (strand-aware)."""
    if i.strand == "+":
        return Interval(i.contig, i.start + trim_start, i.end - trim_stop, i.strand)
    return Interval(i.contig, i.start + trim_stop, i.end - trim_start, i.strand)


def test_sim_truth_planted_features(sim: Sim) -> None:
    gm = sim.model
    n_sorf = n_pseudo = n_rna = 0
    for g in gm.genes:
        if g.feature_type == "sorf":
            w = _interior(g.interval, 3, 3)  # without start and stop codons
            assert classify(w, gm, sim.genome) == ("sorf", g.locus_tag)
            n_sorf += 1
        elif g.feature_type == "pseudogene":
            for strand_off in (0, 1, 2):  # any frame on the same strand (Q3)
                s = g.interval.start + 30 + strand_off
                w = Interval(g.interval.contig, s, s + 45, g.interval.strand)
                assert classify(w, gm, sim.genome) == ("pseudogene", g.locus_tag)
            n_pseudo += 1
    for f in gm.features:
        if f.type in RNA and f.interval.length_bp >= 60:
            for strand in ("+", "-"):
                s = f.interval.start + 10
                w = Interval(f.interval.contig, s, s + 30, strand)  # type: ignore[arg-type]
                assert classify(w, gm, sim.genome) == ("rna_overlap", None)
            n_rna += 1
    assert n_sorf >= 4 and n_pseudo >= 4 and n_rna >= 3


def test_sim_truth_nterm_extensions(sim: Sim) -> None:
    rows = [r for r in sim.truth("features.tsv") if r["true_start"] not in ("", "NA", ".")]
    assert len(rows) >= 4
    for r in rows:
        g = sim.model.by_locus(r["feature_id"])
        ts = int(r["true_start"])  # 1-based first base of the true start codon
        iv = g.interval
        if iv.strand == "+":
            w = Interval(iv.contig, ts - 1, iv.start, "+")
        else:
            w = Interval(iv.contig, iv.end, ts, "-")
        assert w.length_bp % 3 == 0 and w.length_bp >= 15
        assert classify(w, sim.model, sim.genome) == ("nterm_extension", g.locus_tag)
        # also reaching into the annotated CDS
        into = (
            Interval(iv.contig, ts - 1, iv.start + 30, "+")
            if iv.strand == "+"
            else Interval(iv.contig, iv.end - 30, ts, "-")
        )
        assert classify(into, sim.model, sim.genome) == ("nterm_extension", g.locus_tag)


def test_sim_truth_novel_orfs(sim: Sim) -> None:
    rows = sim.truth("novel_orfs.tsv")
    seen: dict[str, int] = {}
    for r in rows:
        orf = Interval(r["contig"], int(r["start"]) - 1, int(r["end"]), r["strand"])  # type: ignore[arg-type]
        w = _interior(orf, 0, 3)  # the translated peptide region: start codon kept, stop dropped
        hosts = [x for x in sim.model.overlapping(w)]
        if r["class"] != "intergenic" and len(hosts) != 1:
            continue  # a second feature may legitimately outrank the truth class
        assert classify(w, sim.model, sim.genome) == (r["class"], None), r
        seen[r["class"]] = seen.get(r["class"], 0) + 1
    assert seen.get("intergenic") == 3
    assert seen.get("cds_outframe", 0) >= 2 and seen.get("antisense", 0) >= 2


def test_sim_origin_crossing_gene(sim: Sim) -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    chrom = next(iter(sim.reps))
    rep = sim.reps[chrom]
    crossing = [g for g in sim.model.genes if g.interval.crosses_origin(rep)]
    assert len(crossing) == 1
    g = crossing[0]
    # an in-frame window that spans the origin
    off = (rep.length - g.interval.start) // 3 - 3 if g.interval.strand == "+" else None
    if off is None:
        off = (g.interval.end - rep.length) // 3 - 3
    w = window(g, rep, off, 8)
    assert len(split_at_origin(w, rep)) == 2
    assert classify(w, sim.model, sim.genome) == ("cds_inframe", g.locus_tag)
    row = {c: None for c in HIT_COLUMNS}
    row.update(
        peptide="ORIGINPEP", query="ORLGLNPEP", search="target", contig=chrom, frame_id=0,
        aa_start=0, aa_end=8, start=w.start, end=w.end, strand=w.strand, frame=1,
        route="exact", n_diff=0, diff_types="", source="database", placement_q=0.0, flags="",
    )  # fmt: skip
    acc = pd.DataFrame([row])
    peps = pd.DataFrame(
        {"peptide": ["ORIGINPEP"], "protein_ids": [None], "search_db": ["other"],
         "source": ["database"], "is_decoy": [False], "is_contaminant": [False]}
    )  # fmt: skip
    out = annotate_peptide_hits(acc, sim.model, sim.genome, peps, "related")
    r = out.iloc[0]
    assert r["gene"] == g.locus_tag and r["uniqueness"] == "unique"
    assert "origin_spanning" in str(r["flags"]).split(";")


# --- annotate_peptide_hits: row-order invariance on sim hits ---------------------------------


@settings(max_examples=40, deadline=None)
@given(data=st.data())
def test_annotate_row_order_invariant(sim: Sim, data: st.DataObject) -> None:
    from gippa.mapping.loci import annotate_peptide_hits

    hits = [_near_gene_hit(data, sim) for _ in range(data.draw(st.integers(1, 12)))]
    names = [data.draw(st.sampled_from(["PA", "PB", "PC"])) for _ in hits]
    rows = []
    seen = set()
    for name, h in zip(names, hits, strict=True):
        if (name, h) in seen:
            continue
        seen.add((name, h))
        row = {c: None for c in HIT_COLUMNS}
        row.update(
            peptide=name, query=name, search="target", contig=h.contig, frame_id=0, aa_start=0,
            aa_end=h.length_bp // 3, start=h.start, end=h.end, strand=h.strand, frame=1,
            route=data.draw(st.sampled_from(["exact", "tolerant"])), n_diff=0, diff_types="",
            source="database", placement_q=0.0, flags="",
        )  # fmt: skip
        rows.append(row)
    peps = pd.DataFrame(
        {"peptide": ["PA", "PB", "PC"], "protein_ids": [None] * 3, "search_db": ["other"] * 3,
         "source": ["database"] * 3, "is_decoy": [False] * 3, "is_contaminant": [False] * 3}
    )  # fmt: skip
    perm = data.draw(st.permutations(range(len(rows))))

    def canon(df: pd.DataFrame) -> list[tuple[Any, ...]]:
        return sorted(
            (
                r["peptide"],
                r["contig"],
                int(r["start"]),
                int(r["end"]),
                r["strand"],
                r["locus_class"],
                None if pd.isna(r["gene"]) else r["gene"],
                r["uniqueness"],
                tuple(sorted(f for f in str(r["flags"] or "").split(";") if f)),
            )
            for r in df.to_dict("records")
        )

    a = annotate_peptide_hits(pd.DataFrame(rows), sim.model, sim.genome, peps, "related")
    b = annotate_peptide_hits(
        pd.DataFrame([rows[i] for i in perm]), sim.model, sim.genome, peps, "related"
    )
    assert canon(a) == canon(b)
