"""Hypothesis properties for the gene model (issue M3-01; GIAP_SPEC.md §6.5, §7.1.2, §11.2).

- For giapsim genomes (random seeds) ``overlapping(g.interval)`` contains ``g`` and ``nt``
  translated with table 11 (initiator read as M) equals ``aa + "*"``.
- ``overlapping`` equals a brute-force origin-aware overlap test on random circular and linear
  replicons, for both strands, and is deterministic.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.model import Interval, Replicon, split_at_origin
from giap.model.genes import Gene, GeneModel, NcFeature, gene_model_from_bakta_json

BASES = "TCAG"
AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE_11 = {
    a + b + c: AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}


def translate11(nt: str) -> str:
    aa = "".join(CODE_11[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))
    return "M" + aa[1:]


@settings(max_examples=8, deadline=None)
@given(seed=st.integers(min_value=0, max_value=2**31 - 1))
def test_giapsim_genes_self_overlap_and_translate(seed: int) -> None:
    import giapsim

    with tempfile.TemporaryDirectory() as d:
        res = giapsim.simulate(giapsim.SimConfig(genome_length=20_000), seed=seed, outdir=d)
        reps = {"chr1": Replicon("chr1", 20_000, "circular", "chromosome")}
        m = gene_model_from_bakta_json(res.bakta_json, reps)
    assert len(m.genes) == res.n_genes
    for g in m.genes:
        assert g in m.overlapping(g.interval)
        assert m.by_locus(g.locus_tag) == g
        assert g.nt is not None
        assert translate11(g.nt) == g.aa + "*"


def _parts(iv: Interval, rep: Replicon) -> tuple[Interval, ...]:
    return split_at_origin(iv, rep)


def _brute(a: Interval, b: Interval, rep: Replicon) -> bool:
    return any(x.start < y.end and y.start < x.end for x in _parts(a, rep) for y in _parts(b, rep))


@st.composite
def interval_on(draw: st.DrawFn, rep: Replicon, strands: str = "+-") -> Interval:
    start = draw(st.integers(min_value=0, max_value=rep.length - 1))
    max_end = start + rep.length if rep.topology == "circular" else rep.length
    end = draw(st.integers(min_value=start + 1, max_value=min(max_end, start + 400)))
    return Interval(rep.id, start, end, draw(st.sampled_from(strands)))


@st.composite
def model_and_query(draw: st.DrawFn) -> tuple[GeneModel, Interval]:
    length = draw(st.integers(min_value=50, max_value=3_000))
    topology = draw(st.sampled_from(["circular", "linear"]))
    rep = Replicon("c1", length, topology)
    other = Replicon("c2", 1_000, "circular")
    ivs = draw(st.lists(interval_on(rep), max_size=25))
    n_nc = draw(st.integers(min_value=0, max_value=len(ivs)))
    order = {"c1": 0, "c2": 1}
    genes = sorted(
        (Gene(f"L{i:04d}", "cds", iv, "M", None, None) for i, iv in enumerate(ivs[n_nc:])),
        key=lambda g: (order[g.interval.contig], g.interval.start, g.interval.end,
                       g.interval.strand, g.locus_tag),
    )  # fmt: skip
    feats = tuple(NcFeature(f"f{i}", "trna", iv) for i, iv in enumerate(ivs[:n_nc]))
    decoy = Gene("OTHER", "cds", Interval("c2", 0, 1_000, "+"), "M", None, None)
    model = GeneModel(replicons={"c1": rep, "c2": other}, genes=(*genes, decoy), features=feats)
    return model, draw(interval_on(rep, "+-."))


@settings(max_examples=300, deadline=None)
@given(mq=model_and_query())
def test_overlapping_equals_brute_force(mq: tuple[GeneModel, Interval]) -> None:
    model, q = mq
    rep = model.replicons["c1"]
    items: list[Gene | NcFeature] = [*model.genes, *model.features]
    expected = {x for x in items if x.interval.contig == "c1" and _brute(x.interval, q, rep)}
    got = model.overlapping(q)
    assert len(got) == len(set(got))  # no duplicates, even when both parts of a split hit
    assert set(got) == expected
    assert model.overlapping(q) == got


@settings(max_examples=20, deadline=None)
@given(perm=st.permutations(list(range(6))))
def test_reader_order_independent(perm: list[int]) -> None:
    rep = Replicon("chr1", 10_000, "circular")
    base = [
        {"type": "cds", "contig": "chr1", "start": s, "stop": e, "strand": st_, "locus": f"G{i}",
         "aa": "MK", "nt": "ATGAAATAA", "start_type": "ATG"}
        for i, (s, e, st_) in enumerate(
            [(101, 400, "+"), (101, 400, "-"), (9901, 50, "+"), (51, 400, "+"),
             (101, 700, "+"), (5001, 5300, "-")]
        )
    ]  # fmt: skip
    with tempfile.TemporaryDirectory() as d:
        p0 = Path(d) / "a.json"
        p1 = Path(d) / "b.json"
        p0.write_text(json.dumps({"features": base}))
        p1.write_text(json.dumps({"features": [base[i] for i in perm]}))
        a = gene_model_from_bakta_json(p0, {"chr1": rep})
        b = gene_model_from_bakta_json(p1, {"chr1": rep})
    assert a.genes == b.genes
    assert [g.locus_tag for g in a.genes] == ["G3", "G0", "G1", "G4", "G5", "G2"]
