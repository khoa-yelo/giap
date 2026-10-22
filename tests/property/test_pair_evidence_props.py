"""Hypothesis property guarding the M6-04 speed-up (PR #143; docs/specs/M6-04-pair-table.md
decision 5, GIAP_SPEC.md §7.4.1).

``pair_frame`` indexes the layers once (``promoter_tracks``, per-strand Tx arrays, a terminator
pre-filter around each gap). On small linear and circular contigs with dense layers anchored at the
gene boundaries (origin-crossing genes and terminators, overlaps, 2- and 3-gene rings, empty
tables), every OP/NO/SW row's evidence equals the still-public per-pair helpers run on the full
layers (``promoter_between``, ``tx_rate_between``, ``terminator_between(..., yes=0.0)``), CV/DV
evidence is NA, and ``rna_between`` matches a brute-force base-set overlap of the gap(s).
"""

from __future__ import annotations

import math
from collections import Counter

import pandas as pd
from helpers.pairs_data import evidence, frame, gene, rna_frame
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.operons import adjacent_pairs, intergenic_distance
from giap.annotate.promoter_calculator import tx_rate_between
from giap.annotate.promoters import promoter_between
from giap.annotate.terminators import Terminator, terminator_between
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

SOURCES = ("TerminatorNet", "TerminatorNet", "TransTermHP", "BacTermFinder")


@st.composite
def worlds(draw: st.DrawFn) -> tuple:
    length = draw(st.integers(120, 500))
    circ = draw(st.booleans())
    rep = Replicon("c1", length, "circular" if circ else "linear")
    genes: list[Gene] = []
    for i in range(draw(st.sampled_from([0, 1, 2, 2, 3, 3, 4, 5, 6]))):
        s = draw(st.integers(0, length - 2))
        e = s + draw(st.integers(1, length // 3))
        genes.append(gene(f"g{i}", s, e if circ else min(e, length), draw(st.sampled_from("+-"))))
    anchors = [x for g in genes for x in (g.interval.start, g.interval.end)]
    anchors += [0, length - 1] if circ else [0]  # calls at and across the origin

    def pos() -> int:  # near a gene boundary, so every edge and the pre-filter margin get hit
        p = draw(st.sampled_from(anchors)) + draw(st.integers(-70, 70))
        return p % length if circ else min(max(p, 0), length - 1)

    strand = st.sampled_from("+-")
    conf = st.sampled_from([0.1, 0.35, 0.6, 0.9])
    terms = []
    for _ in range(draw(st.integers(0, 25))):
        s = pos()
        e = s + draw(st.integers(1, 40))
        iv = Interval("c1", s, e if circ else min(e, length), draw(strand))
        terms.append(Terminator(iv, draw(conf), draw(st.sampled_from(SOURCES))))
    k, m, r = (draw(st.integers(0, 25)) for _ in range(3))
    windows = pd.DataFrame({"contig": ["c1"] * k, "strand": [draw(strand) for _ in range(k)],
                            "end3p": [pos() + 1 for _ in range(k)],
                            "score": [draw(conf) for _ in range(k)]})  # fmt: skip
    tx = pd.DataFrame({"contig": ["c1"] * m, "strand": [draw(strand) for _ in range(m)],
                       "tss": [pos() for _ in range(m)],
                       "tx_rate": [draw(conf) for _ in range(m)]})  # fmt: skip
    tx = tx.astype({"tss": "int64", "tx_rate": "float64"})
    rna = []
    for j in range(r):
        s = pos()
        rna.append((f"r{j}", "trna", s, s + draw(st.integers(1, 30)), draw(strand)))
    if not circ:
        rna = [(f, t, s, min(e, length), x) for f, t, s, e, x in rna if s < length - 1]
    return genes, rep, windows, terms, tx, rna_frame(rna)


def _gap(left: Gene, right: Gene, rep: Replicon) -> set[int]:
    lo, n = left.interval.end, rep.length
    return {(lo + k) % n for k in range(max(intergenic_distance(left, right, rep), 0))}


def _rna_hit(gaps: list[tuple[Gene, Gene]], rna: pd.DataFrame, rep: Replicon) -> bool:
    bases = {p % rep.length for s, e in zip(rna["start"], rna["end"], strict=True)
             for p in range(s, e)}  # fmt: skip
    return any(_gap(a, b, rep) & bases for a, b in gaps)


def _same(got: float, want: float | None) -> bool:
    return (want is None and math.isnan(got)) or (want is not None and got == want)


@settings(max_examples=150, deadline=None)
@given(worlds())
def test_pair_frame_evidence_equals_the_per_pair_helpers(world: tuple) -> None:
    from giap.operons.pairs import adjacent_chain, sandwich_pairs

    genes, rep, windows, terms, tx, rna = world
    replicons = {"c1": rep}
    probs = {(a.locus_tag, b.locus_tag): 0.6 for a, b in adjacent_pairs(genes, replicons)}
    df = frame(genes, replicons, probs, ev=evidence(windows, terms, tx, rna))
    by = {g.locus_tag: g for g in genes}
    cols = ["class", "gene_a", "gene_b", "promoter_between", "tx_rate_between",
            "terminator_between", "flags"]  # fmt: skip
    got = Counter()
    for cls, a, b, prom, txr, term, flags in df[cols].itertuples(index=False):
        if cls in ("CV", "DV"):
            assert all(math.isnan(v) for v in (prom, txr, term)), (cls, a, b)
        else:
            u, d = by[a], by[b]
            assert _same(prom, promoter_between(u, d, windows, rep)), ("promoter", cls, a, b)
            assert _same(txr, tx_rate_between(u, d, tx, rep)), ("tx", cls, a, b)
            want = terminator_between(u, d, terms, rep, yes=0.0)
            assert _same(term, want), ("terminator", cls, a, b, term, want)
        got[(cls == "SW", frozenset((a, b)), "rna_between" in flags)] += 1
    want_rna = Counter(
        (False, frozenset((x.locus_tag, y.locus_tag)), _rna_hit([(x, y)], rna, rep))
        for x, y in adjacent_chain(genes, replicons)
    )
    for x, y, z in sandwich_pairs(genes, replicons):
        key = frozenset((x.locus_tag, z.locus_tag))
        want_rna[(True, key, _rna_hit([(x, y), (y, z)], rna, rep))] += 1
    assert got == want_rna
