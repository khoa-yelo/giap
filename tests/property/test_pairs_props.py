"""Hypothesis properties of the M6-04 pair table (docs/specs/M6-04-pair-table.md, acceptance
criteria; GIAP_SPEC.md §7.4.1).

* the OP/NO/CV/DV rows are exactly ``adjacent_chain``, each pair once, one class per pair, the
  class agreeing with the strands; SW rows are the brute-force sandwiches and never adjacent pairs;
* mirroring (coordinates x -> L - x, every strand flipped) keeps every class and every distance.
  OP/NO/SW keep ``gene_a``/``gene_b`` and Δ (orientation is strand-aware). For CV/DV the left gene
  becomes the right one, so a/b swap and Δ changes sign (test author's reading of the criterion
  "negates no Δ", which only holds for the transcription-oriented classes);
* ``pair_id`` does not move when ``pair_threshold`` goes 0.3 -> 0.7, while pairs with
  0.3 <= uniop_p < 0.7 go OP -> NO.
"""

from __future__ import annotations

from collections import Counter

import pandas as pd
from helpers.pairs_data import frame, gene, gene_tbl, reps
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.operons import adjacent_pairs
from giap.model.genes import Gene

P_VALUES = (0.05, 0.3, 0.45, 0.5, 0.69, 0.7, 0.95)


@st.composite
def layouts(draw: st.DrawFn) -> tuple[list[Gene], dict, dict, pd.DataFrame]:
    n = draw(st.integers(1, 10))
    pos = draw(st.integers(0, 100))
    genes = []
    for i in range(n):
        length = draw(st.integers(60, 300))
        genes.append(gene(f"g{i}", pos, pos + length, draw(st.sampled_from("+-"))))
        pos += length + draw(st.integers(-20, 200))
    length = max(g.interval.end for g in genes) + draw(st.integers(1, 400))
    replicons = reps(("c1", length, draw(st.sampled_from(["linear", "circular"]))))
    probs = {
        (a.locus_tag, b.locus_tag): draw(st.sampled_from(P_VALUES))
        for a, b in adjacent_pairs(genes, replicons)
    }
    ab = {g.locus_tag: draw(st.sampled_from([float("nan"), 12.0, 15.5, 20.25])) for g in genes}
    return genes, replicons, probs, gene_tbl(genes, {k: v for k, v in ab.items() if v == v})


def _key(a: str, b: str) -> frozenset[str]:
    return frozenset((a, b))


@settings(max_examples=60, deadline=None)
@given(layouts())
def test_adjacent_rows_are_the_chain_and_sw_rows_the_sandwiches(lay: tuple) -> None:
    from giap.operons.pairs import adjacent_chain

    genes, replicons, probs, tbl = lay
    df = frame(genes, replicons, probs, tbl=tbl)
    strand = {g.locus_tag: g.interval.strand for g in genes}
    adj = df[df["class"] != "SW"]
    chain = Counter(_key(a.locus_tag, b.locus_tag) for a, b in adjacent_chain(genes, replicons))
    got = Counter(_key(a, b) for a, b in zip(adj["gene_a"], adj["gene_b"], strict=True))
    assert got == chain and set(got.values()) <= {1}
    for cls, a, b in zip(adj["class"], adj["gene_a"], adj["gene_b"], strict=True):
        same = strand[a] == strand[b]
        assert (cls in ("OP", "NO")) == same
        if not same:
            assert (cls, strand[a]) in (("CV", "+"), ("DV", "-"))
    n, circ = len(genes), replicons["c1"].topology == "circular"
    idx = range(n) if circ and n >= 3 else range(max(n - 2, 0))
    tags = [g.locus_tag for g in genes]
    want = Counter(
        _key(tags[i], tags[(i + 2) % n])
        for i in idx
        if strand[tags[i]] == strand[tags[(i + 2) % n]] != strand[tags[(i + 1) % n]]
    )
    sw = df[df["class"] == "SW"]
    if not (circ and n == 3):
        assert Counter(_key(a, b) for a, b in zip(sw["gene_a"], sw["gene_b"], strict=True)) == want
        assert not set(want) & set(chain)


def _mirror(genes: list[Gene], length: int) -> list[Gene]:
    flip = {"+": "-", "-": "+"}
    return [
        gene(g.locus_tag, length - g.interval.end, length - g.interval.start,
             flip[g.interval.strand])
        for g in genes
    ]  # fmt: skip


def _rows(df: pd.DataFrame) -> dict[tuple[str, str, str], tuple[int, float, str]]:
    cols = ["class", "gene_a", "gene_b", "intergenic_bp", "delta_log2", "strand"]
    return {(c, a, b): (int(d), x, s) for c, a, b, d, x, s in df[cols].itertuples(index=False)}


@settings(max_examples=60, deadline=None)
@given(layouts())
def test_mirroring_keeps_classes_orientation_and_delta(lay: tuple) -> None:
    genes, replicons, probs, tbl = lay
    length = replicons["c1"].length
    orig = _rows(frame(genes, replicons, probs, tbl=tbl))
    mprobs = {(b, a): p for (a, b), p in probs.items()}
    mirr = _rows(frame(_mirror(genes, length), replicons, mprobs, tbl=tbl))
    assert len(orig) == len(mirr)
    flip = {"+": "-", "-": "+", ".": "."}
    for (cls, a, b), (d, delta, s) in orig.items():
        swap = cls in ("CV", "DV")
        md, mdelta, ms = mirr[(cls, b, a) if swap else (cls, a, b)]
        assert md == d and ms == flip[s]
        assert (pd.isna(delta) and pd.isna(mdelta)) or mdelta == (-delta if swap else delta)


@settings(max_examples=60, deadline=None)
@given(layouts())
def test_pair_ids_do_not_depend_on_the_threshold(lay: tuple) -> None:
    genes, replicons, probs, tbl = lay
    lo = frame(genes, replicons, probs, tbl=tbl, threshold=0.3)
    hi = frame(genes, replicons, probs, tbl=tbl, threshold=0.7)

    def ids(df: pd.DataFrame) -> dict[str, tuple[str, str]]:
        return dict(zip(df["pair_id"], zip(df["gene_a"], df["gene_b"], strict=True), strict=True))

    assert ids(lo) == ids(hi)
    m = lo.merge(hi, on="pair_id", suffixes=("_lo", "_hi"))
    moved = m["class_lo"] != m["class_hi"]
    band = m["uniop_p_lo"].between(0.3, 0.7, inclusive="left")
    assert (moved == band).all()
    assert (m.loc[moved, "class_lo"] == "OP").all() and (m.loc[moved, "class_hi"] == "NO").all()
