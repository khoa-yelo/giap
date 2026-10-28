"""The M6-04 pair table: classes, orientation, distances, evidence (docs/specs/M6-04-pair-table.md;
GIPPA_SPEC.md §7.4.1, §7.1.3-§7.1.5, §9.3).

Examples 1-11 and the degenerate-case table on hand-built frames (``helpers.pairs_data``). The
``layer_missing:<layer>`` run flag of example 11 is emitted by the ``operons`` step (M6-08) and is
tested there; here a missing layer only means an all-NA evidence column.
"""

from __future__ import annotations

import pandas as pd
import pytest
from helpers.pairs_data import (
    BASE,
    PAIR_COLS,
    TEST_COLS,
    L,
    evidence,
    frame,
    gene,
    gene_tbl,
    genes_of,
    reps,
    ring,
    rna_frame,
    row,
)

from gippa.annotate.promoter_calculator import tx_rate_between
from gippa.annotate.promoters import promoter_between
from gippa.annotate.terminators import Terminator, terminator_between
from gippa.model.coords import Interval

G = genes_of(BASE)
BY = {g.locus_tag: g for g in G}
P_BASE = {("G1", "G2"): 0.8, ("G3", "G4"): 0.9, ("G5", "G1"): 0.1}
EVIDENCE_COLS = ("promoter_between", "tx_rate_between", "terminator_between")


def _base(probs: dict[tuple[str, str], float] | None = None, **kw: object) -> pd.DataFrame:
    return frame(G, ring(), P_BASE if probs is None else probs, **kw)  # type: ignore[arg-type]


def _flags(r: pd.Series) -> object:
    f = r["flags"]
    return () if f is None or (isinstance(f, float) and pd.isna(f)) else f


# ------------------------------------------------------------------ adjacency (decision 1)


def test_adjacent_pairs_same_strand_false_keeps_mixed_pairs() -> None:
    from gippa.annotate.operons import adjacent_pairs

    tags = lambda ps: [(a.locus_tag, b.locus_tag) for a, b in ps]  # noqa: E731
    assert tags(adjacent_pairs(G, ring())) == [("G1", "G2"), ("G3", "G4"), ("G5", "G1")]
    assert tags(adjacent_pairs(G, ring(), same_strand=False)) == [
        ("G1", "G2"), ("G2", "G3"), ("G3", "G4"), ("G4", "G5"), ("G5", "G1"),
    ]  # fmt: skip


def test_adjacent_chain_is_genome_order_with_wrap() -> None:
    from gippa.operons.pairs import adjacent_chain, sandwich_pairs

    chain = [(a.locus_tag, b.locus_tag) for a, b in adjacent_chain(G, ring())]
    assert chain == [("G1", "G2"), ("G2", "G3"), ("G3", "G4"), ("G4", "G5"), ("G5", "G1")]
    assert sandwich_pairs(G, ring()) == []  # example 6: no SW in the base layout


# ------------------------------------------------------------------ examples 1-8, 10


def test_columns_and_one_row_per_adjacent_pair() -> None:
    df = _base()
    assert set(PAIR_COLS) <= set(df.columns)
    assert len(df) == 5 and set(df["class"]) <= {"OP", "NO", "CV", "DV"}


def test_example_1_op_and_example_2_no() -> None:
    r = row(_base(), "G1", "G2")
    assert (r["class"], r["strand"], r["intergenic_bp"], r["uniop_p"]) == ("OP", "+", 50, 0.8)
    r = row(_base({**P_BASE, ("G1", "G2"): 0.2}), "G1", "G2")
    assert (r["class"], r["strand"], r["intergenic_bp"], r["uniop_p"]) == ("NO", "+", 50, 0.2)


def test_threshold_is_inclusive() -> None:
    assert row(_base({**P_BASE, ("G1", "G2"): 0.5}), "G1", "G2")["class"] == "OP"


def test_example_3_convergent() -> None:
    r = row(_base(), "G2", "G3")
    assert (r["class"], r["strand"], r["intergenic_bp"]) == ("CV", ".", 50)
    assert pd.isna(r["uniop_p"]) and all(pd.isna(r[c]) for c in EVIDENCE_COLS)


def test_example_4_minus_strand_op_is_oriented_upstream_first() -> None:
    ab = {"G3": 18.0, "G4": 16.5}
    df = _base(tbl=gene_tbl(G, ab))
    r = row(df, "G4", "G3")
    assert (r["class"], r["strand"], r["intergenic_bp"], r["uniop_p"]) == ("OP", "-", 100, 0.9)
    assert r["delta_log2"] == pytest.approx(1.5)  # example 10, second half
    assert not ((df["gene_a"] == "G3") & (df["gene_b"] == "G4")).any()


def test_example_5_divergent_a_is_left() -> None:
    r = row(_base(), "G4", "G5")
    assert (r["class"], r["strand"], r["intergenic_bp"]) == ("DV", ".", 100)


def test_example_7_sandwich_decoy() -> None:
    rows = [*BASE[:2], ("G3", 950, 1400, "+"), *BASE[3:]]
    gs = genes_of(rows)
    probs = {("G1", "G2"): 0.8, ("G2", "G3"): 0.7, ("G5", "G1"): 0.1}
    df = frame(gs, ring(), probs)
    sw = df[df["class"] == "SW"]
    assert len(sw) == 1
    r = sw.iloc[0]
    assert (r["gene_a"], r["gene_b"], r["strand"], r["intergenic_bp"]) == ("G3", "G5", "+", 200)
    assert pd.isna(r["uniop_p"])
    assert len(df) == 6  # five adjacent rows + the SW row


def test_sw_on_minus_strand_is_oriented_upstream_first() -> None:
    gs = genes_of([("A", 100, 400, "-"), ("B", 500, 800, "+"), ("C", 900, 1200, "-")])
    df = frame(gs, reps(("c1", L, "linear")), {})
    r = row(df, "C", "A", "SW")
    assert (r["strand"], r["intergenic_bp"]) == ("-", 200)


def test_example_8_wrap_pair() -> None:
    r = row(_base(), "G5", "G1")
    assert (r["class"], r["strand"], r["intergenic_bp"]) == ("NO", "+", 10_000 - 2400 + 100)


def test_example_10_delta_and_per_gene_columns() -> None:
    ab = {"G1": 18.0, "G2": 16.5, "G3": 18.0}
    tbl = gene_tbl(G, ab, p_detect={"G1": 0.95, "G2": 0.6})
    df = _base(tbl=tbl)
    r = row(df, "G1", "G2")
    assert r["delta_log2"] == pytest.approx(-1.5)
    assert (r["abundance_log2_a"], r["abundance_log2_b"]) == (18.0, 16.5)
    assert (bool(r["detected_a"]), bool(r["detected_b"])) == (True, True)
    assert (r["p_detect_a"], r["p_detect_b"]) == (0.95, 0.6)
    assert row(df, "G2", "G3")["delta_log2"] == pytest.approx(1.5)  # CV: b - a, a = left
    assert pd.isna(row(df, "G4", "G5")["delta_log2"])  # G4 has no abundance


def test_test_columns_are_na_and_none() -> None:
    df = _base()
    assert (df["test"] == "none").all()
    assert all(df[c].isna().all() for c in TEST_COLS)


def test_pair_ids_follow_contig_start_a_start_b_order() -> None:
    df = _base()
    got = dict(zip(df["pair_id"], zip(df["gene_a"], df["gene_b"], strict=True), strict=True))
    assert got == {
        "PAIR_1": ("G1", "G2"), "PAIR_2": ("G2", "G3"), "PAIR_3": ("G4", "G3"),
        "PAIR_4": ("G4", "G5"), "PAIR_5": ("G5", "G1"),
    }  # fmt: skip


def test_pair_ids_follow_c1_contig_order_not_name_order() -> None:
    gs = genes_of([("B1", 100, 400, "+"), ("B2", 500, 800, "-")], "chrB")
    gs += genes_of([("A1", 100, 400, "+"), ("A2", 500, 800, "-")], "chrA")
    df = frame(gs, reps(("chrB", 5000, "linear"), ("chrA", 5000, "linear")), {})
    assert row(df, "B1", "B2")["pair_id"] == "PAIR_1"
    assert row(df, "A1", "A2")["pair_id"] == "PAIR_2"


# ------------------------------------------------------------------ evidence, examples 9 and 11


def _layers() -> tuple[pd.DataFrame, list[Terminator], pd.DataFrame]:
    windows = pd.DataFrame(
        {"contig": "c1", "strand": ["+", "+", "-"], "end3p": [380, 420, 1450],
         "score": [0.95, 0.7, 0.6]}
    )  # fmt: skip
    terms = [
        Terminator(Interval("c1", 405, 440, "+"), 0.3, "TerminatorNet"),
        Terminator(Interval("c1", 910, 940, "+"), 0.9, "TerminatorNet"),  # in the CV gap
        Terminator(Interval("c1", 1410, 1450, "-"), 0.8, "TerminatorNet"),
    ]
    tx = pd.DataFrame({"contig": "c1", "strand": ["+", "-"], "tss": [420, 1450],
                       "tx_rate": [5.0, 2.0]})  # fmt: skip
    return windows, terms, tx


def test_evidence_columns_use_the_m2_helpers_with_raw_terminator_maximum() -> None:
    windows, terms, tx = _layers()
    df = _base(ev=evidence(windows, terms, tx))
    rep = ring()["c1"]
    for up, down in (("G1", "G2"), ("G4", "G3")):
        r = row(df, up, down)
        u, d = BY[up], BY[down]
        exp = (
            promoter_between(u, d, windows, rep),
            tx_rate_between(u, d, tx, rep),
            terminator_between(u, d, terms, rep, yes=0.0),
        )
        for col, e in zip(EVIDENCE_COLS, exp, strict=True):
            assert (pd.isna(r[col]) and e is None) or r[col] == pytest.approx(e), (up, col)
    r = row(df, "G1", "G2")
    assert (r["promoter_between"], r["tx_rate_between"]) == (0.7, 5.0)
    assert r["terminator_between"] == pytest.approx(0.3)  # below the default yes = 0.5 (dec. 5)


def test_cv_dv_evidence_is_na_even_with_layers() -> None:
    df = _base(ev=evidence(*_layers()))  # terminator_between would raise on mixed strands
    for a, b in (("G2", "G3"), ("G4", "G5")):
        assert all(pd.isna(row(df, a, b)[c]) for c in EVIDENCE_COLS)


@pytest.mark.parametrize("missing", [0, 1, 2])
def test_example_11_missing_layer_gives_all_na_column(missing: int) -> None:
    layers: list[object] = list(_layers())
    layers[missing] = None
    df = _base(ev=evidence(*layers))  # type: ignore[arg-type]
    order = ("promoter_between", "terminator_between", "tx_rate_between")
    assert df[order[missing]].isna().all()
    assert all(df[c].notna().any() for c in order if c != order[missing])


def test_example_9_rna_between() -> None:
    rna = rna_frame([("t1", "trna", 410, 440, "+"), ("r1", "rrna", 3000, 4000, "+")])
    df = _base(ev=evidence(rna=rna))
    assert "rna_between" in _flags(row(df, "G1", "G2"))  # type: ignore[operator]
    assert row(df, "G1", "G2")["class"] == "OP"  # otherwise ordinary
    for a, b in (("G2", "G3"), ("G4", "G3"), ("G4", "G5")):
        assert "rna_between" not in _flags(row(df, a, b))  # type: ignore[operator]
    assert "rna_between" in _flags(row(df, "G5", "G1"))  # type: ignore[operator]  # 3000-4000


def test_rna_between_on_sw_uses_the_union_of_both_gaps() -> None:
    gs = genes_of([("A", 100, 400, "+"), ("B", 500, 800, "-"), ("C", 900, 1200, "+")])
    rna = rna_frame([("n1", "ncrna", 820, 880, "-")])
    df = frame(gs, reps(("c1", L, "linear")), {}, ev=evidence(rna=rna))
    assert "rna_between" in _flags(row(df, "A", "C", "SW"))  # type: ignore[operator]


# ------------------------------------------------------------------ degenerate cases


def test_one_gene_replicon_gives_no_pairs() -> None:
    gs = [gene("A", 100, 400, "+")]
    for topo in ("linear", "circular"):
        df = frame(gs, reps(("c1", L, topo)), {})
        assert len(df) == 0 and set(PAIR_COLS) <= set(df.columns)


def test_circular_two_genes_one_pair_no_sw() -> None:
    gs = genes_of([("A", 100, 400, "+"), ("B", 500, 800, "+")])
    df = frame(gs, ring(), {("A", "B"): 0.9, ("B", "A"): 0.1})  # M2-12 writes both
    assert len(df) == 1
    r = df.iloc[0]
    assert (r["gene_a"], r["gene_b"], r["intergenic_bp"]) == ("A", "B", 100)


def test_linear_no_wrap_pair_and_no_sw_across_ends() -> None:
    rows = [("A", 100, 400, "+"), ("B", 500, 800, "+"), ("C", 900, 1200, "+"),
            ("D", 1300, 1600, "-")]  # fmt: skip
    gs = genes_of(rows)
    probs = {("A", "B"): 0.9, ("B", "C"): 0.9}
    lin = frame(gs, reps(("c1", L, "linear")), probs)
    assert len(lin) == 3 and "SW" not in set(lin["class"])
    circ = frame(gs, ring(), probs)
    assert row(circ, "D", "A", "DV")["intergenic_bp"] == L - 1600 + 100
    assert row(circ, "C", "A", "SW")["intergenic_bp"] == 100 + (L - 1600 + 100)


def test_partial_gene_keeps_its_row_and_gene_table_values() -> None:
    gs = [gene("A", 100, 400, "+", partial=True), gene("B", 450, 900, "+")]
    tbl = gene_tbl(gs, detected={"A": False}, p_detect={"A": float("nan")})
    r = row(frame(gs, reps(("c1", L, "linear")), {("A", "B"): 0.9}, tbl=tbl), "A", "B")
    assert (r["class"], r["intergenic_bp"], bool(r["detected_a"])) == ("OP", 50, False)
    assert pd.isna(r["p_detect_a"])


def test_overlapping_genes_negative_distance() -> None:
    gs = genes_of([("A", 100, 400, "+"), ("B", 396, 900, "+")])
    r = row(frame(gs, reps(("c1", L, "linear")), {("A", "B"): 0.9}), "A", "B")
    assert r["intergenic_bp"] == -4


def test_gene_crossing_the_origin() -> None:
    gs = genes_of([("A", 150, 400, "+"), ("M", 5000, 5300, "+"), ("Z", 9800, 10_100, "+")])
    probs = {("A", "M"): 0.1, ("M", "Z"): 0.1, ("Z", "A"): 0.9}
    r = row(frame(gs, ring(), probs), "Z", "A")
    assert r["intergenic_bp"] == 50


def test_same_strand_pair_missing_from_pairs_raw_raises() -> None:
    with pytest.raises(ValueError, match=r"G3.*G4|G4.*G3"):
        _base({("G1", "G2"): 0.8, ("G5", "G1"): 0.1})


def test_gene_absent_from_gene_table_raises() -> None:
    tbl = gene_tbl(G).iloc[1:]
    with pytest.raises((ValueError, KeyError), match="G1"):
        _base(tbl=tbl)
