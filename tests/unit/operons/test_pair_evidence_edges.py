"""Edge cases of the M6-04 pair table that the #143 speed-up could break and the random property
reaches only by luck (docs/specs/M6-04-pair-table.md rules 4 and 7, decision 5; test author).

A 3-gene ring keeps its wrapping sandwich triples; ``rna_between`` sees RNA after the origin; the
Tx-rate gap is half-open ``[left.end, right.start)`` on both strands (M2-14: strictly between
``up``'s 3' end and ``down``'s start codon); a promoter window far upstream of ``down``'s own
k-range is not ``down``'s window (M2-13); a terminator in the 50 bp upstream slack still counts.
"""

from __future__ import annotations

import pandas as pd
import pytest
from helpers.pairs_data import evidence, frame, genes_of, reps, rna_frame, row

from gippa.annotate.operons import adjacent_pairs
from gippa.annotate.promoter_calculator import tx_rate_between
from gippa.annotate.promoters import promoter_between
from gippa.annotate.terminators import Terminator
from gippa.model.coords import Interval

RING = reps(("c1", 1000, "circular"))
LINE = reps(("c1", 5000, "linear"))


@pytest.mark.parametrize(
    ("layout", "a", "b", "d"),
    [
        ((("A", 100, 300, "+"), ("B", 400, 600, "-"), ("C", 700, 900, "+")), "A", "C", 200),
        ((("A", 100, 300, "-"), ("B", 400, 600, "+"), ("C", 700, 900, "-")), "C", "A", 200),
        # the triple (B, C, A) wraps: d = 100 + (1000 - 900 + 100)
        ((("A", 100, 300, "+"), ("B", 400, 600, "+"), ("C", 700, 900, "-")), "B", "A", 300),
    ],
)
def test_three_gene_ring_keeps_its_sandwich(layout: tuple, a: str, b: str, d: int) -> None:
    gs = genes_of(layout)
    probs = {(x, y): 0.6 for x, y in (("A", "B"), ("B", "C"), ("C", "A"))}
    df = frame(gs, RING, probs)
    sw = df[df["class"] == "SW"]
    assert len(sw) == 1 and row(df, a, b, "SW")["intergenic_bp"] == d


@pytest.mark.parametrize(
    ("c_end", "rna"),
    [(950, (10, 40)), (1050, (60, 70)), (950, (990, 1010))],  # after / across the origin
)
def test_rna_between_on_the_wrap_gap(c_end: int, rna: tuple[int, int]) -> None:
    gs = genes_of([("A", 100, 300, "+"), ("B", 400, 600, "+"), ("C", 800, c_end, "+")])
    probs = {("A", "B"): 0.6, ("B", "C"): 0.6, ("C", "A"): 0.6}
    df = frame(gs, RING, probs, ev=evidence(rna=rna_frame([("t1", "trna", *rna, "-")])))
    assert "rna_between" in row(df, "C", "A")["flags"]
    assert "rna_between" not in row(df, "A", "B")["flags"]
    assert "rna_between" not in row(df, "B", "C")["flags"]


TX_CASES = [  # (layout, replicons, TSS inside the gap, TSS just outside on either side)
    ((("U", 100, 200, "+"), ("D", 300, 400, "+")), LINE, (200, 299), (199, 300)),
    ((("D", 100, 200, "-"), ("U", 300, 400, "-")), LINE, (200, 299), (199, 300)),
    # rings: the wrap gap [990, 1020); a third gene M makes the wrap pair a row of its own
    ((("D", 20, 100, "+"), ("M", 400, 600, "+"), ("U", 900, 990, "+")), RING, (990, 19), (989, 20)),
    ((("U", 20, 100, "-"), ("M", 400, 600, "-"), ("D", 900, 990, "-")), RING, (990, 19), (989, 20)),
]


@pytest.mark.parametrize(("layout", "replicons", "inside", "outside"), TX_CASES)
def test_tx_rate_gap_is_half_open_on_both_strands(
    layout: tuple, replicons: dict, inside: tuple[int, int], outside: tuple[int, int]
) -> None:
    gs = genes_of(layout)
    by = {g.locus_tag: g for g in gs}
    strand = gs[0].interval.strand
    tx = pd.DataFrame(
        {"contig": "c1", "strand": strand, "tss": [*inside, *outside], "tx_rate": [1, 2, 8, 9]}
    ).astype({"tss": "int64", "tx_rate": "float64"})
    probs = {(a.locus_tag, b.locus_tag): 0.6 for a, b in adjacent_pairs(gs, replicons)}
    assert tx_rate_between(by["U"], by["D"], tx, replicons["c1"]) == 2.0
    df = frame(gs, replicons, probs, ev=evidence(tx=tx))
    assert row(df, "U", "D")["tx_rate_between"] == 2.0


def test_promoter_window_far_upstream_of_down_is_not_its_window() -> None:
    gs = genes_of([("U", 100, 200, "+"), ("D", 1000, 1300, "+")])
    # down's windows: end3p in [1000 - 200, 1000 + 20]; 500 is after U's end but not D's window
    w = pd.DataFrame({"contig": "c1", "strand": "+", "end3p": [500, 900], "score": [0.9, 0.3]})
    assert promoter_between(gs[0], gs[1], w, LINE["c1"]) == 0.3
    df = frame(gs, LINE, {("U", "D"): 0.6}, ev=evidence(w))
    assert row(df, "U", "D")["promoter_between"] == 0.3


@pytest.mark.parametrize("minus", [False, True])
def test_terminator_in_the_upstream_slack_counts(minus: bool) -> None:
    # + : window [300 - 50, 400); the call starts 45 bp before U's end
    # - : mirrored, window [600, 700 + 50); the call ends 45 bp after U's start
    layout = [("U", 100, 300, "+"), ("D", 400, 600, "+")]
    iv = Interval("c1", 255, 290, "+")
    if minus:
        layout = [("D", 400, 600, "-"), ("U", 700, 900, "-")]
        iv = Interval("c1", 710, 745, "-")
    gs = genes_of(layout)
    term = Terminator(iv, 0.7, "TerminatorNet")
    probs = {(gs[0].locus_tag, gs[1].locus_tag): 0.6}
    df = frame(gs, LINE, probs, ev=evidence(terminators=[term]))
    assert row(df, "U", "D")["terminator_between"] == 0.7


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(("start", "end"), [(990, 1010), (999, 1020), (10, 30)])
def test_terminator_prefilter_across_the_origin(strand: str, start: int, end: int) -> None:
    # the wrap pair of A [100, 300), B [700, 900) on a 1 kb ring (M on the other strand): the
    # gap [900, 1100) is longer than twice the pre-filter margin; calls at, on and after the origin
    other = "-" if strand == "+" else "+"
    gs = genes_of([("A", 100, 300, strand), ("M", 400, 500, other), ("B", 700, 900, strand)])
    term = Terminator(Interval("c1", start, end, strand), 0.7, "TerminatorNet")
    df = frame(gs, RING, {("B", "A"): 0.6}, ev=evidence(terminators=[term]))
    up, down = ("B", "A") if strand == "+" else ("A", "B")
    assert row(df, up, down)["terminator_between"] == 0.7
