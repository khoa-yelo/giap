"""TerminatorNet windows, matching and the coincidence null (issue M6-12; spec section 4,
examples 6-7; methods docs/methods/operons.md "K-12 cut-off calibration (M6-12)" > TerminatorNet,
B2).

Interface pinned here (``giap.stats.cutoffs``, pure, 0-based half-open coordinates):

- ``Window``: a named tuple ``(gene, strand, lo, hi)``, unrolled across the origin on a circular
  replicon; ``lo >= hi`` is an empty window;
- ``downstream_windows(genes, replicon, slack_bp=50) -> list[Window]``: one window per gene, in
  genome order whatever the input order, via ``giap.annotate.terminators._window`` with the
  genome-order neighbour on **either** strand; a linear-contig gene without a neighbour gets
  ``DOWNSTREAM_NO_NEXT_BP`` (300) past its stop, clipped to the contig;
- ``window_of(iv, windows, replicon) -> str | None``: the locus tag of the first same-strand
  window wholly containing ``iv`` for a shift in {0, +-length} (circular), else None;
- ``interval_gap(a, b, replicon) -> int``: min over the shifts of
  ``max(0, b.start - a.end, a.start - b.end)``; overlap is 0;
- ``terminator_match(iv, truth, replicon, tol_bp=50) -> int | None``: the index of the nearest
  same-contig, same-strand truth interval with gap <= ``tol_bp`` (lowest index on ties);
- ``coincidence_rate(calls, confidence, windows, truth, replicon, threshold, *, rng, R=200,
  tol_bp=50) -> float``: B2 placement null; each eligible call with confidence >= threshold moves
  to a uniformly drawn same-strand window long enough to hold it, keeping its length and its
  offset from the window's 3' side (``lo`` for +, ``hi`` for -); the mean over R placements of the
  matched share. NaN when no call is at or above the threshold;
- ``precision_lower_bound(obs, c) -> float``: ``max(0, (obs - c)/(1 - c))``; NaN when c >= 1.
"""

from __future__ import annotations

import math
import random

import numpy as np
import pytest

from giap.annotate.operons import intergenic_distance
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats.cutoffs import (
    DEFAULT_SEED,
    STEP,
    Window,
    choose_cutoff,
    coincidence_rate,
    downstream_windows,
    interval_gap,
    precision_lower_bound,
    terminator_match,
    threshold_curve,
    window_of,
)

C = "chr1"
LIN = Replicon(C, 10_000, "linear")
CIRC = Replicon(C, 10_000, "circular")


def gene(tag: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(C, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def iv(start: int, end: int, strand: str) -> Interval:
    return Interval(C, start, end, strand)  # type: ignore[arg-type]


def rng() -> np.random.Generator:
    return for_step(STEP, DEFAULT_SEED)


# --- example 6: linear 10 kb contig ----------------------------------------------------------

EX6_GENES = [
    gene("g1", 100, 1000, "+"),
    gene("g2", 1500, 2500, "+"),
    gene("g3", 3000, 4000, "-"),
    gene("g4", 5000, 6000, "-"),
]
EX6_CALLS = {
    "c1": (iv(1100, 1140, "+"), 0.95),
    "c2": (iv(2460, 2490, "+"), 0.90),
    "c3": (iv(2600, 2640, "-"), 0.85),
    "c4": (iv(3500, 3540, "+"), 0.99),
    "c5": (iv(4100, 4130, "-"), 0.60),
    "c6": (iv(1460, 1520, "+"), 0.99),
}
EX6_TRUTH = [iv(1120, 1150, "+"), iv(2560, 2580, "+"), iv(2680, 2700, "-"), iv(4100, 4130, "+")]


def test_example6_windows() -> None:
    w = downstream_windows(EX6_GENES, LIN)
    assert w == [
        Window("g1", "+", 950, 1500),
        Window("g2", "+", 2450, 3000),  # ends at the opposite-strand neighbour g3
        Window("g3", "-", 2500, 3050),  # minus: [neighbour.end, g.start + 50)
        Window("g4", "-", 4000, 5050),
    ]
    assert all(isinstance(x, Window) for x in w)
    assert (w[0].gene, w[0].strand, w[0].lo, w[0].hi) == ("g1", "+", 950, 1500)


def test_windows_genome_order_whatever_input_order() -> None:
    shuffled = EX6_GENES[:]
    random.Random(4).shuffle(shuffled)
    assert downstream_windows(shuffled, LIN) == downstream_windows(EX6_GENES, LIN)


def test_example6_eligibility() -> None:
    w = downstream_windows(EX6_GENES, LIN)
    got = {k: window_of(c, w, LIN) for k, (c, _) in EX6_CALLS.items()}
    assert got == {"c1": "g1", "c2": "g2", "c3": "g3", "c4": None, "c5": "g4", "c6": None}


def test_example6_matches_gaps_and_strand() -> None:
    calls = {k: c for k, (c, _) in EX6_CALLS.items()}
    assert interval_gap(calls["c1"], EX6_TRUTH[0], LIN) == 0  # overlap
    assert interval_gap(calls["c2"], EX6_TRUTH[1], LIN) == 70
    assert interval_gap(calls["c3"], EX6_TRUTH[2], LIN) == 40
    assert terminator_match(calls["c1"], EX6_TRUTH, LIN) == 0
    assert terminator_match(calls["c2"], EX6_TRUTH, LIN) is None  # 70 > 50
    assert terminator_match(calls["c3"], EX6_TRUTH, LIN) == 2  # minus strand
    assert terminator_match(calls["c5"], EX6_TRUTH, LIN) is None  # T4 is on the other strand
    assert interval_gap(calls["c5"], EX6_TRUTH[3], LIN) == 0  # same place, wrong strand


def test_example6_cutoff_precision_and_recall() -> None:
    w = downstream_windows(EX6_GENES, LIN)
    elig = [(c, p) for c, p in EX6_CALLS.values() if window_of(c, w, LIN) is not None]
    conf = [p for _, p in elig]
    y = [int(terminator_match(c, EX6_TRUTH, LIN) is not None) for c, _ in elig]
    assert sorted(zip(conf, y, strict=True), reverse=True) == [
        (0.95, 1),
        (0.90, 0),
        (0.85, 1),
        (0.60, 0),
    ]
    curve = threshold_curve(conf, y)
    assert curve["tp"].tolist() == [1, 1, 2, 2]
    assert curve["fp"].tolist() == [0, 1, 1, 2]
    cut = choose_cutoff(curve)
    assert cut.threshold == 0.95
    assert cut.rule == "precision_0.90"
    assert cut.precision == 1.0
    truth_elig = [t for t in EX6_TRUTH if window_of(t, w, LIN) is not None]
    assert truth_elig == EX6_TRUTH[:3]  # T4 lies in no + window
    called = [c for c, p in elig if p >= cut.threshold]
    hit = [t for t in truth_elig if any(terminator_match(c, [t], LIN) is not None for c in called)]
    assert len(hit) / len(truth_elig) == pytest.approx(1 / 3)
    assert len(hit) / len(EX6_TRUTH) == pytest.approx(1 / 4)  # recall_all


# --- example 7: the circular origin wrap ------------------------------------------------------


def test_example7_circular_wrap_window_eligibility_match_and_distance() -> None:
    ga, gb = gene("gA", 200, 800, "+"), gene("gB", 9500, 9900, "+")
    w = downstream_windows([ga, gb], CIRC)
    # Neighbour wraps to gA, unrolled. The spec's example 7 prints [9450, 10200); gB's stop is
    # 9900, so the window starts at 9900 - 50 = 9850 (the spec's lo is a typo; reported).
    assert Window("gB", "+", 9850, 10200) in w
    assert Window("gA", "+", 750, 9500) in w
    call = iv(20, 60, "+")
    assert window_of(call, w, CIRC) == "gB"  # shifted by +10,000: [10020, 10060)
    truth = [iv(80, 100, "+")]
    assert interval_gap(call, truth[0], CIRC) == 20
    assert terminator_match(call, truth, CIRC) == 0
    assert intergenic_distance(gb, ga, CIRC) == 300


def test_circular_minus_strand_wrap() -> None:
    first, last = gene("m1", 100, 500, "-"), gene("z9", 9000, 9600, "-")
    w = downstream_windows([last, first], CIRC)
    assert Window("m1", "-", 9600, 10150) in w  # previous gene wraps to z9
    assert window_of(iv(9700, 9750, "-"), w, CIRC) == "m1"
    assert window_of(iv(40, 90, "-"), w, CIRC) == "m1"  # [10040, 10090) unrolled
    assert window_of(iv(40, 90, "+"), w, CIRC) is None  # strand must agree
    assert window_of(iv(120, 160, "-"), w, CIRC) is None  # past g.start + 50


def test_gap_across_the_origin() -> None:
    a, b = iv(9980, 9995, "+"), iv(10, 30, "+")
    assert interval_gap(a, b, CIRC) == 15
    assert interval_gap(b, a, CIRC) == 15
    assert terminator_match(a, [b], CIRC) == 0
    assert interval_gap(a, b, LIN) == 9950  # no wrap on a linear contig
    assert terminator_match(a, [b], LIN) is None


# --- window edge cases ------------------------------------------------------------------------


def test_linear_contig_ends_get_300_bp_clipped() -> None:
    w = downstream_windows([gene("a", 200, 900, "-"), gene("b", 9000, 9800, "+")], LIN)
    assert w == [Window("a", "-", 0, 250), Window("b", "+", 9750, 10_000)]
    w2 = downstream_windows([gene("a", 1000, 1900, "-"), gene("b", 5000, 5600, "+")], LIN)
    assert w2 == [Window("a", "-", 700, 1050), Window("b", "+", 5550, 5900)]


def test_overlapping_neighbour_gives_empty_window() -> None:
    w = downstream_windows([gene("a", 100, 1000, "+"), gene("b", 900, 2000, "+")], LIN)
    assert w[0] == Window("a", "+", 950, 900)
    assert window_of(iv(950, 990, "+"), w, LIN) is None
    assert window_of(iv(900, 905, "+"), w, LIN) is None


def test_touching_neighbour_leaves_only_the_50bp_tail() -> None:
    w = downstream_windows([gene("a", 100, 1000, "+"), gene("b", 1000, 2000, "+")], LIN)
    assert w[0] == Window("a", "+", 950, 1000)
    assert window_of(iv(960, 1000, "+"), w, LIN) == "a"
    assert window_of(iv(960, 1001, "+"), w, LIN) is None  # crosses into b


def test_gap_boundary_is_inclusive_at_50() -> None:
    c = iv(1000, 1030, "+")
    assert terminator_match(c, [iv(1080, 1100, "+")], LIN) == 0  # gap 50
    assert terminator_match(c, [iv(1081, 1100, "+")], LIN) is None  # gap 51
    assert terminator_match(c, [iv(940, 950, "+")], LIN) == 0  # gap 50 on the left
    assert terminator_match(c, [iv(1081, 1100, "+")], LIN, tol_bp=51) == 0


def test_match_picks_the_nearest_truth() -> None:
    c = iv(1000, 1030, "+")
    truth = [iv(1070, 1090, "+"), iv(990, 1010, "-"), iv(1035, 1050, "+"), iv(1031, 1040, "+")]
    assert terminator_match(c, truth, LIN) == 3  # gap 1 beats 5 and 40; the minus one is skipped
    assert terminator_match(c, [iv(1040, 1050, "+"), iv(980, 990, "+")], LIN) == 0  # tie: first


def test_other_contig_never_matches() -> None:
    assert terminator_match(iv(10, 40, "+"), [Interval("chr2", 10, 40, "+")], LIN) is None


# --- the coincidence null (B2) ----------------------------------------------------------------

# + windows: gA [950, 1500) (550 bp), gB [2350, 2950) (600 bp).
# - windows: gC [2400, 3000), h1 [4000, 5050), h2 [6000, 6550).
NULL_GENES = [
    gene("gA", 100, 1000, "+"),
    gene("gB", 1500, 2400, "+"),
    gene("gC", 2950, 4000, "-"),
    gene("h1", 5000, 6000, "-"),
    gene("h2", 6500, 7000, "-"),
]
R_NULL = 600


def _windows() -> list[Window]:
    w = downstream_windows(NULL_GENES, LIN)
    assert [(x.gene, x.lo, x.hi) for x in w] == [
        ("gA", 950, 1500),
        ("gB", 2350, 2950),
        ("gC", 2400, 3000),
        ("h1", 4000, 5050),
        ("h2", 6000, 6550),
    ]
    return w


def test_coincidence_keeps_the_offset_from_lo_on_plus() -> None:
    """The call sits on the truth at offset 50 in gA; placed into gA (p = 1/2) it matches, into gB
    ([2400, 2430)) it does not. Uniform-in-window placement would give far less than 1/2."""
    w = _windows()
    c = coincidence_rate(
        [iv(1000, 1030, "+")], [0.9], w, [iv(1000, 1030, "+")], LIN, 0.5, rng=rng(), R=R_NULL
    )
    assert c == pytest.approx(0.5, abs=0.08)


def test_coincidence_keeps_the_offset_from_hi_on_minus() -> None:
    """Offset 30 from h1's hi; placed into gC it lands at [2950, 2970), into h2 at [6500, 6520):
    only h1 matches, so c = 1/3. Measuring from lo (offset 1000) would leave only h1 long enough
    and give 1."""
    w = _windows()
    truth = [iv(5000, 5020, "-")]
    c = coincidence_rate([iv(5000, 5020, "-")], [0.7], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == pytest.approx(1 / 3, abs=0.08)


def test_coincidence_skips_windows_too_short_for_the_call() -> None:
    """Offset 560 + length 30 = 590 > 550: gA cannot hold the call, so it always stays in gB."""
    w = _windows()
    truth = [iv(2920, 2925, "+")]
    c = coincidence_rate([iv(2910, 2940, "+")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=50)
    assert c == 1.0


def test_coincidence_threshold_filters_calls() -> None:
    w = _windows()
    calls = [iv(1000, 1030, "+"), iv(2910, 2940, "+")]
    truth = [iv(2920, 2925, "+")]
    assert coincidence_rate(calls, [0.4, 0.9], w, truth, LIN, 0.5, rng=rng(), R=50) == 1.0
    assert math.isnan(coincidence_rate(calls, [0.4, 0.45], w, truth, LIN, 0.5, rng=rng(), R=50))


def test_coincidence_deterministic_for_the_seed() -> None:
    w = _windows()
    calls = [iv(1000, 1030, "+"), iv(5000, 5020, "-")]
    truth = [iv(1000, 1030, "+"), iv(5000, 5020, "-")]
    a = coincidence_rate(calls, [0.9, 0.8], w, truth, LIN, 0.0, rng=rng(), R=100)
    b = coincidence_rate(calls, [0.9, 0.8], w, truth, LIN, 0.0, rng=rng(), R=100)
    assert a == b


def test_coincidence_default_r_is_200() -> None:
    w = _windows()
    g1, g2 = rng(), rng()
    a = coincidence_rate([iv(1000, 1030, "+")], [0.9], w, [iv(1000, 1030, "+")], LIN, 0.5, rng=g1)
    b = coincidence_rate(
        [iv(1000, 1030, "+")], [0.9], w, [iv(1000, 1030, "+")], LIN, 0.5, rng=g2, R=200
    )
    assert a == b


@pytest.mark.parametrize(
    ("obs", "c", "want"),
    [(0.226, 0.06, (0.226 - 0.06) / 0.94), (1.0, 0.3, 1.0), (0.05, 0.1, 0.0), (0.5, 0.0, 0.5)],
)
def test_precision_lower_bound(obs: float, c: float, want: float) -> None:
    assert precision_lower_bound(obs, c) == pytest.approx(want, abs=1e-15)


def test_precision_lower_bound_undefined_at_c_one() -> None:
    assert math.isnan(precision_lower_bound(0.5, 1.0))
