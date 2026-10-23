"""Review gaps of PR #263 (M6-12 part A, ``giap.stats.cutoffs``): blockers B1-B4 of the code
review and item 1 of the statistician's review.

Spec docs/specs/M6-12-uniop-terminator-cutoffs.md (lead decisions 6-8); methods
docs/methods/operons.md "K-12 cut-off calibration (M6-12)" > TerminatorNet, B2, and R1.

- B1: the coincidence null keeps each call's offset from the window's 3' side; the offsets here
  (170-180 bp) are far above ``tol_bp`` + the call length, so a dropped or sign-flipped offset
  moves the placed call out of match range.
- B2: calls eligible only through a +-length shift, and origin-crossing intervals (end > length)
  in ``downstream_windows``, ``window_of`` and ``interval_gap``.
- B3: the lone gene on a circular replicon is its own neighbour (decision 8); a called terminator
  outside any window raises in ``coincidence_rate`` (decision 8).
- B4: containment is inclusive at the window's ``lo``; a window exactly as long as offset +
  length still holds the call in the null.
- Statistician item 1: ``precision_cv`` against a reference rebuilt from the same rng stream
  (folds ``rng.permutation(n_clusters) % k`` over clusters, held-out TP/FP pooled).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

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
    precision_cv,
    terminator_match,
    threshold_curve,
    window_of,
)

C = "chr1"
LIN = Replicon(C, 10_000, "linear")
CIRC = Replicon(C, 10_000, "circular")
R_NULL = 600


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


# + windows: gA [950, 1500) (550 bp), gB [2350, 2950) (600 bp).
# - windows: gC [2400, 3000) (600 bp), h1 [4000, 5050) (1050 bp), h2 [6000, 6550) (550 bp).
NULL_GENES = [
    gene("gA", 100, 1000, "+"),
    gene("gB", 1500, 2400, "+"),
    gene("gC", 2950, 4000, "-"),
    gene("h1", 5000, 6000, "-"),
    gene("h2", 6500, 7000, "-"),
]


def null_windows() -> list[Window]:
    w = downstream_windows(NULL_GENES, LIN)
    assert [(x.gene, x.strand, x.lo, x.hi) for x in w] == [
        ("gA", "+", 950, 1500),
        ("gB", "+", 2350, 2950),
        ("gC", "-", 2400, 3000),
        ("h1", "-", 4000, 5050),
        ("h2", "-", 6000, 6550),
    ]
    return w


# --- B1: offsets well above tol_bp on both strands -------------------------------------------


def test_b1_plus_offset_180_from_lo_is_kept_exactly() -> None:
    """Call [1130, 1160) sits 180 bp past gA's lo. Placed into gA it is [1130, 1160), into gB
    [2530, 2560); truth sits at both, so c = 1 exactly. Without the offset the call lands flush at
    lo ([950, 980) / [2350, 2380)), 150 bp from either truth: c = 0."""
    w = null_windows()
    truth = [iv(1130, 1160, "+"), iv(2530, 2560, "+")]
    c = coincidence_rate([iv(1130, 1160, "+")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == 1.0


def test_b1_minus_offset_180_from_hi_is_kept_exactly() -> None:
    """Call [4840, 4870) ends 180 bp before h1's hi (5050). Placed into gC it is [2790, 2820),
    into h1 [4840, 4870), into h2 [6340, 6370); truth sits at all three, so c = 1 exactly.
    Dropping the offset puts it flush at hi (150 bp from each truth) and flipping its sign puts it
    330 bp past hi: c = 0 either way."""
    w = null_windows()
    truth = [iv(2790, 2820, "-"), iv(4840, 4870, "-"), iv(6340, 6370, "-")]
    c = coincidence_rate([iv(4840, 4870, "-")], [0.7], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == 1.0


def test_b1_minus_offset_only_one_window_matches() -> None:
    """Same call, truth only at gC's placement [2790, 2820): c is the share of draws into gC, 1/3
    (all three - windows hold 180 + 30 bp)."""
    w = null_windows()
    truth = [iv(2790, 2820, "-")]
    c = coincidence_rate([iv(4840, 4870, "-")], [0.7], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == pytest.approx(1 / 3, abs=0.08)


def test_b1_plus_offset_only_one_window_matches() -> None:
    """Call 180 bp past gA's lo, truth only at gB's placement [2530, 2560): c = 1/2."""
    w = null_windows()
    truth = [iv(2530, 2560, "+")]
    c = coincidence_rate([iv(1130, 1160, "+")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == pytest.approx(0.5, abs=0.08)


# --- B2: wrap-shifted calls in the null and origin-crossing intervals -------------------------


def test_b2_null_minus_call_eligible_only_through_the_shift() -> None:
    """Genes m1 - [300, 700) and z9 - [9000, 9600) on CIRC: m1's window wraps to z9, unrolled
    [9600, 10350); z9's is [700, 9050). The call [150, 180) is in m1 only as [10150, 10180),
    170 bp before hi. Placements: m1 [10150, 10180) (= [150, 180)), z9 [8850, 8880). Truth at
    both: c = 1 exactly. Ignoring the shift gives an offset of 10170, which no window holds."""
    genes = [gene("m1", 300, 700, "-"), gene("z9", 9000, 9600, "-")]
    w = downstream_windows(genes, CIRC)
    assert w == [Window("m1", "-", 9600, 10350), Window("z9", "-", 700, 9050)]
    call = iv(150, 180, "-")
    assert window_of(call, w, CIRC) == "m1"
    truth = [iv(150, 180, "-"), iv(8850, 8880, "-")]
    c = coincidence_rate([call], [0.9], w, truth, CIRC, 0.5, rng=rng(), R=R_NULL)
    assert c == 1.0


def test_b2_null_plus_call_eligible_only_through_the_shift() -> None:
    """Genes ga + [200, 800), gq + [900, 1500), gb + [9500, 9900) on CIRC. Windows: ga [750, 900)
    (150 bp), gq [1450, 9500), gb [9850, 10200) (unrolled). The call [20, 60) is in gb only as
    [10020, 10060), 170 bp past lo, so offset + length = 210 and ga is too short. Placements:
    gq [1620, 1660), gb [10020, 10060) (= [20, 60)); truth at both gives c = 1 exactly. Ignoring
    the shift makes the offset negative, so ga also qualifies and its draws never match."""
    genes = [gene("ga", 200, 800, "+"), gene("gq", 900, 1500, "+"), gene("gb", 9500, 9900, "+")]
    w = downstream_windows(genes, CIRC)
    assert w == [
        Window("ga", "+", 750, 900),
        Window("gq", "+", 1450, 9500),
        Window("gb", "+", 9850, 10200),
    ]
    call = iv(20, 60, "+")
    assert window_of(call, w, CIRC) == "gb"
    truth = [iv(20, 60, "+"), iv(1620, 1660, "+")]
    c = coincidence_rate([call], [0.9], w, truth, CIRC, 0.5, rng=rng(), R=R_NULL)
    assert c == 1.0


def test_b2_null_plus_shifted_call_only_its_own_window_matches() -> None:
    """Same genome, truth only at [20, 60): of the two windows that hold the call, gb matches,
    so c = 1/2."""
    genes = [gene("ga", 200, 800, "+"), gene("gq", 900, 1500, "+"), gene("gb", 9500, 9900, "+")]
    w = downstream_windows(genes, CIRC)
    truth = [iv(20, 60, "+")]
    c = coincidence_rate([iv(20, 60, "+")], [0.9], w, truth, CIRC, 0.5, rng=rng(), R=R_NULL)
    assert c == pytest.approx(0.5, abs=0.08)


def test_b2_origin_crossing_gene_window() -> None:
    """x + [9800, 10200) crosses the origin; its genome-order neighbour wraps to a + [500, 900),
    so its window is [10150, 10500); a's window ends at x's start."""
    genes = [gene("x", 9800, 10_200, "+"), gene("a", 500, 900, "+")]
    w = downstream_windows(genes, CIRC)
    assert w == [Window("a", "+", 850, 9800), Window("x", "+", 10_150, 10_500)]


def test_b2_origin_crossing_minus_gene_window() -> None:
    """y - [9900, 10300) crosses the origin; its previous gene b - [9000, 9500) gives y the window
    [9500, 9950). b's previous gene wraps to y, so b's window runs from y's end to b's start + 50,
    unrolled one turn: [10300, 19050)."""
    genes = [gene("b", 9000, 9500, "-"), gene("y", 9900, 10_300, "-")]
    w = downstream_windows(genes, CIRC)
    assert w == [Window("b", "-", 10_300, 19_050), Window("y", "-", 9500, 9950)]
    assert window_of(iv(9600, 9650, "-"), w, CIRC) == "y"
    assert window_of(iv(400, 450, "-"), w, CIRC) == "b"  # [10400, 10450) after +10,000
    assert window_of(iv(250, 290, "-"), w, CIRC) is None  # inside y's body past the origin


def test_b2_window_of_origin_crossing_and_shifted_calls() -> None:
    genes = [gene("x", 9800, 10_200, "+"), gene("a", 500, 900, "+")]
    w = downstream_windows(genes, CIRC)
    assert window_of(iv(300, 340, "+"), w, CIRC) == "x"  # [10300, 10340) after +10,000
    assert window_of(iv(10_160, 10_190, "+"), w, CIRC) == "x"  # end > length, no shift
    assert window_of(iv(10_480, 10_520, "+"), w, CIRC) is None  # crosses a's start (10,500)
    assert window_of(iv(470, 500, "+"), w, CIRC) == "x"  # flush at hi after the shift
    # A minus call crossing the origin inside m1's unrolled window [9600, 10150).
    m = downstream_windows([gene("m1", 100, 500, "-"), gene("z9", 9000, 9600, "-")], CIRC)
    assert window_of(iv(9990, 10_030, "-"), m, CIRC) == "m1"
    assert window_of(iv(9990, 10_030, "+"), m, CIRC) is None


def test_b2_interval_gap_with_origin_crossing_intervals() -> None:
    a = iv(9990, 10_010, "+")
    assert interval_gap(a, iv(30, 40, "+"), CIRC) == 20
    assert interval_gap(iv(30, 40, "+"), a, CIRC) == 20
    assert interval_gap(a, iv(5, 15, "+"), CIRC) == 0  # [10005, 10015) overlaps
    assert interval_gap(a, iv(9900, 9950, "+"), CIRC) == 40
    assert interval_gap(iv(10_100, 10_150, "+"), iv(20, 30, "+"), CIRC) == 70
    assert terminator_match(a, [iv(70, 80, "+"), iv(30, 40, "+")], CIRC) == 1  # 20 beats 60


# --- B3: decision 8 readings ------------------------------------------------------------------


def test_b3_lone_plus_gene_on_a_circle_is_its_own_neighbour() -> None:
    w = downstream_windows([gene("g", 1000, 2000, "+")], CIRC)
    assert w == [Window("g", "+", 1950, 11_000)]
    assert window_of(iv(100, 130, "+"), w, CIRC) == "g"  # [10100, 10130) unrolled


def test_b3_lone_minus_gene_on_a_circle_is_its_own_neighbour() -> None:
    w = downstream_windows([gene("g", 1000, 2000, "-")], CIRC)
    assert w == [Window("g", "-", 2000, 11_050)]
    assert window_of(iv(1020, 1050, "-"), w, CIRC) == "g"  # [11020, 11050), flush at hi


def test_b3_lone_gene_on_a_linear_contig_is_not_wrapped() -> None:
    w = downstream_windows([gene("g", 1000, 2000, "+")], LIN)
    assert w == [Window("g", "+", 1950, 2300)]


def test_b3_coincidence_raises_on_a_call_in_no_window() -> None:
    """c4-like call: + [3500, 3540) lies in neither + window (gA, gB)."""
    w = null_windows()
    bad = iv(3500, 3540, "+")
    assert window_of(bad, w, LIN) is None
    with pytest.raises(ValueError, match="no window"):
        coincidence_rate([bad], [0.9], w, [iv(3500, 3540, "+")], LIN, 0.5, rng=rng(), R=10)
    with pytest.raises(ValueError, match="no window"):
        coincidence_rate(
            [iv(1000, 1030, "+"), bad],
            [0.9, 0.95],
            w,
            [iv(1000, 1030, "+")],
            LIN,
            0.5,
            rng=rng(),
            R=10,
        )


def test_b3_coincidence_raises_on_a_minus_call_in_no_window() -> None:
    w = null_windows()
    bad = iv(3500, 3540, "-")  # inside gC's body, past gC's hi (3000)
    with pytest.raises(ValueError, match="no window"):
        coincidence_rate([bad], [0.9], w, [bad], LIN, 0.5, rng=rng(), R=10)


# --- B4: half-open boundaries -----------------------------------------------------------------


def test_b4_containment_is_inclusive_at_lo() -> None:
    w = null_windows()
    assert window_of(iv(950, 980, "+"), w, LIN) == "gA"  # start == lo
    assert window_of(iv(949, 980, "+"), w, LIN) is None
    assert window_of(iv(2400, 2430, "-"), w, LIN) == "gC"  # minus window lo
    assert window_of(iv(2399, 2430, "-"), w, LIN) is None


def test_b4_containment_inclusive_at_lo_after_a_shift() -> None:
    m = downstream_windows([gene("m1", 100, 500, "-"), gene("z9", 9000, 9600, "-")], CIRC)
    assert Window("m1", "-", 9600, 10_150) in m
    assert window_of(iv(9600, 9640, "-"), m, CIRC) == "m1"


def test_b4_null_counts_a_window_that_fits_exactly() -> None:
    """Call [2870, 2900) sits 520 bp past gB's lo; offset + length = 550 = gA's length, so gA
    holds it exactly ([1470, 1500), flush at hi). Truth only at that placement: c = 1/2. Excluding
    the exact fit would leave only gB and give 0."""
    w = null_windows()
    truth = [iv(1470, 1500, "+")]
    c = coincidence_rate([iv(2870, 2900, "+")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=R_NULL)
    assert c == pytest.approx(0.5, abs=0.08)


def test_b4_null_drops_a_window_one_bp_too_short() -> None:
    """One bp further (offset 521): gA no longer holds it, so the call always stays in gB."""
    w = null_windows()
    truth = [iv(2871, 2901, "+")]
    c = coincidence_rate([iv(2871, 2901, "+")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=50)
    assert c == 1.0


def test_b4_null_minus_exact_fit() -> None:
    """Minus: call [4000, 4030) is flush at h1's lo, offset 1020 from hi; + 30 = 1050 = h1's
    length, which only h1 has. It must stay eligible and always land in h1."""
    w = null_windows()
    truth = [iv(4000, 4030, "-")]
    c = coincidence_rate([iv(4000, 4030, "-")], [0.9], w, truth, LIN, 0.5, rng=rng(), R=50)
    assert c == 1.0


# --- statistician item 1: precision_cv against a reference ------------------------------------


def _cv_data(n_dir: int = 60, seed: int = 11) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Directons of 1-6 rows with a shared directon effect, so held-out and in-sample precision
    differ. Clusters are contiguous integer labels (first appearance == sorted order)."""
    g = np.random.default_rng(seed)
    s: list[float] = []
    y: list[int] = []
    cl: list[int] = []
    for c in range(n_dir):
        size = int(g.integers(1, 7))
        u = g.normal(scale=1.0)
        p = 1.0 / (1.0 + math.exp(-1.5 * u))
        for _ in range(size):
            v = int(g.random() < p)
            s.append(float(np.round(u + 1.0 * v + g.normal(scale=0.8), 2)))
            y.append(v)
            cl.append(c)
    return np.asarray(s), np.asarray(y), np.asarray(cl)


def _pooled(s: np.ndarray, y: np.ndarray, fold: np.ndarray, k: int, held_out: bool) -> float:
    tp = fp = 0
    for j in range(k):
        train = fold != j
        t = choose_cutoff(threshold_curve(s[train], y[train]), 0.90).threshold
        rows = (~train if held_out else train) & (s >= t)
        tp += int(y[rows].sum())
        fp += int((y[rows] == 0).sum())
    return tp / (tp + fp)


def _reference(s: np.ndarray, y: np.ndarray, cl: np.ndarray, k: int) -> tuple[float, ...]:
    labels, inv = np.unique(cl, return_inverse=True)
    cluster_fold = (rng().permutation(labels.size) % k)[inv]
    row_fold = rng().permutation(s.size) % k
    return (
        _pooled(s, y, cluster_fold, k, held_out=True),
        _pooled(s, y, row_fold, k, held_out=True),
        _pooled(s, y, cluster_fold, k, held_out=False),
    )


@pytest.mark.parametrize("k", [10, 5])
def test_precision_cv_matches_the_cluster_fold_reference(k: int) -> None:
    s, y, cl = _cv_data()
    want, rows, in_sample = _reference(s, y, cl, k)
    # Teeth: row folds and in-sample counting give different numbers on these data.
    assert rows != want
    assert in_sample != want
    got = precision_cv(s, y, cl, rng=rng(), k=k)
    assert got == want


def test_precision_cv_reference_with_shuffled_rows() -> None:
    """Rows of one directon are not contiguous; the folds still follow the cluster labels."""
    s, y, cl = _cv_data(seed=5)
    order = np.random.default_rng(2).permutation(s.size)
    s, y, cl = s[order], y[order], cl[order]
    want, rows, in_sample = _reference(s, y, cl, 10)
    assert rows != want
    assert in_sample != want
    assert precision_cv(s, y, cl, rng=rng()) == want
