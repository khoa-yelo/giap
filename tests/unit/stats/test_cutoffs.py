"""UniOP / TerminatorNet cut-off statistics (issue M6-12 part A; GIPPA_SPEC.md §15.4 steps 3-5,
§7.1.3, §7.4.6; spec docs/specs/M6-12-uniop-terminator-cutoffs.md with lead decisions 1-6;
methods docs/methods/operons.md "K-12 cut-off calibration (M6-12)", the authority where it is
more specific than the spec).

Interface pinned here (``gippa.stats.cutoffs``, pure functions, no I/O):

- constants ``STEP = "m6_12_cutoffs"``, ``DEFAULT_SEED = 20261014``, ``DEFAULT_B = 2000``,
  ``MIN_LABELLED = 500``, ``D0_DEFAULT = 50.0``, ``S_DEFAULT = 20.0``;
- ``require_labelled(y, min_n=MIN_LABELLED) -> (P, N)``: ``ValueError`` below ``min_n`` items or
  with P = 0 or N = 0 (the "< 500 labelled pairs" guard on the full data);
- ``auroc(scores, y) -> float``: Mann-Whitney through average ranks, ties count 1/2; NaN when
  P = 0 or N = 0 (so it composes with ``tu_block_bootstrap``, which drops NaN replicates);
- ``threshold_curve(scores, y) -> pd.DataFrame`` with integer columns ``tp``, ``fp`` and the float
  column ``threshold`` (the distinct observed scores, descending; ties enter together);
- ``choose_cutoff(curve, precision=0.90) -> Cutoff`` with fields ``threshold`` (a Python float,
  the observed score itself, never rounded), ``rule`` (``"precision_0.90"`` | ``"max_f1"``),
  ``precision``, ``recall``, ``f1``, ``tp``, ``fp`` and ``literal_threshold`` (the literal smallest
  candidate with precision >= target, or None if none reaches it). Stable reading (lead
  decision 1): the last candidate before the first failing one, in descending order; if the top
  candidate fails, max-F1 with the smallest t on ties. Integer test ``b*TP >= a*(TP+FP)``;
  F1 compared exactly;
- ``precision_at(scores, y, t)``, ``recall_at(scores, y, t)``: NaN on an empty denominator;
- directon-bootstrap wrappers (each is exactly a fixed sequence of ``tu_block_bootstrap`` calls
  on one DataFrame in input row order, no strata, so the methods' draw order is reproducible):
  ``auroc_ci`` (1 call), ``fixed_cutoff_ci`` (precision call, then recall call),
  ``auroc_diff_ci`` (1 call, AUROC(a) - AUROC(b) on the same replicate), ``distance_fit_ci``
  (d0 call, then s call: separate draws, methods "Bootstrap streams"), ``threshold_boot_q``
  (1 call; 2.5/50/97.5 % linear quantiles of the re-selected t*);
- ``precision_cv(scores, y, cluster, *, rng, k=10, precision=0.90) -> float`` (R1);
- ``wilson_ci(k, n, level=0.95) -> (lo, hi)``; (NaN, NaN) for n = 0;
- ``fit_distance_model(d, y, *, max_iter=100, tol=1e-10) -> DistanceFit`` (B1: exact separation
  check, start at (logit ybar, 0), step-halving Newton; failures never raise) with fields ``d0``,
  ``s``, ``beta0``, ``beta1``, ``converged``, ``n_iter``, ``reason`` and ``calibrated``;
- ``distance_log_loss(d, y, d0, s) -> float`` (mean negative log-likelihood, nats).

The TerminatorNet window / match / coincidence functions are in test_cutoffs_terminators.py.
"""

from __future__ import annotations

import dataclasses
import math
from fractions import Fraction

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import log_loss as sk_log_loss
from sklearn.metrics import roc_auc_score
from statsmodels.stats.proportion import proportion_confint

from gippa.rng import for_step
from gippa.stats.bootstrap import BootCI, tu_block_bootstrap
from gippa.stats.cutoffs import (
    D0_DEFAULT,
    DEFAULT_B,
    DEFAULT_SEED,
    MIN_LABELLED,
    S_DEFAULT,
    STEP,
    Cutoff,
    DistanceFit,
    auroc,
    auroc_ci,
    auroc_diff_ci,
    choose_cutoff,
    distance_fit_ci,
    distance_log_loss,
    fit_distance_model,
    fixed_cutoff_ci,
    precision_at,
    precision_cv,
    recall_at,
    require_labelled,
    threshold_boot_q,
    threshold_curve,
    wilson_ci,
)

NAN = float("nan")


def rng() -> np.random.Generator:
    return for_step(STEP, DEFAULT_SEED)


def split(pairs: list[tuple[float, int]]) -> tuple[list[float], list[int]]:
    return [p for p, _ in pairs], [y for _, y in pairs]


def cut_of(pairs: list[tuple[float, int]]) -> Cutoff:
    s, y = split(pairs)
    return choose_cutoff(threshold_curve(s, y))


def grid(xs: range, n: int, d0: float, s: float) -> tuple[np.ndarray, np.ndarray]:
    """``n`` pairs at each d, ``round(n * sigma)`` of them positive (spec example 5 / B1)."""
    d: list[int] = []
    y: list[int] = []
    for x in xs:
        k = round(n / (1.0 + math.exp((x - d0) / s)))
        d += [x] * n
        y += [1] * k + [0] * (n - k)
    return np.asarray(d), np.asarray(y)


# --- constants -------------------------------------------------------------------------------


def test_constants() -> None:
    assert STEP == "m6_12_cutoffs"
    assert DEFAULT_SEED == 20261014
    assert DEFAULT_B == 2000
    assert MIN_LABELLED == 500
    assert (D0_DEFAULT, S_DEFAULT) == (50.0, 20.0)


# --- AUROC -----------------------------------------------------------------------------------


def test_example1_auroc_tie_counts_half() -> None:
    a = auroc([0.9, 0.8, 0.5, 0.5, 0.2], [1, 1, 1, 0, 0])
    assert a == pytest.approx(5.5 / 6, abs=1e-15)
    assert type(a) is float


def test_example2_auroc() -> None:
    s, y = split(EX2)
    assert auroc(s, y) == pytest.approx(0.8125, abs=1e-15)


def test_auroc_matches_sklearn_on_heavy_ties() -> None:
    g = np.random.default_rng(7)
    s = np.round(g.random(400), 1)  # 11 distinct values
    y = (g.random(400) < 0.4).astype(int)
    assert auroc(s, y) == pytest.approx(roc_auc_score(y, s), abs=1e-12)


def test_auroc_perfect_and_reversed() -> None:
    assert auroc([3.0, 2.0, 1.0, 0.0], [1, 1, 0, 0]) == 1.0
    assert auroc([3.0, 2.0, 1.0, 0.0], [0, 0, 1, 1]) == 0.0


def test_auroc_all_scores_equal_is_half() -> None:
    assert auroc([0.4] * 6, [1, 0, 1, 1, 0, 0]) == 0.5


@pytest.mark.parametrize("y", [[1, 1, 1], [0, 0, 0]])
def test_auroc_one_class_is_nan(y: list[int]) -> None:
    assert math.isnan(auroc([0.1, 0.5, 0.9], y))


@pytest.mark.parametrize(
    ("scores", "y"),
    [
        ([0.1, 0.2], [1]),  # length mismatch
        ([], []),  # empty
        ([0.1, NAN], [1, 0]),  # non-finite score
        ([0.1, math.inf], [1, 0]),
        ([0.1, 0.2], [1, 2]),  # label outside {0, 1}
    ],
)
def test_auroc_and_curve_reject_bad_input(scores: list[float], y: list[int]) -> None:
    with pytest.raises(ValueError):
        auroc(scores, y)
    with pytest.raises(ValueError):
        threshold_curve(scores, y)


# --- the < 500 guard -------------------------------------------------------------------------


def test_require_labelled_returns_counts() -> None:
    y = [1] * 300 + [0] * 200
    assert require_labelled(y) == (300, 200)


@pytest.mark.parametrize("y", [[1] * 300 + [0] * 199, [1] * 600, [0] * 600])
def test_require_labelled_rejects(y: list[int]) -> None:
    with pytest.raises(ValueError):
        require_labelled(y)


def test_require_labelled_min_n_is_a_parameter() -> None:
    assert require_labelled([1, 0, 1], min_n=3) == (2, 1)


# --- threshold curve -------------------------------------------------------------------------

EX2 = [
    (0.99, 1),
    (0.95, 1),
    (0.90, 1),
    (0.90, 0),
    (0.80, 1),
    (0.70, 1),
    (0.60, 0),
    (0.50, 1),
    (0.40, 0),
    (0.10, 0),
]


def test_curve_columns_descending_distinct_with_ties_entering_together() -> None:
    s, y = split(EX2)
    c = threshold_curve(s, y)
    assert {"threshold", "tp", "fp"} <= set(c.columns)
    assert c["threshold"].tolist() == [0.99, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.10]
    assert c["tp"].tolist() == [1, 2, 3, 4, 5, 5, 6, 6, 6]
    assert c["fp"].tolist() == [0, 0, 1, 1, 1, 2, 2, 3, 4]
    assert pd.api.types.is_integer_dtype(c["tp"])
    assert pd.api.types.is_integer_dtype(c["fp"])


def test_curve_is_input_order_independent() -> None:
    s, y = split(EX2)
    a = threshold_curve(s, y)
    b = threshold_curve(s[::-1], y[::-1])
    pd.testing.assert_frame_equal(a.reset_index(drop=True), b.reset_index(drop=True))


# --- the cut-off rule ------------------------------------------------------------------------


def test_example2_precision_rule() -> None:
    c = cut_of(EX2)
    assert c.threshold == 0.95
    assert c.rule == "precision_0.90"
    assert (c.tp, c.fp) == (2, 0)
    assert c.precision == 1.0
    assert c.recall == pytest.approx(2 / 6)
    assert c.f1 == pytest.approx(4 / 8)
    assert c.literal_threshold == 0.95


def test_example3_max_f1_fallback() -> None:
    c = cut_of([(0.9, 0), (0.8, 1), (0.7, 1), (0.6, 0), (0.5, 1)])
    assert c.threshold == 0.5
    assert c.rule == "max_f1"
    assert (c.tp, c.fp) == (3, 2)
    assert c.precision == pytest.approx(0.6)
    assert c.recall == 1.0
    assert c.f1 == pytest.approx(0.75)
    assert c.literal_threshold is None


EX4 = [(t, 1) for t in (0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92)] + [
    (0.91, 0),
    (0.90, 1),
    (0.50, 0),
]


def test_example4_non_monotone_precision_stable_reading() -> None:
    """Lead decision 1: 0.92 (the literal reading gives 0.90, exactly 9/10, reported beside)."""
    curve = threshold_curve(*split(EX4))
    assert curve.loc[curve["threshold"] == 0.90, ["tp", "fp"]].values.tolist() == [[9, 1]]
    c = choose_cutoff(curve)
    assert c.threshold == 0.92
    assert c.rule == "precision_0.90"
    assert (c.tp, c.fp) == (8, 0)
    assert c.precision == 1.0
    assert c.recall == pytest.approx(8 / 9)
    assert c.literal_threshold == 0.90


def test_top_tie_group_below_target_falls_back_to_max_f1() -> None:
    """Methods property 3: 4 positives and 1 negative tied at the top (0.8 < 0.9) send the stable
    rule to max-F1 even though lower thresholds reach 0.90; the literal value is reported."""
    pairs = [(1.0, 1)] * 4 + [(1.0, 0)] + [(0.8, 1)] * 20 + [(0.3, 0)] * 10 + [(0.2, 1)]
    # cumulative precision: 0.8 at 1.0, 24/25 at 0.8, 24/35 at 0.3, 25/36 at 0.2; P = 25
    c = cut_of(pairs)
    assert c.rule == "max_f1"
    # F1 = 2TP/(TP+FP+P): 8/30, 48/50, 48/60, 50/61 -> max at 0.8
    assert c.threshold == 0.8
    assert (c.tp, c.fp) == (24, 1)
    assert c.literal_threshold == 0.8


def test_stable_fallback_differs_from_literal() -> None:
    """Top fails, 0.8 reaches exactly 9/10 (literal), but max-F1 picks 0.6."""
    pairs = [(0.9, 0)] + [(0.8, 1)] * 9 + [(0.7, 0)] + [(0.6, 1)] * 5
    c = cut_of(pairs)
    assert c.rule == "max_f1"
    assert c.threshold == 0.6
    assert (c.tp, c.fp) == (14, 2)
    assert c.literal_threshold == 0.8


def test_precision_exactly_090_counts_as_reached() -> None:
    pairs = [(0.7, 1)] * 9 + [(0.7, 0)] + [(0.3, 0)]
    c = cut_of(pairs)
    assert c.rule == "precision_0.90"
    assert c.threshold == 0.7
    assert (c.tp, c.fp) == (9, 1)
    assert c.precision == pytest.approx(0.9)


def test_equal_max_f1_takes_smallest_t() -> None:
    """F1 = 2/3 exactly at 0.7 (TP 2, FP 1) and at 0.5 (TP 3, FP 3); P = 3."""
    pairs = [(0.9, 0), (0.8, 1), (0.7, 1), (0.6, 0), (0.55, 0), (0.5, 1)]
    curve = threshold_curve(*split(pairs))
    p = 3
    f1 = [Fraction(2 * t, t + f + p) for t, f in zip(curve["tp"], curve["fp"], strict=True)]
    assert f1[2] == f1[5] == Fraction(2, 3) == max(f1)
    c = choose_cutoff(curve)
    assert c.rule == "max_f1"
    assert c.threshold == 0.5
    assert (c.tp, c.fp) == (3, 3)
    assert c.f1 == pytest.approx(2 / 3)


@pytest.mark.parametrize(
    ("y", "rule", "tp", "fp"),
    [([1, 1, 1, 1, 1, 1, 1, 1, 1, 0], "precision_0.90", 9, 1), ([1, 0, 1, 0], "max_f1", 2, 2)],
)
def test_all_scores_equal_one_candidate(y: list[int], rule: str, tp: int, fp: int) -> None:
    c = choose_cutoff(threshold_curve([0.42] * len(y), y))
    assert c.threshold == 0.42
    assert c.rule == rule
    assert (c.tp, c.fp) == (tp, fp)


def test_cutoff_is_the_exact_observed_float_never_rounded() -> None:
    t = 0.8812014060630281  # 2,304 TerminatorNet calls sit at exactly this value
    pairs = [(t, 1)] * 2304 + [(0.9999999999961204, 1)] * 5 + [(0.5, 0)] * 900
    c = cut_of(pairs)
    assert c.rule == "precision_0.90"
    assert type(c.threshold) is float
    assert c.threshold == t
    assert repr(c.threshold) == "0.8812014060630281"
    assert c.tp == 2309


def test_cutoff_keeps_near_one_values_distinct() -> None:
    hi, lo = 0.9999999999961204, 0.9999999999961203
    assert hi != lo
    c = cut_of([(hi, 1)] * 9 + [(lo, 0)] * 5)
    assert c.threshold == hi
    assert repr(c.threshold) == "0.9999999999961204"


def test_cutoff_tp_fp_are_python_ints() -> None:
    c = cut_of(EX2)
    assert type(c.tp) is int and type(c.fp) is int


def test_cutoff_is_a_frozen_dataclass() -> None:
    c = cut_of(EX2)
    assert dataclasses.is_dataclass(c)
    with pytest.raises(dataclasses.FrozenInstanceError):
        c.threshold = 0.1  # type: ignore[misc]


def test_choose_cutoff_other_target_uses_its_own_fraction() -> None:
    """Target 0.75 = 3/4 on example 2: precision 1, 1, 3/4 (exact), 4/5, 5/6, then 5/7 fails at
    0.60, so the stable t* = 0.70; the literal smallest is 0.50 (6/8 = 3/4 exactly)."""
    c = choose_cutoff(threshold_curve(*split(EX2)), precision=0.75)
    assert c.threshold == 0.70
    assert c.rule == "precision_0.75"
    assert (c.tp, c.fp) == (5, 1)
    assert c.literal_threshold == 0.50


# --- precision / recall at a fixed cut-off ---------------------------------------------------


def test_precision_recall_at() -> None:
    s, y = split(EX2)
    assert precision_at(s, y, 0.90) == pytest.approx(3 / 4)
    assert recall_at(s, y, 0.90) == pytest.approx(3 / 6)
    assert math.isnan(precision_at(s, y, 1.5))  # nothing called
    assert recall_at(s, y, 1.5) == 0.0
    assert math.isnan(recall_at([0.2, 0.3], [0, 0], 0.1))  # P = 0


# --- directon bootstrap wrappers -------------------------------------------------------------


def _mixed_clusters(n_dir: int = 40, seed: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Each directon holds at least one positive and one negative, so no replicate is one-class."""
    g = np.random.default_rng(seed)
    s: list[float] = []
    y: list[int] = []
    cl: list[int] = []
    for c in range(n_dir):
        k = int(g.integers(2, 6))
        lab = [1, 0] + [int(v) for v in g.integers(0, 2, k - 2)]
        u = g.normal()
        for v in lab:
            s.append(float(np.round(u + 1.2 * v + g.normal(), 2)))
            y.append(v)
            cl.append(c)
    return np.asarray(s), np.asarray(y), np.asarray(cl)


def _frame(x: np.ndarray, y: np.ndarray, cl: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame({"x": x, "y": y, "directon_id": cl})


def _same(a: BootCI, b: BootCI) -> bool:
    def eq(u: float, v: float) -> bool:
        return (math.isnan(u) and math.isnan(v)) or u == v

    return (
        eq(a.estimate, b.estimate)
        and eq(a.ci_low, b.ci_low)
        and eq(a.ci_high, b.ci_high)
        and a.n_valid == b.n_valid
        and a.B == b.B
    )


def test_auroc_ci_is_one_tu_block_bootstrap_call() -> None:
    s, y, cl = _mixed_clusters()
    got = auroc_ci(s, y, cl, rng=rng(), B=200)
    want = tu_block_bootstrap(
        _frame(s, y, cl), "directon_id", lambda d: auroc(d["x"], d["y"]), 200, rng(), 0.95
    )
    assert isinstance(got, BootCI)
    assert _same(got, want)
    assert got.estimate == auroc(s, y)
    assert got.n_valid == 200
    assert got.ci_low < got.estimate < got.ci_high


def test_auroc_ci_deterministic_and_default_b() -> None:
    s, y, cl = _mixed_clusters(n_dir=15)
    a = auroc_ci(s, y, cl, rng=rng())
    b = auroc_ci(s, y, cl, rng=rng())
    assert a.B == DEFAULT_B
    assert _same(a, b)


def test_auroc_ci_resamples_whole_directons() -> None:
    """Every directon is one positive and one negative tied at the same score. Any multiset of
    whole directons has AUROC exactly 0.5; resampling single pairs would not."""
    v = np.arange(30, dtype=np.float64)
    s = np.repeat(v, 2)
    y = np.tile([1, 0], 30)
    cl = np.repeat(np.arange(30), 2)
    ci = auroc_ci(s, y, cl, rng=rng(), B=300)
    assert (ci.estimate, ci.ci_low, ci.ci_high, ci.n_valid) == (0.5, 0.5, 0.5, 300)
    iid = auroc_ci(s, y, np.arange(60), rng=rng(), B=300)
    assert iid.ci_high - iid.ci_low > 0.05


def test_auroc_ci_drops_one_class_replicates() -> None:
    """Two directons, one all-positive and one all-negative: about half the replicates are
    one-class (NaN) and dropped."""
    s = np.array([0.9, 0.8, 0.2, 0.1])
    y = np.array([1, 1, 0, 0])
    cl = np.array([0, 0, 1, 1])
    ci = auroc_ci(s, y, cl, rng=rng(), B=400)
    assert ci.estimate == 1.0
    assert 120 < ci.n_valid < 280
    want = tu_block_bootstrap(
        _frame(s, y, cl), "directon_id", lambda d: auroc(d["x"], d["y"]), 400, rng(), 0.95
    )
    assert _same(ci, want)


def test_fixed_cutoff_ci_is_precision_then_recall_with_t_held_fixed() -> None:
    s, y, cl = _mixed_clusters()
    t = float(np.median(s))
    p, r = fixed_cutoff_ci(s, y, cl, t, rng=rng(), B=150)
    g = rng()
    df = _frame(s, y, cl)
    wp = tu_block_bootstrap(df, "directon_id", lambda d: precision_at(d["x"], d["y"], t), 150, g)
    wr = tu_block_bootstrap(df, "directon_id", lambda d: recall_at(d["x"], d["y"], t), 150, g)
    assert _same(p, wp) and _same(r, wr)
    assert p.estimate == precision_at(s, y, t)
    assert r.estimate == recall_at(s, y, t)


def test_fixed_cutoff_ci_with_window_clusters() -> None:
    """R2: the same function with window IDs as the cluster gives the window-cluster CI."""
    conf = np.array([0.95, 0.95, 0.90, 0.85, 0.85, 0.6, 0.6, 0.6])
    y = np.array([1, 1, 0, 1, 0, 0, 1, 0])
    win = np.array(["g1", "g1", "g2", "g3", "g3", "g4", "g5", "g5"])
    p, _ = fixed_cutoff_ci(conf, y, win, 0.85, rng=rng(), B=100)
    assert p.estimate == pytest.approx(3 / 5)
    assert p.n_valid >= 90  # NaN only when no called window is drawn


def test_auroc_diff_ci_paired_on_the_same_replicate() -> None:
    s, y, cl = _mixed_clusters()
    g = np.random.default_rng(11)
    other = s + g.normal(0.0, 0.5, s.size)
    got = auroc_diff_ci(s, other, y, cl, rng=rng(), B=150)
    df = pd.DataFrame({"a": s, "b": other, "y": y, "directon_id": cl})
    want = tu_block_bootstrap(
        df, "directon_id", lambda d: auroc(d["a"], d["y"]) - auroc(d["b"], d["y"]), 150, rng()
    )
    assert _same(got, want)
    same = auroc_diff_ci(s, s, y, cl, rng=rng(), B=150)
    assert (same.estimate, same.ci_low, same.ci_high) == (0.0, 0.0, 0.0)


def test_threshold_boot_q_quantiles_of_reselected_cutoff() -> None:
    s, y, cl = _mixed_clusters()
    got = threshold_boot_q(s, y, cl, rng=rng(), B=120)
    vals: list[float] = []

    def stat(d: pd.DataFrame) -> float:
        v = choose_cutoff(threshold_curve(d["x"], d["y"])).threshold
        vals.append(v)
        return v

    tu_block_bootstrap(_frame(s, y, cl), "directon_id", stat, 120, rng())
    q = np.quantile(np.asarray(vals[1:]), [0.025, 0.5, 0.975], method="linear")
    assert got == pytest.approx(tuple(float(v) for v in q), abs=0.0)
    assert got[0] <= got[1] <= got[2]


def test_distance_fit_ci_d0_then_s_separate_draws() -> None:
    g = np.random.default_rng(5)
    n_dir = 120
    cl = np.repeat(np.arange(n_dir), g.integers(1, 5, n_dir))
    d = g.integers(-20, 201, cl.size)
    y = (g.random(cl.size) < 1.0 / (1.0 + np.exp((d - 40) / 15))).astype(int)
    a, b = distance_fit_ci(d, y, cl, rng=rng(), B=60)

    def stat(attr: str):  # type: ignore[no-untyped-def]
        def f(x: pd.DataFrame) -> float:
            fit = fit_distance_model(x["x"].to_numpy(), x["y"].to_numpy())
            return getattr(fit, attr) if fit.converged else NAN

        return f

    gg = rng()
    df = _frame(d.astype(np.float64), y, cl)
    wa = tu_block_bootstrap(df, "directon_id", stat("d0"), 60, gg)
    wb = tu_block_bootstrap(df, "directon_id", stat("s"), 60, gg)
    assert a.estimate == pytest.approx(wa.estimate, rel=1e-9)
    assert (a.ci_low, a.ci_high) == pytest.approx((wa.ci_low, wa.ci_high), rel=1e-9)
    assert (b.ci_low, b.ci_high) == pytest.approx((wb.ci_low, wb.ci_high), rel=1e-9)
    assert (a.n_valid, b.n_valid) == (wa.n_valid, wb.n_valid)


def test_distance_fit_ci_failed_replicates_are_nan() -> None:
    """Two directons whose union overlaps but each alone is separated: one-directon replicates
    fail the separation check (NaN, dropped)."""
    d = np.array([0, 10, 20, 30, 5, 15, 25, 35])
    y = np.array([1, 1, 0, 0, 1, 0, 1, 0])
    cl = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    a, _ = distance_fit_ci(d, y, cl, rng=rng(), B=200)
    assert a.n_valid <= 200


# --- precision_cv (R1) -----------------------------------------------------------------------


def test_precision_cv_perfect_separation_is_one() -> None:
    n = 40
    s = np.r_[np.ones(n), np.zeros(n)]
    y = np.r_[np.ones(n, dtype=int), np.zeros(n, dtype=int)]
    cl = np.r_[np.arange(n), np.arange(n)]  # each directon: one positive, one negative
    assert precision_cv(s, y, cl, rng=rng()) == 1.0


def test_precision_cv_all_scores_equal_is_prevalence() -> None:
    y = np.array([1, 0, 1, 1, 0] * 8)
    cl = np.arange(y.size) // 2
    assert precision_cv(np.full(y.size, 0.3), y, cl, rng=rng()) == pytest.approx(24 / 40)


def test_precision_cv_deterministic_and_bounded() -> None:
    s, y, cl = _mixed_clusters()
    a = precision_cv(s, y, cl, rng=rng())
    assert a == precision_cv(s, y, cl, rng=rng())
    assert 0.0 <= a <= 1.0


def test_precision_cv_needs_k_clusters() -> None:
    with pytest.raises(ValueError):
        precision_cv([0.1, 0.9, 0.5], [0, 1, 1], [0, 1, 2], rng=rng(), k=10)


# --- Wilson ----------------------------------------------------------------------------------


@pytest.mark.parametrize(("k", "n"), [(0, 10), (1, 4), (7, 30), (727, 3217), (30, 30)])
def test_wilson_matches_statsmodels(k: int, n: int) -> None:
    lo, hi = wilson_ci(k, n)
    wlo, whi = proportion_confint(k, n, alpha=0.05, method="wilson")
    assert (lo, hi) == pytest.approx((wlo, whi), abs=1e-12)
    assert 0.0 <= lo <= k / n + 1e-15 and k / n <= hi + 1e-15 and hi <= 1.0


def test_wilson_exact_endpoints_at_zero_and_n() -> None:
    assert wilson_ci(0, 30)[0] == 0.0
    assert wilson_ci(30, 30)[1] == 1.0


def test_wilson_level_and_empty() -> None:
    lo90, hi90 = wilson_ci(7, 30, level=0.90)
    lo95, hi95 = wilson_ci(7, 30)
    assert lo95 < lo90 < hi90 < hi95
    assert all(math.isnan(v) for v in wilson_ci(0, 0))


# --- the d0 / s refit (B1) -------------------------------------------------------------------


def test_example5_refit_recovers_40_15() -> None:
    d, y = grid(range(-20, 201, 10), 1000, 40.0, 15.0)
    f = fit_distance_model(d, y)
    assert isinstance(f, DistanceFit)
    assert f.converged is True and f.calibrated is True and f.reason is None
    assert abs(f.d0 - 40) < 0.2 and abs(f.s - 15) < 0.1
    assert f.d0 == pytest.approx(39.991, abs=1e-3)
    assert f.s == pytest.approx(14.962, abs=1e-3)
    assert f.d0 == pytest.approx(-f.beta0 / f.beta1, rel=1e-12)
    assert f.s == pytest.approx(-1.0 / f.beta1, rel=1e-12)


def test_b1_hard_start_refit_recovers_100_50() -> None:
    """B1: plain Newton from (2.5, -0.05) fails here; the damped fit from (logit ybar, 0) must
    give 100.14 +- 0.05 and 50.06 +- 0.05 (scratch 100.136, 50.055)."""
    d, y = grid(range(-20, 401, 20), 100, 100.0, 50.0)
    assert d.size == 2200
    f = fit_distance_model(d, y)
    assert f.converged is True
    assert f.d0 == pytest.approx(100.14, abs=0.05)
    assert f.s == pytest.approx(50.06, abs=0.05)


def test_refit_matches_statsmodels_mle() -> None:
    import statsmodels.api as sm

    g = np.random.default_rng(21)
    d = g.integers(-60, 560, 2700)
    y = (g.random(d.size) < 1.0 / (1.0 + np.exp((d - 102.7) / 47.2))).astype(int)
    f = fit_distance_model(d, y)
    r = sm.Logit(y, sm.add_constant(d.astype(np.float64))).fit(disp=0, tol=1e-12)
    b0, b1 = (float(v) for v in r.params)
    assert f.converged is True
    assert f.d0 == pytest.approx(-b0 / b1, rel=1e-7)
    assert f.s == pytest.approx(-1.0 / b1, rel=1e-7)


def _assert_failed(f: DistanceFit) -> None:
    assert f.converged is False
    assert f.calibrated is False
    assert (f.d0, f.s) == (50.0, 20.0)
    assert isinstance(f.reason, str) and f.reason


def test_example5_separated_fails_at_the_separation_check() -> None:
    d, _ = grid(range(-20, 201, 10), 1000, 40.0, 15.0)
    y = (d < 60).astype(int)
    f = fit_distance_model(d, y)
    _assert_failed(f)
    assert f.reason == "separation"
    assert f.n_iter == 0  # reported before any iteration


def test_quasi_separated_fails_without_exception() -> None:
    """Three positives added at d = 60: the boundary holds both classes, but no pair inverts."""
    d, _ = grid(range(-20, 201, 10), 1000, 40.0, 15.0)
    y = (d < 60).astype(int)
    f = fit_distance_model(np.r_[d, [60, 60, 60]], np.r_[y, [1, 1, 1]])
    _assert_failed(f)
    assert f.reason == "separation"


def test_reverse_separation_fails() -> None:
    d = np.arange(0, 200, 5)
    _assert_failed(fit_distance_model(d, (d > 100).astype(int)))


def test_beta1_nonnegative_fails() -> None:
    g = np.random.default_rng(1)
    d = g.integers(-20, 200, 500)
    y = (g.random(500) < 1.0 / (1.0 + np.exp(-(d - 100) / 30))).astype(int)  # rising with d
    f = fit_distance_model(d, y)
    _assert_failed(f)
    assert f.beta1 >= 0 or not math.isfinite(f.beta1)


def test_steep_backstop_s_below_one_bp_fails() -> None:
    """An MLE exists (overlap at d = 38..42) but s ~ 0.5 bp: |beta1| > 1 backstop."""
    d, y = grid(range(30, 51), 1000, 40.0, 0.5)
    _assert_failed(fit_distance_model(d, y))


def test_non_convergence_fails() -> None:
    d, y = grid(range(-20, 201, 10), 1000, 40.0, 15.0)
    f = fit_distance_model(d, y, max_iter=1)
    _assert_failed(f)
    assert f.n_iter == 1


@pytest.mark.parametrize(
    ("d", "y"),
    [
        ([0.0, 1e300, 1.0, 2e300], [1, 1, 0, 0]),  # overlap both ways; weights underflow
        ([5.0, 5.0, 5.0, 5.0], [1, 0, 1, 0]),  # one distance: singular design
        ([-1e308, 1e308, 0.0, 1.0], [1, 0, 0, 1]),
    ],
)
def test_hostile_inputs_never_raise(d: list[float], y: list[int]) -> None:
    f = fit_distance_model(d, y)
    if f.converged:
        assert math.isfinite(f.d0) and f.s >= 1.0
    else:
        _assert_failed(f)


@pytest.mark.parametrize(
    ("d", "y"),
    [([1.0, NAN], [1, 0]), ([1.0, 2.0], [1, 3]), ([1.0], [1, 0]), ([], [])],
)
def test_refit_rejects_invalid_input(d: list[float], y: list[int]) -> None:
    with pytest.raises(ValueError):
        fit_distance_model(d, y)


def test_distance_fit_is_frozen() -> None:
    d, y = grid(range(-20, 201, 10), 100, 40.0, 15.0)
    f = fit_distance_model(d, y)
    with pytest.raises(dataclasses.FrozenInstanceError):
        f.d0 = 1.0  # type: ignore[misc]


# --- log-loss --------------------------------------------------------------------------------


def test_log_loss_matches_sklearn_and_refit_beats_placeholder() -> None:
    d, y = grid(range(-20, 401, 20), 100, 100.0, 50.0)
    f = fit_distance_model(d, y)
    p = 1.0 / (1.0 + np.exp((d - 50.0) / 20.0))
    assert distance_log_loss(d, y, 50.0, 20.0) == pytest.approx(sk_log_loss(y, p), rel=1e-10)
    assert distance_log_loss(d, y, f.d0, f.s) < distance_log_loss(d, y, 50.0, 20.0)


def test_log_loss_finite_far_in_the_tails() -> None:
    v = distance_log_loss([-5000, 5000], [1, 0], 50.0, 20.0)
    assert math.isfinite(v) and 0.0 <= v < 1e-50
