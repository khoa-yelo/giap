"""Hypothesis properties of ``giap.stats.cutoffs`` (issue M6-12 part A; spec acceptance criteria
1-2; methods "K-12 cut-off calibration (M6-12)" > "Properties the tests rely on").

Strategies build both classes by construction (separate positive and negative score lists, each
of length >= 1), so no ``assume`` filtering is needed (filter_too_much health check).
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st
from sklearn.metrics import roc_auc_score

from giap.rng import for_step
from giap.stats.cutoffs import (
    STEP,
    auroc,
    auroc_ci,
    choose_cutoff,
    threshold_curve,
)

# Scores on a coarse grid so ties are frequent (UniOP: 153 distinct values over 2,739 pairs).
GRID = st.integers(min_value=0, max_value=12)


@st.composite
def labelled(draw: st.DrawFn, max_side: int = 30) -> tuple[list[float], list[int]]:
    pos = draw(st.lists(GRID, min_size=1, max_size=max_side))
    neg = draw(st.lists(GRID, min_size=1, max_size=max_side))
    items = [(v / 12, 1) for v in pos] + [(v / 12, 0) for v in neg]
    items = draw(st.permutations(items))
    return [s for s, _ in items], [y for _, y in items]


def reference(
    scores: list[float], y: list[int], a: int = 9, b: int = 10
) -> tuple[float, str, float | None]:
    """The methods' rule in Fractions: (t*, rule, literal t)."""
    cands = sorted(set(scores), reverse=True)
    p = sum(y)
    tp = [sum(1 for s, v in zip(scores, y, strict=True) if v == 1 and s >= t) for t in cands]
    fp = [sum(1 for s, v in zip(scores, y, strict=True) if v == 0 and s >= t) for t in cands]
    ok = [Fraction(t, t + f) >= Fraction(a, b) for t, f in zip(tp, fp, strict=True)]
    literal = next((cands[k] for k in range(len(cands) - 1, -1, -1) if ok[k]), None)
    k_star = -1
    for k, good in enumerate(ok):
        if not good:
            break
        k_star = k
    if k_star >= 0:
        return cands[k_star], "precision_0.90", literal
    f1 = [Fraction(2 * t, t + f + p) for t, f in zip(tp, fp, strict=True)]
    best = max(f1)
    k_best = max(k for k, v in enumerate(f1) if v == best)
    return cands[k_best], "max_f1", literal


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_auroc_equals_sklearn_on_tied_inputs(data: tuple[list[float], list[int]]) -> None:
    s, y = data
    assert abs(auroc(s, y) - roc_auc_score(y, s)) <= 1e-12


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_auroc_of_negated_score_is_complement(data: tuple[list[float], list[int]]) -> None:
    s, y = data
    assert abs(auroc(s, y) - (1.0 - auroc([-v for v in s], y))) <= 1e-12


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_auroc_of_flipped_labels_is_complement(data: tuple[list[float], list[int]]) -> None:
    s, y = data
    assert abs(auroc(s, y) - (1.0 - auroc(s, [1 - v for v in y]))) <= 1e-12


@settings(max_examples=200, deadline=None)
@given(labelled())
def test_choose_cutoff_agrees_with_fraction_reference(data: tuple[list[float], list[int]]) -> None:
    s, y = data
    c = choose_cutoff(threshold_curve(s, y))
    t, rule, literal = reference(s, y)
    assert (c.threshold, c.rule, c.literal_threshold) == (t, rule, literal)
    assert c.tp == sum(1 for u, v in zip(s, y, strict=True) if v == 1 and u >= t)
    assert c.fp == sum(1 for u, v in zip(s, y, strict=True) if v == 0 and u >= t)


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_stable_cutoff_properties(data: tuple[list[float], list[int]]) -> None:
    """Every candidate at or above t* passes (precision rule); stable t* >= literal t*."""
    s, y = data
    curve = threshold_curve(s, y)
    c = choose_cutoff(curve)
    if c.rule == "precision_0.90":
        above = curve[curve["threshold"] >= c.threshold]
        assert all(10 * t >= 9 * (t + f) for t, f in zip(above["tp"], above["fp"], strict=True))
        assert c.literal_threshold is not None and c.threshold >= c.literal_threshold
    assert c.threshold in set(s)


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_curve_counts_are_cumulative(data: tuple[list[float], list[int]]) -> None:
    s, y = data
    curve = threshold_curve(s, y)
    tp, fp = curve["tp"].tolist(), curve["fp"].tolist()
    assert tp == sorted(tp) and fp == sorted(fp)
    assert (tp[-1], fp[-1]) == (sum(y), len(y) - sum(y))
    assert curve["threshold"].tolist() == sorted(set(s), reverse=True)


@settings(max_examples=100, deadline=None)
@given(labelled())
def test_cutoff_invariant_to_monotone_transform(data: tuple[list[float], list[int]]) -> None:
    """Rank terms: a strictly increasing map moves t* to f(t*) with the same counts and rule."""
    s, y = data

    def f(v: float) -> float:
        return math.exp(3.0 * v) - 7.0

    a = choose_cutoff(threshold_curve(s, y))
    b = choose_cutoff(threshold_curve([f(v) for v in s], y))
    assert b.threshold == f(a.threshold)
    assert (b.rule, b.tp, b.fp) == (a.rule, a.tp, a.fp)
    assert (b.literal_threshold is None) == (a.literal_threshold is None)
    if a.literal_threshold is not None:
        assert b.literal_threshold == f(a.literal_threshold)
    assert auroc([f(v) for v in s], y) == auroc(s, y)


@settings(max_examples=20, deadline=None)
@given(labelled(max_side=15), st.integers(min_value=0, max_value=2**32 - 1))
def test_auroc_ci_deterministic_for_a_seed(data: tuple[list[float], list[int]], seed: int) -> None:
    s, y = data
    cl = np.arange(len(s)) // 2
    a = auroc_ci(s, y, cl, rng=for_step(STEP, seed), B=40)
    b = auroc_ci(s, y, cl, rng=for_step(STEP, seed), B=40)
    for u, v in ((a.estimate, b.estimate), (a.ci_low, b.ci_low), (a.ci_high, b.ci_high)):
        assert u == v or (math.isnan(u) and math.isnan(v))
    assert a.n_valid == b.n_valid
