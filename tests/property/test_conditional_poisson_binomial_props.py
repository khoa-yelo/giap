"""Hypothesis properties of the conditional Poisson-binomial primitive (issue M8-04a;
docs/specs/M8-04 §1 and its acceptance criteria; ADR-0029).

- T <= 14: ``p_over``, ``p_under``, ``mean`` and the support equal exhaustive enumeration over the
  2^T configurations (exact rationals) to 1e-12 relative;
- T <= 60: the tails equal a ``Fraction`` DP to 1e-12 relative for tails >= 1e-300, and are
  finite floats >= 0 below that;
- p_over + p_under - P(K = k) = 1 +- 1e-12 over the support; p_over non-increasing and p_under
  non-decreasing in k; p_over(lo) = p_under(hi) = 1;
- odds x c (any c > 0) and ``tilt`` to any feasible total leave the law unchanged (1e-12);
- swapping the set and its complement: p_over(C, k) = p_under(complement, N - k) and the means
  add to N; permuting within either side changes nothing;
- monotone likelihood ratio: multiplying the odds of C's units by c >= 1 never lowers p_over or
  the mean, never raises p_under;
- ``poisson_binomial_pmf_trunc`` is non-negative, sums to cdf(kmax) (to 1 when kmax >= n), is the
  prefix of the merged ``poisson_binomial_pmf``;
- ``tilt`` meets the total, preserves the order of q, shifts every logit by one theta, keeps
  q in {0, 1}.
"""

from __future__ import annotations

import itertools
import math
import random
from fractions import Fraction

import numpy as np
from helpers.cpb_oracle import assert_tail_close, cpb_enumerate, cpb_exact, rel_err
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from gippa.stats.poisson_binomial import (
    conditional_poisson_binomial,
    poisson_binomial_cdf,
    poisson_binomial_pmf,
    poisson_binomial_pmf_trunc,
    tilt,
)

RTOL = 1e-12
FAST = settings(max_examples=150, deadline=None)
SLOW = settings(max_examples=40, deadline=None)

# Open probabilities kept away from float underflow of the untilted DP at T <= 60 (the tilt
# handles badly scaled q; see the unit tests), plus exact 0 and 1.
OPEN = st.floats(1e-3, 1 - 1e-3, allow_nan=False)
Q = st.one_of(OPEN, OPEN, OPEN, st.sampled_from([0.0, 1.0]))
# For tests that rescale odds in floating point: rounding of the rescaled q (relative error of
# 1 - q near 1) must stay below the 1e-11 comparison, so q stays in [0.01, 0.99] and c in
# [1/20, 20].
MID = st.floats(0.01, 0.99, allow_nan=False)
QM = st.one_of(MID, MID, MID, st.sampled_from([0.0, 1.0]))


@st.composite
def problems(draw: st.DrawFn, max_t: int, q: st.SearchStrategy[float] = Q) -> tuple:
    """(p_in, p_out, k, N) with S = N of positive probability and k in the support."""
    p_in = draw(st.lists(q, min_size=0, max_size=max_t))
    p_out = draw(st.lists(q, min_size=0, max_size=max_t - len(p_in)))
    arr_in, arr_out = np.asarray(p_in), np.asarray(p_out)
    ones = int((arr_in == 1).sum() + (arr_out == 1).sum())
    pos = int((arr_in > 0).sum() + (arr_out > 0).sum())
    N = draw(st.integers(ones, pos))
    lo = max(int((arr_in == 1).sum()), N - int((arr_out > 0).sum()))
    hi = min(int((arr_in > 0).sum()), N - int((arr_out == 1).sum()))
    assume(lo <= hi)
    k = draw(st.integers(lo, hi))
    return p_in, p_out, k, N


def _odds_scale(q: list[float], c: float) -> list[float]:
    out = []
    for x in q:
        if x in (0.0, 1.0):
            out.append(x)
        else:
            o = x / (1 - x) * c
            out.append(o / (1 + o))
    return out


@FAST
@given(problems(14))
def test_equals_exhaustive_enumeration(prob: tuple) -> None:
    p_in, p_out, k, N = prob
    ex = cpb_enumerate(p_in, p_out, k, N)
    assert ex is not None
    t = conditional_poisson_binomial(p_in, p_out, k, N)
    assert rel_err(t.p_over, ex.p_over) <= RTOL
    assert rel_err(t.p_under, ex.p_under) <= RTOL
    assert rel_err(t.mean, ex.mean) <= RTOL or abs(t.mean - float(ex.mean)) <= 1e-15
    assert (t.lo, t.hi) == (ex.lo, ex.hi)
    assert (t.k, t.N) == (k, N)


@SLOW
@given(problems(60))
def test_equals_fraction_dp_to_the_precision_contract(prob: tuple) -> None:
    p_in, p_out, k, N = prob
    ex = cpb_exact(p_in, p_out, k, N)
    assert ex is not None
    t = conditional_poisson_binomial(p_in, p_out, k, N)
    assert_tail_close(t.p_over, ex.p_over)
    assert_tail_close(t.p_under, ex.p_under)
    assert (t.lo, t.hi) == (ex.lo, ex.hi)


@FAST
@given(problems(14))
def test_tails_overlap_in_exactly_p_k(prob: tuple) -> None:
    p_in, p_out, k, N = prob
    ex = cpb_enumerate(p_in, p_out, k, N)
    assert ex is not None
    t = conditional_poisson_binomial(p_in, p_out, k, N)
    assert abs(t.p_over + t.p_under - float(ex.p_eq) - 1.0) <= RTOL
    assert 0.0 <= t.p_over <= 1.0 + RTOL and 0.0 <= t.p_under <= 1.0 + RTOL
    assert t.lo - 1e-12 <= t.mean <= t.hi + 1e-12


@FAST
@given(problems(40))
def test_tails_are_monotone_in_k(prob: tuple) -> None:
    p_in, p_out, _, N = prob
    t0 = conditional_poisson_binomial(p_in, p_out, _, N)
    over = [conditional_poisson_binomial(p_in, p_out, j, N).p_over for j in range(t0.lo, t0.hi + 1)]
    under = [
        conditional_poisson_binomial(p_in, p_out, j, N).p_under for j in range(t0.lo, t0.hi + 1)
    ]
    assert over[0] == 1.0 or abs(over[0] - 1.0) <= RTOL
    assert abs(under[-1] - 1.0) <= RTOL
    for a, b in itertools.pairwise(over):
        assert b <= a * (1 + 1e-12)
    for a, b in itertools.pairwise(under):
        assert b >= a * (1 - 1e-12)


@FAST
@given(problems(30, QM), st.floats(0.05, 20.0))
def test_odds_scaling_invariance(prob: tuple, c: float) -> None:
    p_in, p_out, k, N = prob
    a = conditional_poisson_binomial(p_in, p_out, k, N)
    b = conditional_poisson_binomial(_odds_scale(p_in, c), _odds_scale(p_out, c), k, N)
    # odds x c moves every q by O(1 ulp) of rounding, i.e. ~T * 1e-16 relative in the law
    assert rel_err(b.p_over, Fraction(a.p_over)) <= 1e-11 or b.p_over == a.p_over
    assert rel_err(b.p_under, Fraction(a.p_under)) <= 1e-11 or b.p_under == a.p_under
    assert abs(b.mean - a.mean) <= 1e-11 * max(1.0, a.mean)
    assert (a.lo, a.hi) == (b.lo, b.hi)


@SLOW
@given(problems(40), st.data())
def test_tilt_to_any_feasible_total_is_invariant(prob: tuple, data: st.DataObject) -> None:
    p_in, p_out, k, N = prob
    q = np.asarray(p_in + p_out, dtype=float)
    ones, pos = int((q == 1).sum()), int((q > 0).sum())
    assume(pos - ones >= 2)  # interior totals: a finite theta exists
    total = data.draw(st.integers(ones + 1, pos - 1))
    qt = tilt(q, total)
    n = len(p_in)
    ex = cpb_exact(p_in, p_out, k, N)  # the exact law of the raw q
    assert ex is not None
    t = conditional_poisson_binomial(qt[:n], qt[n:], k, N)
    assert (t.lo, t.hi) == (ex.lo, ex.hi)
    assert_tail_close(t.p_over, ex.p_over, rtol=1e-11)
    assert_tail_close(t.p_under, ex.p_under, rtol=1e-11)


@FAST
@given(problems(30))
def test_swapping_set_and_complement(prob: tuple) -> None:
    p_in, p_out, k, N = prob
    a = conditional_poisson_binomial(p_in, p_out, k, N)
    b = conditional_poisson_binomial(p_out, p_in, N - k, N)
    assert b.p_under == a.p_over or rel_err(b.p_under, Fraction(a.p_over)) <= 1e-12
    assert b.p_over == a.p_under or rel_err(b.p_over, Fraction(a.p_under)) <= 1e-12
    assert abs(a.mean + b.mean - N) <= 1e-12 * max(1, N)
    assert (b.lo, b.hi) == (N - a.hi, N - a.lo)


@FAST
@given(problems(30), st.randoms(use_true_random=False))
def test_permutation_within_sides_is_invariant(prob: tuple, rnd: random.Random) -> None:
    p_in, p_out, k, N = prob
    a = conditional_poisson_binomial(p_in, p_out, k, N)
    pi, po = list(p_in), list(p_out)
    rnd.shuffle(pi)
    rnd.shuffle(po)
    b = conditional_poisson_binomial(pi, po, k, N)
    assert b.p_over == a.p_over or rel_err(b.p_over, Fraction(a.p_over)) <= 1e-12
    assert b.p_under == a.p_under or rel_err(b.p_under, Fraction(a.p_under)) <= 1e-12


@FAST
@given(problems(30, MID), st.floats(1.0, 20.0))
def test_raising_the_odds_of_c_shifts_k_up(prob: tuple, c: float) -> None:
    p_in, p_out, k, N = prob
    a = conditional_poisson_binomial(p_in, p_out, k, N)
    b = conditional_poisson_binomial(_odds_scale(p_in, c), p_out, k, N)
    assert b.p_over >= a.p_over * (1 - 1e-11)
    assert b.p_under <= a.p_under * (1 + 1e-11)
    assert b.mean >= a.mean - 1e-11 * max(1.0, a.mean)


# ---------------------------------------------------------------- truncated pmf


@FAST
@given(st.lists(Q, min_size=0, max_size=60), st.integers(0, 70))
def test_pmf_trunc_is_the_prefix_and_sums_to_the_cdf(p: list[float], kmax: int) -> None:
    got = poisson_binomial_pmf_trunc(p, kmax)
    assert got.shape == (kmax + 1,)
    assert np.all(got >= 0.0) and np.all(np.isfinite(got))
    full = poisson_binomial_pmf(p)
    m = min(kmax, len(p)) + 1
    np.testing.assert_allclose(got[:m], full[:m], rtol=1e-15, atol=0)
    assert np.all(got[m:] == 0.0)
    assert abs(got.sum() - poisson_binomial_cdf(p, kmax)) <= 1e-12
    if kmax >= len(p):
        assert abs(got.sum() - 1.0) <= 1e-12


# ---------------------------------------------------------------- tilt


@FAST
@given(st.lists(Q, min_size=1, max_size=200), st.data())
def test_tilt_meets_the_total_and_shifts_logits_alike(p: list[float], data: st.DataObject) -> None:
    q = np.asarray(p, dtype=float)
    ones, pos = int((q == 1).sum()), int((q > 0).sum())
    assume(pos - ones >= 2)
    total = data.draw(st.integers(ones + 1, pos - 1))
    qt = tilt(q, total)
    assert qt.shape == q.shape and qt.dtype == np.float64
    assert abs(qt.sum() - total) <= 1e-9 * total
    np.testing.assert_array_equal(qt[q == 0], 0.0)
    np.testing.assert_array_equal(qt[q == 1], 1.0)
    m = (q > 0) & (q < 1)
    assert np.all((qt[m] > 0) & (qt[m] < 1))
    shift = (np.log(qt[m]) - np.log1p(-qt[m])) - (np.log(q[m]) - np.log1p(-q[m]))
    assert np.ptp(shift) <= 1e-9 * max(1.0, float(np.abs(shift).max()))
    # order preserved
    order = np.argsort(q[m], kind="stable")
    assert np.all(np.diff(qt[m][order]) >= 0)


@FAST
@given(st.lists(OPEN, min_size=2, max_size=200), st.floats(-5.0, 5.0), st.data())
def test_tilt_depends_on_q_only_through_the_odds_ratios(
    p: list[float], theta: float, data: st.DataObject
) -> None:
    q = np.asarray(p)
    lg = np.log(q) - np.log1p(-q) + theta
    q2 = 1 / (1 + np.exp(-lg))
    assume(np.all((q2 > 0) & (q2 < 1)))
    total = data.draw(st.integers(1, len(p) - 1))
    a, b = tilt(q, total), tilt(q2, total)
    # the spec tolerance on the total is 1e-9 relative, so the two solutions agree to ~1e-8
    np.testing.assert_allclose(a, b, rtol=1e-7, atol=0)
    assert math.isfinite(float(a.sum()))
