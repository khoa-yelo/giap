"""Conditional Poisson-binomial primitive (issue M8-04a; docs/specs/M8-04 §1, ADR-0029).

Interface pinned here (spec "Public interface", lead decision 7):

- ``poisson_binomial_pmf_trunc(p, kmax) -> np.ndarray``: float64, length ``kmax + 1``,
  P(X = j) for j = 0..kmax (zeros above len(p)); the merged ``_pmf`` recurrence stopped at kmax;
- ``tilt(p, total) -> np.ndarray``: odds x e^theta with sum = total (|sum - total| <= 1e-9 total),
  q in {0, 1} kept fixed;
- ``conditional_poisson_binomial(p_in, p_out, k, N) -> CpbTails``: the law of K = sum_C D given
  S = N; ``p_over`` = P(K >= k), ``p_under`` = P(K <= k) (both include k), ``mean`` = E[K | N],
  ``lo``/``hi`` = the support narrowed to the j of non-zero weight; a k outside it, a zero or
  non-finite support sum, or an invalid probability is a ``ValueError``.

Spec ambiguities resolved here (see the report of the test-author):
- ``poisson_binomial_pmf_trunc`` returns kmax + 1 entries even when kmax > len(p) (zeros above n,
  "j = 0..kmax" in the interface comment); kmax < 0 is a ValueError;
- an empty ``p_out`` (C = U) is a valid primitive call with support {N}; the ORA layer makes the
  row NA (spec degenerate table);
- ``tilt`` with a total outside [#(p == 1), #(p > 0)] is a ValueError (no theta exists); at the
  two ends it returns the limit (open entries -> 0 or -> 1);
- N < 0 or N > len(p_in) + len(p_out) is a ValueError;
- "narrowed to the j with a non-zero weight" means the exact (structural) zeros caused by q in
  {0, 1}: [lo, hi] = [max(#(p_in == 1), N - #(p_out > 0)), min(#(p_in > 0), N - #(p_out == 1))]
  intersected with [max(0, N - |out|), min(|C|, N)]. A k whose weight is positive but underflows
  is inside the support and returns a tail >= 0 (precision contract: "a tail of 0.0 is reachable
  only below ~1e-300"), not a ValueError.
"""

from __future__ import annotations

import dataclasses
import inspect
import math
from fractions import Fraction

import numpy as np
import pytest
from helpers.cpb_oracle import (
    assert_tail_close,
    cpb_enumerate,
    cpb_exact,
    hypergeom_exact,
    pb_pmf_exact,
    rel_err,
)

from giap.stats.poisson_binomial import (
    CpbTails,
    conditional_poisson_binomial,
    poisson_binomial_pmf,
    poisson_binomial_pmf_trunc,
    tilt,
)

RTOL = 1e-12

# Toy universe T6 (spec examples table): q_t = 1 - prod(1 - p_detect) over the TU's genes.
T6_GENES = {
    "TU_1": [0.9, 0.5],
    "TU_2": [0.5],
    "TU_3": [0.2, 0.2],
    "TU_4": [0.5],
    "TU_5": [0.3, 0.3, 0.3],
    "TU_6": [0.1],
}


def _q(ps: list[float]) -> float:
    return 1.0 - math.prod(1.0 - p for p in ps)


def _odds_scale(q: list[float], c: float) -> list[float]:
    o = np.asarray(q) / (1.0 - np.asarray(q)) * c
    return list(o / (1.0 + o))


# ---------------------------------------------------------------- interface


def test_signatures_are_pinned() -> None:
    assert list(inspect.signature(conditional_poisson_binomial).parameters) == [
        "p_in",
        "p_out",
        "k",
        "N",
    ]
    assert list(inspect.signature(poisson_binomial_pmf_trunc).parameters) == ["p", "kmax"]
    assert list(inspect.signature(tilt).parameters) == ["p", "total"]


def test_cpbtails_is_a_frozen_dataclass_with_the_spec_fields() -> None:
    assert dataclasses.is_dataclass(CpbTails)
    names = [f.name for f in dataclasses.fields(CpbTails)]
    assert names == ["k", "N", "p_over", "p_under", "mean", "lo", "hi"]
    t = conditional_poisson_binomial([0.5], [0.5], 1, 1)
    with pytest.raises(dataclasses.FrozenInstanceError):
        t.p_over = 0.0  # type: ignore[misc]


def test_result_types() -> None:
    t = conditional_poisson_binomial([0.3, 0.6], [0.2, 0.4, 0.5], 1, 2)
    assert isinstance(t, CpbTails)
    for v in (t.p_over, t.p_under, t.mean):
        assert isinstance(v, float)
    for v in (t.k, t.N, t.lo, t.hi):
        assert isinstance(v, (int, np.integer))
    assert (t.k, t.N) == (1, 2)


# ---------------------------------------------------------------- spec examples


def test_example1_t6_tu_level() -> None:
    q = {t: _q(g) for t, g in T6_GENES.items()}
    in_c = ["TU_1", "TU_3", "TU_5"]
    p_in = [q[t] for t in in_c]
    p_out = [q[t] for t in q if t not in in_c]
    t = conditional_poisson_binomial(p_in, p_out, 3, 4)
    assert rel_err(t.p_over, Fraction(2134593, 5159698)) <= RTOL
    assert t.p_over == pytest.approx(0.41370503, abs=5e-9)
    assert rel_err(t.p_under, Fraction(1)) <= RTOL
    assert rel_err(t.mean, Fraction(6168059, 2579849)) <= RTOL
    assert t.mean == pytest.approx(2.39086047, abs=5e-9)
    assert 3 - t.mean == pytest.approx(0.60913953, abs=5e-9)  # the spec's effect
    assert (t.lo, t.hi) == (1, 3)  # [max(0, 4 - 3), min(3, 4)]


def test_example1_matches_the_enumeration_oracle() -> None:
    q = {t: _q(g) for t, g in T6_GENES.items()}
    p_in = [q["TU_1"], q["TU_3"], q["TU_5"]]
    p_out = [q["TU_2"], q["TU_4"], q["TU_6"]]
    ex = cpb_enumerate(p_in, p_out, 3, 4)
    assert ex is not None
    t = conditional_poisson_binomial(p_in, p_out, 3, 4)
    assert rel_err(t.p_over, ex.p_over) <= RTOL
    assert rel_err(t.p_under, ex.p_under) <= RTOL
    assert rel_err(t.mean, ex.mean) <= RTOL


def test_example1_gene_level_recorded_values() -> None:
    # gene CPB: C = {g1, g4, g7}, k = 1 (g4), N = 4 detected genes of 10
    p = {f"g{i + 1}": v for i, v in enumerate(x for g in T6_GENES.values() for x in g)}
    c = ["g1", "g4", "g7"]
    t = conditional_poisson_binomial([p[g] for g in c], [p[g] for g in p if g not in c], 1, 4)
    assert t.p_over == pytest.approx(0.97289, abs=5e-6)
    assert t.p_under == pytest.approx(0.56849, abs=5e-6)


def test_example2_tu_packing() -> None:
    t = conditional_poisson_binomial([1 - 0.5**5], [0.5] * 20, 1, 11)
    assert rel_err(t.p_over, Fraction(341, 351)) <= RTOL
    assert t.p_over == pytest.approx(0.97150997, abs=5e-9)
    assert rel_err(t.p_under, Fraction(1)) <= RTOL
    assert rel_err(t.mean, Fraction(341, 351)) <= RTOL
    assert (t.lo, t.hi) == (0, 1)


def test_example2_gene_level_is_the_hypergeometric() -> None:
    t = conditional_poisson_binomial([0.5] * 5, [0.5] * 20, 5, 15)
    assert rel_err(t.p_over, Fraction(13, 230)) <= RTOL
    assert t.p_over == pytest.approx(0.05652174, abs=5e-9)


def test_example4_equal_q_is_the_central_hypergeometric() -> None:
    t = conditional_poisson_binomial([0.4] * 5, [0.4] * 15, 4, 8)
    assert rel_err(t.p_over, Fraction(56, 969)) <= RTOL
    assert rel_err(t.p_under, Fraction(1931, 1938)) <= RTOL
    assert t.p_over == pytest.approx(0.05779154, abs=5e-9)
    assert t.p_under == pytest.approx(0.99638803, abs=5e-9)
    assert rel_err(t.mean, Fraction(2)) <= RTOL  # 8 * 5 / 20
    over = sum((hypergeom_exact(20, 5, 8, j) for j in range(4, 6)), Fraction(0))
    under = sum((hypergeom_exact(20, 5, 8, j) for j in range(5)), Fraction(0))
    assert rel_err(t.p_over, over) <= RTOL
    assert rel_err(t.p_under, under) <= RTOL


EX5_IN = [0.9, 0.5, 0.36]
EX5_OUT = [0.5, 0.657, 0.1, 0.3, 0.8]


def test_example5_values() -> None:
    t = conditional_poisson_binomial(EX5_IN, EX5_OUT, 2, 4)
    assert rel_err(t.p_over, Fraction(10081485, 15841006)) <= RTOL
    assert rel_err(t.p_under, Fraction(72680561, 79205030)) <= RTOL
    assert rel_err(t.mean, Fraction(13531034, 7920503)) <= RTOL
    assert t.p_over == pytest.approx(0.63641697, abs=5e-9)
    assert t.p_under == pytest.approx(0.91762557, abs=5e-9)


def test_example5_odds_times_three_is_invariant() -> None:
    a = conditional_poisson_binomial(EX5_IN, EX5_OUT, 2, 4)
    b = conditional_poisson_binomial(_odds_scale(EX5_IN, 3.0), _odds_scale(EX5_OUT, 3.0), 2, 4)
    assert b.p_over == pytest.approx(a.p_over, rel=RTOL)
    assert b.p_under == pytest.approx(a.p_under, rel=RTOL)
    assert b.mean == pytest.approx(a.mean, rel=RTOL)


def test_example5_tilted_inputs_are_invariant() -> None:
    q = np.array(EX5_IN + EX5_OUT)
    for total in (1, 3, 4, 6):
        qt = tilt(q, total)
        a = conditional_poisson_binomial(EX5_IN, EX5_OUT, 2, 4)
        b = conditional_poisson_binomial(qt[:3], qt[3:], 2, 4)
        assert b.p_over == pytest.approx(a.p_over, rel=RTOL)
        assert b.p_under == pytest.approx(a.p_under, rel=RTOL)


# ---------------------------------------------------------------- both tails, identity


def test_both_tails_include_k_and_sum_to_one_plus_pk() -> None:
    p_in, p_out = [0.2, 0.7, 0.45, 0.9], [0.1, 0.3, 0.5, 0.6, 0.8]
    for N in range(10):
        ex0 = cpb_exact(p_in, p_out, 0, N)
        assert ex0 is not None
        for k in range(ex0.lo, ex0.hi + 1):
            ex = cpb_exact(p_in, p_out, k, N)
            assert ex is not None
            t = conditional_poisson_binomial(p_in, p_out, k, N)
            assert_tail_close(t.p_over, ex.p_over)
            assert_tail_close(t.p_under, ex.p_under)
            assert abs(t.p_over + t.p_under - float(ex.p_eq) - 1.0) <= 1e-12
            assert (t.lo, t.hi) == (ex.lo, ex.hi)


def test_extreme_tails_keep_relative_precision() -> None:
    # never 1 - cdf: deep upper and lower tails stay accurate
    p_in = [0.001] * 40
    p_out = [0.9] * 40
    N = 40
    top = conditional_poisson_binomial(p_in, p_out, 40, N)
    ex = cpb_exact(p_in, p_out, 40, N)
    assert ex is not None
    assert float(ex.p_over) < 1e-150
    assert_tail_close(top.p_over, ex.p_over)
    assert rel_err(top.p_under, ex.p_under) <= RTOL
    low = conditional_poisson_binomial(p_out, p_in, 0, N)  # C = the 0.9 TUs, none detected
    ex = cpb_exact(p_out, p_in, 0, N)
    assert ex is not None
    assert float(ex.p_under) < 1e-150
    assert_tail_close(low.p_under, ex.p_under)


def test_tails_below_1e300_are_nonnegative_floats() -> None:
    p_in = [1e-4] * 60
    p_out = [1 - 1e-4] * 60
    t = conditional_poisson_binomial(p_in, p_out, 60, 60)
    ex = cpb_exact(p_in, p_out, 60, 60)
    assert ex is not None
    assert ex.p_over < Fraction(1, 10**300)
    assert math.isfinite(t.p_over) and t.p_over >= 0.0
    assert (t.lo, t.hi) == (0, 60)  # structural support: the k = 60 weight is positive
    assert t.p_under == pytest.approx(1.0, rel=RTOL)


def test_does_not_mutate_inputs() -> None:
    p_in = np.array([0.2, 0.4, 0.6])
    p_out = np.array([0.3, 0.5])
    a, b = p_in.copy(), p_out.copy()
    conditional_poisson_binomial(p_in, p_out, 1, 2)
    np.testing.assert_array_equal(p_in, a)
    np.testing.assert_array_equal(p_out, b)


def test_accepts_lists_tuples_and_arrays_identically() -> None:
    a = conditional_poisson_binomial([0.2, 0.4], (0.3, 0.5, 0.7), 1, 2)
    b = conditional_poisson_binomial(np.array([0.2, 0.4]), np.array([0.3, 0.5, 0.7]), 1, 2)
    assert a == b


# ---------------------------------------------------------------- degenerate cases


def test_n_zero_support_is_zero_and_every_p_is_one() -> None:
    t = conditional_poisson_binomial([0.3, 0.6, 0.2], [0.5, 0.4], 0, 0)
    assert (t.lo, t.hi) == (0, 0)
    assert t.p_over == pytest.approx(1.0, rel=RTOL)
    assert t.p_under == pytest.approx(1.0, rel=RTOL)
    assert t.mean == 0.0


def test_n_equals_universe_support_is_n_and_every_p_is_one() -> None:
    t = conditional_poisson_binomial([0.3, 0.6, 0.2], [0.5, 0.4], 3, 5)
    assert (t.lo, t.hi) == (3, 3)
    assert t.p_over == pytest.approx(1.0, rel=RTOL)
    assert t.p_under == pytest.approx(1.0, rel=RTOL)
    assert t.mean == pytest.approx(3.0, rel=RTOL)


def test_empty_complement_has_support_n() -> None:
    t = conditional_poisson_binomial([0.3, 0.6, 0.2], [], 2, 2)
    assert (t.lo, t.hi) == (2, 2)
    assert t.p_over == pytest.approx(1.0, rel=RTOL)
    assert t.p_under == pytest.approx(1.0, rel=RTOL)
    assert t.mean == pytest.approx(2.0, rel=RTOL)


def test_empty_category_has_support_zero() -> None:
    t = conditional_poisson_binomial([], [0.3, 0.6], 0, 1)
    assert (t.lo, t.hi) == (0, 0)
    assert t.p_over == pytest.approx(1.0, rel=RTOL)
    assert t.mean == 0.0


def test_certain_and_impossible_units_narrow_the_support() -> None:
    # PB_C = [0, .5, .5] (one TU certain), PB_out = [.25, .5, .25]; N = 2
    t = conditional_poisson_binomial([1.0, 0.5], [0.5, 0.5], 1, 2)
    ex = cpb_exact([1.0, 0.5], [0.5, 0.5], 1, 2)
    assert ex is not None and (ex.lo, ex.hi) == (1, 2)
    assert (t.lo, t.hi) == (1, 2)
    assert rel_err(t.p_over, ex.p_over) <= RTOL
    assert rel_err(t.p_under, ex.p_under) <= RTOL
    # a q = 0 TU in C can never count
    t = conditional_poisson_binomial([0.0, 0.5, 0.5], [0.5], 1, 2)
    assert (t.lo, t.hi) == (1, 2)


def test_p_in_zero_one_mixture_matches_enumeration() -> None:
    p_in, p_out = [1.0, 0.0, 0.4, 0.7], [0.0, 1.0, 0.3, 0.5, 0.9]
    for N in range(2, 7):
        ex0 = cpb_enumerate(p_in, p_out, 0, N)
        assert ex0 is not None
        for k in range(ex0.lo, ex0.hi + 1):
            ex = cpb_enumerate(p_in, p_out, k, N)
            assert ex is not None
            t = conditional_poisson_binomial(p_in, p_out, k, N)
            assert rel_err(t.p_over, ex.p_over) <= RTOL
            assert rel_err(t.p_under, ex.p_under) <= RTOL
            assert rel_err(t.mean, ex.mean) <= RTOL
            assert (t.lo, t.hi) == (ex.lo, ex.hi)


@pytest.mark.parametrize(
    ("p_in", "p_out", "k", "N"),
    [
        ([0.5, 0.5], [1.0], 0, 0),  # a q = 1 TU is undetected (N = 0 impossible)
        ([0.0], [0.0], 1, 1),  # a q = 0 TU is detected: the support sum is 0
        ([0.0, 0.0], [0.0, 0.5], 0, 2),  # N = 2 needs a q = 0 TU detected
        ([1.0, 1.0], [0.5], 1, 1),  # k = 1 needs a q = 1 TU of C undetected
    ],
)
def test_probability_zero_observation_fails_closed(
    p_in: list[float], p_out: list[float], k: int, N: int
) -> None:
    with pytest.raises(ValueError):
        conditional_poisson_binomial(p_in, p_out, k, N)


@pytest.mark.parametrize(
    ("k", "N"),
    [
        (-1, 2),  # k below 0
        (4, 4),  # k > |C|
        (1, 4),  # k < N - |out| (lo = 2)
        (3, 2),  # k > N
        (0, -1),  # N < 0
        (1, 6),  # N > T
    ],
)
def test_k_or_n_outside_the_support_is_a_value_error(k: int, N: int) -> None:
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.3, 0.6, 0.2], [0.5, 0.4], k, N)


def test_k_outside_the_narrowed_support_is_a_value_error() -> None:
    # [max(0, N - 2), min(2, N)] = [0, 2] but a q = 1 TU of C makes K >= 1
    with pytest.raises(ValueError):
        conditional_poisson_binomial([1.0, 0.5], [0.5, 0.5], 0, 2)


@pytest.mark.parametrize("bad", [-0.1, 1.1, float("nan"), float("inf")])
def test_invalid_probability_is_a_value_error(bad: float) -> None:
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.2, bad], [0.5], 1, 1)
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.2], [0.5, bad], 1, 1)


# ---------------------------------------------------------------- the truncated pmf


def test_pmf_trunc_prefix_of_the_full_pmf() -> None:
    rng = np.random.default_rng(4)
    p = rng.beta(2, 2, 80)
    full = poisson_binomial_pmf(p)
    for kmax in (0, 1, 7, 40, 79, 80):
        got = poisson_binomial_pmf_trunc(p, kmax)
        assert isinstance(got, np.ndarray) and got.dtype == np.float64
        assert got.shape == (kmax + 1,)
        np.testing.assert_allclose(got, full[: kmax + 1], rtol=1e-15, atol=0)


def test_pmf_trunc_kmax_above_n_pads_with_zeros() -> None:
    got = poisson_binomial_pmf_trunc([0.5, 0.5], 4)
    np.testing.assert_allclose(got, [0.25, 0.5, 0.25, 0.0, 0.0], rtol=1e-15, atol=0)


def test_pmf_trunc_kmax_zero_and_empty_input() -> None:
    got = poisson_binomial_pmf_trunc([0.1, 0.2, 0.3], 0)
    np.testing.assert_allclose(got, [0.9 * 0.8 * 0.7], rtol=1e-15)
    np.testing.assert_array_equal(poisson_binomial_pmf_trunc([], 0), [1.0])
    np.testing.assert_array_equal(poisson_binomial_pmf_trunc([], 2), [1.0, 0.0, 0.0])


def test_pmf_trunc_rejects_bad_input() -> None:
    with pytest.raises(ValueError):
        poisson_binomial_pmf_trunc([0.5], -1)
    for bad in (-0.1, 1.5, float("nan")):
        with pytest.raises(ValueError):
            poisson_binomial_pmf_trunc([0.5, bad], 1)


def test_pmf_trunc_spec_precision_check_beta_1_4() -> None:
    # spec §1 step 1: 120 Beta(1, 4) probabilities against exact rationals over the full support
    q = np.random.default_rng(1).beta(1, 4, 120)
    got = poisson_binomial_pmf_trunc(q, 120)
    ex = pb_pmf_exact(list(q))
    worst = 0.0
    for j, e in enumerate(ex):
        assert_tail_close(float(got[j]), e)
        if e >= Fraction(1, 10**300):
            worst = max(worst, rel_err(float(got[j]), e))
    assert worst <= 1e-12


# ---------------------------------------------------------------- tilt


def _logit(x: np.ndarray) -> np.ndarray:
    return np.log(x) - np.log1p(-x)


def test_tilt_hits_the_total_and_scales_every_odds_alike() -> None:
    q = np.array([0.95, 0.5, 0.36, 0.5, 0.657, 0.1])
    for total in (1, 2, 3, 4, 5):
        qt = tilt(q, total)
        assert isinstance(qt, np.ndarray) and qt.dtype == np.float64 and qt.shape == q.shape
        assert abs(qt.sum() - total) <= 1e-9 * total
        theta = _logit(qt) - _logit(q)
        np.testing.assert_allclose(theta, theta[0], rtol=0, atol=1e-9)


def test_tilt_is_identity_when_the_total_is_already_met() -> None:
    q = np.array([0.2, 0.3, 0.5])
    qt = tilt(q, 1.0)
    np.testing.assert_allclose(qt, q, rtol=1e-9)


def test_tilt_keeps_zero_and_one_fixed() -> None:
    q = np.array([0.0, 1.0, 0.2, 0.4, 0.6, 1.0, 0.0])
    qt = tilt(q, 4)
    assert qt[0] == 0.0 and qt[6] == 0.0
    assert qt[1] == 1.0 and qt[5] == 1.0
    assert abs(qt.sum() - 4) <= 1e-9 * 4
    assert np.all((qt[2:5] > 0) & (qt[2:5] < 1))


def test_tilt_badly_scaled_q() -> None:
    # q = 1e-6 everywhere, total = T / 2: every tilted q is 1/2
    qt = tilt(np.full(1000, 1e-6), 500)
    np.testing.assert_allclose(qt, 0.5, rtol=1e-9)


def test_tilt_limits_at_the_ends_of_the_feasible_range() -> None:
    q = np.array([1.0, 0.3, 0.6, 0.0])
    lo = tilt(q, 1)  # only the q = 1 TU: the open entries go to 0
    assert np.all(np.isfinite(lo)) and abs(lo.sum() - 1) <= 1e-9
    assert lo[0] == 1.0 and lo[3] == 0.0
    hi = tilt(q, 3)  # every q > 0 TU: the open entries go to 1
    assert np.all(np.isfinite(hi)) and abs(hi.sum() - 3) <= 3e-9
    assert hi[3] == 0.0
    zero = tilt(np.array([0.2, 0.5]), 0)
    assert np.all(zero >= 0) and zero.sum() <= 1e-12


@pytest.mark.parametrize("total", [-1, 0, 4, 5])
def test_tilt_infeasible_total_is_a_value_error(total: float) -> None:
    # sum must lie in [#(q == 1), #(q > 0)] = [1, 3]
    with pytest.raises(ValueError):
        tilt(np.array([1.0, 0.3, 0.6, 0.0]), total)


def test_tilt_rejects_invalid_probabilities() -> None:
    for bad in (-0.1, 1.5, float("nan")):
        with pytest.raises(ValueError):
            tilt(np.array([0.5, bad]), 1)


def test_tilt_does_not_mutate_its_input() -> None:
    q = np.array([0.1, 0.2, 0.3])
    c = q.copy()
    tilt(q, 2)
    np.testing.assert_array_equal(q, c)


def test_tilt_rescues_underflow_equal_q_gives_exact_hypergeometric() -> None:
    # spec acceptance: q = 1e-6 with N = T / 2 returns finite tails (the tilt)
    T, n, N, k = 1000, 100, 500, 60
    qt = tilt(np.full(T, 1e-6), N)
    t = conditional_poisson_binomial(qt[:n], qt[n:], k, N)
    over = sum((hypergeom_exact(T, n, N, j) for j in range(k, n + 1)), Fraction(0))
    under = sum((hypergeom_exact(T, n, N, j) for j in range(k + 1)), Fraction(0))
    assert rel_err(t.p_over, over) <= RTOL
    assert rel_err(t.p_under, under) <= RTOL
    assert t.mean == pytest.approx(N * n / T, rel=RTOL)


def test_tilt_rescues_underflow_unequal_q_matches_raw_rationals() -> None:
    # raw q ~ 1e-6: PB_out(N - j) ~ 1e-600 underflows in float; the conditional law does not
    rng = np.random.default_rng(29)
    q = 1e-6 * rng.uniform(0.2, 5.0, 160)
    n, N = 30, 80
    qt = tilt(q, N)
    assert np.all(np.isfinite(qt)) and np.all((qt > 0) & (qt < 1))
    for k in (5, 15, 25):
        t = conditional_poisson_binomial(qt[:n], qt[n:], k, N)
        ex = cpb_exact(list(q[:n]), list(q[n:]), k, N)  # exact law of the raw q
        assert ex is not None
        assert rel_err(t.p_over, ex.p_over) <= 1e-11
        assert rel_err(t.p_under, ex.p_under) <= 1e-11
        assert rel_err(t.mean, ex.mean) <= 1e-11
