"""The tilt keeps complements (M8-04a review B1, PR #172; docs/specs/M8-04 lead decision 11).

Interface pinned here (decision 11):

- ``tilt_pair(p, total) -> tuple[np.ndarray, np.ndarray]``: the tilted ``(p, comp)``, both
  float64 of the input's shape, ``comp = expit(-(logit p + theta))`` computed directly (never as
  ``1 - p``). Exact 0s and 1s of the input stay exact with comp 1 or 0; every open entry has
  ``p > 0`` and ``comp > 0``, clamped below at the smallest positive normal double. At the two
  ends of the feasible range the limit is returned (open entries exactly 0 with comp 1, or exactly
  1 with comp 0), so that ``tilt == tilt_pair[0]`` there as well (decision 10 keeps the limit);
- ``tilt(p, total)``: unchanged signature, returns ``tilt_pair(p, total)[0]`` and raises
  ``ValueError`` when an interior tilt rounds an open entry to exactly 1.0, or underflows it to
  exactly 0.0 before the clamp;
- ``conditional_poisson_binomial(p_in, p_out, k, N, *, comp_in=None, comp_out=None)``: the
  complements are keyword-only and default to ``1 - p`` (bit-identical to the merged behaviour);
  the support is narrowed only by the exact 0s of p (impossible unit) and of comp (certain unit);
  ``k`` and ``N`` go through ``operator.index`` (a float is a ``TypeError``, N5);
- contract: on ``tilt_pair(q, total)`` split into C and its complement, the call equals the
  UNTILTED exact law of the raw q to 1e-11 relative on p_over, p_under and mean, and never raises
  when that law has a non-empty support.

Resolutions (see the test-author report): a comp array of the wrong length, or a comp outside
[0, 1] or NaN, is a ``ValueError``; the N5 rule applies to N as well as k; numpy integers stay
valid.
"""

from __future__ import annotations

import inspect
import math
from fractions import Fraction

import numpy as np
import pytest
from helpers.cpb_oracle import assert_tail_close, cpb_exact, rel_err
from helpers.cpb_pair_oracle import (
    cpb_grouped_exact,
    cpb_pairs_exact,
    dec_rel_err,
    tilted_reference,
)

from giap.stats.poisson_binomial import (
    conditional_poisson_binomial,
    tilt,
    tilt_pair,
)

RTOL = 1e-12
CONTRACT = 1e-11
TINY = float(np.finfo(np.float64).tiny)
EPS = float(np.finfo(np.float64).eps)


def _rev1() -> np.ndarray:
    """Reviewer reproduction 1: q[0] = 1 - 1e-15 among 999 q = 1e-6, total 500."""
    q = np.full(1000, 1e-6)
    q[0] = 1 - 1e-15
    return q


def _rev2() -> np.ndarray:
    """Reviewer reproduction 2: q[0] = 1 - 1e-9 among 199 q = 0.01, total 100."""
    q = np.full(200, 0.01)
    q[0] = 1 - 1e-9
    return q


def _call(p: np.ndarray, c: np.ndarray, n_in: int, k: int, N: int):  # type: ignore[no-untyped-def]
    return conditional_poisson_binomial(
        p[:n_in], p[n_in:], k, N, comp_in=c[:n_in], comp_out=c[n_in:]
    )


def _assert_law(t, ex) -> None:  # type: ignore[no-untyped-def]
    assert ex is not None
    assert (t.lo, t.hi) == (ex.lo, ex.hi)
    assert_tail_close(t.p_over, ex.p_over, rtol=CONTRACT)
    assert_tail_close(t.p_under, ex.p_under, rtol=CONTRACT)
    assert rel_err(t.mean, ex.mean) <= CONTRACT, (t.mean, float(ex.mean))


# ---------------------------------------------------------------- interface


def test_tilt_pair_signature_and_return_types() -> None:
    # M8-04b decision 12 adds a keyword-only comp=None after the two positional parameters
    assert list(inspect.signature(tilt_pair).parameters)[:2] == ["p", "total"]
    q = np.array([0.2, 0.5, 0.7])
    out = tilt_pair(q, 1)
    assert isinstance(out, tuple) and len(out) == 2
    for a in out:
        assert isinstance(a, np.ndarray) and a.dtype == np.float64 and a.shape == q.shape


def test_tilt_keeps_its_signature() -> None:
    assert list(inspect.signature(tilt).parameters) == ["p", "total"]


def test_conditional_takes_keyword_only_complements() -> None:
    params = inspect.signature(conditional_poisson_binomial).parameters
    assert list(params) == ["p_in", "p_out", "k", "N", "comp_in", "comp_out"]
    for name in ("comp_in", "comp_out"):
        assert params[name].kind is inspect.Parameter.KEYWORD_ONLY
        assert params[name].default is None
    with pytest.raises(TypeError):
        conditional_poisson_binomial([0.5], [0.5], 1, 1, [0.5], [0.5])  # type: ignore[misc]


# ---------------------------------------------------------------- tilt_pair


def test_tilt_pair_keeps_exact_zeros_and_ones() -> None:
    q = np.array([0.0, 1.0, 0.2, 0.4, 0.6, 1.0, 0.0])
    p, c = tilt_pair(q, 4)
    for i in (0, 6):
        assert p[i] == 0.0 and c[i] == 1.0
    for i in (1, 5):
        assert p[i] == 1.0 and c[i] == 0.0
    assert np.all((p[2:5] > 0) & (p[2:5] < 1) & (c[2:5] > 0) & (c[2:5] < 1))
    assert abs(p.sum() - 4) <= 1e-9 * 4


def test_tilt_pair_returns_the_limit_at_the_ends() -> None:
    q = np.array([1.0, 0.3, 0.6, 0.0])
    p, c = tilt_pair(q, 1)  # only the q = 1 unit: the open entries are exactly 0
    np.testing.assert_array_equal(p, [1.0, 0.0, 0.0, 0.0])
    np.testing.assert_array_equal(c, [0.0, 1.0, 1.0, 1.0])
    p, c = tilt_pair(q, 3)  # every q > 0 unit: the open entries are exactly 1
    np.testing.assert_array_equal(p, [1.0, 1.0, 1.0, 0.0])
    np.testing.assert_array_equal(c, [0.0, 0.0, 0.0, 1.0])


@pytest.mark.parametrize("total", [-1, 0, 4, 5])
def test_tilt_pair_infeasible_total_is_a_value_error(total: int) -> None:
    with pytest.raises(ValueError):
        tilt_pair(np.array([1.0, 0.3, 0.6, 0.0]), total)


def test_tilt_pair_rejects_invalid_probabilities() -> None:
    for bad in (-0.1, 1.5, float("nan")):
        with pytest.raises(ValueError):
            tilt_pair(np.array([0.5, bad]), 1)


def test_tilt_pair_does_not_mutate_its_input() -> None:
    q = _rev1()
    c0 = q.copy()
    tilt_pair(q, 500)
    np.testing.assert_array_equal(q, c0)


def test_reviewer_case1_saturated_entry_keeps_a_positive_complement() -> None:
    # logit q[0] + theta ~ 34.5 + 13.8 > 40: p rounds to 1.0 in floating point, comp ~ 1e-21
    q = _rev1()
    p, c = tilt_pair(q, 500)
    assert c[0] > 0.0
    assert abs(p.sum() - 500) <= 1e-9 * 500
    _, ref_c = tilted_reference(q, p, c)
    assert ref_c[0] is not None
    assert dec_rel_err(float(c[0]), ref_c[0]) <= RTOL
    assert float(c[0]) == pytest.approx(1e-21, rel=0.01)


def test_logit_plus_theta_below_minus_40_stays_positive() -> None:
    # theta ~ -30.6 (twenty units near 1, total 1): logit(1e-15) + theta ~ -65
    q = np.array([1e-15] + [1 - 1e-12] * 20)
    p, c = tilt_pair(q, 1)
    assert p[0] > 0.0 and c[0] > 0.0
    ref_p, ref_c = tilted_reference(q, p, c)
    assert ref_p[0] is not None and ref_c[0] is not None
    assert dec_rel_err(float(p[0]), ref_p[0]) <= RTOL
    assert dec_rel_err(float(c[0]), ref_c[0]) <= RTOL


@pytest.mark.parametrize(
    ("q", "total", "i", "side"),
    [
        # p underflows (exp(-721) is subnormal): clamped up to the smallest normal
        (np.array([1e-300] + [1 - 1e-12] * 20), 1, 0, "p"),
        # comp underflows to 0 (logit + theta ~ 780): clamped up to the smallest normal
        (np.array([1 - 2.0**-53, 5e-324, 5e-324, 5e-324]), 2, 0, "comp"),
    ],
)
def test_open_entries_are_clamped_to_the_smallest_normal(
    q: np.ndarray, total: int, i: int, side: str
) -> None:
    p, c = tilt_pair(q, total)
    v = p[i] if side == "p" else c[i]
    assert TINY <= v <= 1e-300
    m = (q > 0) & (q < 1)
    assert np.all(p[m] >= TINY) and np.all(c[m] >= TINY)
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(c))


def _cases() -> list[tuple[str, np.ndarray, int]]:
    rng = np.random.default_rng(1172)
    mixed = np.concatenate(
        [
            rng.uniform(0.01, 0.99, 40),
            1 - 10.0 ** -rng.uniform(9, 15, 6),
            10.0 ** -rng.uniform(9, 15, 6),
            [0.0, 1.0],
        ]
    )
    return [
        ("rev1", _rev1(), 500),
        ("rev2", _rev2(), 100),
        ("near1_pos_theta", np.array([1 - 1e-15, 1 - 1e-9, 0.5] + [1e-4] * 300), 150),
        ("near0_neg_theta", np.array([1e-15, 1e-9, 0.2] + [1 - 1e-12] * 30 + [0.5] * 10), 5),
        ("mixed_low", mixed, 10),
        ("mixed_high", mixed, 45),
    ]


@pytest.mark.parametrize(("name", "q", "total"), _cases(), ids=[c[0] for c in _cases()])
def test_tilt_pair_matches_the_exact_tilt(name: str, q: np.ndarray, total: int) -> None:
    # reference expit(+-(logit q + theta)) from the exact rational q in 60 digits, not 1 - p
    p, c = tilt_pair(q, total)
    assert abs(p.sum() - total) <= 1e-9 * total
    ref_p, ref_c = tilted_reference(q, p, c)
    worst = 0.0
    for i, (rp, rc) in enumerate(zip(ref_p, ref_c, strict=True)):
        if rp is None or rc is None:
            continue
        assert p[i] > 0.0 and c[i] > 0.0
        for got, ref in ((p[i], rp), (c[i], rc)):
            if ref >= 1e-300:
                worst = max(worst, dec_rel_err(float(got), ref))
            else:
                assert TINY <= got <= 1e-300
    assert worst <= RTOL, (name, worst)


def test_tilt_pair_p_and_comp_sum_to_one_within_rounding() -> None:
    for _, q, total in _cases():
        p, c = tilt_pair(q, total)
        assert np.all(np.abs(p + c - 1.0) <= 4 * EPS)


@pytest.mark.parametrize(
    ("q", "total"),
    [
        (np.array([0.95, 0.5, 0.36, 0.5, 0.657, 0.1]), 1),
        (np.array([0.95, 0.5, 0.36, 0.5, 0.657, 0.1]), 4),
        (np.array([0.0, 1.0, 0.2, 0.4, 0.6, 1.0, 0.0]), 4),
        (np.array([1.0, 0.3, 0.6, 0.0]), 1),  # the ends: limits
        (np.array([1.0, 0.3, 0.6, 0.0]), 3),
        (np.full(1000, 1e-6), 500),
        (_rev2(), 100),  # logit + theta ~ 25: close to 1 but not saturated
    ],
)
def test_tilt_is_the_first_array_of_tilt_pair_without_saturation(q: np.ndarray, total: int) -> None:
    np.testing.assert_array_equal(tilt(q, total), tilt_pair(q, total)[0])


def test_tilt_raises_when_an_open_entry_rounds_to_one() -> None:
    with pytest.raises(ValueError):
        tilt(_rev1(), 500)  # reviewer reproduction 1: tilted q[0] == 1.0


def test_tilt_raises_when_an_open_entry_underflows_to_zero() -> None:
    # logit(5e-324) + logit(1/40) ~ -748: expit underflows to exactly 0.0
    q = np.array([5e-324] + [0.5] * 40)
    with pytest.raises(ValueError):
        tilt(q, 1)
    p, c = tilt_pair(q, 1)  # tilt_pair keeps it open
    assert p[0] >= TINY and c[0] == pytest.approx(1.0)


# ---------------------------------------------------------------- the contract


def test_contract_reviewer_case1_both_k() -> None:
    q = _rev1()
    p, c = tilt_pair(q, 500)
    for k in (0, 1):
        t = _call(p, c, 1, k, 500)  # never raises: the support of the raw law is [0, 1]
        _assert_law(t, cpb_grouped_exact(q[:1], [(1e-6, 999)], k, 500))
    t0 = _call(p, c, 1, 0, 500)
    assert (t0.lo, t0.hi) == (0, 1)
    assert 0.0 < t0.p_under < 1e-20


def test_contract_reviewer_case2_p_under_at_zero() -> None:
    q = _rev2()
    p, c = tilt_pair(q, 100)
    for k in (0, 1):
        _assert_law(_call(p, c, 1, k, 100), cpb_grouped_exact(q[:1], [(0.01, 199)], k, 100))
    ex = cpb_grouped_exact(q[:1], [(0.01, 199)], 0, 100)
    assert ex is not None
    assert float(ex.p_under) == pytest.approx(1.0101009825e-11, rel=1e-9)  # the review's number
    assert rel_err(_call(p, c, 1, 0, 100).p_under, ex.p_under) <= CONTRACT


def test_contract_near_one_in_c_with_positive_theta() -> None:
    # within 1e-15 and 1e-9 of 1, theta ~ +9: logit + theta ~ 43.7 and 30
    p_in = [1 - 1e-15, 1 - 1e-9, 0.5]
    q = np.array(p_in + [1e-4] * 300)
    N = 150
    p, c = tilt_pair(q, N)
    for k in range(4):
        _assert_law(_call(p, c, 3, k, N), cpb_grouped_exact(p_in, [(1e-4, 300)], k, N))


def test_contract_near_one_in_the_complement_with_positive_theta() -> None:
    p_in = [0.3] * 10
    groups = [(1 - 1e-15, 2), (1 - 1e-9, 1), (1e-5, 400)]
    q = np.array(p_in + [g for g, n in groups for _ in range(n)])
    N = 200
    p, c = tilt_pair(q, N)
    for k in range(11):
        _assert_law(_call(p, c, 10, k, N), cpb_grouped_exact(p_in, groups, k, N))


def test_contract_near_zero_with_negative_theta() -> None:
    p_in = [1e-15, 1e-9, 0.2]
    groups = [(1 - 1e-12, 30), (0.5, 10)]
    q = np.array(p_in + [g for g, n in groups for _ in range(n)])
    N = 5
    p, c = tilt_pair(q, N)
    for k in range(4):
        _assert_law(_call(p, c, 3, k, N), cpb_grouped_exact(p_in, groups, k, N))


def test_contract_small_mixed_universe_against_the_rational_dp() -> None:
    rng = np.random.default_rng(4172)
    q = np.concatenate([rng.uniform(0.05, 0.95, 12), [1 - 1e-15, 1 - 1e-9, 1e-15, 1e-9, 0.0, 1.0]])
    rng.shuffle(q)
    n_in = 7
    for N in (2, 6, 10, 15):
        p, c = tilt_pair(q, N)
        ex0 = cpb_exact(list(q[:n_in]), list(q[n_in:]), 0, N)
        assert ex0 is not None
        for k in range(ex0.lo, ex0.hi + 1):
            ex = cpb_exact(list(q[:n_in]), list(q[n_in:]), k, N)
            _assert_law(_call(p, c, n_in, k, N), ex)


# ---------------------------------------------------------------- support narrowing


def test_float_one_with_positive_comp_is_not_certain_in_c() -> None:
    p_in, c_in, p_out, c_out = [1.0], [1e-20], [0.5, 0.5], [0.5, 0.5]
    t = conditional_poisson_binomial(p_in, p_out, 0, 1, comp_in=c_in, comp_out=c_out)
    assert (t.lo, t.hi) == (0, 1)
    ex = cpb_pairs_exact(p_in, c_in, p_out, c_out, 0, 1)
    assert ex is not None
    assert rel_err(t.p_under, ex.p_under) <= RTOL
    assert t.p_under == pytest.approx(2e-20, rel=1e-6)
    assert rel_err(t.mean, ex.mean) <= RTOL


def test_float_one_with_positive_comp_is_not_certain_in_the_complement() -> None:
    p_in, c_in, p_out, c_out = [0.5], [0.5], [1.0], [1e-20]
    t = conditional_poisson_binomial(p_in, p_out, 1, 1, comp_in=c_in, comp_out=c_out)
    assert (t.lo, t.hi) == (0, 1)
    ex = cpb_pairs_exact(p_in, c_in, p_out, c_out, 1, 1)
    assert ex is not None
    assert rel_err(t.p_over, ex.p_over) <= RTOL
    assert t.p_over == pytest.approx(1e-20, rel=1e-6)


def test_exact_zero_comp_is_certain_and_exact_zero_p_is_impossible() -> None:
    # unit 0 of C is certain (comp 0), unit 1 of C impossible (p 0): K = 1 + D_2
    p_in, c_in = [1.0, 0.0, 0.4], [0.0, 1.0, 0.6]
    p_out, c_out = [0.3, 1.0], [0.7, 0.0]  # the second unit of the complement is certain
    t = conditional_poisson_binomial(p_in, p_out, 1, 2, comp_in=c_in, comp_out=c_out)
    assert (t.lo, t.hi) == (1, 1)
    assert t.p_over == 1.0 and t.p_under == 1.0 and t.mean == 1.0
    t = conditional_poisson_binomial(p_in, p_out, 2, 3, comp_in=c_in, comp_out=c_out)
    assert (t.lo, t.hi) == (1, 2)
    with pytest.raises(ValueError):
        conditional_poisson_binomial(p_in, p_out, 0, 3, comp_in=c_in, comp_out=c_out)


def test_complements_of_the_wrong_length_are_a_value_error() -> None:
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.5, 0.5], [0.5], 1, 1, comp_in=[0.5], comp_out=[0.5])
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.5], [0.5], 1, 1, comp_in=[0.5], comp_out=[0.5, 0.5])


@pytest.mark.parametrize("bad", [-0.1, 1.5, float("nan")])
def test_invalid_complements_are_a_value_error(bad: float) -> None:
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.5], [0.5], 1, 1, comp_in=[bad], comp_out=[0.5])
    with pytest.raises(ValueError):
        conditional_poisson_binomial([0.5], [0.5], 1, 1, comp_in=[0.5], comp_out=[bad])


# ---------------------------------------------------------------- old behaviour


@pytest.mark.parametrize(
    ("p_in", "p_out", "k", "N"),
    [
        ([0.3, 0.6], [0.2, 0.4, 0.5], 1, 2),
        ([0.95, 0.36, 0.657], [0.5, 0.5, 0.1], 3, 4),
        ([1.0, 0.0, 0.4], [0.3, 1.0, 0.7], 2, 3),
        ([1e-6] * 30, [1e-6] * 50, 5, 20),
        ([], [0.2, 0.7], 0, 1),
    ],
)
def test_without_complements_the_result_is_unchanged(
    p_in: list[float], p_out: list[float], k: int, N: int
) -> None:
    a = conditional_poisson_binomial(p_in, p_out, k, N)
    b = conditional_poisson_binomial(p_in, p_out, k, N, comp_in=None, comp_out=None)
    c = conditional_poisson_binomial(
        p_in,
        p_out,
        k,
        N,
        comp_in=1.0 - np.asarray(p_in, dtype=float),
        comp_out=1.0 - np.asarray(p_out, dtype=float),
    )
    assert a == b == c  # bit-identical: the default complement is 1 - p


def test_without_complements_float_one_is_still_certain() -> None:
    t = conditional_poisson_binomial([1.0, 0.3], [0.5], 1, 1)
    assert (t.lo, t.hi) == (1, 1)
    with pytest.raises(ValueError):
        conditional_poisson_binomial([1.0, 0.3], [0.5], 0, 1)


# ---------------------------------------------------------------- review nits


def test_n1_untilted_underflow_hits_the_support_sum_guard() -> None:
    # untilted q = 1e-6 (T = 1000, C = 100, N = 500): every weight underflows to 0.0; the caller
    # must tilt (decision 10)
    q = np.full(1000, 1e-6)
    with pytest.raises(ValueError, match="sum"):
        conditional_poisson_binomial(q[:100], q[100:], 60, 500)


@pytest.mark.parametrize(("k", "N"), [(3.0, 4), (3.5, 4), (1, 4.0), (np.float64(2.0), 4)])
def test_n5_float_k_or_N_is_a_type_error(k: float, N: float) -> None:
    with pytest.raises(TypeError):
        conditional_poisson_binomial([0.3, 0.6, 0.5], [0.2, 0.4, 0.5], k, N)  # type: ignore[arg-type]


def test_n5_numpy_integers_are_accepted() -> None:
    a = conditional_poisson_binomial([0.3, 0.6, 0.5], [0.2, 0.4, 0.5], np.int64(2), np.int32(4))
    b = conditional_poisson_binomial([0.3, 0.6, 0.5], [0.2, 0.4, 0.5], 2, 4)
    assert a == b
    assert type(a.k) is int and type(a.N) is int
    assert math.isfinite(a.p_over)
    assert Fraction(a.p_over) >= 0
