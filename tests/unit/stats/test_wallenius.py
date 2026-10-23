"""Wallenius ORA with p_detect odds (issue M8-01; GIAP_SPEC.md §7.7.3, §8.2, §8.3).

Interface pinned here (docs/specs/M8-01-stats-wallenius-elim.md, module ``giap.stats.wallenius``):

- ``OraResult``: frozen dataclass ``(M, n, N, k, odds, p_over, p_under)``;
- ``category_odds(p_detect, in_cat) -> float``: mean(p in C) / mean(p outside C); NaN when C or
  its complement is empty or the outside mean is 0; 0.0 allowed, never inf;
- ``wallenius_tails(k, M, n, N, odds) -> (p_over, p_under)``: P(X >= k), P(X <= k), both direct
  sums of the Wallenius pmf (M8-index lead decision 5: never ``.sf(k - 1)`` / ``.cdf(k)``);
  k outside ``[max(0, N - (M - n)), min(n, N)]`` is a ``ValueError``;
- ``wallenius_ora(detected, p_detect, in_cat) -> OraResult`` over the universe it is given;
  NaN p_detect raises naming the indices.

Oracles are written from the definition only (goseq and BiasedUrn are never read or imported):
Wallenius' distribution is the number of red balls among N drawn one at a time without
replacement, each remaining ball drawn with probability proportional to its weight (red ω,
white 1). Tiny cases enumerate every ordered draw in exact rationals; larger cases use the
forward DP over (reds drawn, whites drawn), which sums positive terms only (no cancellation).
"""

from __future__ import annotations

import copy
import dataclasses
import inspect
import itertools
import math
from fractions import Fraction
from functools import cache

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.stats import hypergeom, nchypergeom_wallenius

from giap.stats.wallenius import OraResult, category_odds, wallenius_ora, wallenius_tails

# ----------------------------------------------------------------------------- oracles


def _ordered_enumeration(M: int, n: int, N: int, odds: Fraction) -> dict[int, Fraction]:
    """Exact pmf by summing over every ordered sequence of N distinct balls (balls 0..n-1 red)."""
    weight = [odds if i < n else Fraction(1) for i in range(M)]
    total = sum(weight)
    pmf: dict[int, Fraction] = {}
    for seq in itertools.permutations(range(M), N):
        prob = Fraction(1)
        left = total
        for ball in seq:
            prob *= weight[ball] / left
            left -= weight[ball]
        reds = sum(1 for ball in seq if ball < n)
        pmf[reds] = pmf.get(reds, Fraction(0)) + prob
    return pmf


def _exact_dp(M: int, n: int, N: int, odds: Fraction) -> dict[int, Fraction]:
    """Exact pmf by the urn recursion in rationals (state: reds left, whites left, draws left)."""

    @cache
    def go(r: int, w: int, d: int) -> tuple[tuple[int, Fraction], ...]:
        if d == 0:
            return ((0, Fraction(1)),)
        out: dict[int, Fraction] = {}
        tot = odds * r + w
        if r > 0:
            for x, v in go(r - 1, w, d - 1):
                out[x + 1] = out.get(x + 1, Fraction(0)) + v * odds * r / tot
        if w > 0:
            for x, v in go(r, w - 1, d - 1):
                out[x] = out.get(x, Fraction(0)) + v * w / tot
        return tuple(sorted(out.items()))

    return dict(go(n, M - n, N))


def _float_dp(M: int, n: int, N: int, odds: float) -> np.ndarray:
    """Float pmf over x = 0..n by the forward urn DP (positive terms only)."""
    b = M - n
    P = np.zeros((n + 1, N + 1))
    P[0, 0] = 1.0
    for i in range(n + 1):
        for j in range(min(N - i, b) + 1):
            if i == 0 and j == 0:
                continue
            v = 0.0
            if i > 0:
                v += P[i - 1, j] * odds * (n - i + 1) / (odds * (n - i + 1) + (b - j))
            if j > 0:
                v += P[i, j - 1] * (b - j + 1) / (odds * (n - i) + (b - j + 1))
            P[i, j] = v
    return np.array([P[i, N - i] if 0 <= N - i <= b else 0.0 for i in range(n + 1)])


def _support(M: int, n: int, N: int) -> tuple[int, int]:
    return max(0, N - (M - n)), min(n, N)


def _close(a: float, b: float, rel: float, abs_: float = 0.0) -> bool:
    return abs(a - b) <= max(rel * max(abs(a), abs(b)), abs_)


# scipy's Wallenius pmf (Fog's method) is not normalised to 1e-12: its sum over the support is
# off by up to ~1.0e-12 at M <= 500 (e.g. M, n, N, ω = 79, 16, 67, 2.0) and ~1.3e-11 at
# M = 4000, n = 1000, N = 2000, ω = 0.01. Identities that sum the whole support (p_over + p_under
# - pmf(k) = 1, a tail of 1 at the far end) therefore hold to this tolerance, not 1e-12.
PMF_SUM_TOL = 1e-10

TINY = [
    (4, 2, 2, Fraction(3)),
    (5, 2, 3, Fraction(3, 2)),
    (6, 3, 3, Fraction(1, 4)),
    (6, 1, 4, Fraction(5)),
    (5, 4, 2, Fraction(2, 7)),
]

MEDIUM = [
    (12, 5, 7, Fraction(1, 4)),
    (30, 10, 12, Fraction(7, 3)),
    (40, 8, 20, Fraction(5)),
    (25, 20, 15, Fraction(1, 3)),
]

# ----------------------------------------------------------------------------- interface


def test_signatures_are_pinned() -> None:
    assert list(inspect.signature(category_odds).parameters) == ["p_detect", "in_cat"]
    assert list(inspect.signature(wallenius_tails).parameters) == ["k", "M", "n", "N", "odds"]
    assert list(inspect.signature(wallenius_ora).parameters) == ["detected", "p_detect", "in_cat"]
    fields = [f.name for f in dataclasses.fields(OraResult)]
    assert fields == ["M", "n", "N", "k", "odds", "p_over", "p_under"]


def test_ora_result_is_frozen() -> None:
    r = OraResult(M=4, n=2, N=2, k=1, odds=1.0, p_over=0.5, p_under=0.5)
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.p_over = 0.1  # type: ignore[misc]


def test_tails_return_python_floats() -> None:
    out = wallenius_tails(4, 20, 5, 8, 1.7)
    assert isinstance(out, tuple) and len(out) == 2
    assert all(isinstance(v, float) for v in out)


# ----------------------------------------------------------------------------- oracle checks


def test_oracle_self_check_enumeration_equals_dp() -> None:
    for M, n, N, w in TINY:
        enum = _ordered_enumeration(M, n, N, w)
        dp = _exact_dp(M, n, N, w)
        assert enum == dp
        assert sum(enum.values()) == 1
        dp_f = _float_dp(M, n, N, float(w))
        for x, v in enum.items():
            assert dp_f[x] == pytest.approx(float(v), rel=1e-13, abs=0)


def test_hand_computed_pmf_smallest_case() -> None:
    # M = 4 (2 red of weight 3, 2 white of weight 1), N = 2.
    # P(X = 2) = (6/8)(3/5) = 9/20; P(X = 0) = (2/8)(1/7) = 1/28;
    # P(X = 1) = 1 - 9/20 - 1/28 = 18/35.
    pmf = _ordered_enumeration(4, 2, 2, Fraction(3))
    assert pmf == {0: Fraction(1, 28), 1: Fraction(18, 35), 2: Fraction(9, 20)}
    p_over, p_under = wallenius_tails(2, 4, 2, 2, 3.0)
    assert p_over == pytest.approx(9 / 20, rel=1e-12, abs=0)
    assert p_under == pytest.approx(1.0, abs=1e-12)
    p_over, p_under = wallenius_tails(1, 4, 2, 2, 3.0)
    assert p_over == pytest.approx(18 / 35 + 9 / 20, rel=1e-12, abs=0)
    assert p_under == pytest.approx(1 / 28 + 18 / 35, rel=1e-12, abs=0)
    p_over, p_under = wallenius_tails(0, 4, 2, 2, 3.0)
    assert p_over == pytest.approx(1.0, abs=1e-12)
    assert p_under == pytest.approx(1 / 28, rel=1e-12, abs=0)


@pytest.mark.parametrize(("M", "n", "N", "w"), TINY + MEDIUM)
def test_both_tails_match_exact_urn_oracle(M: int, n: int, N: int, w: Fraction) -> None:
    pmf = _exact_dp(M, n, N, w)
    lo, hi = _support(M, n, N)
    assert set(pmf) == set(range(lo, hi + 1))
    for k in range(lo, hi + 1):
        over = float(sum(v for x, v in pmf.items() if x >= k))
        under = float(sum(v for x, v in pmf.items() if x <= k))
        p_over, p_under = wallenius_tails(k, M, n, N, float(w))
        assert p_over == pytest.approx(over, rel=1e-10, abs=0), k
        assert p_under == pytest.approx(under, rel=1e-10, abs=0), k


# ----------------------------------------------------------------------------- spec examples


def test_example1_scipy_sf_underflow_regression() -> None:
    # sf(54) = 1 - cdf(54) cancels to 0.0 on scipy 1.18.1; the true tail is ~3.2e-14.
    M, n, N, k, w = 4000, 60, 1500, 55, 1.3
    p_over, _ = wallenius_tails(k, M, n, N, w)
    assert p_over > 0.0
    assert p_over == pytest.approx(3.2049528951e-14, rel=1e-9, abs=0)
    direct = float(nchypergeom_wallenius(M, n, N, w).pmf(np.arange(n, k - 1, -1)).sum())
    assert p_over == pytest.approx(direct, rel=1e-6, abs=0)
    oracle = float(_float_dp(M, n, N, w)[k:].sum())
    assert p_over == pytest.approx(oracle, rel=1e-9, abs=0)


def test_example1_scipy_sf_still_buggy_here() -> None:
    """Records why the direct sum exists: the regression case still breaks ``.sf(k - 1)``."""
    sf = float(nchypergeom_wallenius(4000, 60, 1500, 1.3).sf(54))
    if sf != 0.0:
        pytest.skip(f"scipy sf(54) no longer underflows ({sf!r}); example 1 still pins our value")
    assert sf == 0.0


def test_example2_moderate_tail_agrees_with_sf() -> None:
    M, n, N, w = 4000, 60, 1500, 1.3
    p_over, _ = wallenius_tails(40, M, n, N, w)
    assert p_over == pytest.approx(7.356215112e-4, rel=1e-9, abs=0)
    assert p_over == pytest.approx(float(nchypergeom_wallenius(M, n, N, w).sf(39)), rel=1e-8, abs=0)
    assert p_over == pytest.approx(float(_float_dp(M, n, N, w)[40:].sum()), rel=1e-9, abs=0)


def test_tiny_lower_tail_matches_oracle() -> None:
    M, n, N, w = 4000, 60, 1500, 0.7
    oracle = _float_dp(M, n, N, w)
    for k in (0, 1, 2, 5):
        _, p_under = wallenius_tails(k, M, n, N, w)
        assert p_under > 0.0
        assert p_under == pytest.approx(float(oracle[: k + 1].sum()), rel=1e-9, abs=0)
    _, p_under = wallenius_tails(0, M, n, N, w)
    assert p_under == pytest.approx(2.0900264648e-09, rel=1e-8, abs=0)


def test_example3_odds_and_goseq_identity() -> None:
    p = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]
    in_cat = [True, True, False, False, False, False]
    odds = category_odds(p, in_cat)
    assert isinstance(odds, float)
    assert odds == pytest.approx(3.0, rel=1e-12, abs=0)
    avg_w, n_genes, n, total = 0.9, 6, 2, sum(p)
    assert avg_w * (n_genes - n) / (total - n * avg_w) == pytest.approx(3.0, rel=1e-12, abs=0)


@given(st.data())
def test_odds_equals_goseq_form(data: st.DataObject) -> None:
    size = data.draw(st.integers(2, 40))
    p = np.array(data.draw(st.lists(st.floats(0.01, 1.0), min_size=size, max_size=size)))
    in_cat = np.array(data.draw(st.lists(st.booleans(), min_size=size, max_size=size)))
    if in_cat.all() or not in_cat.any():
        in_cat[0] = not in_cat[0]
    n = int(in_cat.sum())
    avg_w = float(p[in_cat].mean())
    goseq = avg_w * (size - n) / (float(p.sum()) - n * avg_w)
    assert category_odds(p, in_cat) == pytest.approx(goseq, rel=1e-12, abs=0)


def test_example4_odds_one_is_central_hypergeometric() -> None:
    p_over, p_under = wallenius_tails(4, 20, 5, 8, 1.0)
    assert p_over == pytest.approx(float(hypergeom(20, 5, 8).sf(3)), abs=1e-12)
    assert p_under == pytest.approx(float(hypergeom(20, 5, 8).cdf(4)), abs=1e-12)


def test_example5_category_is_universe() -> None:
    p = [0.5, 0.4, 0.3, 0.2]
    r = wallenius_ora([True, False, True, False], p, [True] * 4)
    assert (r.M, r.n, r.N, r.k) == (4, 4, 2, 2)
    assert math.isnan(r.odds) and math.isnan(r.p_over) and math.isnan(r.p_under)
    assert math.isnan(category_odds(p, [True] * 4))


def test_empty_category() -> None:
    p = [0.5, 0.4, 0.3, 0.2]
    r = wallenius_ora([True, False, True, False], p, [False] * 4)
    assert (r.M, r.n, r.N, r.k) == (4, 0, 2, 0)
    assert math.isnan(r.odds) and math.isnan(r.p_over) and math.isnan(r.p_under)
    assert math.isnan(category_odds(p, [False] * 4))


def test_outside_mean_zero_is_nan_never_inf() -> None:
    odds = category_odds([0.5, 0.0, 0.0], [True, False, False])
    assert math.isnan(odds)
    r = wallenius_ora([True, False, False], [0.5, 0.0, 0.0], [True, False, False])
    assert math.isnan(r.odds) and math.isnan(r.p_over) and math.isnan(r.p_under)


def test_example6_nothing_detected() -> None:
    p = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]
    in_cat = [True, True, False, False, False, False]
    r = wallenius_ora([False] * 6, p, in_cat)
    assert (r.M, r.n, r.N, r.k) == (6, 2, 0, 0)
    assert r.odds == pytest.approx(3.0, rel=1e-12, abs=0)
    assert r.p_over == 1.0 and r.p_under == 1.0
    assert wallenius_tails(0, 6, 2, 0, 3.0) == (1.0, 1.0)


def test_everything_detected() -> None:
    p = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]
    in_cat = [True, True, False, False, False, False]
    r = wallenius_ora([True] * 6, p, in_cat)
    assert (r.M, r.n, r.N, r.k) == (6, 2, 6, 2)
    assert r.p_over == 1.0 and r.p_under == 1.0
    assert wallenius_tails(2, 6, 2, 6, 3.0) == (1.0, 1.0)


def test_example7_k_at_upper_bound_and_complement_identity() -> None:
    M, n, N, w = 4000, 5, 5, 2.0
    dist = nchypergeom_wallenius(M, n, N, w)
    assert _support(M, n, N) == (0, 5)
    p_over, p_under = wallenius_tails(5, M, n, N, w)
    assert p_over == pytest.approx(float(dist.pmf(5)), rel=1e-12, abs=0)
    assert p_under == pytest.approx(1.0, abs=1e-12)
    for k in range(6):
        o, u = wallenius_tails(k, M, n, N, w)
        assert o + u - float(dist.pmf(k)) == pytest.approx(1.0, abs=1e-12)


def test_k_at_lower_bound() -> None:
    M, n, N, w = 10, 4, 8, 0.6  # lo = 2, hi = 4
    dist = nchypergeom_wallenius(M, n, N, w)
    assert _support(M, n, N) == (2, 4)
    o, u = wallenius_tails(2, M, n, N, w)
    assert u == pytest.approx(float(dist.pmf(2)), rel=1e-12, abs=0)
    assert o == pytest.approx(1.0, abs=1e-12)
    o, u = wallenius_tails(4, M, n, N, w)
    assert o == pytest.approx(float(dist.pmf(4)), rel=1e-12, abs=0)


@pytest.mark.parametrize("k", [1, 5, -1])
def test_k_outside_support_raises(k: int) -> None:
    with pytest.raises(ValueError):
        wallenius_tails(k, 10, 4, 8, 0.6)  # support [2, 4]


def test_odds_zero() -> None:
    # ω = 0: no red ball is ever drawn while a white one is left.
    assert wallenius_tails(0, 10, 4, 5, 0.0) == (1.0, 1.0)
    p_over, p_under = wallenius_tails(1, 10, 4, 5, 0.0)
    assert p_over == 0.0
    assert p_under == 1.0
    p_over, _ = wallenius_tails(4, 10, 4, 5, 0.0)
    assert p_over == 0.0
    p = [0.0, 0.0, 0.5, 0.5, 0.5]
    in_cat = [True, True, False, False, False]
    assert category_odds(p, in_cat) == 0.0
    r = wallenius_ora([False, False, True, True, False], p, in_cat)
    assert r.odds == 0.0 and r.k == 0
    assert r.p_over == 1.0 and r.p_under == 1.0


def test_tail_below_1e300_is_underflow_not_nan_or_negative() -> None:
    M, n, N, k, w = 4000, 1000, 2000, 1000, 0.01
    p_over, p_under = wallenius_tails(k, M, n, N, w)
    assert math.isfinite(p_over) and 0.0 <= p_over < 1e-300
    assert p_under == pytest.approx(1.0, abs=PMF_SUM_TOL)
    p_over, p_under = wallenius_tails(0, 4000, 1000, 2000, 100.0)
    assert p_over == pytest.approx(1.0, abs=PMF_SUM_TOL)
    assert math.isfinite(p_under) and 0.0 <= p_under < 1e-300


# ----------------------------------------------------------------------------- wallenius_ora


def test_ora_counts_and_tails_on_example3_universe() -> None:
    p = [0.9, 0.9, 0.3, 0.3, 0.3, 0.3]
    in_cat = [True, True, False, False, False, False]
    detected = [True, False, True, False, False, False]
    r = wallenius_ora(detected, p, in_cat)
    assert isinstance(r, OraResult)
    assert (r.M, r.n, r.N, r.k) == (6, 2, 2, 1)
    assert r.odds == pytest.approx(3.0, rel=1e-12, abs=0)
    pmf = _exact_dp(6, 2, 2, Fraction(3))
    assert r.p_over == pytest.approx(float(pmf[1] + pmf[2]), rel=1e-12, abs=0)
    assert r.p_under == pytest.approx(float(pmf[0] + pmf[1]), rel=1e-12, abs=0)
    assert (r.p_over, r.p_under) == wallenius_tails(1, 6, 2, 2, r.odds)


def test_ora_nan_p_detect_raises_naming_indices() -> None:
    p = [0.9, 0.9, 0.3, float("nan"), 0.3, float("nan")]
    with pytest.raises(ValueError, match=r"3.*5"):
        wallenius_ora([True] * 3 + [False] * 3, p, [True, True] + [False] * 4)


def test_ora_length_mismatch_raises() -> None:
    with pytest.raises(ValueError):
        wallenius_ora([True, False], [0.5, 0.5, 0.5], [True, False, False])
    with pytest.raises(ValueError):
        wallenius_ora([True, False, True], [0.5, 0.5, 0.5], [True, False])


def test_ora_inputs_not_mutated_and_deterministic() -> None:
    rng = np.random.default_rng(7)
    p = rng.uniform(0.05, 0.95, size=200)
    in_cat = rng.random(200) < 0.15
    detected = rng.random(200) < p
    copies = [copy.deepcopy(x) for x in (detected, p, in_cat)]
    r1 = wallenius_ora(detected, p, in_cat)
    r2 = wallenius_ora(detected, p, in_cat)
    assert r1 == r2
    for orig, cp in zip((detected, p, in_cat), copies, strict=True):
        np.testing.assert_array_equal(orig, cp)
    lists = wallenius_ora(detected.tolist(), p.tolist(), in_cat.tolist())
    assert lists == r1


def test_ora_matches_manual_counts() -> None:
    rng = np.random.default_rng(11)
    p = rng.uniform(0.05, 0.95, size=300)
    in_cat = rng.random(300) < 0.1
    detected = rng.random(300) < p
    r = wallenius_ora(detected, p, in_cat)
    k = int((detected & in_cat).sum())
    assert (r.M, r.n, r.N, r.k) == (300, int(in_cat.sum()), int(detected.sum()), k)
    odds = float(p[in_cat].mean() / p[~in_cat].mean())
    assert r.odds == pytest.approx(odds, rel=1e-12, abs=0)
    over, under = wallenius_tails(k, 300, r.n, r.N, odds)
    assert r.p_over == pytest.approx(over, rel=1e-12, abs=0)
    assert r.p_under == pytest.approx(under, rel=1e-12, abs=0)


# ----------------------------------------------------------------------------- properties


@st.composite
def _urn(draw: st.DrawFn) -> tuple[int, int, int, float]:
    M = draw(st.integers(1, 500))
    n = draw(st.integers(0, M))
    N = draw(st.integers(0, M))
    w = draw(st.floats(0.05, 20.0))
    return M, n, N, w


def _check_complement_identity_and_range(M: int, n: int, N: int, w: float) -> None:
    lo, hi = _support(M, n, N)
    pmf = nchypergeom_wallenius(M, n, N, w).pmf(np.arange(lo, hi + 1))
    # The tails sum SciPy's pmf, so the identity holds only as well as SciPy normalises it.
    tol = PMF_SUM_TOL + abs(float(pmf.sum()) - 1.0)
    prev = math.inf
    for k in range(lo, hi + 1):
        o, u = wallenius_tails(k, M, n, N, w)
        assert 0.0 <= o <= 1.0 and 0.0 <= u <= 1.0
        assert o + u - float(pmf[k - lo]) == pytest.approx(1.0, abs=tol)
        assert o <= prev * (1 + 1e-12) + 1e-300  # p_over non-increasing in k
        prev = o


def test_complement_identity_where_scipy_pmf_sums_past_tolerance() -> None:
    # Hypothesis found this urn: SciPy's pmf sums to 1 + 1.08e-10 (> PMF_SUM_TOL).
    _check_complement_identity_and_range(151, 76, 112, 4.0)


@settings(max_examples=20, deadline=None)
@given(_urn())
def test_property_complement_identity_and_range(urn: tuple[int, int, int, float]) -> None:
    _check_complement_identity_and_range(*urn)


@settings(max_examples=20, deadline=None)
@given(_urn(), st.floats(1.01, 5.0), st.data())
def test_property_p_over_non_decreasing_in_odds(
    urn: tuple[int, int, int, float], factor: float, data: st.DataObject
) -> None:
    M, n, N, w = urn
    lo, hi = _support(M, n, N)
    k = data.draw(st.integers(lo, hi))
    o1, u1 = wallenius_tails(k, M, n, N, w)
    o2, u2 = wallenius_tails(k, M, n, N, w * factor)
    assert o1 <= o2 * (1 + 1e-9) + 1e-15
    assert u2 <= u1 * (1 + 1e-9) + 1e-15


@settings(max_examples=20, deadline=None)
@given(_urn(), st.data())
def test_property_odds_one_is_hypergeometric(
    urn: tuple[int, int, int, float], data: st.DataObject
) -> None:
    M, n, N, _ = urn
    lo, hi = _support(M, n, N)
    k = data.draw(st.integers(lo, hi))
    o, u = wallenius_tails(k, M, n, N, 1.0)
    h = hypergeom(M, n, N)
    assert o == pytest.approx(float(h.pmf(np.arange(k, hi + 1)).sum()), abs=1e-12)
    assert u == pytest.approx(float(h.pmf(np.arange(lo, k + 1)).sum()), abs=1e-12)


@pytest.mark.xfail(
    strict=False,
    reason="scipy's Wallenius pmf loses relative accuracy in some far tails (Fog's precision "
    "constants); the fix is a precision argument or a log-space pmf by ADR (spec: hardening owed)",
)
@pytest.mark.parametrize(
    ("M", "n", "N", "w", "k", "tail"),
    [
        (79, 16, 67, 2.0, 7, "under"),  # true 2.08e-11; scipy-pmf sum is off by ~2.9e-5 rel
        (300, 100, 150, 3.3, 99, "over"),  # true 2.99e-12; scipy-pmf sum is off by ~7.2e-5 rel
        (79, 16, 67, 2.0, 4, "under"),  # true 7.76e-19; scipy pmf returns 0.0
    ],
)
def test_far_tail_accuracy_against_urn_oracle(
    M: int, n: int, N: int, w: float, k: int, tail: str
) -> None:
    pmf = _float_dp(M, n, N, w)
    p_over, p_under = wallenius_tails(k, M, n, N, w)
    if tail == "over":
        assert p_over == pytest.approx(float(pmf[k:].sum()), rel=1e-6, abs=0)
    else:
        assert p_under == pytest.approx(float(pmf[: k + 1].sum()), rel=1e-6, abs=0)
