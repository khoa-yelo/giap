"""``tilt_pair(p, total, *, comp=None)``: the complement input (issue M8-04b; docs/specs/M8-04
lead decision 12, ADR-0035 decision 2).

Interface pinned here:

- ``comp`` is keyword-only and optional; ``comp=None`` is the merged behaviour, bit for bit;
- with ``comp`` the units are the (p, comp) pairs: the tilt solves
  Σ expit(log p - log comp + θ) = total, so 1 - p is never re-derived from a rounded p. An
  entry with ``comp == 0`` is certain (its tilted pair is exactly (1, 0)); an entry with
  ``p == 0`` is impossible ((0, 1)); every other entry is open, *including* p == 1.0 with
  comp > 0 (a TU whose q rounds to 1.0 but whose log complement is finite);
- the feasible range is [#(comp == 0), #(p > 0)], the limits at its ends as in decision 10;
- contract: ``conditional_poisson_binomial`` on the tilted pairs (``comp_in``/``comp_out``)
  equals the exact law of the input pairs (``cpb_pairs_exact``: the law depends on the units
  only through p / comp) to relative 1e-11 on p_over, p_under and mean, for tails >= 1e-300.

Resolutions (reported to the lead): a comp array of another length, NaN or outside [0, 1] is a
``ValueError`` (as ``conditional_poisson_binomial``'s keywords); so is a pair with p == comp == 0
and a pair far from summing to 1 (here 0.3 / 0.3), since it is not a probability and its odds are
not what the caller meant.
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest
from helpers.cpb_oracle import assert_tail_close, rel_err
from helpers.cpb_pair_oracle import cpb_pairs_exact
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from scipy.special import expit

from gippa.stats.poisson_binomial import conditional_poisson_binomial, tilt_pair

CONTRACT = 1e-11
TINY = float(np.finfo(np.float64).tiny)


def _pairs(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(expit(x), expit(-x)): both computed directly, as ``TuUniverse`` keeps them."""
    return expit(x), expit(-x)


def _law(p: np.ndarray, c: np.ndarray, m: np.ndarray, k: int, N: int):  # type: ignore[no-untyped-def]
    pt, ct = tilt_pair(p, N, comp=c)
    return conditional_poisson_binomial(pt[m], pt[~m], k, N, comp_in=ct[m], comp_out=ct[~m])


def _assert_pairs_law(t, p, c, m, k, N) -> None:  # type: ignore[no-untyped-def]
    ex = cpb_pairs_exact(p[m], c[m], p[~m], c[~m], k, N)
    assert ex is not None
    assert (t.lo, t.hi) == (ex.lo, ex.hi)
    assert_tail_close(t.p_over, ex.p_over, rtol=CONTRACT)
    assert_tail_close(t.p_under, ex.p_under, rtol=CONTRACT)
    assert rel_err(t.mean, ex.mean) <= CONTRACT, (t.mean, float(ex.mean))


# ------------------------------------------------------------------------------ interface


def test_comp_is_an_optional_keyword() -> None:
    sig = inspect.signature(tilt_pair)
    assert "comp" in sig.parameters
    par = sig.parameters["comp"]
    assert par.kind is inspect.Parameter.KEYWORD_ONLY
    assert par.default is None


def test_without_comp_the_result_is_bit_identical() -> None:
    q = np.array([0.0, 1.0, 0.3, 0.9, 1e-6, 0.5, 0.999])
    a = tilt_pair(q, 3)
    b = tilt_pair(q, 3, comp=None)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


def test_returns_float64_arrays_of_the_input_shape_and_keeps_inputs() -> None:
    p, c = _pairs(np.linspace(-3.0, 3.0, 9))
    p0, c0 = p.copy(), c.copy()
    pt, ct = tilt_pair(p, 4, comp=c)
    assert pt.dtype == np.float64 and ct.dtype == np.float64
    assert pt.shape == p.shape and ct.shape == c.shape
    np.testing.assert_array_equal(p, p0)
    np.testing.assert_array_equal(c, c0)


def test_tilted_sum_hits_the_total() -> None:
    p, c = _pairs(np.linspace(-4.0, 2.0, 50))
    pt, ct = tilt_pair(p, 20, comp=c)
    assert abs(float(pt.sum()) - 20) <= 1e-9 * 20
    assert np.all(pt > 0) and np.all(ct > 0)
    assert np.max(np.abs(pt + ct - 1.0)) <= 4 * np.finfo(np.float64).eps


def test_matching_complements_agree_with_the_no_comp_path() -> None:
    """For moderate q, comp = 1 - q given explicitly tilts like the merged path (to rounding)."""
    q = np.array([0.1, 0.2, 0.35, 0.5, 0.7, 0.9])
    a, ac = tilt_pair(q, 2)
    b, bc = tilt_pair(q, 2, comp=1.0 - q)
    np.testing.assert_allclose(b, a, rtol=1e-13)
    np.testing.assert_allclose(bc, ac, rtol=1e-13)


# --------------------------------------------------------------- the rounded-q saturation case


def test_a_rounded_one_with_positive_comp_stays_open() -> None:
    """q = 1.0 in float but comp = e^-40: open, so a total excluding it is feasible."""
    p = np.array([1.0, 0.3, 0.3, 0.3])
    c = np.array([math.exp(-40.0), 0.7, 0.7, 0.7])
    pt, ct = tilt_pair(p, 0.5, comp=c)  # below #(comp == 0) + 1: only possible if unit 0 is open
    assert ct[0] > 0.0
    assert pt[0] > 0.0
    assert abs(float(pt.sum()) - 0.5) <= 1e-9


def test_rounded_one_undetected_is_a_valid_observation() -> None:
    """Without comp, a q of 1.0 is certain and k = 0 lies outside the support; with the kept
    complement the observation is possible and its tiny tail matches the exact pair law."""
    p = np.concatenate([[1.0], np.full(199, 0.01)])
    c = np.concatenate([[math.exp(-40.0)], np.full(199, 0.99)])
    m = np.zeros(200, dtype=bool)
    m[0] = True
    N = 3
    t = _law(p, c, m, 0, N)
    _assert_pairs_law(t, p, c, m, 0, N)
    assert 0.0 < t.p_under < 1e-12
    with pytest.raises(ValueError):
        pt, ct = tilt_pair(p, N)
        conditional_poisson_binomial(pt[m], pt[~m], 0, N, comp_in=ct[m], comp_out=ct[~m])


def test_exact_zero_comp_is_certain_and_exact_zero_p_is_impossible() -> None:
    p = np.array([1.0, 0.0, 0.4, 0.6])
    c = np.array([0.0, 1.0, 0.6, 0.4])
    pt, ct = tilt_pair(p, 2, comp=c)
    assert (pt[0], ct[0]) == (1.0, 0.0)
    assert (pt[1], ct[1]) == (0.0, 1.0)
    assert np.all(pt[2:] > 0) and np.all(ct[2:] > 0)


@pytest.mark.parametrize("total", [0.5, 3.5])
def test_feasible_range_counts_certain_units_by_comp(total: float) -> None:
    """[#(comp == 0), #(p > 0)] = [1, 3] for one certain, one impossible and two open units."""
    p = np.array([1.0, 0.0, 0.4, 0.6])
    c = np.array([0.0, 1.0, 0.6, 0.4])
    with pytest.raises(ValueError):
        tilt_pair(p, total, comp=c)


def test_limits_at_the_ends_of_the_range() -> None:
    p = np.array([1.0, 0.4, 1.0])
    c = np.array([0.0, 0.6, 1e-30])
    lo_p, lo_c = tilt_pair(p, 1, comp=c)
    np.testing.assert_array_equal(lo_p, [1.0, 0.0, 0.0])
    np.testing.assert_array_equal(lo_c, [0.0, 1.0, 1.0])
    hi_p, hi_c = tilt_pair(p, 3, comp=c)
    np.testing.assert_array_equal(hi_p, [1.0, 1.0, 1.0])
    np.testing.assert_array_equal(hi_c, [0.0, 0.0, 0.0])


def test_tilt_stress_odds_scaled_by_e_minus_8_with_units_near_one() -> None:
    """The calibration stress design (decision 12) on one fixed draw: every log-odds minus 8, 5 %
    of units at log-odds 36..44 (q rounds to 1.0); the tilted law equals the pair law."""
    rng = np.random.default_rng(20261022)
    T = 40
    b = rng.beta(1.2, 2.5, T)
    x = np.log(b) - np.log1p(-b)
    x[:2] = rng.uniform(36.0, 44.0, 2)
    x -= 8.0
    p, c = _pairs(x)
    m = np.zeros(T, dtype=bool)
    m[[0, 5, 9, 13, 21, 30]] = True
    for N in (2, 6, 15):
        for k in range(0, 7):
            try:
                t = _law(p, c, m, k, N)
            except ValueError:
                ex = cpb_pairs_exact(p[m], c[m], p[~m], c[~m], k, N)
                assert ex is None or not ex.lo <= k <= ex.hi
                continue
            _assert_pairs_law(t, p, c, m, k, N)


# ------------------------------------------------------------------------------- validation


def test_comp_of_the_wrong_length_is_a_value_error() -> None:
    with pytest.raises(ValueError):
        tilt_pair(np.array([0.2, 0.4, 0.6]), 1, comp=np.array([0.8, 0.6]))


@pytest.mark.parametrize("bad", [math.nan, -0.1, 1.5])
def test_invalid_comp_is_a_value_error(bad: float) -> None:
    with pytest.raises(ValueError):
        tilt_pair(np.array([0.2, 0.4, 0.6]), 1, comp=np.array([0.8, bad, 0.4]))


def test_a_pair_with_both_zero_is_a_value_error() -> None:
    with pytest.raises(ValueError):
        tilt_pair(np.array([0.0, 0.4, 0.6]), 1, comp=np.array([0.0, 0.6, 0.4]))


def test_a_pair_far_from_summing_to_one_is_a_value_error() -> None:
    with pytest.raises(ValueError):
        tilt_pair(np.array([0.3, 0.4, 0.6]), 1, comp=np.array([0.3, 0.6, 0.4]))


# --------------------------------------------------------------------- contract (Hypothesis)

_x = st.floats(min_value=-50.0, max_value=50.0, allow_nan=False)


@settings(max_examples=60, deadline=None, suppress_health_check=list(HealthCheck))
@given(
    xs=st.lists(_x, min_size=2, max_size=24),
    n_in=st.integers(min_value=1, max_value=23),
    data=st.data(),
)
def test_contract_against_the_exact_pair_law(xs: list[float], n_in: int, data) -> None:  # type: ignore[no-untyped-def]
    """Any log-odds in [-50, 50] (so some p round to 1.0 or comp to tiny), any split and any
    feasible (k, N): the tilted call equals ``cpb_pairs_exact`` on the raw pairs."""
    x = np.array(xs)
    T = x.size
    n_in = min(n_in, T - 1)
    p, c = _pairs(x)
    m = np.zeros(T, dtype=bool)
    m[:n_in] = True
    N = data.draw(st.integers(min_value=1, max_value=T - 1))
    ex0 = cpb_pairs_exact(p[m], c[m], p[~m], c[~m], 0, N)
    assert ex0 is not None  # every pair is open, so every N in [0, T] is possible
    k = data.draw(st.integers(min_value=ex0.lo, max_value=ex0.hi))
    t = _law(p, c, m, k, N)
    _assert_pairs_law(t, p, c, m, k, N)


def test_log_complement_input_from_a_large_tu() -> None:
    """A TU of 40 genes at p 0.6 has log comp = 40 log 0.4 (q = 1 - 1.2e-16): with comp = exp
    of the log complement it stays open; the law matches the pair oracle."""
    lc = np.array([40 * math.log(0.4), math.log(0.5), math.log(0.7), math.log(0.9), math.log(0.2)])
    p, c = -np.expm1(lc), np.exp(lc)
    m = np.array([True, False, True, False, False])
    for k, N in ((1, 2), (1, 1), (2, 3)):
        _assert_pairs_law(_law(p, c, m, k, N), p, c, m, k, N)
    assert TINY < c[0] < 1e-15
