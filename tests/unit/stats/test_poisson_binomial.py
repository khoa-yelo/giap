"""Exact Poisson-binomial distribution (issue M3-09; GIAP_SPEC.md §7.2.5, §8.3).

Interface pinned here (shared with M4-02, M4-index decision 3; module
``giap.stats.poisson_binomial`` as named in M3-09 and M3-index):

- ``poisson_binomial_pmf(p) -> np.ndarray``: float64, length ``len(p) + 1``, exact DP;
- ``poisson_binomial_cdf(p, k) -> float``: P(X <= k), summed directly from the pmf (never
  ``1 - sf``), so tiny lower tails stay positive;
- ``poisson_binomial_sf(p, k) -> float``: P(X >= k) (§8.3 contract);
- any p outside [0, 1] (or NaN) is a ``ValueError``; empty ``p`` gives pmf ``[1.0]``.
"""

from __future__ import annotations

import inspect
import math

import numpy as np
import pytest

from giap.stats.poisson_binomial import (
    poisson_binomial_cdf,
    poisson_binomial_pmf,
    poisson_binomial_sf,
)


def test_signatures_are_pinned() -> None:
    assert list(inspect.signature(poisson_binomial_pmf).parameters) == ["p"]
    assert list(inspect.signature(poisson_binomial_cdf).parameters) == ["p", "k"]
    assert list(inspect.signature(poisson_binomial_sf).parameters) == ["p", "k"]


def test_example3_pmf_two_fair_coins() -> None:
    pmf = poisson_binomial_pmf([0.5, 0.5])
    assert isinstance(pmf, np.ndarray)
    assert pmf.dtype == np.float64
    assert pmf.shape == (3,)
    np.testing.assert_allclose(pmf, [0.25, 0.5, 0.25], rtol=1e-12, atol=0)


def test_example4_sf_three_coins() -> None:
    # 0.1*0.2*0.7 + 0.1*0.8*0.3 + 0.9*0.2*0.3 + 0.1*0.2*0.3 = 0.098
    sf = poisson_binomial_sf([0.1, 0.2, 0.3], 2)
    assert isinstance(sf, float)
    assert sf == pytest.approx(0.098, rel=1e-9)


def test_pmf_three_coins_full_vector() -> None:
    p = [0.1, 0.2, 0.3]
    expected = [0.9 * 0.8 * 0.7, 0.398, 0.092, 0.006]
    np.testing.assert_allclose(poisson_binomial_pmf(p), expected, rtol=1e-12)


def test_accepts_numpy_and_tuple_input() -> None:
    a = poisson_binomial_pmf(np.array([0.1, 0.2, 0.3]))
    b = poisson_binomial_pmf((0.1, 0.2, 0.3))
    np.testing.assert_array_equal(a, b)


def test_cdf_example_and_type() -> None:
    cdf = poisson_binomial_cdf([0.1, 0.2, 0.3], 1)
    assert isinstance(cdf, float)
    assert cdf == pytest.approx(0.504 + 0.398, rel=1e-12)


def test_cdf_and_sf_complement() -> None:
    p = [0.05, 0.5, 0.9, 0.33, 0.71]
    for k in range(1, len(p) + 1):
        assert poisson_binomial_cdf(p, k - 1) + poisson_binomial_sf(p, k) == pytest.approx(
            1.0, abs=1e-12
        )


def test_cdf_sf_boundaries() -> None:
    p = [0.2, 0.7, 0.4]
    assert poisson_binomial_sf(p, 0) == pytest.approx(1.0, abs=1e-15)
    assert poisson_binomial_cdf(p, len(p)) == pytest.approx(1.0, abs=1e-15)
    assert poisson_binomial_sf(p, len(p)) == pytest.approx(0.2 * 0.7 * 0.4, rel=1e-12)
    assert poisson_binomial_cdf(p, 0) == pytest.approx(0.8 * 0.3 * 0.6, rel=1e-12)


def test_empty_p() -> None:
    pmf = poisson_binomial_pmf([])
    np.testing.assert_array_equal(pmf, [1.0])
    assert poisson_binomial_cdf([], 0) == 1.0
    assert poisson_binomial_sf([], 0) == 1.0


def test_degenerate_probabilities() -> None:
    np.testing.assert_array_equal(poisson_binomial_pmf([0.0, 1.0, 1.0]), [0.0, 0.0, 1.0, 0.0])


@pytest.mark.parametrize("bad", [[-0.01], [1.01], [0.5, float("nan")], [0.2, float("inf")]])
def test_p_outside_unit_interval_is_value_error(bad: list[float]) -> None:
    with pytest.raises(ValueError):
        poisson_binomial_pmf(bad)
    with pytest.raises(ValueError):
        poisson_binomial_cdf(bad, 0)
    with pytest.raises(ValueError):
        poisson_binomial_sf(bad, 1)


def test_tiny_lower_tail_stays_positive_and_exact() -> None:
    """cdf is summed from the pmf, not 1 - sf: P(X = 0) = 0.01**40 = 1e-80."""
    p = [0.99] * 40
    cdf0 = poisson_binomial_cdf(p, 0)
    assert cdf0 > 0.0 and math.isfinite(cdf0)
    assert cdf0 == pytest.approx(0.01**40, rel=1e-9)
    cdf1 = poisson_binomial_cdf(p, 1)
    assert cdf1 == pytest.approx(0.01**40 + 40 * 0.99 * 0.01**39, rel=1e-9)


def test_tiny_upper_tail_stays_positive_and_exact() -> None:
    p = [1e-3] * 30
    sf = poisson_binomial_sf(p, 30)
    assert sf > 0.0
    assert sf == pytest.approx(1e-90, rel=1e-9)


def test_long_vector_is_exact_against_binomial() -> None:
    n, q = 500, 0.069
    pmf = poisson_binomial_pmf([q] * n)
    for j in (0, 10, 34, 35, 100):
        ref = math.comb(n, j) * q**j * (1 - q) ** (n - j)
        assert pmf[j] == pytest.approx(ref, rel=1e-9)
    assert pmf.sum() == pytest.approx(1.0, abs=1e-12)


def test_inputs_not_mutated() -> None:
    p = np.array([0.1, 0.9, 0.5])
    before = p.copy()
    poisson_binomial_pmf(p)
    poisson_binomial_cdf(p, 1)
    poisson_binomial_sf(p, 1)
    np.testing.assert_array_equal(p, before)
