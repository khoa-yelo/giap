"""Exact Poisson-binomial distribution by dynamic programming (GIAP_SPEC.md §7.2.5, §8.3).

Shared by the chance-hit model (M3-09) and the statistics engine (M4-02). Every recurrence step
is a convex combination of non-negative numbers, so tiny tails stay positive and exact; the cdf
is summed directly from the pmf (never ``1 - sf``) and vice versa.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _validate(p: ArrayLike) -> np.ndarray:
    arr = np.asarray(p, dtype=np.float64).ravel()
    if not np.all((arr >= 0.0) & (arr <= 1.0)):  # NaN fails both comparisons
        raise ValueError("Poisson-binomial probabilities must lie in [0, 1]")
    return arr


def _pmf(arr: np.ndarray) -> np.ndarray:
    pmf = np.zeros(arr.size + 1, dtype=np.float64)
    pmf[0] = 1.0
    for i, q in enumerate(arr):
        # pmf[j] <- pmf[j] (1 - q) + pmf[j - 1] q over the support 0..i+1
        pmf[1 : i + 2] = pmf[1 : i + 2] * (1.0 - q) + pmf[: i + 1] * q
        pmf[0] *= 1.0 - q
    return pmf


def poisson_binomial_pmf(p: ArrayLike) -> np.ndarray:
    """P(X = j) for j = 0..n, X = sum of independent Bernoulli(p_i); exact DP, O(n^2)."""
    return _pmf(_validate(p))


def poisson_binomial_cdf(p: ArrayLike, k: int) -> float:
    """P(X <= k), summed directly from the pmf."""
    pmf = _pmf(_validate(p))
    if k < 0:
        return 0.0
    return float(pmf[: k + 1].sum())


def poisson_binomial_sf(p: ArrayLike, k: int) -> float:
    """P(X >= k), summed directly from the pmf (§8.3 contract)."""
    pmf = _pmf(_validate(p))
    return float(pmf[max(k, 0) :].sum())
