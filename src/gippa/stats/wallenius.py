"""Wallenius over-representation analysis with p_detect odds (GIPPA_SPEC.md §7.7.3, §8.2, §8.3).

Per category C over the eligible universe: M genes, n in C, N detected, k detected in C, and the
odds ω = mean(p_detect in C) / mean(p_detect outside C) (goseq's weighting applied to
detectability). Both tails are direct sums of scipy's Wallenius pmf over the tail (``math.fsum``,
correctly rounded, so at least as accurate as summing from the far end): scipy's ``sf`` is
1 - cdf and cancels to 0.0 below ~1e-16 (M8-index lead decision 5). ω = 1 uses the central
hypergeometric pmf and ω = 0 is exact (no red ball is drawn while a white one is left).
Written from the definition; goseq and BiasedUrn are never used.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.stats import hypergeom, nchypergeom_wallenius

NAN = float("nan")
#: odds where scipy's Wallenius pmf is finite and sums to 1 (measured on scipy 1.18.1, #153 review);
#: outside it the pmf turns NaN and, at ~1e±300, scipy's C++ aborts the process
ODDS_RANGE = (1e-6, 1e6)
_PMF_SUM_TOL = 1e-6


@dataclass(frozen=True, slots=True)
class OraResult:
    """Counts, odds and both one-sided p-values of one category; NaN where undefined."""

    M: int
    n: int
    N: int
    k: int
    odds: float
    p_over: float
    p_under: float


def category_odds(p_detect: ArrayLike, in_cat: ArrayLike) -> float:
    """mean(p in C) / mean(p outside C) over finite p; NaN if a side is empty or outside mean 0."""
    p = np.asarray(p_detect, dtype=np.float64).ravel()
    cat = np.asarray(in_cat, dtype=bool).ravel()
    if p.shape != cat.shape:
        raise ValueError(f"p_detect has {p.size} entries but in_cat has {cat.size}")
    finite = np.isfinite(p)
    inside, outside = p[cat & finite], p[~cat & finite]
    if inside.size == 0 or outside.size == 0:
        return NAN
    mean_out = float(outside.mean())
    if mean_out == 0.0:
        return NAN
    return float(inside.mean()) / mean_out


def _tail_sum(pmf: np.ndarray) -> float:
    """Sum of non-negative pmf terms, clipped to [0, 1] (``fsum`` is order-independent)."""
    return min(max(math.fsum(pmf.tolist()), 0.0), 1.0)


def wallenius_tails(k: int, M: int, n: int, N: int, odds: float) -> tuple[float, float]:
    """(P(X >= k), P(X <= k)) for X ~ Wallenius(M, n, N, odds), both summed from the pmf."""
    if not 0 <= n <= M or not 0 <= N <= M:
        raise ValueError(f"need 0 <= n, N <= M; got M={M}, n={n}, N={N}")
    if not (math.isfinite(odds) and odds >= 0.0):
        raise ValueError(f"odds must be finite and >= 0; got {odds!r}")
    lo, hi = max(0, N - (M - n)), min(n, N)
    if not lo <= k <= hi:
        raise ValueError(f"k={k} is outside the support [{lo}, {hi}]")
    if lo == hi:
        return 1.0, 1.0
    if odds == 0.0:  # X = lo almost surely
        return (1.0 if k == lo else 0.0), 1.0
    if not ODDS_RANGE[0] <= odds <= ODDS_RANGE[1]:
        raise ValueError(
            f"odds={odds!r} is outside {ODDS_RANGE} where the Wallenius pmf is reliable "
            f"(M={M}, n={n}, N={N})"
        )
    xs = np.arange(lo, hi + 1)
    if odds == 1.0:
        pmf = np.asarray(hypergeom.pmf(xs, M, n, N), dtype=np.float64)
    else:
        pmf = np.asarray(nchypergeom_wallenius.pmf(xs, M, n, N, odds), dtype=np.float64)
    if not np.isfinite(pmf).all() or abs(math.fsum(pmf.tolist()) - 1.0) > _PMF_SUM_TOL:
        raise ValueError(
            f"Wallenius pmf is not a distribution at M={M}, n={n}, N={N}, odds={odds!r}"
        )
    j = k - lo
    return _tail_sum(pmf[j:]), _tail_sum(pmf[: j + 1])


def wallenius_ora(detected: ArrayLike, p_detect: ArrayLike, in_cat: ArrayLike) -> OraResult:
    """Wallenius ORA of one category over the universe given (no size filter, no BH)."""
    det = np.asarray(detected, dtype=bool).ravel()
    p = np.asarray(p_detect, dtype=np.float64).ravel()
    cat = np.asarray(in_cat, dtype=bool).ravel()
    if not det.size == p.size == cat.size:
        raise ValueError(
            f"lengths differ: detected {det.size}, p_detect {p.size}, in_cat {cat.size}"
        )
    bad = np.flatnonzero(~np.isfinite(p))
    if bad.size:
        raise ValueError(f"p_detect is not finite at indices {bad.tolist()}")
    M, n, N, k = p.size, int(cat.sum()), int(det.sum()), int((det & cat).sum())
    odds = category_odds(p, cat)
    if math.isnan(odds):
        return OraResult(M, n, N, k, odds, NAN, NAN)
    p_over, p_under = wallenius_tails(k, M, n, N, odds)
    return OraResult(M, n, N, k, odds, p_over, p_under)
