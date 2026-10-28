"""Hodges-Lehmann two-sample shift (GIPPA_SPEC.md §7.3.5(b), §7.3.6, §8.3).

``hodges_lehmann(x, y)`` is the median of all pairwise differences ``x_i - y_j`` (mean of the two
middle values for an even count), computed exactly: by direct enumeration for m*n <= 5e4, and
above by bisection on the value, counting pairs per row by binary search on sorted y. Rounding
is monotone, so for fixed x_i the computed differences ``x_i - y_j`` are non-increasing over
ascending y; counts therefore use the same float differences as direct enumeration, and the
result equals ``np.median(np.subtract.outer(x, y))`` bit for bit.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

__all__ = ["hodges_lehmann"]

_DIRECT_MAX = 50_000  # lead decision (M4-04 review): bisection is exact and ~6x faster at 9e5 pairs
_FALLBACK_MAX = 20_000_000  # largest m*n enumerated directly when bisection cannot be used
_MAX_ITER = 4_096  # bisection halves a finite float64 interval; ~2100 steps reach adjacent floats

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]


def _first_le(x: F64, ys: F64, t: float) -> I64:
    """Per row i, the first j in [0, n] with ``x_i - ys[j] <= t`` (ys ascending)."""
    n = ys.size
    lo = np.zeros(x.size, dtype=np.int64)
    hi = np.full(x.size, n, dtype=np.int64)
    for _ in range(64):  # halves [0, n]; n < 2**63
        active = lo < hi
        if not active.any():
            return lo
        mid = (lo + hi) // 2
        cond = (x - ys[np.minimum(mid, n - 1)]) <= t
        hi = np.where(active & cond, mid, hi)
        lo = np.where(active & ~cond, mid + 1, lo)
    raise RuntimeError("unreachable: binary search did not converge")  # pragma: no cover


class _BisectionFailed(Exception):
    """Bisection cannot run (overflowing bounds) or did not converge within ``_MAX_ITER``."""


def _middle(x: F64, ys: F64, k1: int, k2: int) -> tuple[float, float]:
    """The k1-th and k2-th smallest (0-based, k1 <= k2) of the computed ``x_i - ys_j``.

    Raises:
        _BisectionFailed: if the bounds are not finite (differences overflow) or no convergence.
    """
    n = ys.size
    budget = 4 * (x.size + n)
    with np.errstate(over="ignore", invalid="ignore"):
        lo = float(np.nextafter(x.min() - ys[-1], -np.inf))  # count(d <= lo) = 0
        hi = float(x.max() - ys[0])  # count(d <= hi) = m*n
    if not (math.isfinite(lo) and math.isfinite(hi)):  # finite bounds give finite d
        raise _BisectionFailed
    j_lo = np.full(x.size, n, dtype=np.int64)  # per row, first j with d <= lo
    j_hi = np.zeros(x.size, dtype=np.int64)  # per row, first j with d <= hi
    c_lo = 0
    # invariant: count(d <= lo) = c_lo <= k1 <= k2 < count(d <= hi)
    for _ in range(_MAX_ITER):
        if int((j_lo - j_hi).sum()) <= budget:
            break
        mid = 0.5 * lo + 0.5 * hi
        if not math.isfinite(mid):  # pragma: no cover - finite lo, hi give a finite midpoint
            raise _BisectionFailed
        if mid <= lo or mid >= hi:  # adjacent floats: every value in (lo, hi] equals hi
            return hi, hi
        j_mid = _first_le(x, ys, mid)
        c_mid = int((n - j_mid).sum())
        if c_mid <= k1:
            lo, j_lo, c_lo = mid, j_mid, c_mid
        elif c_mid > k2:
            hi, j_hi = mid, j_mid
        else:  # mid separates k1 (<= mid) from k2 (> mid): select each one on its own
            return _middle(x, ys, k1, k1)[0], _middle(x, ys, k2, k2)[0]
    else:  # pragma: no cover - interval halves each step; adjacent floats return above
        raise _BisectionFailed
    cnt = j_lo - j_hi  # candidates in (lo, hi]: per row the j range [j_hi, j_lo)
    rows = np.repeat(np.arange(x.size), cnt)
    first = np.repeat(j_hi - (np.cumsum(cnt) - cnt), cnt)
    vals = x[rows] - ys[first + np.arange(rows.size)]
    r1, r2 = k1 - c_lo, k2 - c_lo
    part = np.partition(vals, [r1, r2])
    return float(part[r1]), float(part[r2])


def hodges_lehmann(x: ArrayLike, y: ArrayLike) -> float:
    """Median of the pairwise differences ``x_i - y_j``; NaN if either sample is empty.

    Differences that overflow float64 (inputs near +-1.8e308) are handled by direct enumeration,
    matching ``np.median(np.subtract.outer(x, y))``, when m*n <= ``_FALLBACK_MAX``.

    Raises:
        ValueError: if either input contains NaN or +-inf, or if the differences overflow and
            m*n exceeds ``_FALLBACK_MAX``.
    """
    xa = np.asarray(x, dtype=np.float64).ravel()
    ya = np.asarray(y, dtype=np.float64).ravel()
    if not (np.isfinite(xa).all() and np.isfinite(ya).all()):
        raise ValueError("hodges_lehmann: inputs contain NaN or inf; drop non-finite values first")
    m, n = xa.size, ya.size
    if m == 0 or n == 0:
        return math.nan
    total = m * n
    if total <= _DIRECT_MAX:
        return _direct(xa, ya)
    ys = np.sort(ya)
    k1, k2 = (total - 1) // 2, total // 2
    try:
        a, b = _middle(xa, ys, k1, k2)
    except _BisectionFailed:
        if total > _FALLBACK_MAX:
            raise ValueError(
                f"hodges_lehmann: pairwise differences overflow float64 and m*n = {total} is too "
                "large to enumerate; rescale the inputs"
            ) from None
        return _direct(xa, ya)
    return a if k1 == k2 else float(np.mean(np.array([a, b])))


def _direct(xa: F64, ya: F64) -> float:
    """Median of all pairwise differences by enumeration (overflow to +-inf as in numpy)."""
    with np.errstate(over="ignore", invalid="ignore"):
        return float(np.median(np.subtract.outer(xa, ya)))
