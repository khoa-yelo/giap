"""Permutation p-values, circular shifts and Besag-Clifford sequential p-values.

GIAP_SPEC.md §7.3.5(b), §8.2, §8.3, decision 14, Appendix B.4. All randomness comes from the
``np.random.Generator`` passed in by the caller (obtained via ``giap.rng.for_step``).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Literal

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

__all__ = [
    "Alternative",
    "apply_shift",
    "besag_clifford",
    "choose_shifts",
    "circular_shifts",
    "perm_p",
]

Alternative = Literal["two-sided", "greater", "less"]

_TIE_RTOL = 1e-12


def _tol(*scale: float) -> float:
    """Float-tie tolerance ``1e-12 * max(1, |x|)``."""
    return _TIE_RTOL * max(1.0, *(abs(s) for s in scale))


def perm_p(obs: float, null: ArrayLike, alternative: Alternative = "two-sided") -> float:
    """Permutation p-value with the +1 rule, p = (1 + b)/(1 + B).

    NaN entries of ``null`` are dropped (B counts finite values); NaN ``obs`` or B = 0 gives NaN.
    The two-sided p is centred on the null mean. Values within ``1e-12 * max(1, |x|)`` count as
    ties, i.e. as exceedances (conservative).

    Raises:
        ValueError: on an unknown ``alternative``.
    """
    o = float(obs)
    z = np.asarray(null, dtype=np.float64).ravel()
    z = z[np.isfinite(z)]
    B = z.size
    if np.isnan(o) or B == 0:
        return float("nan")
    if alternative == "greater":
        b = int(np.count_nonzero(z >= o - _tol(o)))
    elif alternative == "less":
        b = int(np.count_nonzero(z <= o + _tol(o)))
    elif alternative == "two-sided":
        centre = float(z.mean())
        d_obs = abs(o - centre)
        b = int(np.count_nonzero(np.abs(z - centre) >= d_obs - _tol(d_obs, centre)))
    else:
        raise ValueError(f"unknown alternative {alternative!r}")
    return (1.0 + b) / (1.0 + B)


def _check_idx(n: int, idx: ArrayLike) -> npt.NDArray[np.int64]:
    a = np.asarray(idx, dtype=np.int64).ravel()
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if a.size == 0:
        raise ValueError("idx must not be empty")
    if a.min() < 0 or a.max() >= n:
        raise ValueError(f"idx must lie in [0, {n})")
    if np.unique(a).size != a.size:
        raise ValueError("idx must not contain duplicates")
    return a


def circular_shifts(n: int, idx: ArrayLike, circular: bool) -> npt.NDArray[np.int64]:
    """All valid shifts k of the gene-index set ``idx`` on a replicon of ``n`` genes, sorted.

    Circular: k = 1 ... n - 1. Linear: only k for which (idx + k) mod n does not wrap, i.e.
    ``max(idx) + k < n`` (rightward) or ``min(idx) + k >= n`` (leftward translation).

    Raises:
        ValueError: if ``idx`` is empty, has duplicates or leaves [0, n), or ``n < 1``.
    """
    a = _check_idx(n, idx)
    k = np.arange(1, n, dtype=np.int64)
    if circular:
        return k
    keep = (int(a.max()) + k < n) | (int(a.min()) + k >= n)
    return k[keep]


def apply_shift(idx: ArrayLike, k: int, n: int) -> npt.NDArray[np.int64]:
    """Shifted gene-index set ``(idx + k) mod n``, element-wise and in input order."""
    return np.mod(np.asarray(idx, dtype=np.int64) + int(k), int(n)).astype(np.int64)


def choose_shifts(
    valid: ArrayLike, n_max_exhaustive: int, n: int, B: int, rng: np.random.Generator
) -> tuple[npt.NDArray[np.int64], bool]:
    """Shifts to evaluate and whether they are exhaustive.

    All of ``valid`` when ``n <= n_max_exhaustive``; otherwise ``B`` shifts drawn without
    replacement from ``valid`` (all of them when ``B >= len(valid)``). Returned sorted.
    """
    v = np.unique(np.asarray(valid, dtype=np.int64).ravel())
    if n <= n_max_exhaustive:
        return v, True
    if v.size == 0 or v.size <= B:
        return v, False
    picked = rng.choice(v, size=int(B), replace=False)
    return np.sort(picked).astype(np.int64), False


def besag_clifford(
    stat_fn: Callable[[np.random.Generator], float],
    obs: float,
    h: int,
    max_B: int,
    rng: np.random.Generator,
) -> tuple[float, int]:
    """Sequential Monte Carlo p-value (Besag & Clifford 1991), alternative "greater".

    Draws null statistics ``stat_fn(rng)`` one at a time and stops at the h-th draw >= ``obs``
    (ties within ``1e-12 * max(1, |obs|)`` count), returning p = h / L with L the draws so far.
    If ``max_B`` draws give b < h exceedances, p = (b + 1)/(max_B + 1). For a two-sided test the
    caller passes centred absolute statistics.

    Returns:
        (p, B_used).

    Raises:
        ValueError: if ``h < 1`` or ``max_B < h``.
    """
    if h < 1:
        raise ValueError(f"h must be >= 1, got {h}")
    if max_B < h:
        raise ValueError(f"max_B must be >= h, got max_B={max_B}, h={h}")
    o = float(obs)
    thresh = o - _tol(o)
    b = 0
    for L in range(1, max_B + 1):
        if float(stat_fn(rng)) >= thresh:
            b += 1
            if b == h:
                return h / L, L
    return (b + 1.0) / (max_B + 1.0), max_B
