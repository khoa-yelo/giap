"""Exact Poisson-binomial distribution by dynamic programming (GIPPA_SPEC.md §7.2.5, §8.3).

Shared by the chance-hit model (M3-09), the statistics engine (M4-02) and detection ORA (M8-04,
ADR-0029). Every recurrence step is a convex combination of non-negative numbers, so tiny tails
stay positive and exact; the cdf is summed directly from the pmf (never ``1 - sf``) and vice
versa.

``conditional_poisson_binomial`` is the law of K = Σ_C D_t given S = Σ_U D_t = N for independent
units (the conditional-Bernoulli null of docs/specs/M8-04 §1): two truncated DPs, multiplied and
normalised by their support sum. ``tilt_pair`` rescales every odds by one factor e^θ so that
Σ q = N, which leaves that law unchanged and keeps the DPs clear of underflow (the caller tilts
once per run). It returns each complement 1 - q computed directly, because a tilted q close to 1
keeps too little of 1 - q (decision 11); the conditional takes those complements as keywords.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass

import numpy as np
from numba import njit
from numpy.typing import ArrayLike
from scipy.special import expit

_TILT_RTOL = 1e-13  # |Σ q(θ) - total| target; the spec contract is 1e-9 · total
_TILT_MAXITER = 200
_TINY = float(np.finfo(np.float64).tiny)  # the clamp that keeps a tilted open entry open
_PAIR_ATOL = 1e-9  # |p + comp - 1| allowed for a caller's (p, comp) pair


@dataclass(frozen=True, slots=True)
class CpbTails:
    """Both tails, the mean and the support of K given S = N (both tails include k)."""

    k: int
    N: int
    p_over: float  # P(K >= k | S = N)
    p_under: float  # P(K <= k | S = N)
    mean: float  # E[K | S = N]
    lo: int  # support of K given N: the j of non-zero weight
    hi: int


def _validate(p: ArrayLike) -> np.ndarray:
    arr = np.asarray(p, dtype=np.float64).ravel()
    if not np.all((arr >= 0.0) & (arr <= 1.0)):  # NaN fails both comparisons
        raise ValueError("Poisson-binomial probabilities must lie in [0, 1]")
    return arr


def _validate_comp(comp: ArrayLike, p: np.ndarray) -> np.ndarray:
    c = _validate(comp)
    if c.shape != p.shape:
        raise ValueError(f"{c.size} complements for {p.size} probabilities")
    return c


@njit(cache=True, nogil=True)
def _pmf_band(
    arr: np.ndarray, comp: np.ndarray, kmax: int, kmin: int
) -> np.ndarray:  # pragma: no cover - numba-compiled; via _pmf_trunc and the conditional
    """P(X = j) for kmin <= j <= kmax (entries below kmin are 0 on return).

    ``comp`` holds the complements 1 - q (``1.0 - arr`` unless the caller kept them, decision 11).
    The recurrence pmf[j] <- pmf[j](1 - q) + pmf[j - 1] q, stopped at kmax (spec M8-04 §1 step
    1). After unit i, R more units follow, so an entry j < kmin - R can never reach the returned
    band and is skipped; every entry that is computed is the same as in the full DP. Two buffers
    alternate and the inner loop runs over non-overlapping views, so LLVM vectorises it; there is
    no reduction and no FMA contraction (no fastmath), so each entry is bit-identical to the
    in-place recurrence of the merged numpy ``_pmf``.
    """
    n = arr.shape[0]
    x = np.zeros(kmax + 1, dtype=np.float64)
    y = np.zeros(kmax + 1, dtype=np.float64)
    x[0] = 1.0
    for i in range(n):
        q = arr[i]
        r = comp[i]
        src = x if i % 2 == 0 else y
        dst = y if i % 2 == 0 else x
        top = min(i + 1, kmax)  # the support after unit i is 0..i+1
        low = max(kmin - (n - 1 - i), 0)
        if low == 0:
            dst[0] = src[0] * r
            low = 1
        below, same, out = src[low - 1 : top], src[low : top + 1], dst[low : top + 1]
        for j in range(top + 1 - low):
            out[j] = same[j] * r + below[j] * q
    res = y if n % 2 == 1 else x
    res[:kmin] = 0.0
    return res


def _pmf_trunc(arr: np.ndarray, kmax: int) -> np.ndarray:
    return _pmf_band(arr, 1.0 - arr, kmax, 0)


def _pmf(arr: np.ndarray) -> np.ndarray:
    return _pmf_trunc(arr, arr.size)


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


def poisson_binomial_pmf_trunc(p: ArrayLike, kmax: int) -> np.ndarray:
    """P(X = j) for j = 0..kmax (length kmax + 1, zeros above n); the DP stopped at kmax."""
    if kmax < 0:
        raise ValueError("kmax must be >= 0")
    return _pmf_trunc(_validate(p), int(kmax))


def tilt_pair(
    p: ArrayLike, total: float, *, comp: ArrayLike | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """(q, 1 - q) with every odds x e^θ and Σ q = total (spec M8-04 §1 step 2, decision 11).

    θ solves Σ expit(logit q + θ) = total by Newton with a bisection guard; the complement is
    expit(-(logit q + θ)), never 1 - q. q ∈ {0, 1} stay exact (complement 1 or 0); every open
    entry keeps both values >= the smallest normal double, so it is never taken as certain.
    ``total`` must lie in [#(q = 1), #(q > 0)]; at the two ends the limit θ → ∓∞ is returned
    (open entries exactly 0 or 1, decision 10).

    With ``comp`` the units are the (p, comp) pairs (decision 12, ADR-0035): the logit is
    log p - log comp, so a p that rounds to 1.0 with comp > 0 stays open; comp == 0 is certain,
    p == 0 impossible, and the feasible range is [#(comp == 0), #(p > 0)]. A pair with
    p == comp == 0 or with |p + comp - 1| > 1e-9 is a ``ValueError``.
    """
    q = _validate(p)
    if comp is None:
        out, c, _ = _tilt(q, total)
        return out, c
    c = _validate_comp(comp, q)
    if np.any((q == 0.0) & (c == 0.0)) or np.any(np.abs(q + c - 1.0) > _PAIR_ATOL):
        raise ValueError("(p, comp) pairs must be probabilities summing to 1")
    out, c, _ = _tilt(q, total, c)
    return out, c


def tilt(p: ArrayLike, total: float) -> np.ndarray:
    """``tilt_pair(p, total)[0]``; a ``ValueError`` when one array cannot keep an entry open.

    That is an open entry whose tilted value rounds to 1.0, or underflows to 0.0 before the clamp;
    use ``tilt_pair`` and pass the complements to ``conditional_poisson_binomial`` then.
    """
    out, _, saturated = _tilt(_validate(p), total)
    if saturated:
        raise ValueError("an open entry rounds to 0 or 1 after the tilt; use tilt_pair")
    return out


def _tilt(
    q: np.ndarray, total: float, c: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray, bool]:
    """The tilted (q, 1 - q) and whether an interior tilt rounded an open entry to 0 or 1.

    ``c`` is the caller's complement (``tilt_pair(..., comp=)``); None means 1 - q.
    """
    if c is None:
        certain = q == 1.0
        open_ = (q > 0.0) & (q < 1.0)
        out = q.copy()
        comp = 1.0 - q  # exact on q ∈ {0, 1}; the open entries are overwritten
    else:
        certain = c == 0.0
        open_ = (q > 0.0) & (c > 0.0)
        out = np.where(certain, 1.0, 0.0)  # the closed pairs exactly (1, 0) and (0, 1)
        comp = 1.0 - out
    n_one = int(np.count_nonzero(certain))
    n_open = int(np.count_nonzero(open_))
    if not n_one <= total <= n_one + n_open:
        raise ValueError(f"tilt total {total} outside [{n_one}, {n_one + n_open}]")
    target = float(total) - n_one
    if target == 0.0:
        out[open_], comp[open_] = 0.0, 1.0
        return out, comp, False
    if target == n_open:
        out[open_], comp[open_] = 1.0, 0.0
        return out, comp, False
    tail = np.log1p(-q[open_]) if c is None else np.log(c[open_])
    logit = np.log(q[open_]) - tail
    s = logit + _solve_theta(logit, target)
    raw = expit(s)
    saturated = bool(np.any((raw == 0.0) | (raw == 1.0)))
    out[open_] = np.maximum(raw, _TINY)
    comp[open_] = np.maximum(expit(-s), _TINY)
    return out, comp, saturated


def _solve_theta(logit: np.ndarray, target: float) -> float:
    """θ with Σ expit(logit + θ) = target, 0 < target < len(logit); the sum is increasing in θ."""

    def f(theta: float) -> tuple[float, float]:
        s = expit(logit + theta)
        return float(s.sum()) - target, float((s * (1.0 - s)).sum())

    a, b = -1.0, 1.0
    while f(a)[0] > 0.0:
        a *= 2.0
    while f(b)[0] < 0.0:
        b *= 2.0
    theta = 0.0
    for _ in range(_TILT_MAXITER):
        g, dg = f(theta)
        if abs(g) <= _TILT_RTOL * target:
            break
        if g < 0.0:
            a = theta
        else:
            b = theta
        step = theta - g / dg if dg > 0.0 else float("nan")  # nan: bisect
        theta = step if a < step < b else 0.5 * (a + b)
        if b - a <= 4.0 * np.finfo(np.float64).eps * max(1.0, abs(theta)):
            break
    return theta


def conditional_poisson_binomial(
    p_in: ArrayLike,
    p_out: ArrayLike,
    k: int,
    N: int,
    *,
    comp_in: ArrayLike | None = None,
    comp_out: ArrayLike | None = None,
) -> CpbTails:
    """Tails of K = Σ_{t∈C} D_t given S = N (spec M8-04 §1; ADR-0029); exact, no approximation.

    P(K = j | N) ∝ PB_C(j) · PB_Ū(N - j), normalised by the sum over the support. ``p_over`` and
    ``p_under`` both include k, are sums of non-negative terms and are clamped to [0, 1] (the
    numerator and the support sum round apart by an ulp, decision 13). The support [lo, hi] is
    [max(0, N - |Ū|), min(|C|, N)] narrowed by the exact zeros of q (impossible units) and of
    1 - q (certain units); a k outside it, an N outside [0, |C| + |Ū|], or a zero or non-finite
    support sum is a ``ValueError``. Callers pass tilted q and its complements (``tilt_pair``) so
    that no weight underflows; ``comp_in``/``comp_out`` default to 1 - q. k and N must be
    integers (``operator.index``).
    """
    a_in, a_out = _validate(p_in), _validate(p_out)
    c_in = 1.0 - a_in if comp_in is None else _validate_comp(comp_in, a_in)
    c_out = 1.0 - a_out if comp_out is None else _validate_comp(comp_out, a_out)
    k, N = operator.index(k), operator.index(N)
    if not 0 <= N <= a_in.size + a_out.size:
        raise ValueError(f"N = {N} outside [0, {a_in.size + a_out.size}]")
    lo = max(int(np.count_nonzero(c_in == 0.0)), N - int(np.count_nonzero(a_out > 0.0)), 0)
    hi = min(int(np.count_nonzero(a_in > 0.0)), N - int(np.count_nonzero(c_out == 0.0)), N)
    if lo > hi:
        raise ValueError(f"S = {N} has probability 0 under q")
    if not lo <= k <= hi:
        raise ValueError(f"k = {k} outside the support [{lo}, {hi}] given N = {N}")
    pb_in = _pmf_band(a_in, c_in, hi, lo)[lo:]
    pb_out = _pmf_band(a_out, c_out, N - lo, N - hi)[N - hi :][::-1]  # PB_Ū(N - j), j = lo..hi
    w = pb_in * pb_out
    total = float(w.sum())
    if not (total > 0.0 and np.isfinite(total)):
        raise ValueError(f"the support weights given N = {N} sum to {total}")
    i = k - lo
    return CpbTails(
        k=k,
        N=N,
        p_over=min(float(w[i:].sum()) / total, 1.0),  # the two sums round apart: 1 + ulp
        p_under=min(float(w[: i + 1].sum()) / total, 1.0),
        mean=float((np.arange(lo, hi + 1) * w).sum()) / total,
        lo=lo,
        hi=hi,
    )
