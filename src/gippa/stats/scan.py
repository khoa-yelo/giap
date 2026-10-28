"""Window scan, maximum-statistic null, FWER p-values and window merging (GIPPA_SPEC.md §7.3.5(d)).

Issue M5-03 (lead decisions 9-15). Windows of ``wmin .. wmax`` consecutive TUs, capped at n - 1 on
circular replicons (no window covers the whole circle) and at n on linear ones; circular windows
start at every position and wrap. S_w = Σ z / √n_w comes from prefix sums; the same vectorised
kernel serves one series (:func:`scan_max`) and B bootstrap rows (:func:`scan_null_max`), so the
two agree exactly. Windows are ordered by length, then start. Plain numpy: no JIT, so results do
not depend on ``NUMBA_DISABLE_JIT``.
"""

from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

from gippa.stats.permutation import _tol

__all__ = [
    "Cluster",
    "ScanResult",
    "fwer_p",
    "merge_windows",
    "n_windows",
    "scan_max",
    "scan_null_max",
]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

_EMPTY = np.zeros(0, dtype=np.int64)


@dataclass(frozen=True, slots=True)
class ScanResult:
    """All scanned windows of one series.

    Attributes
    ----------
    start : numpy.ndarray of int64
        0-based first position of each window (circular: in [0, n)).
    length : numpy.ndarray of int64
        Number of positions n_w of each window.
    s : numpy.ndarray of float64
        Signed S_w = Σ z / √n_w (positive = high, negative = low).
    max_abs : float
        max |S_w|; NaN when there is no window.
    """

    start: I64
    length: I64
    s: F64
    max_abs: float


@dataclass(frozen=True, slots=True)
class Cluster:
    """Connected component of overlapping significant windows of one direction.

    Attributes
    ----------
    direction : {"high", "low"}
        Sign of S of every member.
    best : int
        Index (into the :class:`ScanResult` arrays) of the representative window.
    members : tuple of int
        Sorted window indices in the cluster.
    """

    direction: Literal["high", "low"]
    best: int
    members: tuple[int, ...]


def _lengths(n: int, wmin: int, wmax: int, circular: bool) -> range:
    if wmin < 1 or wmin > wmax:
        raise ValueError(f"need 1 <= wmin <= wmax, got wmin={wmin}, wmax={wmax}")
    cap = n - 1 if circular else n
    return range(wmin, min(wmax, cap) + 1)


def n_windows(n: int, wmin: int, wmax: int, circular: bool) -> int:
    """Number of windows scanned on a replicon of ``n`` positions.

    Parameters
    ----------
    n : int
        Number of positions (scored TUs).
    wmin, wmax : int
        Window length range, 1 <= wmin <= wmax.
    circular : bool
        Circular replicon (windows wrap, length capped at n - 1).

    Returns
    -------
    int
        The window count (0 when no length fits).

    Raises
    ------
    ValueError
        If wmin < 1 or wmin > wmax.
    """
    return sum(n if circular else n - w + 1 for w in _lengths(n, wmin, wmax, circular))


def _sums(Z: F64, wmin: int, wmax: int, circular: bool) -> Iterator[tuple[int, F64]]:
    """Per window length w, the window sums of every row of ``Z`` (B, n) -> (B, m_w), by start."""
    if not np.isfinite(Z).all():
        raise ValueError("z must be finite (drop TUs without a score before scanning)")
    n = Z.shape[1]
    ws = _lengths(n, wmin, wmax, circular)
    if circular and len(ws):
        Z = np.concatenate([Z, Z[:, : ws[-1]]], axis=1)
    c = np.zeros((Z.shape[0], Z.shape[1] + 1), dtype=np.float64)
    np.cumsum(Z, axis=1, out=c[:, 1:])
    for w in ws:
        m = n if circular else n - w + 1
        yield w, c[:, w : w + m] - c[:, :m]


def _as_matrix(Z: ArrayLike) -> F64:
    a = np.asarray(Z, dtype=np.float64)
    if a.ndim != 2:
        raise ValueError("Z must be 2-D (B, n)")
    return a


def scan_max(z: ArrayLike, wmin: int, wmax: int, circular: bool) -> ScanResult:
    """Scan every window of ``z`` and return the window table and max |S|.

    Parameters
    ----------
    z : array_like of float, shape (n,)
        Detrended TU scores in replicon order; must be finite.
    wmin, wmax : int
        Window length range, 1 <= wmin <= wmax.
    circular : bool
        Circular replicon: windows wrap and are at most n - 1 long.

    Returns
    -------
    ScanResult
        Windows ordered by length, then start; empty with ``max_abs`` NaN when n is too small.

    Raises
    ------
    ValueError
        On a bad window range or non-finite z.
    """
    a = np.asarray(z, dtype=np.float64)
    if a.ndim != 1:
        raise ValueError("z must be 1-D")
    blocks = [(w, d[0] / math.sqrt(w)) for w, d in _sums(a[None, :], wmin, wmax, circular)]
    start = np.concatenate([np.arange(v.size, dtype=np.int64) for _, v in blocks] or [_EMPTY])
    length = np.concatenate([np.full(v.size, w, dtype=np.int64) for w, v in blocks] or [_EMPTY])
    row = np.concatenate([v for _, v in blocks] or [np.zeros(0)])
    max_abs = float(np.abs(row).max()) if row.size else float("nan")
    return ScanResult(start=start, length=length, s=row, max_abs=max_abs)


def scan_null_max(Z: ArrayLike, wmin: int, wmax: int, circular: bool) -> F64:
    """Maximum statistic M_b = max_w |S_w| of each bootstrap row.

    Parameters
    ----------
    Z : array_like of float, shape (B, n)
        Bootstrap replicates of the detrended scores; must be finite.
    wmin, wmax : int
        Window length range, 1 <= wmin <= wmax.
    circular : bool
        Circular replicon.

    Returns
    -------
    numpy.ndarray of float64, shape (B,)
        Equal to ``scan_max(Z[b], ...).max_abs`` row by row (NaN when there is no window).

    Raises
    ------
    ValueError
        On a bad window range, a non-2-D Z or non-finite entries.
    """
    a = _as_matrix(Z)
    out = np.full(a.shape[0], np.nan, dtype=np.float64)
    for w, d in _sums(a, wmin, wmax, circular):
        # correctly rounded division by sqrt(w) > 0 is monotone: max |d| / sqrt(w) == max |S|
        np.fmax(out, np.abs(d, out=d).max(axis=1, initial=0.0) / math.sqrt(w), out=out)
    return out


def fwer_p(s_abs: ArrayLike, null_max: ArrayLike) -> F64:
    """FWER-adjusted p_w = (1 + #{M_b >= |S_w|}) / (1 + B), element-wise ``perm_p(.., "greater")``.

    Parameters
    ----------
    s_abs : array_like of float, shape (k,)
        Observed |S_w|.
    null_max : array_like of float, shape (B,)
        Null maxima M_b; NaN entries are dropped.

    Returns
    -------
    numpy.ndarray of float64, shape (k,)
        p-values (ties within the engine's 1e-12 relative tolerance count); NaN when no M_b is
        valid or |S_w| is not finite.

    Raises
    ------
    ValueError
        If ``s_abs`` or ``null_max`` is not 1-D, or ``s_abs`` has a negative value (decision 22).
    """
    o = np.asarray(s_abs, dtype=np.float64)
    m = np.asarray(null_max, dtype=np.float64)
    if o.ndim != 1 or m.ndim != 1:
        raise ValueError("s_abs and null_max must be 1-D")
    if (np.isfinite(o) & (o < 0)).any():  # -0.0, NaN and ±inf pass (NaN p, decision 19)
        raise ValueError("s_abs must be |S| (got a negative value; pass np.abs(res.s))")
    m = np.sort(m[~np.isnan(m)])
    if m.size == 0:
        return np.full(o.size, np.nan, dtype=np.float64)
    thr = np.array([v - _tol(v) for v in o.tolist()], dtype=np.float64)
    b = m.size - np.searchsorted(m, thr, side="left")
    p: F64 = ((1.0 + b) / (1.0 + m.size)).astype(np.float64)
    p[~np.isfinite(o)] = np.nan  # an infinite |S| cannot come from finite z (decision 19)
    return p


def _find(parent: list[int], i: int) -> int:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


def merge_windows(res: ScanResult, significant: ArrayLike, n: int, circular: bool) -> list[Cluster]:
    """Merge significant same-direction windows sharing >= 1 position into clusters.

    Parameters
    ----------
    res : ScanResult
        The window table.
    significant : array_like of bool, shape (k,)
        Which windows are significant.
    n : int
        Number of positions of the replicon.
    circular : bool
        Positions wrap mod n.

    Returns
    -------
    list of Cluster
        Connected components (transitive); windows with S = 0 have no direction and are dropped.
        Best window = largest |S|, then fewer positions, then smaller start. Ordered by the best
        window's start, "high" before "low" on equal starts.

    Raises
    ------
    ValueError
        If ``significant`` is not 1-D or does not match the window count, or n < 0.
    """
    sig = np.asarray(significant, dtype=bool)
    if sig.ndim != 1:
        raise ValueError("significant must be 1-D")
    if n < 0:
        raise ValueError(f"n must be >= 0, got {n}")
    if sig.size != res.s.size:
        raise ValueError(f"significant has {sig.size} entries for {res.s.size} windows")
    idx = np.flatnonzero(sig & (res.s != 0)).tolist()
    plist = list(range(res.s.size))
    owner: dict[tuple[bool, int], int] = {}  # (high, position) -> first window covering it
    for i in idx:
        high = bool(res.s[i] > 0)
        st, w = int(res.start[i]), int(res.length[i])
        for p in range(st, st + w):
            key = (high, p % n if circular else p)
            j = owner.setdefault(key, i)
            if j != i:
                ri, rj = _find(plist, i), _find(plist, j)
                if ri != rj:
                    plist[max(ri, rj)] = min(ri, rj)
    groups: dict[int, list[int]] = {}
    for i in idx:
        groups.setdefault(_find(plist, i), []).append(i)
    out = []
    for members in groups.values():
        best = min(
            members, key=lambda i: (-abs(res.s[i]), int(res.length[i]), int(res.start[i]), i)
        )
        direction: Literal["high", "low"] = "high" if res.s[best] > 0 else "low"
        out.append(Cluster(direction=direction, best=best, members=tuple(sorted(members))))
    out.sort(key=lambda c: (int(res.start[c.best]), c.direction != "high", c.best))
    return out
