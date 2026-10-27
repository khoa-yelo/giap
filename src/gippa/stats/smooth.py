"""Plain local-linear LOESS detrending with circular padding (GIPPA_SPEC.md §7.3.5(c), (d) step 3).

Issue M5-01 (lead decisions 1-8). The fit at each point is a tricube-weighted local linear
regression over its k nearest neighbours by position, with no robustness iterations. It
reimplements, operation for operation, ``statsmodels.nonparametric.smoothers_lowess.lowess(y, x,
frac=k/N, it=0, delta=0)`` (statsmodels 0.15 ``_smoothers_lowess.pyx``): the same sliding
neighbourhood (ties keep the window on the left), radius, tricube weights, the "fewer than two
weights above 1e-12 -> fit = y" rule and the 1e-12 floor on the weighted variance of x.

The fit is linear in y and depends only on the positions, so it is stored as per-point weights
over each window; ``loess_matrix`` scatters them into the n x n smoother (ADR-0020 applies one
matrix to B permuted replicates) and ``loess_detrend`` applies them to one series.

Circular replicons are padded with the last k points at ``pos - period`` in front and the first k
points at ``pos + period`` behind; the padded copies map back to their original columns.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike

__all__ = ["LoessFit", "loess_detrend", "loess_matrix", "loess_span_points"]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]

_WEIGHT_EPS = 1e-12  # statsmodels: weights counted as non-zero above this
_VAR_FLOOR = 1e-12  # statsmodels: floor on the weighted variance of x


@dataclass(frozen=True, slots=True)
class LoessFit:
    """Result of :func:`loess_detrend`.

    Attributes
    ----------
    fit : numpy.ndarray of float64
        Fitted trend at each input point.
    resid : numpy.ndarray of float64
        ``y - fit``.
    k : int
        Neighbours used per local fit (n below ``min_points``; 0 when n = 0).
    """

    fit: F64
    resid: F64
    k: int


def loess_span_points(n: int, frac: float, min_points: int = 3) -> int:
    """Number of neighbours per local fit: min(n, max(min_points, ⌈frac·n - 1e-9⌉)).

    Parameters
    ----------
    n : int
        Number of points.
    frac : float
        Span as a fraction of the points.
    min_points : int, default 3
        Floor on the span.

    Returns
    -------
    int
        The span k, 0 <= k <= n.
    """
    return int(min(n, max(min_points, math.ceil(frac * n - 1e-9))))


def _validate_pos(pos: ArrayLike, frac: float, circular: bool, period: float | None) -> F64:
    p = np.asarray(pos, dtype=np.float64)
    if p.ndim != 1:
        raise ValueError("pos must be 1-D")
    if not 0.0 < frac <= 1.0:
        raise ValueError(f"frac must lie in (0, 1], got {frac}")
    if not np.all(np.isfinite(p)):
        raise ValueError("pos must be finite")
    if p.size > 1 and not np.all(np.diff(p) > 0):
        raise ValueError("pos must be strictly increasing")
    if circular:
        if period is None or not math.isfinite(period):
            raise ValueError("circular smoothing needs a finite period")
        if p.size and (p[0] < 0 or p[-1] >= period):
            raise ValueError("circular pos must lie in [0, period)")
    return p


def _window_weights(pos: F64, k: int, circular: bool, period: float | None) -> tuple[I64, F64]:
    """Window columns (into the original points) and projection weights, both (n, w)."""
    n = pos.size
    m = k if circular else 0
    cols = np.arange(n, dtype=np.int64)
    if circular:
        assert period is not None
        x = np.concatenate([pos[n - m :] - period, pos, pos[:m] + period])
        cols = np.concatenate([cols[n - m :], cols, cols[:m]])
    else:
        x = pos
    big_n = x.size
    w = min(big_n, max(2, k))  # statsmodels clamps its k to [2, N]
    xi = x[m : m + n]
    # statsmodels slides [L, L + w) right while xi > (x[L] + x[L + w]) / 2; the midpoints are
    # increasing, so the final L is the first one whose midpoint is >= xi.
    mid = (x[: big_n - w] + x[w:]) / 2.0
    left = np.searchsorted(mid, xi, side="left").astype(np.int64)
    idx = left[:, None] + np.arange(w, dtype=np.int64)[None, :]
    xw = x[idx]
    radius = np.maximum(xi - xw[:, 0], xw[:, -1] - xi)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.abs(xw - xi[:, None]) / radius[:, None]
        wt = (1.0 - d**3) ** 3
        wt = wt / wt.sum(axis=1, keepdims=True)
    ok = (wt > _WEIGHT_EPS).sum(axis=1) >= 2
    xbar = (wt * xw).sum(axis=1)
    sxx = np.maximum((wt * (xw - xbar[:, None]) ** 2).sum(axis=1), _VAR_FLOOR)
    proj = wt * (1.0 + (xi - xbar)[:, None] * (xw - xbar[:, None]) / sxx[:, None])
    # no regression (statsmodels' reg_ok False): the fit is the point itself
    self_col = np.arange(m, m + n, dtype=np.int64)[:, None] - left[:, None]
    identity = (np.arange(w)[None, :] == self_col).astype(np.float64)
    proj = np.where(ok[:, None], proj, identity)
    return cols[idx], proj


def _smoother(
    pos: ArrayLike, frac: float, circular: bool, period: float | None, min_points: int
) -> tuple[int, I64, F64]:
    p = _validate_pos(pos, frac, circular, period)
    n = p.size
    k = loess_span_points(n, frac, min_points)
    if n == 0:
        return 0, np.empty((0, 0), dtype=np.int64), np.empty((0, 0), dtype=np.float64)
    if n < min_points:  # local fits do not exist: the trend is the mean
        cols = np.broadcast_to(np.arange(n, dtype=np.int64), (n, n))
        return k, cols, np.full((n, n), 1.0 / n)
    cols, proj = _window_weights(p, k, circular, period)
    return k, cols, proj


def loess_matrix(
    pos: ArrayLike,
    frac: float,
    *,
    circular: bool,
    period: float | None = None,
    min_points: int = 3,
) -> F64:
    """The n x n LOESS smoother H with ``loess_detrend(y, pos, ...).fit == H @ y``.

    Parameters
    ----------
    pos : array_like, shape (n,)
        Strictly increasing finite positions (gene indices or TU ordinals).
    frac : float
        Span fraction in (0, 1].
    circular : bool
        Pad circularly with period ``period``.
    period : float, optional
        Replicon length in position units; required when ``circular``; pos must lie in
        [0, period).
    min_points : int, default 3
        Floor on the span; below it every row is the mean (1/n).

    Returns
    -------
    numpy.ndarray of float64, shape (n, n)

    Raises
    ------
    ValueError
        As :func:`loess_detrend` for ``pos``, ``frac`` and ``period``.
    """
    _, cols, proj = _smoother(pos, frac, circular, period, min_points)
    n = cols.shape[0]
    h = np.zeros((n, n), dtype=np.float64)
    rows = np.broadcast_to(np.arange(n)[:, None], cols.shape)
    np.add.at(h, (rows, cols), proj)
    return h


def loess_detrend(
    y: ArrayLike,
    pos: ArrayLike,
    frac: float,
    *,
    circular: bool,
    period: float | None = None,
    min_points: int = 3,
) -> LoessFit:
    """Detrend y with plain local-linear LOESS along ``pos`` (tricube weights, ``it = 0``).

    Parameters
    ----------
    y : array_like, shape (n,)
        Finite values.
    pos : array_like, shape (n,)
        Strictly increasing finite positions (gene indices or TU ordinals, gaps allowed).
    frac : float
        Span fraction in (0, 1]; k = :func:`loess_span_points` (n, frac, min_points).
    circular : bool
        Pad circularly: the last k points at ``pos - period`` in front, the first k at
        ``pos + period`` behind.
    period : float, optional
        Replicon length in position units; required when ``circular``.
    min_points : int, default 3
        For 1 <= n < min_points the fit is mean(y) everywhere and k = n.

    Returns
    -------
    LoessFit
        ``fit``, ``resid = y - fit`` and ``k``.

    Raises
    ------
    ValueError
        If y and pos differ in length, y is not finite, pos is not strictly increasing or not
        finite, frac is outside (0, 1], or circular without a period or with pos outside
        [0, period).

    Notes
    -----
    With k = 3 the fit equals y for any strictly increasing positions (the farthest neighbour has
    tricube weight 0), so every residual is 0; callers flag replicons with k < 4.
    """
    ya = np.asarray(y, dtype=np.float64)
    if ya.ndim != 1 or ya.shape != np.shape(pos):
        raise ValueError("y and pos must be 1-D arrays of equal length")
    if not np.all(np.isfinite(ya)):
        raise ValueError("y must be finite")
    k, cols, proj = _smoother(pos, frac, circular, period, min_points)
    fit = (proj * ya[cols]).sum(axis=1) if ya.size else np.empty(0, dtype=np.float64)
    fit = np.asarray(fit, dtype=np.float64)
    return LoessFit(fit=fit, resid=ya - fit, k=k)
