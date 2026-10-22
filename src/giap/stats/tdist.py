"""Student-t maximum likelihood and the two-sided p (GIAP_SPEC.md §7.4.2, §7.4.3, §8.3).

The within-operon distribution of the abundance difference Delta over OP pairs is modelled as a
Student-t with free location, scale and degrees of freedom, fitted by maximum likelihood with df
confined to ``df_bounds`` (the spec's [2, 50]). The optimiser is Nelder-Mead on
``(loc, log scale, df)`` with the df bound given to the optimiser, so a df sitting on a bound is a
legitimate optimum and is returned exactly at the bound with ``converged = True``.

Conventions: NaN (and +-inf) values are dropped before the fit; ``n`` counts the finite values; a
fit that cannot be made is returned as ``StudentT(nan, nan, nan, n, False)`` rather than raised, and
every downstream value computed from it is NaN.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy import optimize, stats

__all__ = ["StudentT", "fit_student_t", "student_t_two_sided"]

F64 = npt.NDArray[np.float64]

#: fewer finite values than this and the three-parameter fit is not attempted (spec, degenerate
#: cases table)
MIN_FIT_N = 5
_DF_START = 8.0
_SCALE_FLOOR = 1e-6  # zero MAD: the spec's starting scale
#: a fitted scale below this fraction of max(1, range of x) is a collapse onto tied values (the
#: likelihood has no maximum there), not an estimate
_COLLAPSE_FRAC = 1e-8
_MAXITER = 2000


@dataclass(frozen=True, slots=True)
class StudentT:
    """A fitted Student-t.

    Attributes:
        loc: location (NaN when the fit was not made).
        scale: scale, strictly positive (NaN when the fit was not made).
        df: degrees of freedom inside ``df_bounds`` (NaN when the fit was not made).
        n: number of finite values the fit used.
        converged: whether the optimiser stopped on its tolerances rather than on iterations.
    """

    loc: float
    scale: float
    df: float
    n: int
    converged: bool


def _finite(x: ArrayLike, name: str) -> F64:
    """1-D float copy of ``x`` with the non-finite values dropped."""
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {a.shape}")
    return a[np.isfinite(a)]


def _unestimated(n: int) -> StudentT:
    return StudentT(loc=math.nan, scale=math.nan, df=math.nan, n=n, converged=False)


def _neg_log_likelihood(theta: F64, x: F64) -> float:
    """Negative log-likelihood of ``x`` under t(loc, exp(log_scale), df); inf when degenerate."""
    loc, log_scale, df = float(theta[0]), float(theta[1]), float(theta[2])
    if not (math.isfinite(loc) and math.isfinite(log_scale) and df > 0.0):
        return math.inf
    scale = math.exp(log_scale) if log_scale > -700.0 else 0.0
    if scale <= 0.0:
        return math.inf
    with np.errstate(divide="ignore", invalid="ignore"):
        total = float(np.sum(stats.t.logpdf(x, df, loc=loc, scale=scale)))
    return -total if math.isfinite(total) else math.inf


def _snap_to_bound(theta: F64, x: F64, lo: float, hi: float) -> float:
    """The df of ``theta``, replaced by a bound when the bound does not fit worse.

    A simplex stops just short of a bound it is pressed against, which leaves a df like 49.9997
    where the optimum is the bound itself. Both bounds are therefore scored at the fitted location
    and scale and the best of the three degrees of freedom is kept.
    """
    free = _neg_log_likelihood(theta, x)
    best_df, best = float(theta[2]), free
    for bound in (lo, hi):
        value = _neg_log_likelihood(np.array([theta[0], theta[1], bound]), x)
        if value <= best:
            best_df, best = bound, value
    return best_df


def fit_student_t(x: ArrayLike, df_bounds: tuple[float, float] = (2.0, 50.0)) -> StudentT:
    """Maximum likelihood Student-t fit of ``x`` with df confined to ``df_bounds``.

    Nelder-Mead on ``(loc, log scale, df)`` from ``loc = median(x)``,
    ``scale = 1.4826 * MAD(x)`` (floored at 1e-6 when the MAD is zero) and ``df = 8``, with
    ``maxiter = 2000``. The df bound is enforced by the optimiser, so a fit that wants a lighter or
    heavier tail than the bounds allow is returned exactly at the bound and counted as converged.

    Args:
        x: values; NaN and +-inf are dropped.
        df_bounds: ``(low, high)`` for the degrees of freedom, both > 0 and low < high.

    Returns:
        The fit. ``StudentT(nan, nan, nan, n, False)`` when fewer than ``MIN_FIT_N`` finite values
        are given or the optimiser leaves the parameters degenerate (a zero, non-finite or collapsed
        scale: many tied values make the likelihood unbounded), never an exception.

    Raises:
        ValueError: if ``x`` is not one-dimensional or ``df_bounds`` is not a valid interval.
    """
    lo, hi = float(df_bounds[0]), float(df_bounds[1])
    if not (0.0 < lo < hi) or not math.isfinite(hi):
        raise ValueError(f"df_bounds must satisfy 0 < low < high < inf, got {df_bounds}")
    a = _finite(x, "x")
    n = int(a.size)
    if n < MIN_FIT_N:
        return _unestimated(n)
    loc0 = float(np.median(a))
    mad = float(np.median(np.abs(a - loc0)))
    scale0 = max(1.4826 * mad, _SCALE_FLOOR)
    start = np.array([loc0, math.log(scale0), min(max(_DF_START, lo), hi)], dtype=np.float64)
    result = optimize.minimize(
        _neg_log_likelihood,
        start,
        args=(a,),
        method="Nelder-Mead",
        bounds=[(None, None), (None, None), (lo, hi)],
        options={"maxiter": _MAXITER, "maxfev": 10 * _MAXITER},
    )
    loc, log_scale, df = (float(v) for v in result.x)
    df = _snap_to_bound(np.array([loc, log_scale, df]), a, lo, hi)
    scale = math.exp(log_scale) if -700.0 < log_scale < 700.0 else 0.0
    collapsed = scale < _COLLAPSE_FRAC * max(1.0, float(np.ptp(a)))
    if collapsed or not (scale > 0.0 and math.isfinite(loc) and lo <= df <= hi):
        return _unestimated(n)
    return StudentT(loc=loc, scale=scale, df=df, n=n, converged=bool(result.success))


def student_t_two_sided(x: ArrayLike, fit: StudentT) -> F64:
    """Two-sided p of each value of ``x`` under ``fit``: ``2 * min(F, 1 - F)``.

    ``F = F_t((x - loc) / scale, df)``, so the p is exactly 1.0 at ``x = loc`` and asymmetric
    around 0 whenever ``loc != 0``.

    Args:
        x: values to test.
        fit: the fitted t; an unestimated fit (NaN parameters) gives all-NaN output.

    Returns:
        The p-values in input order, in [0, 1], NaN where ``x`` is NaN.

    Raises:
        ValueError: if ``x`` is not one-dimensional.
    """
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"x must be one-dimensional, got shape {a.shape}")
    usable = math.isfinite(fit.loc) and fit.scale > 0.0 and fit.df > 0.0
    if not usable or a.size == 0:
        return np.full(a.shape, np.nan, dtype=np.float64)
    cdf: F64 = np.asarray(stats.t.cdf((a - fit.loc) / fit.scale, fit.df), dtype=np.float64).reshape(
        a.shape
    )
    return np.clip(2.0 * np.minimum(cdf, 1.0 - cdf), 0.0, 1.0)
