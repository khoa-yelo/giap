"""Student-t maximum likelihood and the two-sided p (GIAP_SPEC.md §7.4.2, §7.4.3, §8.3).

The within-operon distribution of the abundance difference Delta over OP pairs is modelled as a
Student-t with free location, scale and degrees of freedom, fitted by maximum likelihood with df
confined to ``df_bounds`` (the spec's [2, 50]). The optimiser is Nelder-Mead on
``(loc, log scale, df)`` with the df bound given to the optimiser, so a df sitting on a bound is a
legitimate optimum and is returned exactly at the bound with ``converged = True``.

Conventions: NaN (and +-inf) values are dropped before the fit; ``n`` counts the finite values; a
fit that cannot be made is returned as ``StudentT(nan, nan, nan, n, False)`` rather than raised, and
every downstream value computed from it is NaN. M6-13: :func:`fit_student_t_scaled`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy import optimize, stats

__all__ = ["DF_BOUNDS", "GAMMA_BOUNDS", "ScaledStudentT", "StudentT", "fit_student_t",
           "fit_student_t_scaled", "scale_at", "scaled_t_two_sided",
           "student_t_two_sided"]  # fmt: skip

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
#: bounds of the degrees of freedom nu (both fits)
DF_BOUNDS = (2.0, 50.0)
#: the covariate SD below which the scale is held constant (M6-13 S-5: a relevance guard)
MIN_ABUNDANCE_SD_LOG2 = 0.1
#: bounds of gamma, the change in log sigma per SD of the covariate (M6-13 S-4)
GAMMA_BOUNDS = (-2.0, 2.0)
_MAXITER_SCALED = 4000
_STEP = 0.1  # the explicit initial simplex (M6-13 S-6)


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


def fit_student_t(x: ArrayLike, df_bounds: tuple[float, float] = DF_BOUNDS) -> StudentT:
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


@dataclass(frozen=True, slots=True)
class ScaledStudentT:
    """x ~ loc + scale0 exp(gamma z) t_df, z = (covariate - center) / sd (M6-13).

    ``gamma`` is per SD of the covariate (0.0 for the constant null), NaN parameters when
    unestimated (``model = "none"``). ``reason`` says why a non-loglinear null is in force.
    z is clamped to [z_min, z_max] (its range over the fit pairs) when evaluated. ``lr_stat`` and
    ``lr_p`` are reported, never used.
    """

    loc: float
    scale0: float
    gamma: float
    df: float
    n: int
    converged: bool
    model: Literal["loglinear", "constant", "none"]
    reason: str | None
    center: float
    sd: float
    z_min: float
    z_max: float
    nll: float
    nll_constant: float
    lr_stat: float
    lr_p: float


def _scaled_nll(theta: F64, x: F64, z: F64) -> float:
    """Negative log-likelihood of x under mu + sigma0 exp(gamma z) t_nu; inf when degenerate."""
    loc, log_s0, gamma, df = (float(v) for v in theta)
    if not (math.isfinite(loc) and math.isfinite(log_s0) and math.isfinite(gamma) and df > 0.0):
        return math.inf
    log_sigma = log_s0 + gamma * z
    if not np.all(np.abs(log_sigma) < 700.0):
        return math.inf
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        total = float(np.sum(stats.t.logpdf((x - loc) / np.exp(log_sigma), df) - log_sigma))
    return -total if math.isfinite(total) else math.inf


def _snap(theta: F64, x: F64, z: F64, i: int, lo: float, hi: float) -> F64:
    """``theta`` with entry ``i`` replaced by a bound when the bound does not fit worse."""
    best, out = _scaled_nll(theta, x, z), theta
    for bound in (lo, hi):
        t = theta.copy()
        t[i] = bound
        value = _scaled_nll(t, x, z)
        if value <= best:
            best, out = value, t
    return out


def fit_student_t_scaled(
    x: ArrayLike,
    covariate: ArrayLike,
    *,
    df_bounds: tuple[float, float] = DF_BOUNDS,
    gamma_bounds: tuple[float, float] = GAMMA_BOUNDS,
    min_sd: float = MIN_ABUNDANCE_SD_LOG2,
    min_n: int = MIN_FIT_N,
) -> ScaledStudentT:
    """Maximum likelihood fit of :class:`ScaledStudentT` (docs/specs/M6-13 step 3).

    Pairs with a non-finite x or covariate are dropped together, the rest sorted on
    (x, covariate), so the fit is order-free. The constant :func:`fit_student_t` comes first:
    unestimated, not converged or fewer than ``min_n`` pairs is ``model = "none"``. A covariate
    SD below ``min_sd`` is the constant null (``flat_abundance``); a failed four-parameter fit is
    the constant null with ``not_converged``, ``collapsed`` or ``worse_than_constant``. Never
    raises on data; ValueError on non-1-D or mismatched inputs or invalid ``gamma_bounds``.
    """
    glo, ghi = float(gamma_bounds[0]), float(gamma_bounds[1])
    if not (math.isfinite(glo) and math.isfinite(ghi) and glo < ghi):
        raise ValueError(f"gamma_bounds must satisfy -inf < low < high < inf, got {gamma_bounds}")
    xa = np.array(x, dtype=np.float64, ndmin=1)
    ma = np.array(covariate, dtype=np.float64, ndmin=1)
    if xa.ndim != 1 or xa.shape != ma.shape:
        raise ValueError(f"x and covariate must be 1-D of one shape, got {xa.shape}, {ma.shape}")
    keep = np.isfinite(xa) & np.isfinite(ma)
    order = np.lexsort((ma[keep], xa[keep]))
    xa, ma = xa[keep][order], ma[keep][order]
    n = int(xa.size)
    center = float(np.mean(ma)) if n else math.nan
    sd = float(np.std(ma, ddof=1)) if n > 1 else math.nan
    z = (ma - center) / sd if sd > 0.0 else np.zeros_like(ma)
    zr = (float(z.min()), float(z.max())) if sd > 0.0 else (math.nan, math.nan)
    nan = math.nan
    const = fit_student_t(xa, df_bounds) if n >= min_n else _unestimated(n)
    if not const.converged:
        return ScaledStudentT(nan, nan, nan, nan, n, False, "none", None, center, sd, *zr,
                              nan, nan, nan, nan)  # fmt: skip
    nll_c = _neg_log_likelihood(np.array([const.loc, math.log(const.scale), const.df]), xa)

    def constant(reason: str) -> ScaledStudentT:
        return ScaledStudentT(const.loc, const.scale, 0.0, const.df, n, True, "constant", reason,
                              center, sd, *zr, nll_c, nll_c, nan, nan)  # fmt: skip

    if not sd >= min_sd:
        return constant("flat_abundance")
    lo, hi = float(df_bounds[0]), float(df_bounds[1])
    start = np.array([const.loc, math.log(const.scale), 0.0, const.df], dtype=np.float64)
    step = np.array([_STEP * const.scale, _STEP, _STEP, _STEP * const.df])
    if start[3] + step[3] > hi:
        step[3] = -step[3]
    result = optimize.minimize(
        _scaled_nll,
        start,
        args=(xa, z),
        method="Nelder-Mead",
        bounds=[(None, None), (None, None), (glo, ghi), (lo, hi)],
        options={"maxiter": _MAXITER_SCALED, "maxfev": 10 * _MAXITER_SCALED,
                 "initial_simplex": np.vstack([start, start + np.diag(step)])},
    )  # fmt: skip
    theta = _snap(
        _snap(np.asarray(result.x, dtype=np.float64), xa, z, 3, lo, hi), xa, z, 2, glo, ghi
    )
    if not (bool(result.success) and np.all(np.isfinite(theta))):
        return constant("not_converged")
    log_sigma = theta[1] + theta[2] * z
    floor = _COLLAPSE_FRAC * max(1.0, float(np.ptp(xa)))
    if float(np.min(log_sigma)) < math.log(floor):
        return constant("collapsed")
    nll = _scaled_nll(theta, xa, z)
    if not nll <= nll_c + 1e-9 * max(1.0, abs(nll_c)):
        return constant("worse_than_constant")
    lr = max(0.0, 2.0 * (nll_c - nll))
    loc, log_s0, gamma, df = (float(v) for v in theta)
    return ScaledStudentT(loc, math.exp(log_s0), gamma, df, n, True, "loglinear", None, center, sd,
                          *zr, nll, nll_c, lr, float(stats.chi2.sf(lr, 1)))  # fmt: skip


def scale_at(covariate: ArrayLike, fit: ScaledStudentT) -> F64:
    """sigma0 exp(gamma clamp(z)) at each raw covariate value, z clamped to [z_min, z_max].

    The constant null (gamma = 0) gives sigma0 at every finite covariate. NaN where the
    covariate is NaN or the fit is unestimated.
    """
    m = np.array(covariate, dtype=np.float64, ndmin=1)
    if m.ndim != 1:
        raise ValueError(f"covariate must be one-dimensional, got shape {m.shape}")
    if not (fit.scale0 > 0.0 and math.isfinite(fit.gamma)):
        return np.full(m.shape, np.nan, dtype=np.float64)
    if fit.gamma == 0.0:
        return np.where(np.isfinite(m), fit.scale0, np.nan)
    z = np.clip((m - fit.center) / fit.sd, fit.z_min, fit.z_max)
    return np.asarray(fit.scale0 * np.exp(fit.gamma * z), dtype=np.float64)


def scaled_t_two_sided(x: ArrayLike, covariate: ArrayLike, fit: ScaledStudentT) -> F64:
    """Two-sided p ``2 min(F, 1 - F)``, ``F = F_t((x - loc) / sigma_i; df)``, sigma_i from
    :func:`scale_at`. Input order kept; NaN where x or the covariate is NaN or the fit is
    unestimated.
    """
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"x must be one-dimensional, got shape {a.shape}")
    sigma = scale_at(covariate, fit)
    if sigma.shape != a.shape:
        raise ValueError(f"x and covariate must have one shape, got {a.shape}, {sigma.shape}")
    if not (math.isfinite(fit.loc) and fit.df > 0.0) or a.size == 0:
        return np.full(a.shape, np.nan, dtype=np.float64)
    cdf: F64 = np.asarray(stats.t.cdf((a - fit.loc) / sigma, fit.df), dtype=np.float64)
    return np.clip(2.0 * np.minimum(cdf, 1.0 - cdf), 0.0, 1.0)
