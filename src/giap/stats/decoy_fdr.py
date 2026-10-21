"""Distance-reweighted decoy FDR, monotone q and local fdr (GIAP_SPEC.md §7.4.4, §8.2, §8.3).

Targets are candidate missed operons and decoys are gene pairs that mimic them (convergent and
sandwich pairs, §8.2). Scores are ``s = -|Delta|``, higher is better as everywhere in
``giap.stats``, and each pair carries a distance bin from ``bin_index`` so that the decoys can be
reweighted to the distance mix of the targets:

    FDR(t) = pi0 * sum_b N_NO,b * (#{decoy in b : s >= t} + 1) / (N_decoy,b + 1)
             / #{target : s >= t}

``N_NO,b`` counts **all** targets in bin b, not only those above t (M6-01 spec step 4), so FDR(t)
is an expected count of false targets over the observed count of calls. A bin with targets and no
decoys contributes its whole target count, which is the maximally conservative value and is kept
unclipped in ``fdr``; ``q`` clips to [0, 1].

Lead decisions of docs/specs/M6-01-lead-decisions.md that this module implements:

1. ``max_decoy_fdr`` takes the element-wise maximum over the decoy sets **first** and monotonises
   afterwards (monotone-of-max >= max-of-monotones, and only the former is an upper bound on the
   FDP of the list it selects).
2. "Monotone" is ``giap.stats.tda.monotone_q``: the running minimum over the looser thresholds.
3. ``DecoyFdr.lfdr`` follows the same "larger wins" rule, the per-bin weights of step 4 feeding the
   null density of each decoy set.
4. At the loosest threshold every decoy set gives FDR = pi0 exactly (the +1 terms cancel), so
   ``which`` is tied there: the first key of ``decoys`` wins.
5. A decoy set with no usable decoy is skipped in ``max_decoy_fdr`` only; ``reweighted_fdr``
   applies the formula as written and returns the conservative value.
7. The KDE bandwidth is scipy's ``gaussian_kde(bw_method="silverman")``, taken from the
   **unweighted** sample.
8. ``_reflected_density`` and ``_lfdr_densities`` expose ``f0`` and ``f_NO`` so that the reflection
   at zero can be tested directly: it cancels in the lfdr ratio and no test of the public return
   value can see it.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy import stats

from giap.stats.tda import monotone_q

__all__ = ["DecoyFdr", "bin_index", "local_fdr", "max_decoy_fdr", "reweighted_fdr"]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
Decoys = Mapping[str, tuple[ArrayLike, ArrayLike]]

#: fewer values than this in either sample and the local fdr is undefined (spec, degenerate cases)
MIN_KDE_N = 5
#: grid steps per bandwidth for the binned KDE evaluation (relative error ~1e-4)
_KDE_STEPS_PER_BW = 40
#: the Gaussian kernel is truncated at this many bandwidths (exp(-50) = 2e-22)
_KDE_TRUNCATION = 10.0
_KDE_MAX_GRID = 1 << 21


@dataclass(frozen=True, slots=True)
class DecoyFdr:
    """The decoy FDR of a target list, one entry per target in input order.

    Attributes:
        q: monotone q-value in [0, 1]; NaN for a target without a score or a bin.
        fdr: the raw step function at the target's score, unclipped (it can exceed 1).
        which: the decoy set that gave the larger FDR at that score ("" where none survives).
        lfdr: local fdr in [0, 1]; NaN where the KDE is not defined.
    """

    q: F64
    fdr: F64
    which: npt.NDArray[np.str_]
    lfdr: F64


def _scores(x: ArrayLike, name: str) -> F64:
    """1-D float copy of ``x``."""
    a = np.array(x, dtype=np.float64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {a.shape}")
    return a


def _bins(x: ArrayLike, name: str) -> I64:
    """1-D int64 copy of ``x``."""
    a = np.array(x, dtype=np.int64, ndmin=1)
    if a.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {a.shape}")
    return a


def _check_pi0(pi0: float) -> None:
    if not 0.0 < float(pi0) <= 1.0:
        raise ValueError(f"pi0 must lie in (0, 1], got {pi0}")


def _check_pair(s: F64, b: I64, what: str) -> None:
    if s.shape != b.shape:
        raise ValueError(f"{what} scores and bins differ in length: {s.size} != {b.size}")


def _usable(s: F64, b: I64, n_bins: int) -> npt.NDArray[np.bool_]:
    """Entries with a score and a bin in ``range(n_bins)``; everything else is excluded."""
    return ~np.isnan(s) & (b >= 0) & (b < n_bins)


def bin_index(distance_bp: ArrayLike, edges: Sequence[float]) -> I64:
    """Half-open distance bins ``[edges[i], edges[i + 1])``, 0-based.

    Args:
        distance_bp: intergenic distances; NaN gives bin -1.
        edges: ``len(edges) - 1`` ascending bin edges, e.g. ``operons.distance_bins_bp``
            ``[-inf, 0, 50, 150, 300, inf]``.

    Returns:
        The bin of each distance, -1 where the distance is NaN or outside every bin.

    Raises:
        ValueError: if fewer than two edges are given or they are not ascending.
    """
    e = np.asarray(edges, dtype=np.float64)
    if e.ndim != 1 or e.size < 2:
        raise ValueError(f"edges must hold at least two values, got shape {e.shape}")
    if np.any(np.diff(e) <= 0.0):
        raise ValueError("edges must be strictly ascending")
    d = _scores(distance_bp, "distance_bp")
    idx = np.searchsorted(e, d, side="right").astype(np.int64) - 1
    outside = np.isnan(d) | (idx < 0) | (idx > e.size - 2)
    return np.where(outside, np.int64(-1), idx)


def reweighted_fdr(
    target_s: ArrayLike,
    target_bin: ArrayLike,
    decoy_s: ArrayLike,
    decoy_bin: ArrayLike,
    n_bins: int,
    pi0: float = 1.0,
) -> F64:
    """The distance-reweighted decoy FDR at each target's score, unclipped.

    ``FDR(t) = pi0 * sum_b N_NO,b * (D_b(t) + 1) / (N_decoy,b + 1) / T(t)`` with ``N_NO,b`` the
    number of targets in bin b, ``D_b(t)`` the decoys of bin b scoring at or above t, and ``T(t)``
    the targets scoring at or above t. Bins without targets contribute nothing; bins without
    decoys contribute ``N_NO,b``, so the result may exceed 1 (spec example 7). An empty decoy set
    is not special-cased here (lead decision 5): every bin then contributes its target count.

    Args:
        target_s: target scores ``-|Delta|``; NaN is excluded.
        target_bin: target distance bins; anything outside ``range(n_bins)`` is excluded.
        decoy_s: decoy scores; NaN is excluded.
        decoy_bin: decoy distance bins; anything outside ``range(n_bins)`` is excluded.
        n_bins: number of distance bins.
        pi0: prior null proportion, in (0, 1].

    Returns:
        The FDR at each target's score, in input order, NaN for an excluded target. Tied scores
        share one value.

    Raises:
        ValueError: if ``pi0`` is outside (0, 1], ``n_bins`` is not positive, or a score and bin
            array differ in length.
    """
    _check_pi0(pi0)
    if n_bins < 1:
        raise ValueError(f"n_bins must be positive, got {n_bins}")
    t_s, t_b = _scores(target_s, "target_s"), _bins(target_bin, "target_bin")
    d_s, d_b = _scores(decoy_s, "decoy_s"), _bins(decoy_bin, "decoy_bin")
    _check_pair(t_s, t_b, "target")
    _check_pair(d_s, d_b, "decoy")
    out = np.full(t_s.shape, np.nan, dtype=np.float64)
    ok_t = _usable(t_s, t_b, n_bins)
    if not ok_t.any():
        return out
    s, b = t_s[ok_t], t_b[ok_t]
    ok_d = _usable(d_s, d_b, n_bins)
    n_no = np.bincount(b, minlength=n_bins).astype(np.float64)
    n_decoy = np.bincount(d_b[ok_d], minlength=n_bins).astype(np.float64)
    n_target_above = s.size - np.searchsorted(np.sort(s), s, side="left")
    expected = np.zeros(s.size, dtype=np.float64)
    for bin_ in range(n_bins):
        if n_no[bin_] == 0.0:
            continue
        in_bin = np.sort(d_s[ok_d & (d_b == bin_)])
        above = in_bin.size - np.searchsorted(in_bin, s, side="left")
        expected += n_no[bin_] * (above + 1.0) / (n_decoy[bin_] + 1.0)
    out[ok_t] = float(pi0) * expected / n_target_above
    return out


def _decoy_weights(
    t_b: I64, ok_t: npt.NDArray[np.bool_], d_b: I64, ok_d: npt.NDArray[np.bool_], n_bins: int
) -> F64:
    """Per-decoy weights ``N_NO,b / N_decoy,b`` of the usable decoys (step 4's reweighting)."""
    n_no = np.bincount(t_b[ok_t], minlength=n_bins).astype(np.float64)
    n_decoy = np.bincount(d_b[ok_d], minlength=n_bins).astype(np.float64)
    ratio = np.divide(n_no, n_decoy, out=np.zeros(n_bins), where=n_decoy > 0.0)
    return ratio[d_b[ok_d]]


def max_decoy_fdr(
    target_s: ArrayLike,
    target_bin: ArrayLike,
    decoys: Decoys,
    n_bins: int,
    pi0: float = 1.0,
) -> DecoyFdr:
    """The larger decoy FDR over ``decoys``, its monotone q and the local fdr.

    ``reweighted_fdr`` is run once per decoy set; the element-wise maximum is taken first and
    monotonised afterwards with ``giap.stats.tda.monotone_q`` (lead decisions 1 and 2), then
    clipped to [0, 1] to give ``q``. ``lfdr`` is the element-wise maximum of ``local_fdr`` over the
    same sets, each with the per-bin weights ``N_NO,b / N_decoy,b`` (lead decision 3); a set whose
    KDE is undefined does not mask one whose KDE is defined.

    At the loosest threshold every set gives exactly ``pi0``, so ``which`` is tied; the tie is
    broken by the order of ``decoys``, the first key winning (lead decision 4).

    Args:
        target_s: target scores ``-|Delta|``; NaN is excluded.
        target_bin: target distance bins; anything outside ``range(n_bins)`` is excluded.
        decoys: decoy sets by name, each a ``(scores, bins)`` pair. A set with no usable decoy is
            skipped (lead decision 5); if no set survives, ``q``, ``fdr`` and ``lfdr`` are NaN.
        n_bins: number of distance bins.
        pi0: prior null proportion, in (0, 1].

    Returns:
        The :class:`DecoyFdr` of the target list, in input order.

    Raises:
        ValueError: if ``pi0`` is outside (0, 1], ``n_bins`` is not positive, or an array pair
            differs in length.
    """
    _check_pi0(pi0)
    if n_bins < 1:
        raise ValueError(f"n_bins must be positive, got {n_bins}")
    t_s, t_b = _scores(target_s, "target_s"), _bins(target_bin, "target_bin")
    _check_pair(t_s, t_b, "target")
    ok_t = _usable(t_s, t_b, n_bins)
    names: list[str] = []
    curves: list[F64] = []
    lfdrs: list[F64] = []
    for name, (d_score, d_bin) in decoys.items():
        d_s, d_b = _scores(d_score, f"{name} decoy_s"), _bins(d_bin, f"{name} decoy_bin")
        _check_pair(d_s, d_b, f"{name} decoy")
        ok_d = _usable(d_s, d_b, n_bins)
        if not ok_d.any():
            continue
        names.append(name)
        curves.append(reweighted_fdr(t_s, t_b, d_s, d_b, n_bins, pi0=pi0))
        weights = _decoy_weights(t_b, ok_t, d_b, ok_d, n_bins)
        one = np.full(t_s.shape, np.nan, dtype=np.float64)
        one[ok_t] = local_fdr(t_s[ok_t], d_s[ok_d], weights, pi0=pi0)
        lfdrs.append(one)
    if not names:
        nan = np.full(t_s.shape, np.nan, dtype=np.float64)
        return DecoyFdr(q=nan, fdr=nan.copy(), which=np.full(t_s.shape, "", dtype="<U1"), lfdr=nan)
    stacked = np.vstack(curves)
    fdr: F64 = stacked.max(axis=0)
    # argmax keeps the first maximum, so the first key of ``decoys`` wins a tie (lead decision 4)
    which = np.where(ok_t, np.asarray(names)[stacked.argmax(axis=0)], "")
    lfdr: F64 = lfdrs[0].copy()
    for other in lfdrs[1:]:
        lfdr = np.fmax(lfdr, other)  # an undefined KDE does not mask a defined one
    q = np.full(t_s.shape, np.nan, dtype=np.float64)
    if ok_t.any():
        order = np.argsort(-t_s[ok_t], kind="stable")  # strictest threshold first
        q_ok = np.empty(int(ok_t.sum()), dtype=np.float64)
        q_ok[order] = monotone_q(fdr[ok_t][order])
        q[ok_t] = q_ok
    return DecoyFdr(q=q, fdr=fdr, which=which, lfdr=lfdr)


def _silverman_bandwidth(sample: F64) -> float:
    """Silverman's bandwidth of the unweighted ``sample``, from ``scipy.stats.gaussian_kde``.

    Lead decision 7: the rule is scipy's ``bw_method="silverman"``, i.e.
    ``(3 n / 4) ** (-1/5) * sd(sample)`` in one dimension, and it is taken from the unweighted
    sample so that the decoy weights move the density and not the smoothing.
    """
    try:
        kde = stats.gaussian_kde(sample, bw_method="silverman")
    except (ValueError, np.linalg.LinAlgError):  # a degenerate (zero-variance) sample
        return math.nan
    h = float(np.sqrt(kde.covariance[0, 0]))
    return h if math.isfinite(h) and h > 0.0 else math.nan


def _reflected_density(sample: F64, weights: F64, x: F64) -> F64:
    """Gaussian KDE of ``sample`` on [0, inf), reflected at zero, evaluated at ``|x|``.

    ``f(x) = sum_j w_j [phi_h(x - sample_j) + phi_h(x + sample_j)]`` with ``w`` normalised to sum
    to 1 and ``h`` from :func:`_silverman_bandwidth`; the reflection makes ``f`` a density on
    [0, inf) instead of leaking mass below the hard boundary at ``|Delta| = 0``, and it is the
    reason this helper is exposed: the reflection cancels in the lfdr ratio (lead decision 8).

    The sum is evaluated by linear binning on a uniform grid of ``h / 40`` and one convolution with
    the kernel truncated at 10 bandwidths, then linearly interpolated at ``|x|`` - exact to about
    1e-4 relative and O(n + G K) instead of O(n m), which the 50,000-target budget needs.

    Args:
        sample: the values the density is built on (``|Delta|``, non-negative).
        weights: non-negative per-value weights; only their ratios matter.
        x: evaluation points; ``f`` is symmetric, so ``x`` and ``-x`` give the same value.

    Returns:
        The density at each point of ``x``, all-NaN when the bandwidth is undefined.
    """
    out = np.full(x.shape, np.nan, dtype=np.float64)
    total = float(weights.sum())
    h = _silverman_bandwidth(sample)
    if not math.isfinite(h) or total <= 0.0:
        return out
    w = weights / total
    eval_at = np.abs(x)
    span = float(max(np.abs(sample).max(), eval_at.max(), h))
    delta = h / _KDE_STEPS_PER_BW
    half = min(math.ceil(span / delta) + 1, _KDE_MAX_GRID)
    delta = span / half if half == _KDE_MAX_GRID else delta
    grid = np.arange(-half, half + 1, dtype=np.float64) * delta  # symmetric about 0
    position = (sample + half * delta) / delta
    left = np.clip(np.floor(position), 0.0, 2.0 * half - 1.0)
    frac = np.clip(position - left, 0.0, 1.0)
    index = left.astype(np.int64)
    binned = np.bincount(index, weights=w * (1.0 - frac), minlength=grid.size)
    binned += np.bincount(index + 1, weights=w * frac, minlength=grid.size)
    binned = binned[: grid.size] + binned[: grid.size][::-1]  # the reflection at zero
    reach = max(math.ceil(_KDE_TRUNCATION * _KDE_STEPS_PER_BW), 1)
    offsets = np.arange(-reach, reach + 1, dtype=np.float64) * delta
    kernel = np.exp(-0.5 * (offsets / h) ** 2) / (h * math.sqrt(2.0 * math.pi))
    density = np.convolve(binned, kernel)[reach : reach + grid.size]  # centred on each node
    return np.interp(eval_at, grid, density)


def _lfdr_densities(target_s: ArrayLike, decoy_s: ArrayLike, decoy_w: ArrayLike) -> tuple[F64, F64]:
    """``(f0, f_NO)`` at the target scores: the weighted decoy and unweighted target densities.

    Exposed for testing the reflection at zero (lead decision 8), which cancels in the ratio
    :func:`local_fdr` returns.
    """
    t_s = _scores(target_s, "target_s")
    d_s, d_w = _scores(decoy_s, "decoy_s"), _scores(decoy_w, "decoy_w")
    if d_s.shape != d_w.shape:
        raise ValueError(f"decoy_s and decoy_w differ in length: {d_s.size} != {d_w.size}")
    nan = np.full(t_s.shape, np.nan, dtype=np.float64)
    ok_t = ~np.isnan(t_s)
    ok_d = ~np.isnan(d_s) & np.isfinite(d_w) & (d_w > 0.0)
    if int(ok_t.sum()) < MIN_KDE_N or int(ok_d.sum()) < MIN_KDE_N:
        return nan, nan.copy()
    target_abs, decoy_abs = -t_s[ok_t], -d_s[ok_d]
    f0, f_no = nan, nan.copy()
    f0[ok_t] = _reflected_density(decoy_abs, d_w[ok_d], target_abs)
    f_no[ok_t] = _reflected_density(target_abs, np.ones_like(target_abs), target_abs)
    return f0, f_no


def local_fdr(target_s: ArrayLike, decoy_s: ArrayLike, decoy_w: ArrayLike, pi0: float = 1.0) -> F64:
    """Local fdr ``min(1, pi0 * f0(s) / f_NO(s))`` at each target score.

    ``f0`` is the Gaussian KDE of the decoy ``|Delta|`` with the per-decoy weights ``decoy_w``
    (the distance reweighting of :func:`reweighted_fdr`, normalised) and ``f_NO`` the unweighted
    Gaussian KDE of the target ``|Delta|``; both are reflected at zero and both use Silverman's
    bandwidth from their own unweighted sample (lead decision 7). Only the ratios of ``decoy_w``
    matter.

    Args:
        target_s: target scores ``-|Delta|``; NaN is excluded.
        decoy_s: decoy scores ``-|Delta|``; NaN is excluded.
        decoy_w: per-decoy weights, non-negative; zero-weight decoys are dropped.
        pi0: prior null proportion, in (0, 1].

    Returns:
        The local fdr in [0, 1] at each target, in input order. All NaN when either sample has
        fewer than ``MIN_KDE_N`` usable values or its bandwidth is undefined, and NaN at a target
        where ``f_NO`` is zero.

    Raises:
        ValueError: if ``pi0`` is outside (0, 1] or ``decoy_s`` and ``decoy_w`` differ in length.
    """
    _check_pi0(pi0)
    f0, f_no = _lfdr_densities(target_s, decoy_s, decoy_w)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(f_no > 0.0, float(pi0) * f0 / f_no, np.nan)
    return np.clip(ratio, 0.0, 1.0)
