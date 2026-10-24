"""TU-block bootstrap, percentile CIs and the circular block bootstrap.

GIAP_SPEC.md §7.3.5(a)-(b), §8.2 ("resample TUs, not genes"), §8.3; decision 14. All randomness
comes from the ``np.random.Generator`` passed in by the caller (obtained via
``giap.rng.for_step``), as in ``giap.stats.permutation``.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.typing import ArrayLike

__all__ = ["BootCI", "circular_block_bootstrap", "percentile_ci", "tu_block_bootstrap"]


@dataclass(frozen=True, slots=True)
class BootCI:
    """Point estimate with a bootstrap percentile CI over ``n_valid`` of ``B`` replicates."""

    estimate: float
    ci_low: float
    ci_high: float
    n_valid: int
    B: int


def _check_level(level: float) -> None:
    if not 0.0 < level < 1.0:
        raise ValueError(f"level must be in (0, 1), got {level}")


def percentile_ci(boot: ArrayLike, level: float = 0.95) -> tuple[float, float]:
    """Percentile CI: linear quantiles at ``(1 - level)/2`` and ``(1 + level)/2``.

    Raises:
        ValueError: if ``level`` is not in (0, 1).
    """
    _check_level(level)
    b = np.asarray(boot, dtype=np.float64).ravel()
    if b.size == 0:
        return (math.nan, math.nan)
    lo, hi = np.quantile(b, [(1.0 - level) / 2.0, (1.0 + level) / 2.0], method="linear")
    return (float(lo), float(hi))


def tu_block_bootstrap(
    df: pd.DataFrame,
    tu_col: str,
    stat_fn: Callable[[pd.DataFrame], float],
    B: int,
    rng: np.random.Generator,
    level: float = 0.95,
    strata_col: str | None = None,
    *,
    columns: Sequence[str] | None = None,
) -> BootCI:
    """Bootstrap ``stat_fn`` by resampling whole TUs with replacement (within strata).

    ``estimate = stat_fn(df)``. Each replicate draws, per stratum (sorted order; one stratum when
    ``strata_col`` is None), as many TU IDs as the stratum has and concatenates their rows in
    their input order. ``stat_fn`` is called exactly ``B + 1`` times; its errors propagate. NaN
    replicate statistics are dropped; fewer than ``B/2`` valid gives CI (NaN, NaN). A TU's
    stratum is that of its first row.

    With ``columns``, ``stat_fn`` receives those columns as arrays instead of a data frame --
    ``stat_fn(*arrays)`` on every row for the estimate, on each replicate's rows (same draws,
    same rows) for the replicates -- so no frame is built per replicate.

    Raises:
        ValueError: on NA TU IDs, ``B < 1`` or ``level`` outside (0, 1).
    """
    _check_level(level)
    if B < 1:
        raise ValueError(f"B must be >= 1, got {B}")
    tu = df[tu_col]
    if bool(tu.isna().any()):
        raise ValueError(f"column {tu_col!r} has NA TU IDs")

    strata = None if strata_col is None else df[strata_col].to_numpy()
    boot = np.empty(B, dtype=np.float64)
    if columns is None:
        estimate = float(stat_fn(df))
        for b, rows in enumerate(_resampled_rows(tu, strata, B, rng)):
            boot[b] = float(stat_fn(df.iloc[rows].reset_index(drop=True)))
    else:
        cols = [df[c].to_numpy() for c in columns]
        estimate = float(stat_fn(*cols))
        for b, rows in enumerate(_resampled_rows(tu, strata, B, rng)):
            boot[b] = float(stat_fn(*(c[rows] for c in cols)))
    return _boot_ci(estimate, boot, B, level)


def _resampled_rows(
    tu: pd.Series, strata: npt.NDArray[Any] | None, B: int, rng: np.random.Generator
) -> Iterator[npt.NDArray[np.int64]]:
    """Row indices of each of the ``B`` replicates (draws made from ``rng`` up front)."""
    # TU codes in order of first appearance; rows grouped by TU, input order within a TU.
    codes, _ = pd.factorize(tu, sort=False)
    codes = np.asarray(codes, dtype=np.int64)
    n_tu = int(codes.max()) + 1 if codes.size else 0
    order = np.argsort(codes, kind="stable")
    lens = np.bincount(codes, minlength=n_tu).astype(np.int64)
    starts = np.zeros(n_tu, dtype=np.int64)
    if n_tu:
        starts[1:] = np.cumsum(lens)[:-1]

    first_row = np.unique(codes, return_index=True)[1]
    if strata is None:
        groups = [np.arange(n_tu, dtype=np.int64)]
    else:
        tu_stratum = strata[first_row]
        keys = sorted(pd.unique(tu_stratum).tolist())
        groups = [np.flatnonzero(tu_stratum == k) for k in keys]

    # Draw order fixed: strata in sorted order, each a (B, n_s) block of draws.
    draws = [s[rng.integers(0, s.size, size=(B, s.size))] for s in groups if s.size]
    picks = np.concatenate(draws, axis=1) if draws else np.zeros((B, 0), dtype=np.int64)

    for b in range(B):
        t = picks[b]
        ln = lens[t]
        total = int(ln.sum())
        offs = np.repeat(starts[t] - (np.cumsum(ln) - ln), ln) + np.arange(total)
        yield order[offs]


def _boot_ci(estimate: float, boot: npt.NDArray[np.float64], B: int, level: float) -> BootCI:
    valid = boot[~np.isnan(boot)]
    n_valid = int(valid.size)
    if n_valid == 0 or n_valid < B / 2:
        lo, hi = math.nan, math.nan
    else:
        lo, hi = percentile_ci(valid, level)
    return BootCI(estimate=estimate, ci_low=lo, ci_high=hi, n_valid=n_valid, B=B)


def circular_block_bootstrap(
    x: ArrayLike, block_len: int, B: int, rng: np.random.Generator
) -> npt.NDArray[np.float64]:
    """Politis-Romano circular block bootstrap of ``x``; shape (B, n).

    Per replicate, ``ceil(n / block_len)`` start positions uniform on [0, n); each block takes
    ``block_len`` consecutive values wrapping modulo n; blocks are concatenated and truncated to
    n. n = 0 gives shape (B, 0).

    Raises:
        ValueError: if ``block_len`` is not in [1, n] (n > 0) or ``B < 0``.
    """
    v = np.asarray(x, dtype=np.float64).ravel()
    n = v.size
    if B < 0:
        raise ValueError(f"B must be >= 0, got {B}")
    if n == 0:
        return np.zeros((B, 0), dtype=np.float64)
    if not 1 <= block_len <= n:
        raise ValueError(f"block_len must be in [1, {n}], got {block_len}")
    k = -(-n // block_len)
    st = np.asarray(rng.integers(0, n, size=(B, k)), dtype=np.int64)
    idx = (st[:, :, None] + np.arange(block_len, dtype=np.int64)) % n
    out: npt.NDArray[np.float64] = v[idx.reshape(B, k * block_len)[:, :n]]
    return out
