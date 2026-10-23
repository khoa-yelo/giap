"""Pairwise Spearman rho with a TU-block bootstrap CI and Bland-Altman limits (M9-22;
GIAP_SPEC.md §7.9.3, §8.2, §8.3). Pure: randomness only from the ``rng`` passed in
(``giap.rng.for_step("compare", seed)``); see docs/methods/compare.md.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from itertools import combinations

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import rankdata

from giap.stats.bootstrap import tu_block_bootstrap

__all__ = ["COLUMNS", "METHOD", "abundance_pairs", "spearman"]

COLUMNS = (
    "run_a", "run_b", "n_shared", "spearman_rho", "ci_low", "ci_high", "ba_bias_log2",
    "ba_loa_low", "ba_loa_high", "method", "n_boot", "flags",
)  # fmt: skip
METHOD = "spearman_tu_bootstrap"
MIN_SHARED = 3
MIN_UNITS = 2
LOA_Z = 1.96
_STATS = ("spearman_rho", "ci_low", "ci_high", "ba_bias_log2", "ba_loa_low", "ba_loa_high")


def spearman(x: npt.ArrayLike, y: npt.ArrayLike) -> float:
    """Spearman rho with average ranks for ties; NaN when either side is constant."""
    rx = rankdata(np.asarray(x, dtype=float)) if np.size(x) else np.zeros(0)
    ry = rankdata(np.asarray(y, dtype=float)) if np.size(y) else np.zeros(0)
    dx, dy = rx - rx.mean(), ry - ry.mean()
    den = math.sqrt(float(dx @ dx) * float(dy @ dy))
    if den == 0.0:
        return math.nan
    return float(np.clip(float(dx @ dy) / den, -1.0, 1.0))


def _measured(frame: pd.DataFrame) -> pd.Series:
    """Detected, not imputed, with a log2 abundance; indexed by locus tag."""
    ok = (
        frame["detected"].fillna(False).astype(bool)
        & ~frame["abundance_imputed"].fillna(False).astype(bool)
        & frame["abundance_log2"].notna()
    )
    return pd.Series(ok.to_numpy(), index=frame["locus_tag"].to_numpy())


def _has_abundance(frame: pd.DataFrame) -> bool:
    kind = frame["abundance_type"].astype(object)
    return bool((frame["abundance_log2"].notna() & kind.notna() & (kind != "none")).any())


def _scale(frame: pd.DataFrame) -> tuple[frozenset[str], frozenset[str]]:
    kinds = frame["abundance_type"].dropna().astype(str)
    return frozenset(kinds[kinds != "none"]), frozenset(frame["abundance_scale"].dropna())


def _rho_frame(d: pd.DataFrame) -> float:
    return spearman(d["a"].to_numpy(), d["b"].to_numpy())


def _pair(
    la: str, lb: str, fa: pd.DataFrame, fb: pd.DataFrame, tus: pd.Series, B: int,
    rng: np.random.Generator,
) -> dict[str, object]:  # fmt: skip
    row: dict[str, object] = dict.fromkeys(_STATS, math.nan)
    row.update(run_a=la, run_b=lb, method=METHOD, n_boot=B)
    flags = [f"missing_in:{lab}" for lab, f in ((la, fa), (lb, fb)) if not _has_abundance(f)]
    ok_b = _measured(fb)
    ok = _measured(fa) & ok_b.reindex(fa["locus_tag"].to_numpy(), fill_value=False).to_numpy()
    tags = fa["locus_tag"].to_numpy()[ok.to_numpy()]
    row["n_shared"] = len(tags)
    if flags:
        row["flags"] = ";".join(flags)
        return row
    if _scale(fa) != _scale(fb):
        flags.append("scale_differs")
    if len(tags) < MIN_SHARED:
        row["flags"] = ";".join(["too_few_shared", *flags])
        return row
    a = fa.set_index("locus_tag").loc[tags, "abundance_log2"].to_numpy(dtype=float)
    b = fb.set_index("locus_tag").loc[tags, "abundance_log2"].to_numpy(dtype=float)
    d = a - b
    bias, half = float(d.mean()), LOA_Z * float(d.std(ddof=1))
    row.update(ba_bias_log2=bias, ba_loa_low=bias - half, ba_loa_high=bias + half)
    units = tus.reindex(tags).astype(object)
    units = units.where(units.notna(), pd.Series(tags, index=units.index, dtype=object))
    df = pd.DataFrame({"tu": units.to_numpy(), "a": a, "b": b})
    row["spearman_rho"] = _rho_frame(df)
    if df["tu"].nunique() < MIN_UNITS:
        flags.append("ci_few_units")
    else:
        ci = tu_block_bootstrap(df, "tu", _rho_frame, B=B, rng=rng)
        row.update(ci_low=ci.ci_low, ci_high=ci.ci_high)
        if not (math.isfinite(ci.ci_low) and math.isfinite(ci.ci_high)):
            flags.append("ci_unestimated")
    row["flags"] = ";".join(flags) if flags else None
    return row


def abundance_pairs(
    runs: Mapping[str, pd.DataFrame], tus: pd.Series, B: int, rng: np.random.Generator
) -> pd.DataFrame:
    """One row per unordered pair (label order) of the ``compare_abundance`` columns.

    Shared genes are detected and measured (``abundance_imputed`` false) in both runs, in the
    first run's row order; ``tus`` (``tu_id`` by locus tag; a gene without one is its own unit)
    are the resampling units of the percentile CI (``tu_block_bootstrap``, B resamples).
    Bland-Altman on log2: bias = mean(a - b), limits = bias +- 1.96 SD (ddof = 1).
    """
    rows = [_pair(la, lb, runs[la], runs[lb], tus, B, rng) for la, lb in combinations(runs, 2)]
    out = pd.DataFrame(rows, columns=list(COLUMNS))
    return out.astype({"n_shared": "Int64", "n_boot": "Int64"})
