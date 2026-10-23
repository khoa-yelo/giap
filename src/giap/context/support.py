"""Plot support of the context_stats step: the region and island nulls (V18;
docs/specs/M9-17-plot-support-analysis.md, ADR-0043). Pure.

``null_regions``: per tested (region, contig, stat) in ``regions_stats`` row order, the used
shifts ascending with the studentised statistic the test compared (``RegionNull``).
``null_regions_grid``: lead decision 2's fallback, the 1st-99th (numpy linear) percentiles of each
null's finite values, third column ``percentile``. ``null_islands``: per scanned replicon,
direction ``abs``, its own max |S_w| per replicate and its largest observed |S_w|; the max over
replicons of a replicate is the genome-wide M_b (decision 5).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.context.islands import IslandScan
from giap.context.regions import RegionNull

__all__ = ["null_islands", "null_regions", "null_regions_grid"]

_KEY = ["region_id", "stat"]


def null_regions(table: pd.DataFrame, nulls: Sequence[RegionNull]) -> pd.DataFrame | None:
    """Exhaustive ``null_regions`` rows; None when no region was tested."""
    by = {(n.region_id, n.contig, n.stat): n for n in nulls}
    tested = table[(table["subset"] == "all") & (table["method"] != "descriptive")]
    keys = zip(tested["region_id"].astype(str), tested["contig"].astype(str),
               tested["stat"].astype(str), strict=True)  # fmt: skip
    parts = []
    for key in dict.fromkeys(k for k in keys if k in by):
        n = by[key]
        o = np.argsort(n.shifts, kind="stable")
        parts.append(pd.DataFrame({"region_id": n.region_id, "stat": n.stat,
                                   "shift_k": n.shifts[o], "null_value": n.null[o],
                                   "observed": n.observed}))  # fmt: skip
    return pd.concat(parts, ignore_index=True) if parts else None


def null_regions_grid(exhaustive: pd.DataFrame) -> pd.DataFrame:
    """The 1st-99th percentile grid of each (region, stat) null of ``null_regions``."""
    pct = np.arange(1, 100)
    parts = []
    for (rid, stat), grp in exhaustive.groupby(_KEY, sort=False):
        v = grp["null_value"].to_numpy(float)
        parts.append(pd.DataFrame({"region_id": rid, "stat": stat, "percentile": pct,
                                   "null_value": np.percentile(v[~np.isnan(v)], pct),
                                   "observed": grp["observed"].iloc[0]}))  # fmt: skip
    return pd.concat(parts, ignore_index=True)


def null_islands(scan: IslandScan, nulls: Sequence[npt.NDArray[np.float64]]) -> pd.DataFrame | None:
    """``null_islands`` rows of the scanned replicons (``nulls``: ``scan_with_nulls``'s
    per-replicon maxima); None when none was scanned."""
    parts = []
    for ts, ok, win, per in zip(scan.series, scan.scanned, scan.windows, nulls, strict=True):
        if not ok:
            continue
        s = np.abs(np.asarray(win.s, dtype=float))
        obs = float(np.nanmax(s)) if s.size else np.nan
        parts.append(pd.DataFrame({"contig": ts.contig, "direction": "abs",
                                   "replicate": np.arange(1, per.size + 1), "max_scan_s": per,
                                   "observed_max": obs}))  # fmt: skip
    return pd.concat(parts, ignore_index=True) if parts else None
