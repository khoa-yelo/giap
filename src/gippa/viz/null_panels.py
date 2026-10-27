"""V18 panel (c): the null distributions against the observed statistic (M9-21; GIPPA_SPEC.md
§10.2 V18, §7.3.5(b, d), §7.8.2 step 4; ADR-0018, ADR-0034).

Pure: the M9-17 support files (``null_regions``, ``null_islands``, ``null_determinants``) and the
tables come from a ``RunView``; nothing is computed beyond binning. Each small multiple is a
``Multiple``: 30 equal bins from the smallest finite null value (the last bin closed), or, for
the ``percentile`` grid form of ``null_regions``, a step CDF. The observed statistic is the one
the support file stores, i.e. the statistic the test compared (regions: studentised O - E,
ADR-0018, not the reported raw O - E).
Lead decision 4: regions are keyed on (region_id, contig); the islands block opens with the
genome-wide maximum-scan null that the FWER p uses.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from gippa.viz.style import OKABE_ITO

if TYPE_CHECKING:
    from gippa.io.run_view import RunView

N_BINS = 30
MAX_REGIONS = 12
MAX_REPLICONS = 7
REGION_STAT = "o_minus_e"
XLABEL = {
    "regions_stats": "test statistic (studentised O \u2212 E)",
    "islands": "max |S_w| over windows",
    "determinants_model": "mean CV R\u00b2 (null: permuted response)",
}
GRID_TEXT = "step CDF of the percentile grid"
GENOME = "genome / abs"
PER_REPLICON = "per-replicon maximum"
NOT_DRAWN = "n_not_drawn"
UNCALIBRATED_NOTE = "(uncalibrated, ADR-0034)"
TOO_FEW = "too_few_rows"


@dataclass(frozen=True)
class Multiple:
    """One small multiple: ``values`` are the finite null draws (histogram) or the grid values
    (``grid`` holds the percentiles / 100, the step CDF); ``extra`` are ``(metric, value)`` data
    rows written after the ``observed`` row."""

    table: str
    family: str
    values: tuple[float, ...]
    observed: float
    caption: str
    grid: tuple[float, ...] | None = None
    extra: tuple[tuple[str, float], ...] = ()


def _num(col: pd.Series) -> np.ndarray:
    out: np.ndarray = pd.to_numeric(col, errors="coerce").to_numpy(dtype=float)
    return out


def _fmt(v: float) -> str:
    return f"{v:.2g}" if math.isfinite(v) else "not available"


def _caption(text: str, observed: float) -> str:
    return text if math.isfinite(observed) else f"{text}\nobserved NA"


def _flags(cell: object) -> set[str]:
    return set() if pd.isna(cell) else set(str(cell).split(";"))


def histogram(values: Sequence[float]) -> tuple[list[float], list[int]]:
    """30 ascending left edges from ``min(values)`` and the counts ``[e_k, e_k+1)``, the last bin
    taking the rest; the width is ``(max - min) / 30``, or ``1 / 30`` for a constant null."""
    arr = np.asarray(values, dtype=float)
    lo, hi = float(arr.min()), float(arr.max())
    width = (hi - lo) / N_BINS if hi > lo else 1.0 / N_BINS
    edges = lo + width * np.arange(N_BINS, dtype=float)
    idx = np.clip(np.searchsorted(edges, arr, side="right") - 1, 0, N_BINS - 1)
    return [float(e) for e in edges], np.bincount(idx, minlength=N_BINS).astype(int).tolist()


def _finite(col: pd.Series) -> tuple[float, ...]:
    v = _num(col)
    return tuple(float(x) for x in v[np.isfinite(v)])


def _first(col: pd.Series) -> float:
    v = _num(col)
    return float(v[0]) if len(v) else math.nan


def _regions(view: RunView) -> list[Multiple]:
    """Up to 12 (region_id, contig) keys with the smallest q; a region on more than one contig
    gets one multiple per contig, family ``"<rid> (<contig>) / o_minus_e"`` (lead decision 4)."""
    nulls, stats = view.support("null_regions"), view.table("regions_stats")
    if nulls is None or stats is None or nulls.empty:
        return []
    nulls = nulls[nulls["stat"] == REGION_STAT]
    s = stats[(stats["stat"] == REGION_STAT) & (stats["subset"] == "all")]
    keys = list(zip(nulls["region_id"].astype(str), nulls["contig"].astype(str), strict=True))
    ids, ctgs = s["region_id"].astype(str).tolist(), s["contig"].astype(str).tolist()
    contigs_of: dict[str, set[str]] = {}
    for rid, ctg in zip(ids, ctgs, strict=True):
        contigs_of.setdefault(rid, set()).add(ctg)
    have = set(keys)
    ranked = sorted(
        (qi, pi if math.isfinite(pi) else math.inf, rid, ctg)
        for qi, pi, rid, ctg in zip(_num(s["q_value"]), _num(s["p_value"]), ids, ctgs, strict=True)
        if math.isfinite(qi) and (rid, ctg) in have
    )[:MAX_REGIONS]
    grid = "percentile" in nulls.columns
    out = []
    for qi, pi, rid, ctg in ranked:
        g = nulls[[k == (rid, ctg) for k in keys]]
        if grid:
            g = g.sort_values("percentile", kind="stable")
        values = tuple(_num(g["null_value"]).tolist()) if grid else _finite(g["null_value"])
        if not values:
            continue
        name = f"{rid} ({ctg})" if len(contigs_of[rid]) > 1 else rid
        obs = _first(g["observed"])
        out.append(
            Multiple(
                "regions_stats",
                f"{name} / {REGION_STAT}",
                values,
                obs,
                _caption(f"p {_fmt(pi)}, q {_fmt(qi)}", obs),
                tuple(float(x) / 100 for x in _num(g["percentile"])) if grid else None,
            )
        )
    return out


def _islands(view: RunView) -> list[Multiple]:
    """The genome-wide ``genome / abs`` null (replicate-wise maximum over the scanned replicons,
    NaN skipped, matched by row position; FWER p of the top island), then at most 7
    per-replicon multiples by ``observed_max`` descending, ties by contig (lead decision 4)."""
    nulls, isl = view.support("null_islands"), view.table("islands")
    if nulls is None or nulls.empty:
        return []
    groups = [
        (str(contig), str(direction), _num(g["max_scan_s"]), _first(g["observed_max"]))
        for (contig, direction), g in nulls.groupby(["contig", "direction"], sort=False)
    ]
    width = max(len(v) for _, _, v, _ in groups)
    stack = np.full((len(groups), width), np.nan)
    for i, (_, _, v, _) in enumerate(groups):
        stack[i, : len(v)] = v
    finite = np.isfinite(stack)
    genome = np.where(finite, stack, -np.inf).max(axis=0)[finite.any(axis=0)]
    if not len(genome):
        return []
    on = np.array([]) if isl is None else _num(isl["p_value"])
    p = float(np.min(on[np.isfinite(on)])) if np.isfinite(on).any() else math.nan
    obs_all = np.array([o for *_, o in groups])
    obs = float(np.max(obs_all[np.isfinite(obs_all)])) if np.isfinite(obs_all).any() else math.nan
    scanned = [g for g in groups if np.isfinite(g[2]).any()]
    scanned.sort(key=lambda g: (-g[3] if math.isfinite(g[3]) else math.inf, g[0], g[1]))
    out = [
        Multiple(
            "islands",
            GENOME,
            tuple(float(x) for x in genome),
            obs,
            _caption(f"FWER p {_fmt(p)}", obs),
            extra=((NOT_DRAWN, float(max(0, len(scanned) - MAX_REPLICONS))),),
        )
    ]
    for contig, direction, v, o in scanned[:MAX_REPLICONS]:
        values = tuple(float(x) for x in v[np.isfinite(v)])
        family = f"{contig} / {direction}"
        out.append(Multiple("islands", family, values, o, _caption(PER_REPLICON, o)))
    return out


def _determinants(view: RunView) -> list[Multiple]:
    nulls, dm = view.support("null_determinants"), view.table("determinants_model")
    if nulls is None or nulls.empty:
        return []
    dm = pd.DataFrame({c: [] for c in ("model", "metric", "value", "flags")}) if dm is None else dm
    out = []
    for model, g in nulls.groupby("model", sort=False):
        rows = dm[dm["model"].astype(str) == str(model)]
        flags: set[str] = set().union(*map(_flags, rows["flags"]))
        values = _finite(g["cv_r2"])
        if TOO_FEW in flags or not values:
            continue
        null_p = _first(rows["value"][rows["metric"] == "null_p"])
        obs = _first(g["observed_cv_r2_mean"])
        caption = _caption(f"null p {_fmt(null_p)} {UNCALIBRATED_NOTE}", obs)
        out.append(Multiple("determinants_model", str(model), values, obs, caption))
    return out


def multiples(view: RunView) -> list[Multiple]:
    """Every panel-c small multiple: regions (ascending q), islands, determinants."""
    return [*_regions(view), *_islands(view), *_determinants(view)]


def rows(m: Multiple) -> list[dict[str, object]]:
    """The data-TSV rows of one small multiple (panel ``c``): the points, the ``observed`` row,
    then the ``extra`` rows."""
    base: dict[str, object] = {"panel": "c", "table": m.table, "family": m.family}
    if m.grid is not None:
        pts = [base | {"x": x, "y": y} for x, y in zip(m.values, m.grid, strict=True)]
    else:
        edges, counts = histogram(m.values)
        pts = [base | {"bin": e, "n": n} for e, n in zip(edges, counts, strict=True)]
    extra = [base | {"metric": k, "value": v} for k, v in m.extra]
    return [*pts, base | {"metric": "observed", "value": m.observed}, *extra]


def draw(ax: Axes, m: Multiple) -> None:
    """The histogram (or step CDF) and the observed line of one small multiple."""
    if m.grid is not None:
        ax.step(m.values, m.grid, where="post", color=OKABE_ITO["sky_blue"], linewidth=1.0)
        ax.text(0.03, 0.95, GRID_TEXT, transform=ax.transAxes, va="top", fontsize=6)
        ax.set_ylabel("cumulative fraction", fontsize=6)
    else:
        edges, counts = histogram(m.values)
        width = edges[1] - edges[0]
        ax.bar(edges, counts, width=width, align="edge", color=OKABE_ITO["sky_blue"])
        ax.set_ylabel("count", fontsize=6)
    if math.isfinite(m.observed):
        ax.axvline(m.observed, color=OKABE_ITO["vermillion"], linewidth=1.0)
    ax.set_title(f"{m.table}\n{m.family}\n{m.caption}", fontsize=6)
    ax.set_xlabel(XLABEL[m.table], fontsize=6)
