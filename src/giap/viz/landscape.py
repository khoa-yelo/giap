"""V2 ``abundance_landscape``: where along each replicon are proteins abundant? (GIAP_SPEC.md §10.2
V2, §7.3.5(d); M9-09 and its lead decisions 1-12). Draws ``genes.tsv``, ``islands.tsv`` and the
M9-03 landmarks; the only computation is the display trend.

One panel per replicon with genes (FASTA order): measured detected genes (detected, not imputed,
not partial, finite ``abundance_log2``) as points, ``detected = false`` genes as a rug coloured
by status (every gene when no gene has an abundance value; the title is "No measured
abundance: detection only" whenever no gene is a point, lead decision 9), a
gene-level LOESS trend at the islands' span (``giap.stats.smooth.loess_detrend``; the legend
names the span and nothing is claimed from it), islands with ``p_value ≤ stats.alpha.context``
shaded (an origin-crossing island is one row, two spans) and oriC / ter lines when
``plots/support/replicon_landmarks.tsv`` calls them.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from giap.stats.smooth import loess_detrend, loess_span_points
from giap.viz import base, style

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = ("contig", "kind", "locus_tag", "position_bp", "start", "end", "value", "category")
KINDS = ("point", "rug", "trend", "island", "oric", "ter")
_INT_COLUMNS = ("position_bp", "start", "end")
NO_ABUNDANCE = "No measured abundance: detection only"
NO_ISLANDS = "islands table absent: no islands shaded"
GREY = "#555555"
_ROW_IN, _HEAD_IN, _FOOT_IN = 1.5, 0.45, 0.3


class Replicon(NamedTuple):
    contig: str
    length: int
    circular: bool


def replicons(view: RunView, genes: pd.DataFrame) -> list[Replicon]:
    """Replicons in FASTA (preflight) order, then any contig only ``genes.tsv`` names (linear,
    length = its largest gene end)."""
    pre = view.preflight() or {}
    out = [
        Replicon(str(r["id"]), int(r["length"]), r.get("topology") == "circular")
        for r in pre.get("replicons", [])
    ]
    known = {r.contig for r in out}
    for contig in sorted(set(genes["contig"].astype(str)) - known):
        ends = pd.to_numeric(genes.loc[genes["contig"].astype(str) == contig, "end"])
        out.append(Replicon(contig, int(ends.max()), False))
    return out


def _bools(col: pd.Series) -> np.ndarray:
    return np.asarray(col.astype("boolean").fillna(False), dtype=bool)


def _mids(frame: pd.DataFrame, length: int) -> np.ndarray:
    """0-based midpoints ``(start0 + end) / 2 mod L`` (unwrapped origin-crossing genes)."""
    s = np.asarray(pd.to_numeric(frame["start"]), dtype=float)
    e = np.asarray(pd.to_numeric(frame["end"]), dtype=float)
    return np.mod((s + e) / 2.0, float(length))


def _row(contig: str, kind: str, **cells: Any) -> dict[str, Any]:
    return {"contig": contig, "kind": kind, **cells}


def _binned_points(
    contig: str, mids: np.ndarray, ys: np.ndarray, length: int
) -> list[dict[str, Any]] | None:
    """Bin medians at the bin centres when a replicon has more than 20,000 points."""
    edges = base.binned(mids.astype(np.int64), mids.astype(np.int64) + 1, length)
    if edges is None:
        return None
    which = np.searchsorted(edges, mids, side="right") - 1
    rows = []
    for b in np.unique(which):
        centre = int((edges[b] + edges[b + 1]) // 2)
        median = float(np.median(ys[which == b]))
        rows.append(_row(contig, "point", position_bp=centre, value=median, category="detected"))
    return rows


def _trend(rep: Replicon, pts: list[dict[str, Any]], frac: float) -> list[dict[str, Any]]:
    """The LOESS fit on the drawn points in position order (mean of tied positions)."""
    x = np.array([p["x"] for p in pts], dtype=float)
    y = np.array([p["value"] for p in pts], dtype=float)
    ux, first, inv = np.unique(x, return_index=True, return_inverse=True)
    if loess_span_points(len(ux), frac) < 3:
        return []
    ym = np.bincount(inv, weights=y) / np.bincount(inv)
    fit = loess_detrend(ym, ux, frac, circular=rep.circular, period=float(rep.length)).fit
    return [
        _row(
            rep.contig,
            "trend",
            locus_tag=pts[i].get("locus_tag"),
            position_bp=pts[i]["position_bp"],
            value=float(f),
            x=float(ux[j]),
        )
        for j, (i, f) in enumerate(zip(first, fit, strict=True))
    ]


def point_mask(genes: pd.DataFrame) -> np.ndarray:
    """Genes drawn as points: detected, measured, not imputed, not partial."""
    y = np.asarray(pd.to_numeric(genes["abundance_log2"]), dtype=float)
    keep = _bools(genes["detected"]) & ~_bools(genes["abundance_imputed"])
    mask: np.ndarray = keep & ~_bools(genes["partial"]) & np.isfinite(y)
    return mask


def _gene_rows(
    rep: Replicon, genes: pd.DataFrame, frac: float, all_rug: bool
) -> list[dict[str, Any]]:
    g = genes[genes["contig"].astype(str) == rep.contig]
    if not len(g):
        return []
    mids = _mids(g, rep.length)
    detected = _bools(g["detected"])
    y = np.asarray(pd.to_numeric(g["abundance_log2"]), dtype=float)
    point = point_mask(g)
    rug = np.ones(len(g), dtype=bool) if all_rug else ~detected
    tags = g["locus_tag"].astype(str).to_numpy()
    status = g["detection_status"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for i in np.flatnonzero(rug):
        rows.append(_row(rep.contig, "rug", locus_tag=tags[i], position_bp=int(mids[i]),
                         x=float(mids[i]), category=status[i]))  # fmt: skip
    idx = np.flatnonzero(point & ~rug)
    idx = idx[np.lexsort((tags[idx], mids[idx]))]
    pts = _binned_points(rep.contig, mids[idx], y[idx], rep.length)
    if pts is None:
        pts = [
            _row(
                rep.contig,
                "point",
                locus_tag=tags[i],
                position_bp=int(mids[i]),
                value=float(y[i]),
                category=status[i],
                x=float(mids[i]),
            )
            for i in idx
        ]
    for p in pts:
        p.setdefault("x", float(p["position_bp"]))
    return [*rows, *pts, *_trend(rep, pts, frac)]


def _island_rows(rep: Replicon, islands: pd.DataFrame, alpha: float) -> list[dict[str, Any]]:
    isl = islands[islands["contig"].astype(str) == rep.contig]
    p = np.asarray(pd.to_numeric(isl["p_value"]), dtype=float)
    rows = []
    for (_, r), pv in zip(isl.iterrows(), p, strict=True):
        if not (np.isfinite(pv) and pv <= alpha):
            continue
        start, end = int(r["start"]), int(r["end"])
        mid = ((start + end) // 2) % rep.length
        rows.append(_row(rep.contig, "island", position_bp=mid, start=start, end=end,
                         category=str(r["direction"])))  # fmt: skip
    return rows


def _landmark_rows(rep: Replicon, landmarks: pd.DataFrame | None) -> list[dict[str, Any]]:
    """oriC / ter rows from ``replicon_landmarks.tsv`` (1-based positions) when ``called``."""
    if landmarks is None or not len(landmarks):
        return []
    hit = landmarks[
        (landmarks["contig"].astype(str) == rep.contig)
        & (landmarks["status"].astype(str) == "called")
    ]
    rows = []
    for _, r in hit.iterrows():
        for kind, col in (("oric", "oric_position"), ("ter", "ter_position")):
            v = pd.to_numeric(pd.Series([r[col]])).iloc[0]
            if pd.notna(v):
                rows.append(_row(rep.contig, kind, position_bp=int(v) - 1, x=float(int(v) - 1)))
    return rows


def _sort_key(r: Mapping[str, Any]) -> tuple[Any, ...]:
    pos = r.get("position_bp")
    return (r["contig"], KINDS.index(r["kind"]), math.inf if pos is None else pos,
            r.get("locus_tag") or "")  # fmt: skip


def data_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    """The V2 data TSV rows (0-based ``position_bp`` and ``start``; see ``Instance.coords``)."""
    ordered = sorted(rows, key=_sort_key)
    cols: dict[str, pd.Series] = {}
    for c in COLUMNS:
        values = [r.get(c) for r in ordered]
        if c in _INT_COLUMNS:
            cols[c] = pd.Series(values, dtype="Int64")
        elif c == "value":
            cols[c] = pd.Series([math.nan if v is None else v for v in values], dtype=float)
        else:
            cols[c] = pd.Series(values, dtype=object)
    return pd.DataFrame(cols, columns=list(COLUMNS))


def _lag1_clause(view: RunView, order: Sequence[str]) -> str:
    """The title clause "; lag-1 Moran's I <I> (q = <q>)" from the gene-label row of the first
    FASTA contig with a finite q; "" otherwise."""
    ac = view.table("autocorrelation")
    if ac is None or not len(ac):
        return ""
    sel = ac[
        (ac["statistic"].astype(str) == "morans_i")
        & (ac["null_model"].astype(str) == "gene_label")
        & (pd.to_numeric(ac["lag"]) == 1)
    ]
    for contig in [*order, *sorted(set(sel["contig"].astype(str)) - set(order))]:
        hit = sel[sel["contig"].astype(str) == contig]
        for v, q in zip(pd.to_numeric(hit["value"]), pd.to_numeric(hit["q_value"]), strict=True):
            if pd.notna(q) and pd.notna(v):
                return f"; lag-1 Moran's I {float(v):g} (q = {float(q):g})"
    return ""


def title(rows: Sequence[Mapping[str, Any]], alpha: float, clause: str) -> str:
    if not any(r["kind"] == "point" for r in rows):
        return NO_ABUNDANCE
    dirs = [r["category"] for r in rows if r["kind"] == "island"]
    n_high, n_low = dirs.count("high"), dirs.count("low")
    return f"{n_high} high and {n_low} low abundance islands at FWER ≤ {alpha:g}{clause}"


def _panel(ax: Axes, rep: Replicon, rows: Sequence[Mapping[str, Any]], frac: float) -> None:
    scale, unit = (1e6, "Mb") if rep.length > 2_000_000 else (1e3, "kb")
    kind: dict[str, list[Mapping[str, Any]]] = {k: [] for k in KINDS}
    for r in rows:
        kind[r["kind"]].append(r)
    for r in kind["island"]:
        colour, hatch = style.ISLAND.get(str(r["category"]), style.ISLAND["high"])
        start, end = int(r["start"]), int(r["end"])
        spans = [(start, min(end, rep.length))]
        if end > rep.length:
            spans.append((0, end - rep.length))
        for x0, x1 in spans:
            ax.axvspan(x0 / scale, x1 / scale, facecolor=colour, edgecolor=colour, hatch=hatch,
                       alpha=0.25, linewidth=0, zorder=0)  # fmt: skip
    pts = kind["point"]
    if pts:
        ax.plot([p["x"] / scale for p in pts], [p["value"] for p in pts], linestyle="none",
                marker="o", markersize=2, color=style.DETECTION["detected"],
                markeredgewidth=0, label="detected (measured)", zorder=2)  # fmt: skip
    if kind["trend"]:
        ax.plot([t["x"] / scale for t in kind["trend"]], [t["value"] for t in kind["trend"]],
                color="#000000", linewidth=1.0, label=f"trend (LOESS, span {frac:g})",
                zorder=3)  # fmt: skip
    for status in sorted({str(r["category"]) for r in kind["rug"]}):
        xs = [r["x"] / scale for r in kind["rug"] if r["category"] == status]
        colour = style.STATUS_COLOUR.get(status, style.DETECTION["not_detected"])
        ax.vlines(xs, 0.0, 0.07, transform=ax.get_xaxis_transform(), colors=colour,
                  linewidth=0.6, label=status.replace("_", " "))  # fmt: skip
    for k, ls, name in (("oric", "-", "oriC"), ("ter", ":", "ter")):
        for r in kind[k]:
            ax.axvline(r["x"] / scale, color=GREY, linestyle=ls, linewidth=0.8, label=name)
    ax.set_xlim(0.0, rep.length / scale)
    ax.set_xlabel(f"position ({unit})")
    ax.set_ylabel("log2 abundance")
    ax.set_title(rep.contig, loc="left")
    if pts:
        ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=6, ncol=4,
                  borderaxespad=0.1)  # fmt: skip


@base.plot(
    "V02", "abundance_landscape", needs=("genes_tsv",), optional=("islands", "autocorrelation")
)
def abundance_landscape(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V2: abundance along each replicon with the display trend, rug, islands and landmarks."""
    genes = view.table("genes_tsv")
    assert genes is not None  # guaranteed by ``needs``
    islands = view.table("islands")
    alpha = float(cfg.stats.alpha.context)
    frac = float(cfg.context.islands.detrend_span_frac)
    reps = [r for r in replicons(view, genes) if (genes["contig"].astype(str) == r.contig).any()]
    if not reps:
        return [base.empty_state("V02", "abundance_landscape", "No genes", COLUMNS, "double")]
    measured = bool(point_mask(genes).any())
    raw = np.asarray(pd.to_numeric(genes["abundance_log2"]), dtype=float)
    no_abundance = not np.isfinite(raw).any()  # no value at all: every gene is a rug tick
    landmarks = view.support("replicon_landmarks")
    per: dict[str, list[dict[str, Any]]] = {}
    for rep in reps:
        rows = _gene_rows(rep, genes, frac, all_rug=no_abundance)
        if islands is not None:
            rows += _island_rows(rep, islands, alpha)
        per[rep.contig] = rows + _landmark_rows(rep, landmarks)
    every = [r for rows in per.values() for r in rows]
    text = title(every, alpha, _lag1_clause(view, [r.contig for r in reps]) if measured else "")
    height = _HEAD_IN + _FOOT_IN + _ROW_IN * len(reps)
    fig = style.new_figure("double", height)
    fig.set_label("V02_abundance_landscape")
    grid = fig.add_gridspec(len(reps), 1, left=0.09, right=0.98, bottom=_FOOT_IN / height + 0.06,
                            top=1 - _HEAD_IN / height, hspace=0.75)  # fmt: skip
    for i, rep in enumerate(reps):
        _panel(fig.add_subplot(grid[i, 0], label=rep.contig), rep, per[rep.contig], frac)
    fig.suptitle(text, y=1 - 0.12 / height, va="top")
    if islands is None:
        fig.text(0.01, 0.01, NO_ISLANDS, fontsize=6, color=GREY, va="bottom")
    data = data_frame(every)
    return [base.Instance(None, fig, data, text, coords=("position_bp", "start"))]
