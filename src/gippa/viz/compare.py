# ruff: noqa: RUF001, RUF002, RUF003  (ρ is the spec's own character)
"""V16 ``compare``: UpSet of detected genes across runs and the pairwise abundance concordance
matrix (GIPPA_SPEC.md §10.2 V16, §7.9.3; M9-23 and its lead decisions 1-3; ADR-0046).

Written only by ``gippa compare``: the registered ``(view, cfg)`` function is a skip on an ordinary
run. ``compare_plot`` draws ``compare_detection`` (panel a: exact intersections as bars over a
dot matrix, detected genes per run at the left; the 30 largest drawn, the rest summed as
``other``, the ``none`` count a footnote) and ``compare_abundance`` (panel b: lower triangle
cells coloured by ρ with "ρ r [lo, hi]", upper triangle the shared measured log2 of each pair,
read from each run's ``genes.tsv``, with the identity line; non-empty count bins when the drawn
points would exceed 20,000). Custom matplotlib: UpSetPlot fails on the locked pandas (index
decision 29). Pure: frames in, one ``Instance`` out.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.axes import Axes
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from gippa.viz import base, style
from gippa.viz.interactive import interactive

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

COLUMNS = ("panel", "intersection", "n_genes", "runs", "run_a", "run_b", "spearman_rho",
           "ci_low", "ci_high", "n_shared", "locus_tag", "x", "y")  # fmt: skip
EMPTY = "No gene detected in any run"
SKIP = "input_absent:compare"  # the gippa compare tables: V16 is written only by gippa compare
MAX_BARS = 30
MAX_POINTS = 20_000
BLUE, GREY, TEXT = style.DETECTION["detected"], "#BBBBBB", "#555555"
_BINS = 24  # per side, at most; fewer when many pairs share the MAX_POINTS // 2 bin budget
_CMAP = colormaps["cividis"]


def _g(v: float) -> str:
    return f"{v:.2g}"


def _layout(k: int) -> dict[str, float]:
    """Heights in inches, top to bottom: title, panel (a) bars and dots, panel (b) square."""
    dots = max(0.6, 0.16 * k)
    side = min(max(0.55 * k + 0.6, 2.8), 6.2)
    return {"bars_top": 0.65, "bars": 1.3, "dots": dots, "b_top": 0.65 + 1.3 + dots + 1.0,
            "side": side, "height": 0.65 + 1.3 + dots + 1.0 + side + 0.45}  # fmt: skip


def _box(
    fig: Figure, x: float, top: float, w: float, h: float
) -> tuple[float, float, float, float]:
    """Inches from the top-left corner -> a figure-fraction ``add_axes`` rectangle."""
    fw, fh = (float(v) for v in fig.get_size_inches())
    return (x / fw, 1 - (top + h) / fh, w / fw, h / fh)


def _measured(frame: pd.DataFrame) -> pd.Series:
    """log2 of the detected, not imputed, measured genes, indexed by locus tag."""
    ok = (frame["detected"].fillna(False).astype(bool)
          & ~frame["abundance_imputed"].fillna(False).astype(bool)
          & frame["abundance_log2"].notna())  # fmt: skip
    sub = frame[ok.to_numpy()]
    return pd.Series(sub["abundance_log2"].to_numpy(dtype=float), index=sub["locus_tag"])


def _panel_a(detection: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """The drawn intersections (``other`` last when more than 30) and the ``none`` count."""
    inter = detection[detection["section"] == "intersection"]
    names = inter["intersection"].astype(str).to_numpy()
    counts = inter["n_genes"].astype("Int64").fillna(0).astype(int).to_numpy()
    none = int(counts[names == "none"].sum())
    keep = sorted(((n, c) for n, c in zip(names, counts, strict=True) if n != "none" and c > 0),
                  key=lambda nc: (-nc[1], nc[0]))  # fmt: skip
    rows = [("a", n, c, n.replace("&", ";")) for n, c in keep[:MAX_BARS]]
    if len(keep) > MAX_BARS:
        rows.append(("a", "other", sum(c for _, c in keep[MAX_BARS:]), None))
    return pd.DataFrame(rows, columns=["panel", "intersection", "n_genes", "runs"]), none


def _panel_b(abundance: pd.DataFrame, genes: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Per pair (label order): its ρ row, then one row per shared measured gene."""
    parts = []
    for _, r in abundance.iterrows():
        ra, rb = str(r["run_a"]), str(r["run_b"])
        stat = {"run_a": ra, "run_b": rb, "n_shared": r["n_shared"],
                **{c: float(r[c]) for c in ("spearman_rho", "ci_low", "ci_high")}}  # fmt: skip
        parts.append(pd.DataFrame([stat]))
        xa, yb = _measured(genes[ra]), _measured(genes[rb])
        tags = [t for t in xa.index if t in yb.index]
        xy = {"x": xa[tags].to_numpy(), "y": yb[tags].to_numpy()}
        parts.append(pd.DataFrame({"run_a": ra, "run_b": rb, "locus_tag": tags, **xy}))
    out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return out.assign(panel="b")


def _frame(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    out = pd.concat([a, b], ignore_index=True).reindex(columns=list(COLUMNS))
    for c in ("n_genes", "n_shared"):
        out[c] = pd.array(out[c], dtype="Int64")
    for c in ("spearman_rho", "ci_low", "ci_high", "x", "y"):
        out[c] = out[c].astype(float)
    for c in ("panel", "intersection", "runs", "run_a", "run_b", "locus_tag"):
        out[c] = [None if pd.isna(v) else str(v) for v in out[c]]
    return out


def title(labels: list[str], inter: pd.DataFrame, stats: pd.DataFrame) -> str:
    """``inter``: every intersection row (not only the 30 drawn)."""
    n_all = int(inter.loc[inter["intersection"] == "&".join(labels), "n_genes"].sum())
    rho = stats["spearman_rho"].to_numpy(dtype=float)
    rho = rho[np.isfinite(rho)]
    med = _g(float(np.median(rho))) if rho.size else "NA"
    return (f"{base.plural(n_all, 'gene')} detected in all {len(labels)} runs; "
            f"median pairwise ρ {med}")  # fmt: skip


def _upset(fig: Figure, a: pd.DataFrame, labels: list[str], sizes: list[int]) -> None:
    drawn = a[a["intersection"] != "other"]
    k, m = len(labels), len(drawn)
    lay = _layout(k)
    dots_top = lay["bars_top"] + lay["bars"]
    bars = fig.add_axes(_box(fig, 2.2, lay["bars_top"], 4.8, lay["bars"] - 0.05), label="a")
    dots = fig.add_axes(_box(fig, 2.2, dots_top, 4.8, lay["dots"]), label="a_dots", sharex=bars)
    sets = fig.add_axes(_box(fig, 0.35, dots_top, 1.0, lay["dots"]), label="a_sets", sharey=dots)
    xs = np.arange(m)
    bars.bar(xs, drawn["n_genes"].to_numpy(dtype=float), width=0.7, color=BLUE)
    bars.set_ylabel("genes")
    bars.tick_params(axis="x", bottom=False, labelbottom=False)
    bars.set_title("(a) genes detected per exact set of runs", loc="left")
    for side in ("top", "right"):
        bars.spines[side].set_visible(False)
    ys = np.arange(k)
    gx, gy = np.meshgrid(xs, ys)
    dots.scatter(gx.ravel(), gy.ravel(), s=12, color=GREY, linewidths=0, zorder=1)
    for x, runs in zip(xs, drawn["runs"], strict=True):
        hit = [labels.index(r) for r in str(runs).split(";")]
        dots.plot([x, x], [min(hit), max(hit)], color=BLUE, linewidth=1.0, zorder=2)
        dots.scatter([x] * len(hit), hit, s=12, color=BLUE, linewidths=0, zorder=3)
    dots.set_xlim(-0.6, max(m, 1) - 0.4)
    dots.set_ylim(k - 0.5, -0.5)
    dots.set_yticks(ys, labels)
    dots.tick_params(axis="both", length=0, labelbottom=False)
    for spine in dots.spines.values():
        spine.set_visible(False)
    sets.barh(ys, sizes, height=0.6, color=BLUE)
    sets.invert_xaxis()
    sets.set_xlabel("detected genes")
    sets.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "left"):
        sets.spines[side].set_visible(False)


def _cell_text(r: pd.Series[Any], sep: str | None) -> str:
    """ "ρ r [lo, hi]"; no CI when ``sep`` is None (8 or more runs: the CI is in the TSV)."""
    rho, lo, hi = (float(r[c]) for c in ("spearman_rho", "ci_low", "ci_high"))
    if not math.isfinite(rho) or sep is None:
        return f"ρ {_g(rho)}" if math.isfinite(rho) else "ρ NA"
    ci = f" [{_g(lo)}, {_g(hi)}]" if math.isfinite(lo) and math.isfinite(hi) else " [NA]"
    return f"ρ {_g(rho)}{sep}{ci.strip()}"


def _scatter(ax: Axes, row: int, col: int, pts: pd.DataFrame, bins: int) -> None:
    """Points, or (``bins`` > 0) one square per non-empty bin of a ``bins`` x ``bins`` grid."""
    x, y = pts["x"].to_numpy(dtype=float), pts["y"].to_numpy(dtype=float)
    ax.plot([col + 0.08, col + 0.92], [row + 0.92, row + 0.08], color=TEXT, lw=0.5, zorder=2)
    if not len(x):
        return
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    span = hi - lo or 1.0
    u, v = 0.08 + 0.84 * (x - lo) / span, 0.08 + 0.84 * (y - lo) / span
    if not bins:
        ax.scatter(col + u, row + 1 - v, s=2, color=BLUE, linewidths=0, zorder=3)
        return
    counts, vs, us = np.histogram2d(v, u, bins=bins, range=((0.08, 0.92), (0.08, 0.92)))
    i, j = np.nonzero(counts)
    x0, x1, y0, y1 = col + us[j], col + us[j + 1], row + 1 - vs[i], row + 1 - vs[i + 1]
    verts = np.stack([np.column_stack(c) for c in ((x0, y0), (x1, y0), (x1, y1), (x0, y1))], 1)
    ax.add_collection(PolyCollection(list(verts), array=counts[i, j], cmap="viridis",
                                     linewidths=0, zorder=3))  # fmt: skip


def _matrix(fig: Figure, b: pd.DataFrame, abundance: pd.DataFrame, labels: list[str]) -> None:
    k = len(labels)
    lay = _layout(k)
    side, left = lay["side"], (7.2 - lay["side"]) / 2 - 0.3
    ax = fig.add_axes(_box(fig, left, lay["b_top"], side, side), label="b")
    sep = " " if k <= 3 else "\n" if k < 8 else None
    norm = Normalize(0.0, 1.0)
    stats = b[b["locus_tag"].isna()]
    pts = b[b["locus_tag"].notna()]
    bins = min(_BINS, math.isqrt(MAX_POINTS // 2 // max(len(stats), 1)))
    bins = bins if len(pts) > MAX_POINTS else 0
    for i, lab in enumerate(labels):
        ax.text(i + 0.5, i + 0.5, lab, ha="center", va="center", fontsize=7, weight="bold")
    for _, r in stats.iterrows():
        ia, ib = labels.index(str(r["run_a"])), labels.index(str(r["run_b"]))
        rho = float(r["spearman_rho"])
        face = GREY if not math.isfinite(rho) else _CMAP(norm(max(rho, 0.0)))
        hatch = "///" if math.isfinite(rho) and rho < 0 else None
        ax.add_patch(Rectangle((ia, ib), 1, 1, facecolor=face, edgecolor="white", hatch=hatch))
        dark = math.isfinite(rho) and rho < 0.55
        cell = _cell_text(r, sep)
        box = {"facecolor": face, "edgecolor": "none", "pad": 1} if hatch else None
        ax.text(ia + 0.5, ib + 0.5, cell, ha="center", va="center", fontsize=6, bbox=box,
                color="white" if dark else "black", linespacing=1.1)  # fmt: skip
        sub = pts[(pts["run_a"] == r["run_a"]) & (pts["run_b"] == r["run_b"])]
        _scatter(ax, ia, ib, sub, bins)
    ax.set_xlim(0, k)
    ax.set_ylim(k, 0)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("(b) Spearman ρ [95 % CI] (lower), shared log2 (upper)", loc="left")
    cax = fig.add_axes(
        _box(fig, left + side + 0.25, lay["b_top"] + side / 2, 0.1, side / 2), label="b_colourbar"
    )
    steps = np.linspace(0.0, 1.0, 11)  # stepped: < 50 solids stay vector in the SVG
    fig.colorbar(ScalarMappable(norm=norm, cmap=_CMAP), cax=cax, boundaries=steps, ticks=steps[::5],
                 label="ρ (< 0 hatched)")  # fmt: skip


def compare_plot(
    detection: pd.DataFrame, abundance: pd.DataFrame, genes: Mapping[str, pd.DataFrame]
) -> base.Instance:
    """V16 from the ``gippa compare`` tables and each run's ``genes.tsv`` (labels in order)."""
    labels = [str(lab) for lab in genes]
    a, none = _panel_a(detection)
    if a.empty:
        return base.empty_state("V16", "compare", EMPTY, COLUMNS, "double")
    b = _panel_b(abundance, genes)
    data = _frame(a, b)
    inter = detection[detection["section"] == "intersection"]
    text = title(labels, inter, b[b["locus_tag"].isna()] if "locus_tag" in b else b)
    matrix = detection[detection["section"] == "matrix"]
    hit = matrix[matrix["detected"].fillna(False).astype(bool)]
    sizes = [int((hit["run"] == lab).sum()) for lab in labels]
    fig = style.new_figure("double", _layout(len(labels))["height"])
    fig.set_label("V16_compare")
    _upset(fig, a, labels, sizes)
    _matrix(fig, b, abundance, labels)
    notes = [f"{base.plural(none, 'gene')} detected in no run (none)"]
    if (a["intersection"] == "other").any():
        rest = int(a.loc[a["intersection"] == "other", "n_genes"].iloc[0])
        notes.append(f"{base.plural(rest, 'gene')} in smaller intersections: other (data TSV)")
    height = float(fig.get_size_inches()[1])
    fig.text(0.05, 1 - (_layout(len(labels))["b_top"] - 0.3) / height, "; ".join(notes),
             fontsize=6, color=TEXT)  # fmt: skip
    few = abundance["flags"].fillna("").astype(str).str.split(";").map(
        lambda f: "too_few_shared" in f).any()  # fmt: skip
    fig.text(0.05, 0.15 / height, "upper cells: shared measured log2 abundance, x = row run, "
             "y = column run, grey line = identity" + "; ρ NA: too few shared measured genes"
             * bool(few), fontsize=6, color=TEXT)  # fmt: skip
    fig.suptitle(text, y=1 - 0.12 / height, va="top")
    return base.Instance(None, fig, data, text)


@base.plot("V16", "compare")
def compare(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """V16 is written only by ``gippa compare`` (§10.2): a skip on an ordinary run."""
    return base.Skip(SKIP)


@interactive("V16")
def compare_twin(data: pd.DataFrame, text: str) -> Any:
    """The concordance scatter: one trace per pair (the first shown), identity line, ρ in the
    trace name."""
    import plotly.graph_objects as go

    fig = go.Figure()
    b = data[data["panel"].astype(str) == "b"]
    tags = b["locus_tag"]
    stats, pts = b[tags.isna()], b[tags.notna()]
    for i, (_, r) in enumerate(stats.iterrows()):
        sub = pts[(pts["run_a"] == r["run_a"]) & (pts["run_b"] == r["run_b"])]
        rho = float(r["spearman_rho"])
        name = f"{r['run_a']} vs {r['run_b']} (ρ {_g(rho) if math.isfinite(rho) else 'NA'})"
        fig.add_trace(go.Scattergl(
            x=sub["x"].astype(float).tolist(), y=sub["y"].astype(float).tolist(),
            text=sub["locus_tag"].astype(str).tolist(), mode="markers", name=name,
            marker={"size": 4, "color": BLUE}, visible=True if i == 0 else "legendonly",
            hovertemplate="%{text}<br>x %{x:.3g}<br>y %{y:.3g}<extra></extra>",
        ))  # fmt: skip
    vals = pd.concat([pts["x"], pts["y"]]).astype(float)
    if len(vals):
        lo, hi = float(vals.min()), float(vals.max())
        fig.add_trace(go.Scatter(x=[lo, hi], y=[lo, hi], mode="lines", name="identity",
                                 line={"color": TEXT, "width": 1}))  # fmt: skip
    fig.update_layout(title=text, xaxis_title="log2 abundance (run a)",
                      yaxis_title="log2 abundance (run b)", template="simple_white")  # fmt: skip
    return fig
