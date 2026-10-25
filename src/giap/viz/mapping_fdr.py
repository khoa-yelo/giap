"""V6 ``mapping_fdr``: can I trust tolerant hits and de novo peptides? (GIAP_SPEC.md §10.2 V6;
M9-05). Draws the M3 ``mapping_fdr.tsv`` estimates; computes none.

(a) targets vs. decoys per ``mapping.length_classes`` class x ``n_diff`` (log scale, 0 at the
floor); (b) the per-m ``fdr_stratum`` (``length_class = all``) with a dashed line at
``stats.alpha.placement`` and ``n_accepted`` above each point, NA strata as open markers at the
floor; (c) de novo precision per score bin with its CI, only with ``denovo_precision`` rows.
Panels (a) and (b) use the ``reverse`` rows, or the ``shuffle`` rows when only those have targets
(M9-05 lead decision 3). The title is built from the panel (b) rows (index decision 5).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from giap.mapping.placement_fdr import length_class_labels
from giap.viz import base, style

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = (
    "panel",
    "length_class",
    "n_diff",
    "score_bin",
    "metric",
    "value",
    "ci_low",
    "ci_high",
    "flags",
)
TARGET = style.OKABE_ITO["blue"]
DECOY = style.OKABE_ITO["vermillion"]
FOOTNOTE = "c assumed 0.99 (no reference peptides)"
EMPTY = "no peptide placements"
#: the whole-figure empty state (``base.empty_state``) when ``mapping_fdr`` has no rows at all
EMPTY_TITLE = "No peptide placements or de novo precision in this run"
_FLOOR_A = 0.5  # panel (a) log-axis floor: 0 counts are drawn here


def _py(value: object) -> object:
    """A plain Python cell: ``None`` for NA, ``int`` / ``float`` for numbers."""
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _num(value: object) -> float:
    v = _py(value)
    return math.nan if v is None else float(str(v))


def placement_rows(fdr: pd.DataFrame) -> pd.DataFrame:
    """The ``placement`` rows of one decoy type: ``reverse`` unless only ``shuffle`` has targets."""
    rows = fdr[fdr["section"] == "placement"]
    by_type = [rows[rows["decoy_type"] == t] for t in ("reverse", "shuffle")]
    with_targets = [r for r in by_type if r["targets"].fillna(0).sum() > 0]
    present = [r for r in by_type if len(r)]
    return (with_targets or present or [rows.iloc[0:0]])[0]


def _row(
    panel: str,
    lc: object,
    m: object,
    score_bin: object,
    metric: str,
    value: object,
    lo: object = None,
    hi: object = None,
    flags: object = None,
) -> list[object]:
    return [panel, lc, m, score_bin, metric, *(_py(v) for v in (value, lo, hi, flags))]


def data_frame(fdr: pd.DataFrame, classes: Sequence[str]) -> pd.DataFrame:
    """The V06 data TSV: panel (a), (b), (c) rows sorted as in M9-05."""
    rows: list[list[object]] = []
    place = placement_rows(fdr)
    per_class = place[place["length_class"] != "all"]
    extra = [c for c in dict.fromkeys(per_class["length_class"].astype(str)) if c not in classes]
    strata = sorted({int(m) for m in place["n_diff"].dropna()})
    for lc in [*classes, *extra]:
        for m in strata:
            hit = per_class[(per_class["length_class"] == lc) & (per_class["n_diff"] == m)]
            for metric in ("decoys", "targets"):
                value = hit[metric].iloc[0] if len(hit) else None
                rows.append(_row("a", lc, m, None, metric, value))
    totals = place[place["length_class"] == "all"]
    for m in strata:
        hit = totals[totals["n_diff"] == m]
        for metric in ("fdr_stratum", "n_accepted"):
            rows.append(
                _row("b", "all", m, None, metric, hit[metric].iloc[0] if len(hit) else None)
            )
    denovo = fdr[fdr["section"] == "denovo_precision"].sort_values("score_bin", kind="mergesort")
    for r in denovo.itertuples(index=False):
        rows.append(
            _row(
                "c",
                None,
                None,
                _py(r.score_bin),
                "precision",
                r.precision,
                r.ci_low,
                r.ci_high,
                r.flags,
            )
        )
    return pd.DataFrame(rows, columns=list(COLUMNS), dtype=object)


def _strata(b: pd.DataFrame) -> list[tuple[int, float, int]]:
    """``(n_diff, fdr_stratum or NaN, n_accepted)`` per stratum of the panel (b) rows."""
    out = []
    for m in dict.fromkeys(b["n_diff"]):
        sub = b[b["n_diff"] == m].set_index("metric")["value"]
        acc = _num(sub.get("n_accepted"))
        out.append((int(m), _num(sub.get("fdr_stratum")), 0 if math.isnan(acc) else int(acc)))
    return out


def _join(strata: Sequence[int]) -> str:
    return ", ".join(str(m) for m in strata)


def title(data: pd.DataFrame, alpha: float) -> str:
    """From the panel (b) rows: the strata with FDR ≤ alpha, the accepted count, the NA strata."""
    strata = _strata(data[data["panel"] == "b"])
    if not strata:
        bins = int((data["panel"] == "c").sum())
        tail = f"; de novo precision in {base.plural(bins, 'score bin')}" if bins else ""
        return f"Placement FDR unestimated: no peptide placements{tail}"
    ok = [m for m, f, _ in strata if not math.isnan(f) and f <= alpha]
    na = [m for m, f, _ in strata if math.isnan(f)]
    accepted = f"{sum(n for _, _, n in strata)} FDR-tested placements accepted"
    if ok:
        parts = [f"Placement FDR ≤ {alpha:g} in strata {_join(ok)}", accepted]
    elif len(na) < len(strata):
        parts = [f"No stratum with placement FDR ≤ {alpha:g}", accepted]
    else:
        parts = [accepted]
    if na:
        noun = "stratum" if len(na) == 1 else "strata"
        parts.append(f"{noun} {_join(na)} unestimated (no targets)")
    return "; ".join(parts)


def _empty(ax: Axes) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.5, EMPTY, ha="center", va="center", transform=ax.transAxes, color="#555555")


def _panel_a(ax: Axes, a: pd.DataFrame) -> None:
    classes = list(dict.fromkeys(a["length_class"]))
    strata = list(dict.fromkeys(a["n_diff"]))
    width = 0.8 / (2 * len(strata))
    ax.set_yscale("log")
    top = _FLOOR_A
    minor_x, minor_lab = [], []
    for i, lc in enumerate(classes):
        for j, m in enumerate(strata):
            x0 = i - 0.4 + 2 * j * width
            minor_x.append(x0 + width)
            minor_lab.append(str(m))
            for k, (metric, colour) in enumerate((("targets", TARGET), ("decoys", DECOY))):
                hit = a[(a["length_class"] == lc) & (a["n_diff"] == m) & (a["metric"] == metric)]
                v = _num(hit["value"].iloc[0]) if len(hit) else math.nan
                h = max(v, _FLOOR_A) - _FLOOR_A if not math.isnan(v) else 0.0
                top = max(top, _FLOOR_A + h)
                label = metric if (i, j) == (0, 0) else None
                ax.bar(
                    x0 + (k + 0.5) * width,
                    h,
                    width,
                    bottom=_FLOOR_A,
                    color=colour,
                    label=label,
                    linewidth=0,
                )
    ax.set_ylim(_FLOOR_A, max(10.0, top * 2))
    ax.set_xticks(minor_x, minor_lab)
    ax.tick_params(axis="x", length=2)
    for i, lc in enumerate(classes):
        ax.text(i, -0.13, lc, transform=ax.get_xaxis_transform(), ha="center", va="top")
    ax.set_xlabel("differences per length class (aa)", labelpad=12)
    ax.set_ylabel("placements (0 at the floor)")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("(a) targets vs. decoys", loc="left")


def _panel_b(ax: Axes, b: pd.DataFrame, alpha: float) -> None:
    strata = _strata(b)
    finite = [f for _, f, _ in strata if not math.isnan(f)]
    ax.set_xlim(-0.6, len(strata) - 0.4 if strata else 0.6)
    ax.set_ylim(0.0, 1.3 * max([alpha, *finite]))
    ax.axhline(alpha, linestyle="--", color="#555555", linewidth=0.8)
    est = [(m, f) for m, f, _ in strata if not math.isnan(f)]
    if est:
        ax.plot(
            [m for m, _ in est],
            [f for _, f in est],
            marker="o",
            linestyle="none",
            color=TARGET,
            markersize=4,
        )
    na = [m for m, f, _ in strata if math.isnan(f)]
    if na:
        ax.plot(
            na,
            [0.0] * len(na),
            marker="o",
            linestyle="none",
            markerfacecolor="none",
            markeredgecolor="#999999",  # Okabe-Ito grey (§10.1)
            markersize=4,
            clip_on=False,
        )
    for m, f, n in strata:
        y = 0.0 if math.isnan(f) else f
        ax.annotate(
            str(n), (m, y), xytext=(0, 4), textcoords="offset points", ha="center", fontsize=6
        )
        if math.isnan(f):
            ax.annotate(
                "NA (no targets)",
                (m, y),
                xytext=(0, 20),
                textcoords="offset points",
                ha="center",
                fontsize=6,
                color="#555555",
            )
    ax.set_xticks([m for m, _, _ in strata])
    ax.set_xlabel("differences (n_diff)")
    ax.set_ylabel("stratum FDR")
    ax.set_title(f"(b) stratum FDR (dashed: {alpha:g})", loc="left")


def _panel_c(ax: Axes, c: pd.DataFrame) -> bool:
    """Draw panel (c); ``True`` when a row is ``c_assumed`` (the footnote is needed)."""
    x = np.array([_num(v) for v in c["score_bin"]])
    y = np.array([_num(v) for v in c["value"]])
    lo = np.array([_num(v) for v in c["ci_low"]])
    hi = np.array([_num(v) for v in c["ci_high"]])
    assumed = np.array(["c_assumed" in str(f).split(";") for f in c["flags"]], dtype=bool)
    ci = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
    if ci.any():
        ax.errorbar(
            x[ci],
            y[ci],
            yerr=[y[ci] - lo[ci], hi[ci] - y[ci]],
            fmt="none",
            ecolor="#000000",
            elinewidth=0.8,
            capsize=0,
        )
    for hollow in (False, True):
        sel = np.isfinite(y) & (assumed == hollow)
        if sel.any():
            ax.plot(
                x[sel],
                y[sel],
                marker="o",
                linestyle="none",
                markersize=4,
                markeredgecolor=TARGET,
                markerfacecolor="none" if hollow else TARGET,
            )
    vals = np.concatenate([y[np.isfinite(y)], lo[np.isfinite(lo)]])
    ax.set_ylim(max(0.0, float(vals.min()) - 0.05) if len(vals) else 0.0, 1.01)
    ax.set_xticks([int(v) for v in x if math.isfinite(v)])
    ax.set_xlabel("de novo score bin")
    ax.set_ylabel("precision (CI)")
    ax.set_title("(c) de novo precision", loc="left")
    return bool(assumed.any())


@base.plot("V06", "mapping_fdr", needs=("mapping_fdr",))
def mapping_fdr(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V6: placement targets/decoys, stratum FDR against alpha, de novo precision."""
    fdr = view.table("mapping_fdr")
    assert fdr is not None  # guaranteed by ``needs``
    alpha = float(cfg.stats.alpha.placement)
    classes = length_class_labels([list(c) for c in cfg.mapping.length_classes])
    data = data_frame(fdr, classes)
    if not len(data):
        return [base.empty_state("V06", "mapping_fdr", EMPTY_TITLE, COLUMNS, "double")]
    panels = [p for p in ("a", "b", "c") if p in set(data["panel"]) or p != "c"]
    fig = style.new_figure("double", 2.8)
    fig.set_label("V06_mapping_fdr")
    grid = fig.add_gridspec(
        1,
        len(panels),
        width_ratios=[2.6, 1.5, 1.3][: len(panels)],
        left=0.08,
        right=0.98,
        bottom=0.24,
        top=0.8,
        wspace=0.35,
    )
    axes = {p: fig.add_subplot(grid[0, i], label=p) for i, p in enumerate(panels)}
    sub = {p: data[data["panel"] == p] for p in panels}
    if len(sub["a"]) and len(sub["b"]):
        _panel_a(axes["a"], sub["a"])
        _panel_b(axes["b"], sub["b"], alpha)
        if (placement_rows(fdr)["decoy_type"] == "shuffle").any():  # lead decision 3 fallback
            for p in ("a", "b"):
                # on a second line: one line runs into the next panel's title (§14.4)
                axes[p].set_title(axes[p].get_title(loc="left") + "\n(shuffle decoys)", loc="left")
    else:
        _empty(axes["a"])
        _empty(axes["b"])
    if "c" in axes and _panel_c(axes["c"], sub["c"]):
        fig.text(0.98, 0.02, FOOTNOTE, ha="right", va="bottom", fontsize=6, color="#555555")
    text = title(data, alpha)
    fig.suptitle(text, y=0.97)
    return [base.Instance(None, fig, data, text)]
