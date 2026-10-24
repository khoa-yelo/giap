"""V18 ``stat_diagnostics`` pass 1: (a) detectability reliability, (b) p-value histograms per test
family (M9-07; GIAP_SPEC.md §10.2 V18, §7.3.3, §8.1, §8.2, §9.2); (c) the null distributions
against the observed statistic (M9-21, ``giap.viz.null_panels``).

Pure: tables come from a ``RunView``; the figure is built without pyplot. Histogram bins are
exact: a p-value is binned by its shortest decimal form as a rational number, so a reported
``k/20`` (a permutation p with B + 1 = 20, 40, ...) lands in bin ``k``; 1 is in the last bin.
A p outside [0, 1] (including +-inf) is never clamped into an end bar: it is counted per family
and reported as a panel-b row ``metric = n_outside_0_1`` after that family's bin rows. NaN is
missing: neither binned nor counted.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from fractions import Fraction
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from giap.viz import null_panels
from giap.viz.base import Instance, plot
from giap.viz.style import OKABE_ITO, new_figure

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

#: the test families, in the spec's order: (table, key columns); the single source (spec table)
FAMILIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("regions_stats", ("stat",)),
    ("autocorrelation", ("statistic", "null_model")),
    ("pairs", ("test",)),
    ("enrichment_detection", ("collection", "direction")),
    ("enrichment_gsea", ("collection",)),
)
#: tables with a ``p_value`` column that are deliberately not histogrammed: the islands' p is a
#: family-wise maximum-statistic p (panel c, M9-21). ``determinants_model`` has ``null_p``.
NOT_DRAWN: tuple[str, ...] = ("islands",)
COLUMNS = ("panel", "table", "family", "bin", "x", "y", "n", "metric", "value")
N_BINS = 20
MAX_DRAWN = 16
GRID_COLS = 4
CONSTANT_TEXT = "detectability constant: no calibration curve"
ABSENT_TEXT = "detectability not available"
NO_P_TEXT = "no p-values in this run"
NO_NULL_TEXT = "null draws not available"
OUTSIDE_METRIC = "n_outside_0_1"
#: the one explanatory row of an input-free run (``panel`` NA, ``value`` = the panel texts)
MESSAGE_METRIC = "empty_state"
_MAX_AREA = 60.0  # pt², the marker area of the largest bin
_GREY = OKABE_ITO["grey"]
_TEXT = "#555555"

_Row = dict[str, object]


def _row(panel: str, table: str, **cells: object) -> _Row:
    return {c: cells.get(c) for c in COLUMNS} | {"panel": panel, "table": table}


def _g2(v: float) -> str:
    return "NA" if not math.isfinite(v) else f"{v:.2g}"


def in_unit(p: float) -> bool:
    """Whether ``p`` is in the closed interval [0, 1] (False for NaN and +-inf)."""
    return 0.0 <= p <= 1.0


def p_bin(p: float) -> int:
    """The bin ``0..19`` of ``p`` in [0, 1]: ``[k/20, (k+1)/20)``, the last closed.

    Raises ``ValueError`` outside [0, 1]: such values are counted, never clamped.
    """
    if not in_unit(float(p)):
        raise ValueError(f"p-value outside [0, 1]: {p!r}")
    return min(N_BINS - 1, math.floor(Fraction(repr(float(p))) * N_BINS))


# --- panel (a) --------------------------------------------------------------------------------


def _detectability(view: RunView) -> pd.DataFrame | None:
    ds = view.table("detection_summary")
    if ds is None:
        return None
    return ds[ds["section"] == "detectability"]


def _reliability_rows(det: pd.DataFrame) -> tuple[list[_Row], dict[str, float]]:
    """Bin rows in numeric bin order, then the overall metrics sorted by name."""
    bins = det[det["group_type"] == "bin"]
    wide = bins.pivot_table(index="group", columns="metric", values="value", aggfunc="first")
    rows = []
    for g in sorted(wide.index, key=lambda s: int(str(s))):
        rec = wide.loc[g]

        def get(m: str, rec: pd.Series = rec) -> float:
            return float(rec[m]) if m in rec.index and pd.notna(rec[m]) else math.nan

        n = get("n")
        rows.append(
            _row(
                "a",
                "detection_summary",
                bin=int(str(g)),
                x=get("mean_p"),
                y=get("observed_rate"),
                n=int(n) if math.isfinite(n) else None,
            )
        )
    over = det[det["group_type"] == "overall"]
    overall = {str(m): float(v) for m, v in zip(over["metric"], over["value"], strict=True)}
    for m in sorted(overall):
        rows.append(_row("a", "detection_summary", metric=m, value=overall[m]))
    return rows, overall


def _draw_reliability(ax: Axes, rows: Sequence[_Row], overall: dict[str, float]) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("mean predicted p_detect")
    ax.set_ylabel("observed detection rate")
    ax.set_title("a  detectability reliability", fontsize=7, loc="left")
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color=_GREY, linewidth=0.8)
    pts = [r for r in rows if r["metric"] is None]
    xs = [float(str(r["x"])) for r in pts]
    ys = [float(str(r["y"])) for r in pts]
    ns = [float(str(r["n"])) if r["n"] is not None else 0.0 for r in pts]
    top = max(ns) if ns and max(ns) > 0 else 1.0
    colour = OKABE_ITO["blue"]
    ax.plot(xs, ys, color=colour, linewidth=0.8)
    ax.scatter(xs, ys, s=[_MAX_AREA * n / top for n in ns], color=colour, zorder=3)
    lines = [
        f"slope {_g2(overall.get('calibration_slope', math.nan))}",
        f"intercept {_g2(overall.get('calibration_intercept', math.nan))}",
        f"Brier {_g2(overall.get('brier', math.nan))}",
    ]
    ax.text(0.04, 0.96, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=6)


def _message(ax: Axes, text: str) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=_TEXT)


# --- panel (b) --------------------------------------------------------------------------------


_Family = tuple[str, str, list[int], int]


def _family_counts(view: RunView) -> list[_Family]:
    """(table, family, 20 counts, n outside [0, 1]) for every family with a non-NaN p, in table
    then family order."""
    out = []
    for table, keys in FAMILIES:
        frame = view.table(table)
        if frame is None or frame.empty:
            continue
        p = pd.to_numeric(frame["p_value"], errors="coerce").to_numpy(dtype=float)
        keep = ~np.isnan(p)
        if table == "regions_stats":
            keep &= (frame["subset"] == "all").to_numpy(dtype=bool)
        elif table == "pairs":
            keep &= (frame["test"] != "none").to_numpy(dtype=bool)
        labels = [
            " / ".join("NA" if pd.isna(v) else str(v) for v in vals)
            for vals in zip(*(frame[k] for k in keys), strict=True)
        ]
        groups: dict[str, list[int]] = {}
        outside: dict[str, int] = {}
        for label, value, ok in zip(labels, p, keep, strict=True):
            if not ok:
                continue
            counts = groups.setdefault(label, [0] * N_BINS)
            outside.setdefault(label, 0)
            if in_unit(float(value)):
                counts[p_bin(float(value))] += 1
            else:
                outside[label] += 1
        out += [(table, fam, groups[fam], outside[fam]) for fam in sorted(groups)]
    return out


def _draw_histogram(
    ax: Axes, table: str, family: str, counts: Sequence[int], n_outside: int
) -> None:
    edges = [k / N_BINS for k in range(N_BINS)]
    ax.bar(edges, counts, width=1 / N_BINS, align="edge", color=OKABE_ITO["sky_blue"])
    ax.set_xlim(0.0, 1.0)
    extra = f", {n_outside} outside [0, 1]" if n_outside else ""
    ax.set_title(f"{table}\n{family}, n = {sum(counts)}{extra}", fontsize=6)
    ax.set_xlabel("p-value", fontsize=6)
    ax.set_ylabel("count", fontsize=6)


# --- the plot ---------------------------------------------------------------------------------


@plot(
    "V18",
    "stat_diagnostics",
    optional=("detection_summary", *(t for t, _ in FAMILIES), "islands", "determinants_model"),
)
def stat_diagnostics(view: RunView, cfg: GiapConfig) -> list[Instance]:
    det = _detectability(view)
    a_rows, overall = ([], {}) if det is None or det.empty else _reliability_rows(det)
    has_bins = any(r["metric"] is None for r in a_rows)
    if det is None or det.empty:
        state = "absent"
    elif overall.get("detectability_constant") == 1.0:
        state = "constant"
    elif not has_bins:
        state = "absent"
    else:
        state = "fit"

    families = _family_counts(view)
    order = {t: i for i, (t, _) in enumerate(FAMILIES)}
    ranked = sorted(families, key=lambda f: (-sum(f[2]), order[f[0]], f[1]))
    keep = {(t, f) for t, f, _, _ in ranked[:MAX_DRAWN]}
    drawn = [f for f in families if (f[0], f[1]) in keep]

    nulls = null_panels.multiples(view)
    n_rows = math.ceil(len(drawn) / GRID_COLS)
    c_rows = max(1, math.ceil(len(nulls) / GRID_COLS))
    fig = new_figure("double", 2.6 + 1.5 * n_rows + 1.8 * c_rows, layout="constrained")
    fig.set_label("V18_stat_diagnostics")
    ratios = [1.6] + [1.0] * n_rows + [1.2] * c_rows
    gs = fig.add_gridspec(1 + n_rows + c_rows, GRID_COLS, height_ratios=ratios)
    ax_a = fig.add_subplot(gs[0, :2], label="a")
    if state == "fit":
        _draw_reliability(ax_a, a_rows, overall)
    else:
        _message(ax_a, CONSTANT_TEXT if state == "constant" else ABSENT_TEXT)
    if not drawn:
        _message(fig.add_subplot(gs[0, 2:], label="b"), NO_P_TEXT)
    for i, (table, family, counts, n_out) in enumerate(drawn):
        ax = fig.add_subplot(gs[1 + i // GRID_COLS, i % GRID_COLS], label=f"b:{table} / {family}")
        _draw_histogram(ax, table, family, counts, n_out)
    if not nulls:
        _message(fig.add_subplot(gs[1 + n_rows, :], label="c"), NO_NULL_TEXT)
    for i, m in enumerate(nulls):
        cell = gs[1 + n_rows + i // GRID_COLS, i % GRID_COLS]
        null_panels.draw(fig.add_subplot(cell, label=f"c:{m.table} / {m.family}"), m)

    k = len(drawn)
    if state == "fit":
        title = (
            f"Detectability calibration slope {_g2(overall.get('calibration_slope', math.nan))} "
            f"(intercept {_g2(overall.get('calibration_intercept', math.nan))}), "
            f"Brier {_g2(overall.get('brier', math.nan))}"
        )
    elif state == "constant":
        title = f"Detectability constant; {k} test families shown"
    else:
        title = f"Detectability not available; {k} test families shown"
    fig.suptitle(title)

    b_rows: list[_Row] = []
    for table, family, counts, n_out in families:
        b_rows += [
            _row("b", table, family=family, bin=k_ / N_BINS, n=c) for k_, c in enumerate(counts)
        ]
        if n_out:
            b_rows.append(_row("b", table, family=family, metric=OUTSIDE_METRIC, value=n_out))
    null_rows = [r for m in nulls for r in null_panels.rows(m)]
    every = [*a_rows, *b_rows, *null_rows]
    if not every:  # no input at all (decision 7): say why; a header-only TSV explains nothing
        # empty a_rows means no overall row, so the state is "absent" (never "constant")
        note = "; ".join((ABSENT_TEXT, NO_P_TEXT, NO_NULL_TEXT))
        every = [{c: None for c in COLUMNS} | {"metric": MESSAGE_METRIC, "value": note}]
    cells = [[r.get(c) for c in COLUMNS] for r in every]
    data = pd.DataFrame(cells, columns=list(COLUMNS), dtype=object)
    return [Instance(suffix=None, figure=fig, data=data, title=title)]
