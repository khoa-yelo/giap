"""V11 ``operon_pairs``: do operon partners look alike? (GIPPA_SPEC.md §10.2 V11, §7.4, §7.4.5;
M9-11). Draws ``pairs.tsv`` and ``operon_flags.tsv`` values; computes no statistic.

(a) density-normalised histograms of |Δ| = |``delta_log2``| for OP, NO, CV and SW on 30 shared
bins over [0, max |Δ|] (index decision 16; [0, 1] when the maximum is 0); (b) the M6 ``lfdr``
against |Δ| with a dashed line at ``stats.alpha.operon_flags``; (c) upstream vs. downstream
``abundance_log2`` of the OP / NO pairs, coloured by class, the flags with ``q ≤ alpha`` drawn as
larger outlined markers, one shape per flag. Per the M9-11 lead decisions: ``gene_a`` is the
upstream gene on both strands (``pairs.tsv`` is in transcription order), pairs flagged
``abundance_imputed`` are left out of (a) and (c), and the title's ``k`` counts the
``operon_flags.tsv`` rows with ``q ≤ alpha``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from gippa.viz import base, style

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

COLUMNS = (
    "panel",
    "pair_id",
    "class",
    "abs_delta_log2",
    "lfdr",
    "upstream_log2",
    "downstream_log2",
    "flag",
    "q_value",
    "bin",
    "density",
    "n",
)
CLASS_ORDER = ("OP", "NO", "CV", "SW", "DV")
DENSITY_CLASSES = ("OP", "NO", "CV", "SW")
CLASS_COLOUR: Mapping[str, str] = MappingProxyType(
    {
        "OP": style.OKABE_ITO["blue"],
        "NO": style.OKABE_ITO["orange"],
        "CV": style.OKABE_ITO["bluish_green"],
        "SW": style.OKABE_ITO["reddish_purple"],
    }
)
#: the 7 §7.4.5 flags -> marker shapes (panel c)
FLAG_MARKERS: Mapping[str, str] = MappingProxyType(
    {
        "internal_promoter_supported": "^",
        "internal_terminator_supported": "v",
        "unexplained_discordance": "D",
        "conflict_terminator_in_operon": "s",
        "novel_operon_candidate": "*",
        "coabundant_despite_terminator": "P",
        "other": "X",
    }
)
N_BINS = 30
GREY = "#555555"
EMPTY_PAIRS = "No gene pairs in this run"
NOTE_B = "lfdr ranks pairs; flags use q"
THRESHOLD_LABEL = "flag threshold (q)"
_SMALL, _LARGE = 10.0, 36.0  # panel (c) marker areas (pt²)


def _py(value: object) -> object:
    """A plain Python cell: ``None`` for NA, ``int`` / ``float`` for numbers."""
    if value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value)):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return None if math.isnan(float(value)) else float(value)
    return value


def _floats(col: pd.Series) -> np.ndarray:
    return np.asarray(pd.to_numeric(col, errors="coerce"), dtype=float)


def _imputed(flags: pd.Series) -> np.ndarray:
    return np.array(
        ["abundance_imputed" in str(f).split(";") if _py(f) is not None else False for f in flags],
        dtype=bool,
    )


def _flag_lookup(flags: pd.DataFrame | None) -> dict[str, tuple[str, float]]:
    """``pair_id`` -> ``(flag, q_value)`` (the first row of a pair)."""
    out: dict[str, tuple[str, float]] = {}
    if flags is None:
        return out
    for pid, flag, q in zip(
        flags["pair_id"], flags["flag"], _floats(flags["q_value"]), strict=True
    ):
        out.setdefault(str(pid), (str(flag), float(q)))
    return out


def _row(panel: str, **cells: object) -> list[object]:
    return [panel, *(_py(cells.get(c)) for c in COLUMNS[1:])]


def _by_class(rows: Sequence[list[object]]) -> list[list[object]]:
    """Stable sort by class rank: pair rows keep ``pairs.tsv`` order (lead decision 6)."""
    return sorted(rows, key=lambda r: CLASS_ORDER.index(str(r[2])))


def _measured(pairs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(class, |Δ|, kept)``: kept = finite ``delta_log2`` and not ``abundance_imputed``."""
    cls = pairs["class"].astype(str).to_numpy()
    absd = np.abs(_floats(pairs["delta_log2"]))
    return cls, absd, _kept(pairs) & np.isfinite(absd)


def _kept(pairs: pd.DataFrame) -> np.ndarray:
    if "flags" not in pairs.columns:
        return np.ones(len(pairs), dtype=bool)
    return ~_imputed(pairs["flags"])


def edges(abs_delta: np.ndarray) -> np.ndarray:
    """The 31 shared bin edges over [0, max |Δ|], or [0, 1] when the maximum is 0."""
    top = float(abs_delta.max()) if len(abs_delta) else 0.0
    return np.linspace(0.0, top if top > 0 else 1.0, N_BINS + 1)


def data_frame(pairs: pd.DataFrame, flags: pd.DataFrame | None) -> pd.DataFrame:
    """The V11 data TSV (M9-11 "Data TSV" and lead decisions 3, 5, 6)."""
    cls, absd, measured = _measured(pairs)
    pid = pairs["pair_id"].astype(str).to_numpy()
    up, down = _floats(pairs["abundance_log2_a"]), _floats(pairs["abundance_log2_b"])
    lfdr = _floats(pairs["lfdr"])
    lookup = _flag_lookup(flags)
    bins = edges(absd[measured & np.isin(cls, DENSITY_CLASSES)])
    rows: list[list[object]] = []
    for c in CLASS_ORDER:
        sel = measured & (cls == c)
        if c != "DV" and sel.any():
            dens, _ = np.histogram(absd[sel], bins=bins, density=True)
            rows += [_row("a", **{"class": c, "bin": i, "density": d}) for i, d in enumerate(dens)]
        rows.append(_row("a", **{"class": c, "n": int(sel.sum())}))
    for i in np.flatnonzero(np.isfinite(lfdr)):
        flag, q = lookup.get(pid[i], (None, None))
        rows.append(
            _row(
                "b",
                pair_id=pid[i],
                **{"class": cls[i]},
                abs_delta_log2=absd[i],
                lfdr=lfdr[i],
                flag=flag,
                q_value=q,
            )
        )
    same = _kept(pairs) & np.isin(cls, ("OP", "NO")) & np.isfinite(up) & np.isfinite(down)
    c_rows = []
    for i in np.flatnonzero(same):
        flag, q = lookup.get(pid[i], (None, None))
        c_rows.append(
            _row(
                "c",
                pair_id=pid[i],
                **{"class": cls[i]},
                abs_delta_log2=absd[i],
                lfdr=lfdr[i],
                upstream_log2=up[i],
                downstream_log2=down[i],
                flag=flag,
                q_value=q,
            )
        )
    b_rows = [r for r in rows if r[0] == "b"]
    a_rows = [r for r in rows if r[0] == "a"]
    ordered = [*a_rows, *_by_class(b_rows), *_by_class(c_rows)]
    return pd.DataFrame(ordered, columns=list(COLUMNS), dtype=object)


def _num(value: object) -> float:
    v = _py(value)
    return math.nan if v is None else float(str(v))


def medians(pairs: pd.DataFrame) -> dict[str, float]:
    """Median |Δ| of the panel (a) OP and NO pairs (NaN without one)."""
    cls, absd, measured = _measured(pairs)
    return {
        c: float(np.median(absd[measured & (cls == c)]))
        if (measured & (cls == c)).any()
        else math.nan
        for c in ("OP", "NO")
    }


def title(meds: Mapping[str, float], k: int, alpha: float, n: int | None = None) -> str:
    """Median |Δ| of OP and NO (``%.2g``) and the flag count; the "too few" form without one.

    ``k`` counts the ``operon_flags.tsv`` rows with ``q ≤ alpha`` and ``n`` all its rows. When
    some rows miss ``alpha`` (``n > k``) the clause names both, "0 of 1 candidate flags at
    q ≤ 0.1", so a shown flag row with a large q is not read as a flagged pair."""
    if n is not None and n > k:
        flagged = f"{k} of {n} candidate flags at q ≤ {alpha:g}"
    else:
        flagged = f"{k} pairs flagged at q ≤ {alpha:g}"
    if any(math.isnan(meds.get(c, math.nan)) for c in ("OP", "NO")):
        return f"{flagged}; too few measured pairs for densities"
    return f"Median |Δ| OP {meds['OP']:.2g} vs NO {meds['NO']:.2g} log2; {flagged}"


def _empty(ax: Axes, message: str) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, color=GREY)


def _panel_a(ax: Axes, a: pd.DataFrame, bins: np.ndarray) -> None:
    counts = {str(r["class"]): int(_num(r["n"])) for _, r in a[a["bin"].isna()].iterrows()}
    hist = a[a["bin"].notna()]
    if not len(hist):
        _empty(ax, "no measured pairs")
        return
    for c in DENSITY_CLASSES:
        dens = _floats(hist.loc[hist["class"] == c, "density"])
        if len(dens):
            label = f"{c} (n = {counts.get(c, 0)})"
            ax.stairs(dens, bins, color=CLASS_COLOUR[c], linewidth=1.0, label=label)
    ax.set_xlim(0.0, float(bins[-1]))
    ax.set_xlabel("|Δ| partner log2 abundance")
    ax.set_ylabel("density")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title("(a) |Δ| by pair class", loc="left")


def _panel_b(ax: Axes, b: pd.DataFrame, alpha: float) -> None:
    if not len(b):
        _empty(ax, "no lfdr values")
        return
    x, y = _floats(b["abs_delta_log2"]), _floats(b["lfdr"])
    order = np.argsort(x, kind="stable")
    ax.plot(x[order], y[order], marker="o", markersize=3, linewidth=0.8, color=CLASS_COLOUR["NO"])
    ax.axhline(alpha, linestyle="--", color=GREY, linewidth=0.8, label=THRESHOLD_LABEL)
    ax.text(0.02, 0.98, NOTE_B, transform=ax.transAxes, fontsize=6, color=GREY, va="top")
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("|Δ| partner log2 abundance")
    ax.set_ylabel("local fdr")
    ax.legend(frameon=False, loc="center right")  # clear of the alpha line and the note
    ax.set_title("(b) local fdr", loc="left")


def _panel_c(ax: Axes, c: pd.DataFrame, alpha: float) -> None:
    if not len(c):
        _empty(ax, "no measured same-strand pairs")
        return
    x, y = _floats(c["upstream_log2"]), _floats(c["downstream_log2"])
    q = _floats(c["q_value"])
    flag = np.array([str(f) if _py(f) is not None else "" for f in c["flag"]], dtype=object)
    cls = c["class"].astype(str).to_numpy()
    hit = (flag != "") & np.isfinite(q) & (q <= alpha)
    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    if hi <= lo:
        lo, hi = lo - 1.0, hi + 1.0
    ax.plot([lo, hi], [lo, hi], color=style.OKABE_ITO["grey"], linewidth=0.6, linestyle=":")
    for k in ("OP", "NO"):
        sel = (cls == k) & ~hit
        if sel.any():
            ax.scatter(x[sel], y[sel], s=_SMALL, c=CLASS_COLOUR[k], linewidths=0, label=k, zorder=2)
    for f, marker in FLAG_MARKERS.items():
        sel = hit & (flag == f)
        if sel.any():
            ax.scatter(
                x[sel],
                y[sel],
                s=_LARGE,
                marker=marker,
                c=[CLASS_COLOUR[k] for k in cls[sel]],
                edgecolors="#000000",
                linewidths=0.8,
                label=f.replace("_", " "),
                zorder=3,
            )
    ax.set_xlabel("upstream partner log2 abundance")
    ax.set_ylabel("downstream partner log2 abundance")
    handles, labels = ax.get_legend_handles_labels()
    ax.get_figure(root=True).legend(  # type: ignore[union-attr]
        handles,
        [f"(c) {t}" for t in labels],
        frameon=False,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=min(len(labels), 5),
        fontsize=6,
    )
    ax.set_title("(c) partners (OP, NO)", loc="left")


@base.plot("V11", "operon_pairs", needs=("pairs",), optional=("operon_flags",))
def operon_pairs(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """V11: |Δ| densities per pair class, the lfdr curve and the partner scatter."""
    pairs = view.table("pairs")
    assert pairs is not None  # guaranteed by ``needs``
    if not len(pairs):
        return [base.empty_state("V11", "operon_pairs", EMPTY_PAIRS, COLUMNS, "double")]
    flags = view.table("operon_flags")
    alpha = float(cfg.stats.alpha.operon_flags)
    data = data_frame(pairs, flags)
    sub = {p: data[data["panel"] == p] for p in ("a", "b", "c")}
    cls, absd, measured = _measured(pairs)
    bins = edges(absd[measured & np.isin(cls, DENSITY_CLASSES)])
    k = 0 if flags is None else int((_floats(flags["q_value"]) <= alpha).sum())
    fig = style.new_figure("double", 3.0)
    fig.set_label("V11_operon_pairs")
    grid = fig.add_gridspec(1, 3, left=0.07, right=0.98, bottom=0.25, top=0.84, wspace=0.32)
    axes = {p: fig.add_subplot(grid[0, i], label=p) for i, p in enumerate(("a", "b", "c"))}
    _panel_a(axes["a"], sub["a"], bins)
    _panel_b(axes["b"], sub["b"], alpha)
    _panel_c(axes["c"], sub["c"], alpha)
    n = 0 if flags is None else len(flags)
    text = title(medians(pairs), k, alpha, n)
    fig.suptitle(text, y=0.97)
    return [base.Instance(None, fig, data, text)]
