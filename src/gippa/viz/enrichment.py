"""V14 ``enrichment``: which functions are high or low? (GIPPA_SPEC.md §10.2 V14, §7.7.3-§7.7.5,
§9.3 ``enrichment_gsea.tsv`` / ``enrichment_detection.tsv``; ADR-0025/0026/0027/0029/0035; M9-19
and its lead decisions 1-4).

(a) The fgsea running sums of the six sets in ``plots/support/gsea_running`` (M9-17's order),
2 x 3, with the set-member ticks and the leading edge shaded up to the peak; a set flagged
``genomic_cluster`` / ``detection_confounded`` gets a "not calibrated" banner and a dashed outline.
The support's signed max deviation must match the table's ``es`` (to the ``%.6g`` cell rounding,
relative 1e-5, on top of 1e-9): a mismatch is a bug and raises. (b) NES against -log10 q, one
marker per finite row, coloured by collection (fixed hues), flagged rows hollow, redundant rows
small and grey, no error bars (ADR-0027), the ten most significant non-redundant rows labelled in
a column at the right; a q of 0 sits at the cap (``neg_log10_q``), noted in the footnote. (c)
Detection ORA: the rows at q <= alpha, else the 20 smallest q (decision 1), O/E = ``stat`` /
(``stat`` - ``effect``) on a log2 axis coloured by direction, coverage at the right, -log10 q in a
strip with the alpha line; no bars, "CI unestimated" (decision 2). Pure: tables from the RunView.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from gippa.viz import base, style
from gippa.viz.interactive import interactive

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

COLUMNS = ("panel", "collection", "term", "name", "direction", "rank_position", "running_es",
           "in_set", "nes", "q_value", "o_over_e", "coverage", "flags", "redundant")  # fmt: skip
COLLECTIONS: tuple[str, ...] = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
#: fixed per collection (the Okabe-Ito cycle without yellow, weak on white)
_HUES = ("orange", "sky_blue", "bluish_green", "blue", "vermillion", "reddish_purple")
COLLECTION_COLOUR: Mapping[str, str] = {
    c: style.OKABE_ITO[h] for c, h in zip(COLLECTIONS, _HUES, strict=True)
}
DIRECTION_COLOUR: Mapping[str, str] = {"over": style.OKABE_ITO["blue"],
                                       "under": style.OKABE_ITO["vermillion"]}  # fmt: skip
GREY, TEXT = style.OKABE_ITO["grey"], "#555555"
UNCALIBRATED = ("genomic_cluster", "detection_confounded")  # ADR-0025, ADR-0026
TOP_SETS, N_LABELS, N_SMALLEST, N_DRAWN_C = 6, 10, 20, 30
TWIN_ALPHA = 0.05  # the twin sees only the TSV: its alpha line is at the default (decision 6)
FOOT_GSEA = "GSEA is conditional on detection, with no p_detect adjustment"
FOOT_COVERAGE = "coverage = share of eligible genes with ≥ 1 term"
XLABEL_C = "observed / expected detected TUs"
NO_RESOURCES = "No annotation resources installed (gippa setup --go --kegg --cog)"
EMPTY = "No GSEA or detection ORA rows in this run"
NO_ABUNDANCE, NO_SUPPORT = "no abundance input", "running scores not available"
NO_TESTABLE = "no testable sets"
OE_UNDEFINED = "O/E undefined (E = 0)"
LOG10_Q = "\N{MINUS SIGN}log10 q"
_REL = 1e-5  # the %.6g rounding of the support and table cells
_GSEA = ("collection", "term", "name", "es", "nes", "redundant", "coverage", "q_value", "flags")
_ORA = ("collection", "term", "name", "direction", "stat", "effect", "coverage", "ci_low",
        "ci_high", "q_value", "flags")  # fmt: skip


def _flags(v: object) -> list[str]:
    return str(v).split(";") if isinstance(v, str) and v else []


def _f(v: object) -> float:
    return math.nan if v is None or v is pd.NA or pd.isna(v) else float(str(v))


def _bool(v: object) -> bool:
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    return str(v).strip().lower() in ("true", "1", "1.0")


def _order(coll: object) -> int:
    return COLLECTIONS.index(str(coll)) if str(coll) in COLLECTIONS else len(COLLECTIONS)


def _short(text: str, n: int) -> str:
    return text if len(text) <= n else text[: n - 1] + "…"


def _g(v: float, spec: str = ".2g") -> str:
    return format(v, spec) if math.isfinite(v) else "NA"


def _row(panel: str, r: Mapping[str, Any], **kw: Any) -> dict[str, Any]:
    out = dict.fromkeys(COLUMNS)
    out.update(panel=panel, collection=r["collection"], term=r["term"], name=r.get("name"),
               flags=r.get("flags"), **kw)  # fmt: skip
    if panel != "c":
        out["redundant"] = _bool(r.get("redundant"))
    return out


def _key(r: Mapping[str, Any]) -> tuple[Any, ...]:
    rank = r["rank_position"]
    return (r["panel"], _order(r["collection"]), str(r["term"]), str(r["direction"] or ""),
            -1 if rank is None else rank)  # fmt: skip


def max_deviation(running: Sequence[float]) -> float:
    """fgsea's ES: the signed maximum deviation of the running sum from 0."""
    lo, hi = min(0.0, *running), max(0.0, *running)
    return hi if hi > -lo else lo


def running_sets(
    support: pd.DataFrame, gsea: Mapping[tuple[str, str], Mapping[str, Any]]
) -> list[tuple[Mapping[str, Any], pd.DataFrame]]:
    """The first six support sets with their table rows; a set absent from the table or whose
    max deviation differs from its ``es`` raises ``ValueError`` (examples 4, 5)."""
    out = []
    pairs = zip(support["collection"].astype(str), support["term"].astype(str), strict=True)
    keys = dict.fromkeys(pairs)
    for coll, term in list(keys)[:TOP_SETS]:
        row = gsea.get((coll, term))
        if row is None:
            raise ValueError(f"gsea_running set {coll}:{term} is not in enrichment_gsea")
        sub = support[(support["collection"].astype(str) == coll)
                      & (support["term"].astype(str) == term)]  # fmt: skip
        sub = sub.sort_values("rank_position")
        dev, es = max_deviation([float(v) for v in sub["running_es"]]), _f(row["es"])
        if not abs(dev - es) <= 1e-9 + _REL * abs(es):
            raise ValueError(f"gsea_running {coll}:{term}: max deviation {dev:g} != es {es:g}")
        out.append((row, sub))
    return out


def _by_q(ora: pd.DataFrame) -> pd.DataFrame:
    """``ora`` by q, ties broken on collection (enum order), term and direction (N6)."""
    keys = pd.DataFrame({"q": ora["q_value"].astype(float),
                         "c": [_order(c) for c in ora["collection"]],
                         "t": ora["term"].astype(str),
                         "d": ora["direction"].astype(str)}, index=ora.index)  # fmt: skip
    return ora.loc[keys.sort_values(["q", "c", "t", "d"], kind="stable").index]


def ora_selection(ora: pd.DataFrame, level: float) -> pd.DataFrame:
    """The rows at q <= ``level``, else the 20 smallest q (lead decision 1)."""
    q = ora["q_value"].astype(float)
    sig = ora[q <= level]
    return sig if len(sig) else _by_q(ora[q.notna()]).iloc[:N_SMALLEST]


def drawn_c(sel: pd.DataFrame) -> pd.DataFrame:
    """The panel-(c) rows drawn: the 30 smallest q (N5), O/E-descending top to bottom, ties on
    collection (enum order), term and direction; NaN O/E last."""
    top = _by_q(sel).iloc[:N_DRAWN_C]
    keys = pd.DataFrame({"o": top["o_over_e"].astype(float),
                         "c": [_order(c) for c in top["collection"]],
                         "t": top["term"].astype(str),
                         "d": top["direction"].astype(str)}, index=top.index)  # fmt: skip
    idx = keys.sort_values(["o", "c", "t", "d"], ascending=[False, True, True, True],
                           kind="stable", na_position="last").index  # fmt: skip
    return top.loc[idx].reset_index(drop=True)


def o_over_e(stat: float, effect: float) -> float:
    """``stat / (stat - effect)`` (ADR-0035: ``effect`` = k - E[K | N]); NaN when E <= 0."""
    expected = stat - effect
    return stat / expected if expected > 0 else math.nan


def title(gsea: pd.DataFrame, ora: pd.DataFrame, level: float) -> str:
    red = gsea["redundant"].fillna(False).astype(bool)
    n = int(((gsea["q_value"].astype(float) <= level) & ~red).sum())
    m = int((ora["q_value"].astype(float) <= level).sum())
    if not n and not m:
        return f"No set or category passes q ≤ {level:g}"
    sets = "set" if n == 1 else "sets"
    cats = "category" if m == 1 else "categories"
    return (f"{n} {sets} enriched by abundance at q ≤ {level:g} (conditional on detection); "
            f"{m} {cats} detected more or less than expected")  # fmt: skip


def neg_log10_q(q: np.ndarray, level: float) -> tuple[np.ndarray, float]:
    """-log10 q with q = 0 at the cap: the smallest non-zero q, or alpha when that is smaller."""
    pos = q[q > 0]
    cap = -math.log10(min(float(pos.min()) if len(pos) else level, level))
    return np.where(q > 0, -np.log10(np.where(q > 0, q, 1.0)), cap), cap


def _note(fig: Figure, rect: tuple[float, float, float, float], label: str, text: str) -> None:
    ax = fig.add_axes(rect, label=label)
    ax.set_axis_off()
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=TEXT)


def _draw_set(ax: Axes, row: Mapping[str, Any], sub: pd.DataFrame) -> None:
    x = sub["rank_position"].astype(float).to_numpy()
    y = sub["running_es"].astype(float).to_numpy()
    colour = COLLECTION_COLOUR.get(str(row["collection"]), GREY)
    ax.plot(x, y, color=colour, lw=0.8)
    peak = int(np.argmax(np.abs(y)))
    edge = slice(0, peak + 1) if y[peak] > 0 else slice(peak, len(y))
    ax.fill_between(x[edge], y[edge], 0, color=colour, alpha=0.3, lw=0)
    ax.axhline(0, color=GREY, lw=0.5)
    lo, hi = min(0.0, float(y.min())), max(0.0, float(y.max()))
    tick = lo - 0.12 * (hi - lo or 1.0)
    members = x[sub["in_set"].astype(str).str.lower().isin(("true", "1"))]
    ax.plot(members, np.full(len(members), tick), ls="none", marker="|", ms=4, color="#000000")
    ax.set_ylim(tick - 0.08 * (hi - lo or 1.0), hi + 0.08 * (hi - lo or 1.0))
    ax.set_xlim(0.5, len(x) + 0.5)
    name = _short(f"{row['term']} {row['name'] or ''}".strip(), 34)
    bad = [f for f in UNCALIBRATED if f in _flags(row["flags"])]
    ax.set_title(f"{name}\n({row['collection']}) NES {_g(_f(row['nes']))}, q "
                 f"{_g(_f(row['q_value']))}", fontsize=6, pad=13 if bad else 4)  # fmt: skip
    if bad:  # N4: between the title and the plot area
        ax.annotate(f"not calibrated ({', '.join(bad)})", xy=(0.5, 1.0), xycoords="axes fraction",
                    xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
                    fontsize=6, color="#000000",
                    bbox={"facecolor": "#DDDDDD", "edgecolor": "none", "pad": 1})  # fmt: skip
        for spine in ax.spines.values():
            spine.set_linestyle("--")


def _scatter_q(ax: Axes, x: np.ndarray, y: np.ndarray, q: np.ndarray, **kw: Any) -> None:
    """One scatter for q > 0 (``o``) and one for q = 0 (``^``, at the cap; N2). Per-point
    sequences in ``kw`` are split with the points."""
    for mask, marker in ((q != 0, "o"), (q == 0, "^")):
        if not mask.any():
            continue
        part = {k: np.asarray(v, dtype=object)[mask].tolist()
                if isinstance(v, (list, np.ndarray)) else v for k, v in kw.items()}  # fmt: skip
        ax.scatter(x[mask], y[mask], marker=marker, **part)


LABEL_X = 1.03  # (b)'s label column, in axes fraction


def _draw_b(ax: Axes, rows: pd.DataFrame, level: float, height_pt: float) -> float | None:
    """The dot plot; returns the q = 0 cap when one is drawn."""
    q = rows["q_value"].astype(float).to_numpy()
    y, cap = neg_log10_q(q, level)
    x = rows["nes"].astype(float).to_numpy()
    red = np.array([_bool(v) for v in rows["redundant"].fillna(False)], dtype=bool)
    bad = np.array([any(f in UNCALIBRATED for f in _flags(v)) for v in rows["flags"]], dtype=bool)
    colour = [GREY if r else COLLECTION_COLOUR.get(str(c), GREY)
              for c, r in zip(rows["collection"], red, strict=True)]  # fmt: skip
    face = ["none" if b else c for c, b in zip(colour, bad, strict=True)]
    _scatter_q(ax, x, y, q, s=np.where(red, 6.0, 18.0), facecolors=face, edgecolors=colour,
               lw=0.8)  # fmt: skip
    ax.axhline(-math.log10(level), ls="--", color=GREY, lw=0.6)
    span = float(x.max() - x.min()) if len(x) else 1.0
    xlim = (float(x.min()) - 0.1 * span - 0.2, float(x.max()) + 0.1 * span + 0.2)
    ylim = (0.0, max(float(y.max()), -math.log10(level)) * 1.1 + 0.1)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    order = np.lexsort((rows["term"].astype(str), -np.abs(x), q))
    top = sorted([i for i in order if not red[i]][:N_LABELS], key=lambda i: y[i])
    gap = (ylim[1] - ylim[0]) * 8.0 / height_pt
    placed: list[float] = []
    for i in top:
        placed.append(max(y[i], placed[-1] + gap) if placed else y[i])
    shift = max(0.0, (placed[-1] if placed else 0.0) - ylim[1] + gap / 2)
    # N3: leaders in data coordinates (they stay on the label column when the axes is resized):
    # label start, an elbow at the axes edge, the point
    x_label = xlim[0] + LABEL_X * (xlim[1] - xlim[0])
    for i, yl in zip(top, placed, strict=True):
        r = rows.iloc[i]
        ax.plot([x_label, xlim[1], x[i]], [yl - shift, yl - shift, y[i]], color=GREY, lw=0.4,
                clip_on=False)  # fmt: skip
        ax.text(LABEL_X, yl - shift, _short(f"{r['term']} {r['name'] or ''}".strip(), 26),
                transform=ax.get_yaxis_transform(), va="center", ha="left",
                fontsize=6)  # fmt: skip
    present = [c for c in COLLECTIONS if c in set(rows["collection"].astype(str)[~red])]
    handles = [Line2D([], [], ls="none", marker="o", ms=4, color=COLLECTION_COLOUR[c], label=c)
               for c in present]  # fmt: skip
    if bad.any():
        handles.append(Line2D([], [], ls="none", marker="o", ms=4, mfc="none", color=TEXT,
                              label="not calibrated"))  # fmt: skip
    if red.any():
        handles.append(Line2D([], [], ls="none", marker="o", ms=2.5, color=GREY, label="redundant"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=6, handletextpad=0.2)
    ax.set_xlabel("NES")
    ax.set_ylabel(LOG10_Q)
    return cap if (q == 0).any() else None


def _draw_c(ax: Axes, ax_q: Axes, rows: pd.DataFrame, level: float) -> list[str]:
    """O/E per row (log2) and its -log10 q strip; returns the footnotes it needs."""
    notes = []
    oe = rows["o_over_e"].astype(float).to_numpy()
    y = np.arange(len(rows), dtype=float)
    pos = oe[oe > 0]
    floor = float(pos.min()) / 2 if len(pos) else 0.5
    xs = np.where(oe > 0, oe, np.where(oe == 0, floor, np.nan))
    colour = [DIRECTION_COLOUR.get(str(d), GREY) for d in rows["direction"]]
    ax.scatter(xs, y, s=18, c=colour, lw=0)
    ax.axvline(1.0, color=GREY, lw=0.5)
    ax.set_xscale("log", base=2)
    fin = xs[np.isfinite(xs)]
    lo, hi = (float(fin.min()), float(fin.max())) if len(fin) else (0.5, 2.0)
    ax.set_xlim(min(lo, 1.0) / 1.5, max(hi, 1.0) * 1.5)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    labels = [_short(f"{r['term']} {r['name'] or ''}".strip(), 22) + f" ({r['collection']})"
              for _, r in rows.iterrows()]  # fmt: skip
    ax.set_yticks(y, labels)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlabel(XLABEL_C, loc="right")  # clear of the q strip's label
    for yi, cov in zip(y, rows["coverage"].astype(float), strict=True):
        ax.text(1.04, yi, _g(cov, ".2f"), transform=ax.get_yaxis_transform(), va="center",
                fontsize=6)  # fmt: skip
    ax.text(1.04, 1.01, "coverage", transform=ax.transAxes, va="bottom", fontsize=6)
    q = rows["q_value"].astype(float).to_numpy()
    qy, cap = neg_log10_q(q, level)
    _scatter_q(ax_q, qy, y, q, s=10, c=colour, lw=0)
    ax_q.axvline(-math.log10(level), ls="--", color=GREY, lw=0.6)
    ax_q.set_ylim(len(rows) - 0.5, -0.5)
    ax_q.set_yticks([])
    ax_q.set_xlim(0.0, max(cap, -math.log10(level)) * 1.15)
    ax_q.set_xlabel(LOG10_Q)
    if (q == 0).any():  # N1
        notes.append(f"(c) q = 0 drawn at the cap \N{MINUS SIGN}log10 q = {cap:.3g}")
    if (oe == 0).any():
        notes.append("O/E = 0 drawn at the axis floor")
    n_nan = int(np.isnan(oe).sum())
    if n_nan:  # N11
        notes.append(f"{OE_UNDEFINED}: {n_nan} not drawn")
    return notes


def _frames(view: RunView) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame | None]:
    """The two tables (absent: zero rows) and the support file (absent: ``None``)."""
    gsea, ora = view.table("enrichment_gsea"), view.table("enrichment_detection")
    return (pd.DataFrame(columns=_GSEA) if gsea is None else gsea,
            pd.DataFrame(columns=_ORA) if ora is None else ora,
            view.support("gsea_running"))  # fmt: skip


@base.plot("V14", "enrichment", optional=("enrichment_gsea", "enrichment_detection"),
           module="function")  # fmt: skip
def enrichment(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """V14: GSEA running scores, NES vs -log10 q and detection ORA O/E."""
    if not view.has("enrichment_gsea") and not view.has("enrichment_detection"):
        return base.Skip("input_absent:enrichment_gsea")
    gsea, ora, support = _frames(view)
    level = float(cfg.stats.alpha.enrichment)
    if not len(gsea) and not len(ora):
        missing = all(f"layer_missing:{x}" in view.run_flags() for x in ("go", "kegg", "cog"))
        msg = NO_RESOURCES if missing else EMPTY
        return [base.empty_state("V14", "enrichment", msg, COLUMNS, "double")]
    by_key = {(str(r["collection"]), str(r["term"])): r for r in gsea.to_dict("records")}
    sets = running_sets(support, by_key) if len(gsea) and support is not None else []
    nes, q = gsea["nes"].astype(float), gsea["q_value"].astype(float)
    finite = gsea[np.isfinite(nes) & np.isfinite(q)]
    finite = finite.sort_values(["collection", "term"], kind="stable").reset_index(drop=True)
    sel = ora_selection(ora, level).copy()
    sel["o_over_e"] = [o_over_e(_f(a), _f(e)) for a, e in zip(sel["stat"], sel["effect"],
                                                              strict=True)]  # fmt: skip
    text = title(gsea, ora, level)
    no_b = NO_TESTABLE if len(gsea) else NO_ABUNDANCE
    fig = _figure(sets, finite, sel, level, NO_SUPPORT if len(gsea) else NO_ABUNDANCE, no_b)
    fig.suptitle(text.replace("; ", ";\n"), y=1 - 0.08 / fig.get_figheight(), va="top")
    return [base.Instance(None, fig, data_frame(sets, finite, sel), text)]


def data_frame(
    sets: Sequence[tuple[Mapping[str, Any], pd.DataFrame]], finite: pd.DataFrame, sel: pd.DataFrame
) -> pd.DataFrame:
    """``V14_enrichment.tsv``: (a) the support rows, (b) the drawn GSEA rows, (c) every selected
    ORA row (the figure draws the 30 smallest q of them; decision 6); sorted by panel, collection
    (enum order), term, direction, rank; ``redundant`` for the GSEA rows."""
    rows: list[dict[str, Any]] = []
    for r, sub in sets:
        for p in sub.to_dict("records"):
            rows.append(_row("a", r, rank_position=int(p["rank_position"]),
                             running_es=float(p["running_es"]),
                             in_set=str(p["in_set"]).lower() in ("true", "1"),
                             nes=_f(r["nes"]), q_value=_f(r["q_value"])))  # fmt: skip
    for r in finite.to_dict("records"):
        rows.append(_row("b", r, nes=_f(r["nes"]), q_value=_f(r["q_value"]),
                         coverage=_f(r["coverage"])))  # fmt: skip
    for r in sel.to_dict("records"):
        rows.append(_row("c", r, direction=r["direction"], q_value=_f(r["q_value"]),
                         o_over_e=r["o_over_e"], coverage=_f(r["coverage"])))  # fmt: skip
    return pd.DataFrame(sorted(rows, key=_key), columns=list(COLUMNS))


def _footnote(notes: Sequence[str], width: int = 140) -> str:
    """The footnotes, two fixed ones on the first line, the rest wrapped at ``width`` chars."""
    lines = ["; ".join(notes[:2])]
    cur = ""
    for n in notes[2:]:
        if cur and len(cur) + 2 + len(n) > width:
            lines.append(cur)
            cur = n
        else:
            cur = f"{cur}; {n}" if cur else n
    return "\n".join([*lines, *([cur] if cur else [])])


#: inches: (b) from 0.55 (1.9 wide at most, 1.2 at least); (c) ends at 6.15, >= 0.8 wide
B_LEFT, B_WIDTH, B_MIN, C_RIGHT, C_MIN, C_LEFT, GAP = 0.55, 1.9, 1.2, 6.15, 0.8, 5.1, 0.1
#: text widths move by up to ~9 % with the output dpi (glyph hinting, 72-600 dpi): measured
#: widths are padded so the layout holds at any dpi the writer uses
HINTING = 1.1


def _label_widths(fig: Figure, ax_b: Axes | None, ax_c: Axes) -> tuple[float, float]:
    """B1, measured with the Agg renderer: the width (inches) of (b)'s label texts past their
    anchor at ``LABEL_X``, and of (c)'s tick labels left of its axes."""
    canvas = FigureCanvasAgg(fig)
    canvas.draw()  # type: ignore[no-untyped-call]
    renderer = canvas.get_renderer()  # type: ignore[no-untyped-call]
    text_w = 0.0
    if ax_b is not None:
        box = ax_b.get_window_extent(renderer)
        anchor = box.x0 + LABEL_X * box.width
        ends = [t.get_window_extent(renderer).x1 for t in ax_b.texts if t.get_text()]
        text_w = (max([anchor, *ends]) - anchor) / fig.dpi
    x0 = ax_c.get_window_extent(renderer).x0
    starts = [t.get_window_extent(renderer).x0 for t in ax_c.get_yticklabels()
              if t.get_text() and t.get_visible()]  # fmt: skip
    return HINTING * text_w, HINTING * (x0 - min([x0, *starts])) / fig.dpi


def _figure(
    sets: Sequence[tuple[Mapping[str, Any], pd.DataFrame]],
    finite: pd.DataFrame,
    sel: pd.DataFrame,
    level: float,
    no_sets: str,
    no_b: str,
) -> Figure:
    """The three panels (inches from the top: (a) 1.1-3.75, (b) and (c) from 4.45); (c)'s left
    edge follows the measured width of (b)'s label column and of its own (B1)."""
    drawn = drawn_c(sel) if len(sel) else sel
    top_bc, h_bc = 4.45, max(1.7, 0.13 * len(drawn) + 0.1)
    height = top_bc + h_bc + 0.85
    fig = style.new_figure("double", height)
    fig.set_label("V14_enrichment")

    def rect(left: float, top: float, w: float, h: float) -> tuple[float, float, float, float]:
        return (left / style.DOUBLE_IN, 1 - (top + h) / height, w / style.DOUBLE_IN, h / height)

    fig.text(0.07 / style.DOUBLE_IN, 1 - 0.55 / height, "(a) GSEA running scores, top sets",
             fontsize=7, va="top")  # fmt: skip
    if not sets:
        _note(fig, rect(0.55, 1.1, 6.5, 2.65), "a_note", no_sets)
    for i, (r, sub) in enumerate(sets):
        ax = fig.add_axes(rect(0.55 + 2.35 * (i % 3), 1.1 + 1.65 * (i // 3), 1.85, 1.0), label="a")
        _draw_set(ax, r, sub)
        if i // 3 == 1 or i + 3 >= len(sets):
            ax.set_xlabel("rank position (TUs)")
        if i % 3 == 0:
            ax.set_ylabel("running ES (enrichment score)")
    notes = [FOOT_GSEA, FOOT_COVERAGE]
    ax_b = None
    if len(finite):
        ax_b = fig.add_axes(rect(B_LEFT, top_bc, B_WIDTH, h_bc), label="b")
        ax_b.set_title("(b) NES vs \N{MINUS SIGN}log10 q", loc="left", fontsize=7)
        cap = _draw_b(ax_b, finite, level, h_bc * 72)
        if cap is not None:
            notes.append(f"q = 0 drawn at the cap \N{MINUS SIGN}log10 q = {cap:.3g}")
    else:
        _note(fig, rect(0.55, top_bc, 2.9, h_bc), "b", no_b)
    if len(drawn):
        ax_c = fig.add_axes(rect(C_LEFT, top_bc, C_RIGHT - C_LEFT, h_bc), label="c")
        ax_q = fig.add_axes(rect(6.62, top_bc, 0.45, h_bc), label="c_q")
        c_title = ax_c.set_title("(c) detection ORA", loc="left", fontsize=7)
        notes += _draw_c(ax_c, ax_q, drawn, level)
        if drawn["ci_low"].isna().any() or drawn["ci_high"].isna().any():
            notes.append("ORA points: CI unestimated (no bars)")
        if len(sel) > len(drawn):  # N5
            notes.append(f"(c) shows the {len(drawn)} smallest q; "
                         f"+{len(sel) - len(drawn)} more in the TSV")  # fmt: skip
        text_w, tick_w = _label_widths(fig, ax_b, ax_c)
        b_w = B_WIDTH
        if ax_b is not None:
            c_left = B_LEFT + LABEL_X * b_w + text_w + GAP + tick_w
            if C_RIGHT - c_left < C_MIN:  # narrow (b) before (c)
                b_w = max(B_MIN, b_w - (C_MIN - (C_RIGHT - c_left)) / LABEL_X)
                ax_b.set_position(rect(B_LEFT, top_bc, b_w, h_bc))
            c_left = B_LEFT + LABEL_X * b_w + text_w + GAP + tick_w
        else:
            c_left = max(C_LEFT, B_LEFT + tick_w)
        c_w = max(C_RIGHT - c_left, 0.3)
        ax_c.set_position(rect(c_left, top_bc, c_w, h_bc))
        c_title.set_x(-tick_w / c_w)  # over the label column, clear of "coverage"
    else:
        _note(fig, rect(4.0, top_bc, 3.0, h_bc), "c", "no detection ORA rows")
    fig.text(0.07 / style.DOUBLE_IN, 0.42 / height, _footnote(notes), fontsize=6, color=TEXT,
             va="top")  # fmt: skip
    return fig


@interactive("V14", COLUMNS)
def enrichment_twin(data: pd.DataFrame, text: str) -> Any:
    """(a) the running sums, (b) NES vs -log10 q by collection (q = 0 at the cap, the alpha line
    at the default 0.05, redundant rows small and grey), (c) O/E (log) by direction, O/E order."""
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(rows=1, cols=3, subplot_titles=["(a) running scores",
                        "(b) NES vs \N{MINUS SIGN}log10 q", "(c) detection ORA (O/E)"])  # fmt: skip
    panel = {p: data[data["panel"].astype(str) == p] for p in "abc"}
    for (coll, term), sub in panel["a"].groupby(["collection", "term"], sort=False):
        fig.add_trace(go.Scatter(x=sub["rank_position"].astype(float).tolist(),
                                 y=sub["running_es"].astype(float).tolist(), mode="lines",
                                 name=f"{term} ({coll})",
                                 line={"color": COLLECTION_COLOUR.get(str(coll), GREY)}),
                      row=1, col=1)  # fmt: skip
    b = panel["b"]
    q_all = b["q_value"].astype(float).to_numpy()
    y_all = neg_log10_q(q_all, TWIN_ALPHA)[0] if len(b) else q_all
    for coll in COLLECTIONS:
        m = (b["collection"].astype(str) == coll).to_numpy()
        if not m.any():
            continue
        sub, q, y = b[m], q_all[m], y_all[m]
        red = [_bool(v) for v in sub["redundant"].fillna(False)]
        hover = [f"{t} {n} · NES {_g(v, '.3g')} · q {_g(w, '.3g')} · {';'.join(_flags(f))}"
                 for t, n, v, w, f in zip(sub["term"], sub["name"], sub["nes"].astype(float), q,
                                          sub["flags"], strict=True)]  # fmt: skip
        sym = [("triangle-up" if w == 0 else "circle")
               + ("-open" if any(x in UNCALIBRATED for x in _flags(f)) else "")
               for w, f in zip(q, sub["flags"], strict=True)]  # fmt: skip
        fig.add_trace(go.Scatter(x=sub["nes"].astype(float).tolist(), y=y.tolist(),
                                 mode="markers", name=coll, hovertext=hover, hoverinfo="text",
                                 marker={"color": [GREY if r else COLLECTION_COLOUR[coll]
                                                   for r in red],
                                         "size": [5 if r else 9 for r in red],
                                         "symbol": sym}),
                      row=1, col=2)  # fmt: skip
    if len(b):
        fig.add_hline(y=-math.log10(TWIN_ALPHA), line_dash="dash", line_color=GREY,
                      line_width=1, row=1, col=2)  # fmt: skip
    c = panel["c"]
    if len(c):
        c = c.assign(_oe=c["o_over_e"].astype(float))  # bottom to top: O/E ascending
        c = c.sort_values("_oe", kind="stable", na_position="first")
        labels = [f"{t} ({k}, {d})" for t, k, d in zip(c["term"], c["collection"], c["direction"],
                                                      strict=True)]  # fmt: skip
        hover = [f"{lab} · O/E {_g(o, '.3g')} · coverage {_g(v, '.2f')} · CI unestimated"
                 for lab, o, v in zip(labels, c["_oe"], c["coverage"].astype(float),
                                      strict=True)]  # fmt: skip
        fig.add_trace(go.Scatter(x=c["_oe"].tolist(), y=labels,
                                 mode="markers", name="detection ORA", hovertext=hover,
                                 hoverinfo="text", marker={"color": [DIRECTION_COLOUR.get(
                                     str(d), GREY) for d in c["direction"]]}),
                      row=1, col=3)  # fmt: skip
        fig.update_yaxes(categoryorder="array", categoryarray=labels, row=1, col=3)
    fig.update_xaxes(title_text="rank position (TUs)", row=1, col=1)
    fig.update_xaxes(title_text="NES", row=1, col=2)
    fig.update_yaxes(title_text=LOG10_Q, row=1, col=2)
    fig.update_xaxes(title_text=XLABEL_C, type="log", row=1, col=3)
    fig.update_layout(template=style.plotly_template(), title={"text": text}, height=420)
    return fig
