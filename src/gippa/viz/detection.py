"""V7 ``obs_vs_exp`` (the O - E forest) and V8 ``cog_abundance`` (M9-06; GIPPA_SPEC.md §10.2).

Pure: tables come from a ``RunView`` and the plots return ``Instance`` objects; nothing is written
here. V7 draws the M4 estimates and CIs as they are tabled (``detection_summary`` ``basic`` rows,
``regions_stats`` ``o_minus_e``/``all`` rows; ADR-0018: the raw O - E with its CI, the q value
from the studentised test). V8 draws measured log2 abundances only (M9-index decision 13).
"""

from __future__ import annotations

import math
import statistics
import textwrap
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import pandas as pd
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from gippa.io.cog_categories import cog_descriptions
from gippa.viz.base import Instance, Skip, empty_state, plot, plural
from gippa.viz.style import COG_SUPER, COG_SUPER_COLOUR, REGION, cog_colour, new_figure

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

V07_COLUMNS = (
    "block", "group", "region_type", "n", "observed", "expected", "o_minus_e",
    "ci_low", "ci_high", "q_value", "flags",
)  # fmt: skip
V08_COLUMNS = (
    "category", "super_category", "locus_tag", "abundance_log2", "detected",
    "median_log2", "n_genes", "n_detected", "n_imputed",
)  # fmt: skip
V07_XLABEL = "observed \N{MINUS SIGN} expected detected genes"
V08_YLABEL = "log2 abundance"
NO_ABUNDANCE = "No measured abundance; detected fraction per COG category"
BASIC_BLOCKS = ("replicon", "feature_class", "cog_category")
BLACK = "#000000"
NONE_GREY = "#999999"
MAX_STRIP = 50  # a violin above this many measured values (§10.2 V8)
MIN_HIGHEST = 5  # measured genes a category needs before the title calls its median the highest
_WRAP = {"single": 50, "double": 105}  # characters per title line at 8 pt


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _num(v: object) -> float | None:
    return None if _na(v) else float(str(v))


def _true(v: object) -> bool:
    return not _na(v) and bool(v)


def _flags(v: object) -> list[str]:
    return [] if _na(v) or not str(v) else str(v).split(";")


def _frame(rows: Sequence[Sequence[object]], columns: Sequence[str]) -> pd.DataFrame:
    cols = {c: pd.Series([r[i] for r in rows], dtype=object) for i, c in enumerate(columns)}
    return pd.DataFrame(cols, columns=list(columns))


def _pct(x: float | None, d: int | None) -> str:
    if x is None or not d:
        return "NA"
    return f"{int(100 * x / d + 0.5)} %"


def _figure(width: Literal["single", "double"], height: float, label: str, title: str) -> Figure:
    fig = new_figure(width, height, layout="constrained")
    fig.set_label(label)
    fig.suptitle(textwrap.fill(title, _WRAP[width]), fontsize=8)
    return fig


# --- V7 ---------------------------------------------------------------------------------------


def _basic_rows(summary: pd.DataFrame, no_cog: bool) -> list[list[object]]:
    b = summary[summary["section"].astype(str) == "basic"]
    value = {
        (str(r["group_type"]), str(r["group"]), str(r["metric"])): _num(r["value"])
        for r in b.to_dict("records")
    }
    rows: list[list[object]] = []
    for gtype in BASIC_BLOCKS:
        if gtype == "cog_category" and no_cog:
            continue
        sub = b[(b["group_type"].astype(str) == gtype) & (b["metric"].astype(str) == "o_minus_e")]
        for r in sub.to_dict("records"):
            group = str(r["group"])
            n_genes = value.get((gtype, group, "n_genes"))
            n = int(n_genes) if n_genes is not None else _num(r["n"])
            obs, exp = (
                value.get((gtype, group, "n_detected")),
                value.get((gtype, group, "expected")),
            )
            ome = _num(r["value"])
            if ome is None and exp == 0:
                ome = obs if obs is not None else 0.0
            flags = None if _na(r["flags"]) or not str(r["flags"]) else str(r["flags"])
            rows.append([gtype, group, None, None if n is None else int(n), obs, exp, ome,
                         _num(r["ci_low"]), _num(r["ci_high"]), None, flags])  # fmt: skip
    return rows


def _region_rows(
    regions: pd.DataFrame, types: Sequence[str]
) -> tuple[list[list[object]], list[str | None]]:
    keep = (regions["stat"].astype(str) == "o_minus_e") & (regions["subset"].astype(str) == "all")
    recs = regions[keep].to_dict("records")

    def key(r: Mapping[str, object]) -> tuple[int, str, bool, float, str]:
        t, q = str(r["region_type"]), _num(r["q_value"])
        idx = types.index(t) if t in types else len(types)
        return (idx, t, q is None, 0.0 if q is None else q, str(r["region_id"]))

    rows: list[list[object]] = []
    labels: list[str | None] = []
    for r in sorted(recs, key=key):
        obs, exp, ome, q = (_num(r[c]) for c in ("observed", "expected", "effect", "q_value"))
        if ome is None and exp == 0:
            ome = obs if obs is not None else 0.0
        n = _num(r["n_eligible"])
        flags = None if _na(r["flags"]) or not str(r["flags"]) else str(r["flags"])
        rows.append(["region", str(r["region_id"]), str(r["region_type"]),
                     None if n is None else int(n), obs, exp, ome,
                     _num(r["ci_low"]), _num(r["ci_high"]), q, flags])  # fmt: skip
        label = None if _na(r["label"]) or not str(r["label"]) else str(r["label"])
        labels.append(f"q = {q:.2g}" if q is not None else f"{label} (untested)" if label else None)
    return rows, labels


def v07_title(data: pd.DataFrame, has_regions: bool, alpha: float) -> str:
    """The V7 title from its data TSV (index decision 5; M9-06 lead decisions 1, 2)."""
    recs = data.to_dict("records")
    if has_regions:
        tested = [r for r in recs if r["block"] == "region" and r["q_value"] is not None]
        best = None
        for r in tested:
            if r["q_value"] <= alpha and (best is None or r["q_value"] < best["q_value"]):
                best = r
        if best is not None:
            return (
                f"{best['group']} ({best['region_type']}): {_pct(best['observed'], best['n'])} "
                f"detected vs {_pct(best['expected'], best['n'])} expected "
                f"(q = {best['q_value']:.2g})"
            )
        le = "\N{LESS-THAN OR EQUAL TO}"
        k = len(tested)
        return f"No region departs from detectability at q {le} {alpha:g}; {k} regions tested"
    top = None
    for r in recs:
        lo, hi, ome = r["ci_low"], r["ci_high"], r["o_minus_e"]
        if ome is None or _ci_unestimated(r) or lo <= 0 <= hi:
            continue  # only a row that carries its error estimate is headlined (decision 5)
        if top is None or abs(ome) > abs(top["o_minus_e"]):
            top = r
    if top is None:
        return "Detection matches detectability in every group with a CI"
    return (
        f"{top['group']}: O \N{MINUS SIGN} E = {top['o_minus_e']:.1f} "
        f"(95 % CI {top['ci_low']:.1f} to {top['ci_high']:.1f})"
    )


def _ci_unestimated(r: Mapping[str, object]) -> bool:
    lo, hi = r["ci_low"], r["ci_high"]
    return lo is None or hi is None or lo == hi or "ci_unestimated" in _flags(r["flags"])


def _xlim(recs: Sequence[Mapping[str, object]]) -> tuple[float, float]:
    """The estimates, CI ends and 0, padded 8 %, at least 2 genes wide."""
    xs = [0.0] + [float(v) for r in recs for c in ("o_minus_e", "ci_low", "ci_high")
                  if isinstance(v := r[c], float | int) and math.isfinite(v)]  # fmt: skip
    lo, hi = min(xs), max(xs)
    pad = max(0.08 * (hi - lo), (2.0 - (hi - lo)) / 2, 0.0)
    return lo - pad, hi + pad


def _forest(data: pd.DataFrame, notes: Sequence[str | None], title: str) -> Figure:
    recs = data.to_dict("records")
    # lead decision 6: double above 20 rows; also when an "(untested)" label would leave a
    # single-width forest no room (decision recorded in the M9-06 report)
    wide = len(recs) > 20 or any(n is not None and n.endswith("(untested)") for n in notes)
    width: Literal["single", "double"] = "double" if wide else "single"
    ys: list[float] = []
    y, prev = 0.0, None
    for r in recs:
        if prev is not None and r["block"] != prev:
            y -= 0.6
        ys.append(y)
        y, prev = y - 1.0, r["block"]
    fig = _figure(width, 1.2 + 0.17 * (1.0 - ys[-1]), "V07_obs_vs_exp", title)
    ax = fig.add_subplot()
    ax.axvline(0.0, color="#777777", lw=0.6, zorder=0)
    unest = [_ci_unestimated(r) for r in recs]
    note_dx = 62 if any(unest) else 4  # the right-hand columns: "CI unestimated", then q / label
    for r, yv, note, no_ci in zip(recs, ys, notes, unest, strict=True):
        flags = _flags(r["flags"])
        colour = REGION.get(str(r["region_type"]), BLACK) if r["block"] == "region" else BLACK
        lo, hi, ome = r["ci_low"], r["ci_high"], r["o_minus_e"]
        if lo is not None and hi is not None:
            ax.plot([lo, hi], [yv, yv], color=colour, lw=0.9, solid_capstyle="butt")
        if ome is not None:
            face = "none" if "small_region" in flags else colour
            ax.plot([ome], [yv], marker="o", ms=3.5, ls="none", mfc=face, mec=colour, mew=0.8)
        right = [("CI unestimated", 4)] if no_ci else []
        right += [(note, note_dx)] if note is not None else []
        for text, dx in right:
            ax.annotate(text, (1.0, yv), xycoords=("axes fraction", "data"), xytext=(dx, 0),
                        textcoords="offset points", va="center", fontsize=6)  # fmt: skip
    labels = [f"{r['group']} (n = {'NA' if r['n'] is None else r['n']})" for r in recs]
    ax.set_yticks(ys, labels=labels)
    ax.set_ylim(ys[-1] - 0.7, 0.7)
    ax.set_xlim(*_xlim(recs))
    ax.set_xlabel(V07_XLABEL)
    ax.spines[["top", "right"]].set_visible(False)
    return fig


@plot("V07", "obs_vs_exp", needs=("detection_summary",), optional=("regions_stats",))
def obs_vs_exp(view: RunView, cfg: GippaConfig) -> list[Instance] | Skip:
    """§10.2 V7: the O - E forest of the basic groups and the regions."""
    summary = view.table("detection_summary")
    if summary is None:
        return Skip("input_absent:detection_summary")
    regions = view.table("regions_stats")
    rows = _basic_rows(summary, "layer_missing:cog" in view.run_flags())
    notes: list[str | None] = [None] * len(rows)
    if regions is not None:
        region_rows, labels = _region_rows(regions, list(cfg.context.region_types))
        rows, notes = rows + region_rows, notes + labels
    if not rows:
        msg = "No detection summary or region rows in this run"
        return [empty_state("V07", "obs_vs_exp", msg, V07_COLUMNS, "single")]
    data = _frame(rows, V07_COLUMNS)
    title = v07_title(data, regions is not None, cfg.stats.alpha.context)
    return [Instance(suffix=None, figure=_forest(data, notes, title), data=data, title=title)]


# --- V8 ---------------------------------------------------------------------------------------


@dataclass
class _Cat:
    letter: str
    n_genes: int = 0
    n_detected: int = 0
    n_imputed: int = 0
    genes: list[tuple[str, float | None, bool]] = field(default_factory=list)

    @property
    def values(self) -> list[float]:
        return [v for _, v, _ in self.genes if v is not None]

    @property
    def median(self) -> float | None:
        vals = self.values
        return float(statistics.median(vals)) if vals else None

    @property
    def super_category(self) -> str:
        return "none" if self.letter == "none" else COG_SUPER.get(self.letter, "poorly")

    @property
    def colour(self) -> str:
        return NONE_GREY if self.letter == "none" else cog_colour(self.letter)


def _letters(cog: object) -> list[str]:
    letters = sorted({c for c in ("" if _na(cog) else str(cog).upper()) if c.isalpha()})
    return letters if letters and str(cog).lower() != "none" else ["none"]


def _categories(genes: pd.DataFrame) -> list[_Cat]:
    """Per-category counts and measured values, in drawing order (M9-06 V8 rows)."""
    cats: dict[str, _Cat] = {}
    for r in genes.to_dict("records"):
        if _true(r["partial"]):
            continue
        det, imp = _true(r["detected"]), _true(r["abundance_imputed"])
        v = _num(r["abundance_log2"])
        drawn = v if det and not imp and v is not None and math.isfinite(v) else None
        for letter in _letters(r["cog_category"]):
            c = cats.setdefault(letter, _Cat(letter))
            c.n_genes, c.n_detected, c.n_imputed = (
                c.n_genes + 1,
                c.n_detected + det,
                c.n_imputed + imp,
            )
            c.genes.append((str(r["locus_tag"]), drawn, det))

    def key(c: _Cat) -> tuple[bool, bool, float, str]:
        m = c.median
        return (c.letter == "none", m is None, 0.0 if m is None else -m, c.letter)

    return sorted(cats.values(), key=key)


def _v08_data(cats: Sequence[_Cat]) -> pd.DataFrame:
    rows: list[list[object]] = []
    for c in cats:
        sup = c.super_category
        rows.append([c.letter, sup, None, None, None, c.median, c.n_genes, c.n_detected,
                     c.n_imputed])  # fmt: skip
        for tag, v, det in sorted(c.genes, key=lambda g: g[0]):
            rows.append([c.letter, sup, tag, v, det, None, None, None, None])
    return _frame(rows, V08_COLUMNS)


def v08_title(cats: Sequence[_Cat], flags: frozenset[str]) -> str:
    """The V8 title (M9-06 lead decisions 4, 5)."""
    no_cog = all(c.letter == "none" for c in cats)
    tail = "no COG annotation" + (" (layer_missing:cog)" if "layer_missing:cog" in flags else "")
    top = next((c for c in cats if c.values), None)
    if top is None:
        return f"{NO_ABUNDANCE}; {tail}" if no_cog else NO_ABUNDANCE
    m, n = top.median, len(top.values)
    if top.letter == "none":
        lead = f"Median abundance {m:.1f} log2 ({plural(n, 'gene')})"
        return f"{lead}; {tail}" if no_cog else f"{lead} in genes without a COG letter only"
    best = next((c for c in cats if c.letter != "none" and len(c.values) >= MIN_HIGHEST), None)
    if best is None:  # K-12 survey: no "highest" claim on fewer than MIN_HIGHEST genes
        most = max((c for c in cats if c.letter != "none"), key=lambda c: len(c.values))
        return (f"No COG category has {MIN_HIGHEST} or more measured genes (most: "
                f"{most.letter}, {plural(len(most.values), 'gene')})")  # fmt: skip
    m, n = best.median, len(best.values)
    desc = cog_descriptions().get(best.letter, best.letter)
    return (f"{best.letter} ({desc}) has the highest median abundance among categories "
            f"with ≥ {MIN_HIGHEST} genes ({m:.1f} log2, "
            f"{plural(n, 'gene')})")  # fmt: skip


def _jitter(n: int) -> list[float]:
    """Deterministic strip offsets: a golden-ratio sequence over the locus-tag rank."""
    return [0.0] if n == 1 else [((i * 0.6180339887) % 1.0 - 0.5) * 0.5 for i in range(n)]


def _abundance_figure(cats: Sequence[_Cat], title: str) -> Figure:
    fig = _figure("double", 3.0, "V08_cog_abundance", title)
    ax = fig.add_subplot()
    measured = any(c.values for c in cats)
    for x, c in enumerate(cats):
        if not measured:
            frac = c.n_detected / c.n_genes if c.n_genes else 0.0
            ax.bar([x], [frac], width=0.7, color=c.colour)
            continue
        genes = sorted((g for g in c.genes if g[1] is not None), key=lambda g: g[0])
        vals = [float(v) for _, v, _ in genes if v is not None]
        if len(vals) > MAX_STRIP:
            parts = ax.violinplot([vals], positions=[x], widths=0.7, showextrema=False,
                                  showmedians=True)  # fmt: skip
            bodies: list[PolyCollection] = list(parts["bodies"])  # type: ignore[call-overload]
            for body in bodies:
                body.set_facecolor(c.colour)
                body.set_edgecolor(c.colour)
                body.set_alpha(0.6)
            parts["cmedians"].set_color(BLACK)
        elif vals:
            xs = [x + d for d in _jitter(len(vals))]
            ax.scatter(xs, vals, s=7, color=c.colour, linewidths=0)
    labels = [f"{c.n_detected}/{c.n_genes}" for c in cats]
    # rotate the n/N labels when they would not fit side by side (6 pt digits ~ 3.6 pt wide)
    slot_pt = 0.8 * fig.get_figwidth() * 72 / max(len(cats), 1)
    turn = 90 if max(len(t) for t in labels) * 3.6 + 2 > slot_pt else 0
    for x, text in enumerate(labels):
        ax.annotate(text, (x, 1.0), xycoords=("data", "axes fraction"), xytext=(0, 2),
                    textcoords="offset points", ha="center", va="bottom", fontsize=6,
                    rotation=turn)  # fmt: skip
    ax.set_xticks(range(len(cats)), labels=[c.letter for c in cats])
    ax.set_xlim(-0.6, len(cats) - 0.4)
    ax.set_xlabel("COG category")
    if measured:
        ax.set_ylabel(V08_YLABEL)
    else:
        ax.set_ylabel("detected fraction")
        ax.set_ylim(0.0, 1.0)
    sups = list(dict.fromkeys(c.super_category for c in cats))
    handles = [
        Patch(color=NONE_GREY if s == "none" else COG_SUPER_COLOUR[s], label=s) for s in sups
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False,
              fontsize=6)  # fmt: skip
    ax.spines[["top", "right"]].set_visible(False)
    return fig


@plot("V08", "cog_abundance", needs=("genes_tsv",))
def cog_abundance(view: RunView, cfg: GippaConfig) -> list[Instance] | Skip:
    """§10.2 V8: log2 abundance per COG category, sorted by median, with detected fractions."""
    genes = view.table("genes_tsv")
    if genes is None:
        return Skip("input_absent:genes_tsv")
    cats = _categories(genes)
    if not cats:
        msg = "No non-partial genes in this run"
        return [empty_state("V08", "cog_abundance", msg, V08_COLUMNS, "double")]
    title = v08_title(cats, view.run_flags())
    fig = _abundance_figure(cats, title)
    return [Instance(suffix=None, figure=fig, data=_v08_data(cats), title=title)]
