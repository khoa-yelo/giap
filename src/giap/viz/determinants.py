# ruff: noqa: RUF001, RUF003  (Δ, ², ± and – are the spec's own characters)
"""V15 ``determinants``: what explains abundance? (GIAP_SPEC.md §10.2 V15, §7.8.2, §9.3; M9-20;
ADR-0031..0034). Draws the M8 ``determinants_*`` tables and the M9-17 support files
``determinants_pd`` / ``determinants_shap``; computes no statistic.

(a) drop-one-block ΔR² with its CI, elastic net (primary, ADR-0033) beside LightGBM; (b) partial
dependence of ``cai``, ``tai``, ``tir_ostir`` and ``promotech_max`` ("promoter score", lead
decision 2), the ``pd_low`` to ``pd_high`` spread shaded, a rug of the model rows' values; (c) the
SHAP beeswarm of the top 20 features, jittered and coloured by rank within the feature; (d)
observed vs. LightGBM out-of-fold predicted, unexplained genes and ribosomal references marked
(model genes binned on a log count scale above 20,000 rows). The title names the elastic net's
largest block whose repeat t-interval has both ends finite and lies wholly above 0, ties in
BLOCKS order (decisions 1, 4, 5 and 6); the interval is named uncalibrated (ADR-0031).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.patches import Patch

from giap.determinants import CAVEATS
from giap.viz import base, style

if TYPE_CHECKING:
    from giap.config import GiapConfig
    from giap.io.run_view import RunView

COLUMNS = ("panel", "model", "block", "feature", "x", "y", "ci_low", "ci_high", "locus_tag",
           "category", "flags")  # fmt: skip
BLOCKS = ("transcription", "translation", "codon", "protein", "context")
MODELS = ("elastic_net", "lightgbm")
LEGEND = {"elastic_net": "elastic net (primary, ADR-0033)",
          "lightgbm": "LightGBM (secondary check)"}  # fmt: skip
MARKER = {"elastic_net": ("o", True), "lightgbm": ("D", False)}  # (marker, filled)
MODEL_COLOUR = {"elastic_net": style.OKABE_ITO["blue"],
                "lightgbm": style.OKABE_ITO["bluish_green"]}  # fmt: skip
PD_FEATURES = ("cai", "tai", "tir_ostir", "promotech_max")
PD_XLABEL = {"cai": "CAI (index)", "tai": "tAI (index)", "tir_ostir": "OSTIR TIR (au)",
             "promotech_max": "promoter score (promotech_max)"}  # fmt: skip
#: the layer a PD feature comes from (decision 4: fallback when no run flag names it)
PD_LAYER = {"cai": "codon", "tai": "codon", "tir_ostir": "ostir", "promotech_max": "promoters"}
SHAP_TOP = 20
MAX_DRAWN = 20_000  # the M9-01 element rule
_C_EXTRA = 1  # panel (c)'s zero line counts towards the rule too
UNEXPLAINED = "abundance_not_explained_by_sequence"
CAT_COLOUR = {"model": style.OKABE_ITO["blue"], "unexplained": style.OKABE_ITO["vermillion"],
              "ribosomal_reference": style.OKABE_ITO["orange"]}  # fmt: skip
CAT_LABEL = {"model": "model genes",
             "unexplained": "unexplained (> 2 SD of LightGBM residuals)",
             "ribosomal_reference": "ribosomal reference (not in R²)"}  # fmt: skip
GREY = "#555555"
UNI_FOOT = "univariate CIs uncalibrated (ADR-0032)"
DELTA_FOOT = "ΔR² intervals uncalibrated: CV repeats share data (ADR-0031)"
A_HEADER = "(a) block ΔR² (t-interval over CV repeats, uncalibrated)"
D_HEADER = "(d) observed vs. LightGBM out-of-fold"
HEX_BAR, HEX_PROXY = "genes per bin (log count)", "model genes (binned)"
MINUS = "\u2212"


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _py(v: object) -> object:
    if _na(v):
        return None
    if isinstance(v, np.generic):
        v = v.item()
    return None if isinstance(v, float) and math.isnan(v) else v


def _floats(col: Sequence[object] | pd.Series) -> np.ndarray:
    return np.array([math.nan if _na(v) else float(str(v)) for v in col], dtype=float)


def _flags(v: object) -> list[str]:
    return [] if _na(v) or not str(v) else str(v).split(";")


def _row(panel: str, **cells: object) -> list[object]:
    return [panel, *(_py(cells.get(c)) for c in COLUMNS[1:])]


def _f2(v: float) -> str:
    """Two decimals with a Unicode minus (decision 5, S2), but at least two significant digits
    (``base.sig``) below 0.1 so a small non-zero value never reads 0.00; NaN reads
    "unestimated"."""
    if math.isnan(v):
        return "unestimated"
    out = f"{v:.2f}" if v == 0 or abs(v) >= 0.1 else base.sig(v)
    return out.replace("-", MINUS)


# ---------------------------------------------------------------- data


def _metric(model: pd.DataFrame, name: str, metric: str) -> pd.Series | None:
    sel = model[(model["model"] == name) & (model["metric"] == metric)]
    return None if not len(sel) else sel.iloc[0]


def too_few(model: pd.DataFrame) -> bool:
    return any("too_few_rows" in _flags(f) for f in model["flags"])


def model_rows(features: pd.DataFrame) -> pd.Series:
    """§7.8.2 rule 1: detected, measured, non-imputed abundance, not ribosomal."""

    def flag(c: str) -> pd.Series:
        return features[c].fillna(False).astype(bool)

    measured = features["abundance_log2"].notna()
    return flag("detected") & ~flag("abundance_imputed") & ~flag("ribosomal_protein") & measured


def counts(features: pd.DataFrame | None, attribution: pd.DataFrame | None) -> tuple[int, str]:
    """Model rows (detected, non-imputed, not ribosomal; decision 4) and their distinct TUs."""
    if features is None:
        att = [] if attribution is None else list(attribution["flags"])
        return sum("cai_reference" not in _flags(f) for f in att), "NA"
    keep = model_rows(features)
    return int(keep.sum()), str(features.loc[keep, "tu_id"].nunique())


def cv_text(model: pd.DataFrame, name: str = "elastic_net") -> str:
    """``<mean> ± <sd>`` CV R² of ``name``; "unestimated" (no ± sd) when the mean is NA."""
    rows = (_metric(model, name, m) for m in ("cv_r2_mean", "cv_r2_sd"))
    mean, sd = (math.nan if r is None else float(_floats([r["value"]])[0]) for r in rows)
    if math.isnan(mean):
        return "unestimated"
    return f"{_f2(mean)} ± {_f2(sd)}" if not math.isnan(sd) else f"{_f2(mean)} (SD unestimated)"


def title(model: pd.DataFrame, n_tus: tuple[int, str] | None = None) -> str:
    """The finding: the elastic net's CV R² and its largest block ΔR² whose repeat t-interval is
    finite and lies wholly above 0 (decisions 4-6; ties in BLOCKS order)."""
    if n_tus is not None:
        return f"Too few genes for the determinants model ({n_tus[0]} rows, {n_tus[1]} TUs)"
    head = f"Elastic net CV R² {cv_text(model)}; "
    d = model[(model["model"] == "elastic_net") & (model["metric"] == "delta_r2")]
    best: tuple[float, str, float, float] | None = None
    for b in BLOCKS:
        g = d[d["feature_or_block"] == b]
        for v, lo, hi in zip(_floats(g["value"]), _floats(g["ci_low"]), _floats(g["ci_high"]),
                             strict=True):  # fmt: skip
            above = math.isfinite(lo) and math.isfinite(hi) and lo > 0 and math.isfinite(v)
            if above and (best is None or v > best[0]):
                best = (v, b, lo, hi)
    if best is None:
        return head + "no block's ΔR² interval lies above 0"
    v, b, lo, hi = best
    return head + (f"largest block ΔR²: {b} {_f2(v)} (repeat t-interval {_f2(lo)}–{_f2(hi)}, "
                   "uncalibrated)")  # fmt: skip


def _cap(n: int, k: int) -> np.ndarray:
    """At most ``k`` of ``n`` sorted positions, evenly spaced (deterministic)."""
    if n <= k:
        return np.arange(n)
    return np.unique(np.linspace(0, n - 1, k).round().astype(int))


def _rank(vals: np.ndarray, tags: np.ndarray) -> np.ndarray:
    out = np.empty(len(vals))
    out[np.lexsort((tags, vals))] = np.arange(len(vals))
    return out / max(len(vals) - 1, 1)


def rows_a(model: pd.DataFrame) -> list[list[object]]:
    d = model[model["metric"] == "delta_r2"]
    rows = []
    for m in MODELS:
        for b in BLOCKS:
            for _, r in d[(d["model"] == m) & (d["feature_or_block"] == b)].iterrows():
                rows.append(_row("a", model=m, block=b, x=r["value"], ci_low=r["ci_low"],
                                 ci_high=r["ci_high"], flags=r["flags"]))  # fmt: skip
    return rows


def rows_b(pdf: pd.DataFrame, features: pd.DataFrame | None) -> list[list[object]]:
    rows = []
    for f in PD_FEATURES:
        g = pdf[pdf["feature"] == f].sort_values("grid_value", kind="stable")
        rows += [_row("b", model="lightgbm", feature=f, x=r["grid_value"], y=r["pd_mean"],
                      ci_low=r["pd_low"], ci_high=r["pd_high"], category="pd_mean")
                 for _, r in g.iterrows()]  # fmt: skip
        if len(g) and features is not None and f in features.columns:
            rug = features[model_rows(features) & features[f].notna()]
            rug = rug.iloc[np.lexsort((rug["locus_tag"].astype(str), _floats(rug[f])))]
            rug = rug.iloc[_cap(len(rug), MAX_DRAWN // len(PD_FEATURES))]
            rows += [_row("b", model="lightgbm", feature=f, x=v, locus_tag=t, category="rug")
                     for t, v in zip(rug["locus_tag"], rug[f], strict=True)]  # fmt: skip
    return rows


def top_features(shap: pd.DataFrame) -> list[str]:
    m = shap.assign(a=shap["shap"].abs()).groupby("feature")["a"].mean()
    return sorted(m.index.astype(str), key=lambda f: (-m[f], f))[:SHAP_TOP]


def rows_c(shap: pd.DataFrame) -> list[list[object]]:
    feats = top_features(shap)
    rows = []
    for j, f in enumerate(feats):
        g = shap[shap["feature"].astype(str) == f]
        tags, x = g["locus_tag"].astype(str).to_numpy(), _floats(g["shap"])
        order = np.lexsort((tags, x))[_cap(len(g), (MAX_DRAWN - _C_EXTRA) // len(feats))]
        colour = _rank(_floats(g["feature_value"]), tags)
        for r, i in enumerate(order):  # jitter by rank: the golden-ratio sequence in ±0.35
            y = len(feats) - 1 - j + 0.7 * ((r * 0.6180339887) % 1.0 - 0.5)
            rows.append(_row("c", model="lightgbm", feature=f, x=x[i], y=round(y, 6),
                             ci_low=colour[i], locus_tag=tags[i]))  # fmt: skip
    return rows


def _category(flags: list[str]) -> str:
    if "cai_reference" in flags:
        return "ribosomal_reference"
    return "unexplained" if UNEXPLAINED in flags else "model"


def rows_d(att: pd.DataFrame | None) -> list[list[object]]:
    if att is None:
        return []
    att = att[att["predicted_log2"].notna()]
    att = att.iloc[np.lexsort((att["locus_tag"].astype(str), _floats(att["predicted_log2"])))]
    return [_row("d", model="lightgbm", x=r["predicted_log2"], y=r["abundance_log2"],
                 locus_tag=r["locus_tag"], category=_category(_flags(r["flags"])),
                 flags=r["flags"]) for _, r in att.iterrows()]  # fmt: skip


# ---------------------------------------------------------------- drawing


def _note(ax: Axes, text: str, label: str) -> None:
    ax.set_axis_off()
    ax.set_title(label, loc="left")
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, color=GREY,
            fontsize=6, wrap=True)  # fmt: skip


def _panel_a(ax: Axes, a: pd.DataFrame) -> None:
    blocks = [b for b in BLOCKS if b in set(a["block"])]
    for k, m in enumerate(MODELS):
        d = a[a["model"] == m]
        if not len(d):
            continue
        y = np.array([blocks.index(b) for b in d["block"]]) + (k - 0.5) * 0.3
        x, lo, hi = _floats(d["x"]), _floats(d["ci_low"]), _floats(d["ci_high"])
        ok = np.isfinite(lo) & np.isfinite(hi)
        col, (mk, filled) = MODEL_COLOUR[m], MARKER[m]
        ax.errorbar(x[ok], y[ok], xerr=[x[ok] - lo[ok], hi[ok] - x[ok]], fmt="none", ecolor=col,
                    elinewidth=0.8, capsize=2)  # fmt: skip
        ax.plot(x, y, linestyle="none", marker=mk, markersize=4, color=col,
                markerfacecolor=col if filled else "none", label=LEGEND[m])  # fmt: skip
        for xi, yi in zip(x[~ok], y[~ok], strict=True):
            ax.annotate("CI unestimated", (xi, yi), xytext=(5, 0), textcoords="offset points",
                        fontsize=6, va="center", color=GREY)  # fmt: skip
    ax.axvline(0.0, color=GREY, linewidth=0.6, linestyle=":")
    ax.set_yticks(range(len(blocks)), blocks)
    ax.set_ylim(len(blocks) + 0.5, -0.6)  # room for the legend below the last block
    ax.set_xlabel("ΔR² (drop-one-block, CV)")
    ax.locator_params(axis="x", nbins=4)
    ax.legend(frameon=False, loc="best", fontsize=6)
    ax.set_title(A_HEADER, loc="left", x=-0.3)  # over the block labels: clear of (b)


def _panel_b(fig: Any, cell: Any, pdf: pd.DataFrame | None, b: pd.DataFrame,
             flags: frozenset[str], spread: str) -> None:  # fmt: skip
    if pdf is None:
        text = "not available (no determinants_pd support file)"
        _note(fig.add_subplot(cell, label="b"), text, "(b) partial dependence")
        return
    grid = cell.subgridspec(2, 2, hspace=0.9, wspace=0.45)
    for i, f in enumerate(PD_FEATURES):
        ax = fig.add_subplot(grid[i // 2, i % 2], label=f"b_{f}")
        head = f"(b) partial dependence, LightGBM (secondary); {spread}" if i == 0 else ""
        d = b[b["feature"] == f]
        if not len(d):
            layer = next((fl for fl in sorted(flags) if fl == f"layer_missing:{PD_LAYER[f]}"),
                         f"layer_missing:{PD_LAYER[f]}")  # fmt: skip
            _note(ax, f"{PD_XLABEL[f].split(' (')[0]}:\nnot available ({layer})", head)
            continue
        g, rug = d[d["category"] == "pd_mean"], d[d["category"] == "rug"]
        x = _floats(g["x"])
        ax.fill_between(x, _floats(g["ci_low"]), _floats(g["ci_high"]), color=style.OKABE_ITO[
            "sky_blue"], alpha=0.5, linewidth=0)  # fmt: skip
        ax.plot(x, _floats(g["y"]), color=style.OKABE_ITO["blue"], linewidth=1.0)
        if len(rug):
            ax.plot(_floats(rug["x"]), np.zeros(len(rug)), linestyle="none", marker="|",
                    markersize=3, color="#000000", alpha=0.4,
                    transform=ax.get_xaxis_transform())  # fmt: skip
        ax.locator_params(axis="both", nbins=3)
        ax.set_xlabel(PD_XLABEL[f])
        if i % 2 == 0:
            ax.set_ylabel("PD (log2 abundance)")
        if head:
            ax.set_title(head, loc="left", fontsize=7, wrap=True)  # wraps at the figure edge


def _panel_c(ax: Axes, c: pd.DataFrame) -> None:
    feats = list(dict.fromkeys(c["feature"]))
    sc = ax.scatter(_floats(c["x"]), _floats(c["y"]), c=_floats(c["ci_low"]), cmap=style.CONTINUOUS,
                    vmin=0.0, vmax=1.0, s=2, linewidths=0)  # fmt: skip
    ax.set_yticks([len(feats) - 1 - j for j in range(len(feats))], feats)
    ax.set_ylim(-0.6, len(feats) - 0.4)
    ax.axvline(0.0, color=GREY, linewidth=0.6, linestyle=":")
    ax.set_xlabel("SHAP (log2 abundance)")
    ax.set_title("(c) SHAP, top 20; LightGBM (secondary)", loc="left")
    cax = ax.inset_axes((1.03, 0.0, 0.03, 0.35))
    bar = ax.get_figure(root=True).colorbar(sc, cax=cax, ticks=[0.0, 1.0])  # type: ignore[union-attr]
    bar.ax.set_yticklabels(["low", "high"], fontsize=6)
    bar.set_label("feature value rank", fontsize=6)


def _panel_d(ax: Axes, d: pd.DataFrame, r2: str) -> None:
    """Above 20,000 rows the model genes are binned (log counts, decision 5 B4); the flagged
    genes stay points."""
    x, y = _floats(d["x"]), _floats(d["y"])
    cat = d["category"].astype(str).to_numpy()
    binned = len(d) > MAX_DRAWN
    proxy: list[Any] = []
    if binned:
        sel = cat == "model"
        hb = ax.hexbin(x[sel], y[sel], gridsize=40, bins="log", cmap=style.CONTINUOUS, mincnt=1,
                       linewidths=0)  # fmt: skip
        fig, box = ax.get_figure(root=True), ax.get_position()  # a figure axes: findable
        cax = fig.add_axes((box.x1 + 0.01, box.y0, 0.012, 0.35 * box.height),  # type: ignore[union-attr]
                           label="d_colourbar")  # fmt: skip
        bar = fig.colorbar(hb, cax=cax)  # type: ignore[union-attr]
        bar.ax.tick_params(labelsize=6)
        bar.set_label(HEX_BAR, fontsize=6)
        swatch = mpl.colormaps[style.CONTINUOUS](0.6)
        proxy.append(Patch(facecolor=swatch, linewidth=0, label=HEX_PROXY))
    for k, colour in CAT_COLOUR.items():
        sel = cat == k
        if sel.any() and not (binned and k == "model") and sel.sum() <= MAX_DRAWN:
            cross = k == "ribosomal_reference"
            ax.scatter(
                x[sel],
                y[sel],
                s=14 if cross else 5,
                marker="x" if cross else "o",
                c=colour,
                linewidths=0.8 if cross else 0,
                label=CAT_LABEL[k],
                zorder=3,
            )
    lo, hi = float(np.nanmin([x.min(), y.min()])), float(np.nanmax([x.max(), y.max()]))
    ax.plot([lo, hi], [lo, hi], color=style.OKABE_ITO["grey"], linewidth=0.6, linestyle=":")
    ax.text(1.0, 1.02, r2, transform=ax.transAxes, ha="right", va="bottom", fontsize=6)
    ax.set_xlabel("predicted log2 abundance (LightGBM out-of-fold)")
    ax.set_ylabel("observed log2 abundance")
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=[*proxy, *handles], frameon=False, loc="best", fontsize=6)
    ax.set_title(D_HEADER, loc="left", pad=14)  # above the R² line


def footnotes(model: pd.DataFrame) -> list[str]:
    notes = []
    p = _metric(model, "elastic_net", "null_p")
    if p is not None and not _na(p["value"]):
        cal = " (uncalibrated: valid only for exchangeable TUs, ADR-0034)"
        notes.append(f"null p {float(p['value']):.2g}"
                     + (cal if "null_uncalibrated" in _flags(p["flags"]) else ""))  # fmt: skip
    if (model["metric"] == "delta_r2").any():
        notes.append(DELTA_FOOT)
    if (model["model"] == "univariate").any():
        notes.append(UNI_FOOT)
    return [*notes, CAVEATS[0]]


@base.plot(
    "V15",
    "determinants",
    needs=("determinants_model",),
    optional=("determinants_attribution", "determinants_features"),
    layers=("abundance",),
)
def determinants(view: RunView, cfg: GiapConfig) -> list[base.Instance] | base.Skip:
    """V15: block ΔR², partial dependence, SHAP beeswarm and observed vs. predicted."""
    model = view.table("determinants_model")
    assert model is not None  # guaranteed by ``needs``
    att, features = view.table("determinants_attribution"), view.table("determinants_features")
    if too_few(model):  # M8 then writes NA predictions only: no panel has data
        text = title(model, counts(features, att))
        return [base.empty_state("V15", "determinants", text, COLUMNS, "double")]
    text = title(model)
    pdf, shap = view.support("determinants_pd"), view.support("determinants_shap")
    rows = rows_a(model)
    rows += [] if pdf is None else rows_b(pdf, features)
    rows += [] if shap is None or not len(shap) else rows_c(shap)
    data = pd.DataFrame([*rows, *rows_d(att)], columns=list(COLUMNS), dtype=object)
    sub = {p: data[data["panel"] == p] for p in "abcd"}
    fig = style.new_figure("double", 6.6)
    fig.set_label("V15_determinants")
    outer = fig.add_gridspec(2, 2, left=0.17, right=0.95, bottom=0.15, top=0.91, hspace=0.42,
                             wspace=0.5)  # fmt: skip
    _panel_a(fig.add_subplot(outer[0, 0], label="a"), sub["a"])
    n_fits = cfg.determinants.outer_folds * cfg.determinants.repeats
    spread = f"spread over {n_fits} CV models (not a CI)"
    _panel_b(fig, outer[0, 1], pdf, sub["b"], view.run_flags(), spread)
    ax_c = fig.add_subplot(outer[1, 0], label="c")
    if len(sub["c"]):
        _panel_c(ax_c, sub["c"])
    else:
        _note(ax_c, "not available (no determinants_shap support file)", "(c) SHAP")
    ax_d = fig.add_subplot(outer[1, 1], label="d")
    if len(sub["d"]):
        r2 = (f"CV R² {cv_text(model, 'lightgbm')} (LightGBM, plotted); "
              f"elastic net {cv_text(model)} (primary)")  # fmt: skip
        _panel_d(ax_d, sub["d"], r2)
    else:
        _note(ax_d, "no out-of-fold predictions", D_HEADER)
    for i, note in enumerate(footnotes(model)):
        fig.text(0.01, 0.075 - 0.02 * i, note, fontsize=6, color=GREY)
    fig.suptitle(text, y=0.985)
    return [base.Instance(None, fig, data, text)]
