"""Plotly twin of V15 ``determinants`` (M9-20; GIAP_SPEC.md §10.2 V15): the four panels from the
data TSV, (b) one sub-plot per PD feature. Pure: no file I/O."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from giap.viz.determinants import (
    CAT_COLOUR,
    CAT_LABEL,
    COLUMNS,
    LEGEND,
    MODEL_COLOUR,
    MODELS,
    PD_FEATURES,
    PD_XLABEL,
)
from giap.viz.interactive import interactive
from giap.viz.style import OKABE_ITO
from giap.viz.twins import DOT, ci, finish, floats, fmt, points, rgba, titles

_SPECS = [[{"colspan": 2}, None, {"colspan": 2}, None], [{"colspan": 4}, None, None, None],
          [{}, {}, {}, {}]]  # fmt: skip
_AT = {"a": (1, 1), "d": (1, 3), "c": (2, 1)}


@interactive("V15", COLUMNS)
def v15(data: pd.DataFrame, title: str) -> Any:
    """(a) ΔR² with CI bars per model; (b) PD mean and spread; (c) SHAP points; (d) observed vs.
    predicted by category, the locus tag in the hover text."""
    names = ["(a) block ΔR²", "(d) observed vs. predicted", "(c) SHAP (LightGBM, secondary)",
             *(f"(b) {PD_XLABEL[f]}" for f in PD_FEATURES)]  # fmt: skip
    fig = make_subplots(rows=3, cols=4, specs=_SPECS, subplot_titles=names, vertical_spacing=0.1)
    panel = data["panel"].astype(str)
    a = data[panel == "a"]
    for m in MODELS:
        d = a[a["model"].astype(str) == m]
        x, lo, hi = floats(d["x"]), floats(d["ci_low"]), floats(d["ci_high"])
        hover = [f"{b}{DOT}ΔR² {fmt(v)}{DOT}{ci(lo_, hi_)}"
                 for b, v, lo_, hi_ in zip(d["block"], d["x"], d["ci_low"], d["ci_high"],
                                           strict=True)]  # fmt: skip
        fig.add_trace(points(len(d), x=x, y=list(d["block"]), mode="markers", name=LEGEND[m],
                             marker={"color": MODEL_COLOUR[m]}, text=hover,
                             error_x={"type": "data", "symmetric": False, "array": hi - x,
                                      "arrayminus": x - lo}), *_AT["a"])  # fmt: skip
    titles(fig, {"row": 1, "col": 1}, "ΔR² (drop-one-block, CV)", None)
    b = data[(panel == "b") & data["category"].astype(str).eq("pd_mean")]
    for k, f in enumerate(PD_FEATURES):
        g = b[b["feature"].astype(str) == f]
        x, lo, hi = floats(g["x"]), floats(g["ci_low"]), floats(g["ci_high"])
        fig.add_trace(go.Scatter(x=[*x, *x[::-1]], y=[*lo, *hi[::-1]], fill="toself",
                                 fillcolor=rgba(OKABE_ITO["sky_blue"], 0.5), mode="lines",
                                 line={"width": 0}, showlegend=False, hoverinfo="skip"),
                      3, k + 1)  # fmt: skip
        fig.add_trace(go.Scatter(x=x, y=floats(g["y"]), mode="lines", showlegend=False,
                                 line={"color": OKABE_ITO["blue"]}), 3, k + 1)  # fmt: skip
        ylab = "PD (log2 abundance)" if k == 0 else None
        titles(fig, {"row": 3, "col": k + 1}, PD_XLABEL[f], ylab)
    c = data[panel == "c"]
    hover = [f"{t}{DOT}{f}{DOT}SHAP {fmt(s)}" for t, f, s in zip(c["locus_tag"], c["feature"],
                                                                 c["x"], strict=True)]  # fmt: skip
    dots = {"color": floats(c["ci_low"]), "size": 3, "colorscale": "Viridis"}
    fig.add_trace(points(len(c), x=floats(c["x"]), y=floats(c["y"]), mode="markers", text=hover,
                         showlegend=False, marker=dots), *_AT["c"])  # fmt: skip
    feats = list(dict.fromkeys(c["feature"].astype(str)))
    fig.update_yaxes(tickvals=list(np.arange(len(feats))[::-1]), ticktext=feats, row=2, col=1)
    titles(fig, {"row": 2, "col": 1}, "SHAP (log2 abundance)", None)
    d = data[panel == "d"]
    for cat, colour in CAT_COLOUR.items():
        s = d[d["category"].astype(str) == cat]
        hover = [f"{t}{DOT}predicted {fmt(x)}{DOT}observed {fmt(y)}"
                 for t, x, y in zip(s["locus_tag"], s["x"], s["y"], strict=True)]  # fmt: skip
        fig.add_trace(points(len(s), x=floats(s["x"]), y=floats(s["y"]), mode="markers",
                             name=CAT_LABEL[cat], text=hover, marker={"color": colour, "size": 4}),
                      *_AT["d"])  # fmt: skip
    titles(fig, {"row": 1, "col": 3}, "predicted log2 abundance (out-of-fold)",
           "observed log2 abundance")  # fmt: skip
    return finish(fig, title, 1000)
