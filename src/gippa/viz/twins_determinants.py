"""Plotly twin of V15 ``determinants`` (M9-20; GIPPA_SPEC.md §10.2 V15): the four panels from the
data TSV, (b) one sub-plot per PD feature; (a) dodged numeric y (circle elastic net, open diamond
LightGBM), reversed, "CI unestimated" marked; the static footnotes as an annotation (decision 5,
S3). Pure: no file I/O."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gippa.determinants import CAVEATS
from gippa.viz.determinants import (
    A_HEADER,
    BLOCKS,
    CAT_COLOUR,
    CAT_LABEL,
    COLUMNS,
    D_HEADER,
    DELTA_FOOT,
    GREY,
    LEGEND,
    MODEL_COLOUR,
    MODELS,
    PD_FEATURES,
    PD_XLABEL,
)
from gippa.viz.interactive import interactive
from gippa.viz.style import OKABE_ITO
from gippa.viz.twins import DOT, ci, finish, floats, fmt, points, rgba, titles

_SYMBOL = {"elastic_net": "circle", "lightgbm": "diamond-open"}

_SPECS = [[{"colspan": 2}, None, {"colspan": 2}, None], [{"colspan": 4}, None, None, None],
          [{}, {}, {}, {}]]  # fmt: skip
_AT = {"a": (1, 1), "d": (1, 3), "c": (2, 1)}


@interactive("V15", COLUMNS)
def v15(data: pd.DataFrame, title: str) -> Any:
    """(a) ΔR² with CI bars per model; (b) PD mean and spread; (c) SHAP points; (d) observed vs.
    predicted by category, the locus tag in the hover text."""
    names = [A_HEADER, D_HEADER, "(c) SHAP (LightGBM, secondary)",
             *(f"(b) {PD_XLABEL[f]}" for f in PD_FEATURES)]  # fmt: skip
    fig = make_subplots(rows=3, cols=4, specs=_SPECS, subplot_titles=names, vertical_spacing=0.1)
    panel = data["panel"].astype(str)
    a = data[panel == "a"]
    blocks = [b for b in BLOCKS if b in set(a["block"].astype(str))]
    for k, m in enumerate(MODELS):
        d = a[a["model"].astype(str) == m]
        x, lo, hi = floats(d["x"]), floats(d["ci_low"]), floats(d["ci_high"])
        y = np.array([blocks.index(str(b)) for b in d["block"]], dtype=float) + (k - 0.5) * 0.3
        hover = [f"{b}{DOT}ΔR² {fmt(v)}{DOT}{ci(lo_, hi_)}"
                 for b, v, lo_, hi_ in zip(d["block"], d["x"], d["ci_low"], d["ci_high"],
                                           strict=True)]  # fmt: skip
        fig.add_trace(points(len(d), x=x, y=y, mode="markers", name=LEGEND[m], text=hover,
                             marker={"color": MODEL_COLOUR[m], "symbol": _SYMBOL[m]},
                             error_x={"type": "data", "symmetric": False, "array": hi - x,
                                      "arrayminus": x - lo}), *_AT["a"])  # fmt: skip
        na = ~(np.isfinite(lo) & np.isfinite(hi))
        if na.any():
            fig.add_trace(go.Scatter(x=x[na], y=y[na], mode="text", text=["CI unestimated"]
                                     * int(na.sum()), textposition="middle right",
                                     textfont={"color": GREY, "size": 9}, showlegend=False,
                                     hoverinfo="skip"), *_AT["a"])  # fmt: skip
    fig.update_yaxes(tickvals=list(range(len(blocks))), ticktext=blocks, autorange="reversed",
                     row=1, col=1)  # fmt: skip
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
    notes = [DELTA_FOOT] if len(a) else []
    fig.add_annotation(text="<br>".join([*notes, CAVEATS[0]]), x=0.0, y=-0.05, xref="paper",
                       yref="paper", xanchor="left", yanchor="top", align="left",
                       showarrow=False, font={"color": GREY, "size": 9})  # fmt: skip
    finish(fig, title, 1000)
    fig.update_layout(margin={"t": 70, "l": 60, "r": 20, "b": 110})
    return fig
