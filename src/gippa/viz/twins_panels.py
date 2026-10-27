"""Plotly twins of V10, V11, V17 and V18 (M9-15; GIPPA_SPEC.md §10.2): the multi-panel plots.
Pre-binned rows are bars at the TSV's bins, never Plotly's binning (behaviour 3). Pure."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gippa.viz.correlogram import COLUMNS as V10_COLUMNS
from gippa.viz.correlogram import (
    NA_FLAGS,
    NULL_LABEL,
    NULL_STYLE,
    NULLS,
    STATISTICS,
    YLABEL,
    plotted,
)
from gippa.viz.diagnostics import COLUMNS as V18_COLUMNS
from gippa.viz.diagnostics import FAMILIES, MAX_DRAWN, N_BINS, OUTSIDE_METRIC
from gippa.viz.interactive import interactive
from gippa.viz.operon_pairs import CLASS_COLOUR, DENSITY_CLASSES, NOTE_B
from gippa.viz.operon_pairs import COLUMNS as V11_COLUMNS
from gippa.viz.qc import COLUMNS as V17_COLUMNS
from gippa.viz.qc import IMPUTED, KIND_COLOUR, LEVEL_COLOUR, LEVELS, PLACED, UNMAPPED
from gippa.viz.style import OKABE_ITO
from gippa.viz.twins import (
    BLACK,
    DOT,
    MUTED,
    Meta,
    alpha,
    contigs,
    finish,
    flag_list,
    floats,
    fmt,
    isna,
    note,
    rgba,
    titles,
)

# ---------------------------------------------------------------- V10


def _ref(row: int, col: int, cols: int) -> str:
    n = (row - 1) * cols + col
    return "" if n == 1 else str(n)


@interactive("V10", V10_COLUMNS)
def v10(data: pd.DataFrame, title: str, meta: Meta = None) -> Any:
    """Per replicon: Moran's I and BB joins against the lag with a filled envelope per null;
    a filled marker where the TU-block q is at most alpha (``meta``; else ``significant``)."""
    level, known = alpha(meta, "context")
    order = contigs(data, meta)
    cols, rows = min(3, len(order)), 2 * -(-len(order) // 3)
    names = [""] * (rows * cols)
    for k, contig in enumerate(order):
        for i, panel in enumerate("ab"):
            names[(2 * (k // 3) + i) * cols + k % 3] = f"({panel}) {contig}"
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=names, vertical_spacing=0.12)
    shown = data.copy()  # flagged rows draw nothing, whatever their value (M9-09)
    na = [not NA_FLAGS.isdisjoint(flag_list(f)) for f in shown["flags"]]
    shown.loc[na, ["value", "null_mean", "null_lo", "null_hi", "effect", "q_value"]] = math.nan
    shown["y"] = plotted(list(shown["statistic"]), *(floats(shown[c]) for c in
                                                     ("value", "effect", "null_mean")))  # fmt: skip
    first = True
    for k, contig in enumerate(order):
        for i, stat in enumerate(STATISTICS):
            at = {"row": 2 * (k // 3) + i + 1, "col": k % 3 + 1}
            d = shown[(shown["contig"].astype(str) == contig) & (shown["statistic"] == stat)]
            titles(fig, at, "lag (genes)", YLABEL[stat])
            if not np.isfinite(floats(d["y"])).any():
                flags = ", ".join(dict.fromkeys(f for v in d["flags"] for f in flag_list(v)))
                note(fig, f"not estimable ({flags})", _ref(at["row"], at["col"], cols))
                continue
            for null in NULLS:
                s = d[d["null_model"] == null]
                band, a, line = NULL_STYLE[null]
                x, lo, hi, mean = (floats(s[c]) for c in ("lag", "null_lo", "null_hi", "null_mean"))
                ok = np.isfinite(x) & np.isfinite(lo) & np.isfinite(hi)
                fig.add_trace(go.Scatter(x=[*x[ok], *x[ok][::-1]], y=[*lo[ok], *hi[ok][::-1]],
                                         fill="toself", fillcolor=rgba(band, a), mode="lines",
                                         line={"width": 0, "color": band}, hoverinfo="skip",
                                         name=f"{NULL_LABEL[null]} 95 % envelope",
                                         legendgroup=null, showlegend=first), **at)  # fmt: skip
                ok = np.isfinite(x) & np.isfinite(mean)
                fig.add_trace(go.Scatter(x=x[ok], y=mean[ok], mode="lines", showlegend=False,
                                         line={"dash": "dash", "color": line, "width": 1},
                                         hoverinfo="skip"), **at)  # fmt: skip
            tu, gl = (d[d["null_model"] == n].set_index("lag") for n in ("tu_block", "gene_label"))
            pts: list[tuple[float, float, str, str]] = []
            for lag in sorted({int(v) for v in d["lag"].dropna()}):
                src = tu if lag in tu.index and not isna(tu.at[lag, "y"]) else gl
                if lag not in src.index or isna(src.at[lag, "y"]):
                    continue
                r, q = src.loc[lag], tu.at[lag, "q_value"] if lag in tu.index else math.nan
                own = lag in tu.index and str(tu.at[lag, "significant"]).lower() == "true"
                sig = (not isna(q) and float(q) <= level) if known else own
                raw = (f"raw BB {fmt(r['value'])}",) if stat == "join_count_bb" else ()
                text = DOT.join((f"lag {lag}", f"{YLABEL[stat]} {fmt(r['y'])}", *raw,
                                 f"TU-block envelope {fmt(r['null_lo'])} to {fmt(r['null_hi'])}",
                                 f"q {fmt(q)}"))  # fmt: skip
                sym = "circle" if sig else "circle-open"
                pts.append((float(lag), float(r["y"]), sym, text))
            fig.add_trace(go.Scatter(x=[p[0] for p in pts], y=[p[1] for p in pts],
                                     mode="lines+markers", hovertext=[p[3] for p in pts],
                                     hoverinfo="text", name="observed", showlegend=first,
                                     line={"color": BLACK, "width": 1},
                                     marker={"color": BLACK, "size": 6,
                                             "symbol": [p[2] for p in pts]}),
                          **at)  # fmt: skip
            first = False
    return finish(fig, title, 120 + 260 * rows, legend={"orientation": "h", "y": -0.08})


# ---------------------------------------------------------------- V11

#: the §7.4.5 flags -> Plotly marker symbols (the static's shapes)
FLAG_SYMBOL = {"internal_promoter_supported": "triangle-up", "unexplained_discordance": "diamond",
               "internal_terminator_supported": "triangle-down", "other": "x",
               "conflict_terminator_in_operon": "square", "novel_operon_candidate": "star",
               "coabundant_despite_terminator": "cross"}  # fmt: skip
DX = "|Δ| partner log2 abundance"
PARTNERS = ("upstream partner log2 abundance", "downstream partner log2 abundance")


@interactive("V11", V11_COLUMNS)
def v11(data: pd.DataFrame, title: str, meta: Meta = None) -> Any:
    """(a) the TSV's density bins per class (at the bin index, lead decision 6), (b) the lfdr
    line, (c) the partners, flags with q at most alpha as symbols."""
    level, known = alpha(meta, "operon_flags")
    fig = make_subplots(rows=1, cols=3, subplot_titles=("(a) |Δ| by pair class", "(b) local fdr",
                                                        "(c) partners (OP, NO)"))  # fmt: skip
    a = data[data["panel"] == "a"]
    counts = {str(c): fmt(n) for c, n in zip(a["class"], a["n"], strict=True) if not isna(n)}
    for c in DENSITY_CLASSES:
        h = a[a["bin"].notna() & (a["class"] == c)]
        if len(h):
            fig.add_trace(go.Bar(x=floats(h["bin"]), y=floats(h["density"]), opacity=0.55,
                                 name=f"{c} (n = {counts.get(c, '0')})",
                                 marker={"color": CLASS_COLOUR[c]}), row=1, col=1)  # fmt: skip
    b = data[data["panel"] == "b"]
    b = b.iloc[np.argsort(floats(b["abs_delta_log2"]), kind="stable")]
    cells = zip(b["pair_id"], b["abs_delta_log2"], b["lfdr"], strict=True)
    text = [f"{p}{DOT}|Δ| {fmt(x)}{DOT}lfdr {fmt(y)}" for p, x, y in cells]
    fig.add_trace(go.Scatter(x=floats(b["abs_delta_log2"]), y=floats(b["lfdr"]), hovertext=text,
                             mode="lines+markers", hoverinfo="text", name="lfdr",
                             line={"color": CLASS_COLOUR["NO"]}), row=1, col=2)  # fmt: skip
    fig.add_hline(y=level, line={"dash": "dash", "color": MUTED, "width": 1}, row=1, col=2,
                  annotation_text=f"flag threshold (q) {level:g}" if known else None)  # fmt: skip
    note(fig, NOTE_B, "2", 0.3, 0.97)
    c = data[data["panel"] == "c"]
    q, flag = floats(c["q_value"]), np.array(["" if isna(f) else str(f) for f in c["flag"]])
    hit = (flag != "") & np.isfinite(q) & (q <= level)
    groups = [(k, (c["class"] == k).to_numpy() & ~hit, "circle") for k in ("OP", "NO")]
    groups += [(f.replace("_", " "), hit & (flag == f), s) for f, s in FLAG_SYMBOL.items()]
    for name, sel, symbol in [g for g in groups if g[1].any()]:
        s, big = c[sel], symbol != "circle"
        text = [DOT.join((str(p), str(k), f.replace("_", " ") or "no flag", f"q {fmt(v)}"))
                for p, k, f, v in zip(s["pair_id"], s["class"], flag[sel], s["q_value"],
                                      strict=True)]  # fmt: skip
        colour = [CLASS_COLOUR.get(str(k), MUTED) for k in s["class"]]
        fig.add_trace(go.Scatter(x=floats(s["upstream_log2"]), y=floats(s["downstream_log2"]),
                                 mode="markers", name=name, hovertext=text, hoverinfo="text",
                                 marker={"symbol": symbol, "size": 11 if big else 6,
                                         "color": colour,
                                         "line": {"color": BLACK, "width": 1 if big else 0}}),
                      row=1, col=3)  # fmt: skip
    for col, (xt, yt) in enumerate(((DX, "density"), (DX, "local fdr"), PARTNERS), start=1):
        titles(fig, {"row": 1, "col": col}, xt, yt)
    return finish(fig, title, 400, barmode="overlay", bargap=0,
                  legend={"orientation": "h", "y": -0.25})  # fmt: skip


# ---------------------------------------------------------------- V17


@interactive("V17", V17_COLUMNS)
def v17(data: pd.DataFrame, title: str) -> Any:
    """Five panels as the static (a-e) and the warning list as a table trace."""
    fig = make_subplots(rows=2, cols=3, vertical_spacing=0.18, horizontal_spacing=0.1,
                        specs=[[{}, {}, {}], [{}, {}, {"type": "table"}]],
                        subplot_titles=("(a) preflight checks", "(b) contig lengths",
                                        "(c) peptide lengths", "(d) abundance",
                                        "(e) missingness (non-partial genes)",
                                        "warnings and errors"))  # fmt: skip
    p = {k: data[data["panel"] == k] for k in "abcde"}
    lvl = p["a"][p["a"]["group"] == "level"].set_index("key")["value"]
    fig.add_trace(go.Bar(x=list(LEVELS), y=[float(lvl.get(x, 0)) for x in LEVELS],
                         marker={"color": [LEVEL_COLOUR[x] for x in LEVELS]}, showlegend=False),
                  row=1, col=1)  # fmt: skip
    msgs = [str(m) for m in p["a"].loc[p["a"]["group"] == "message", "key"]] or ["no warnings"]
    fig.add_trace(go.Table(header={"values": ["rule"]}, cells={"values": [msgs]}), row=2, col=3)
    for panel, at in (("b", (1, 2)), ("e", (2, 2))):  # horizontal bars, top to bottom
        h = p[panel]
        x = floats(h["value"]) if panel == "b" else 100 * np.nan_to_num(floats(h["value"]))
        colour = [KIND_COLOUR.get(str(g), KIND_COLOUR["unknown"]) for g in h["group"]]
        text = [DOT.join((str(k), str(g), fmt(v))) for k, g, v in
                zip(h["key"], h["group"], h["value"], strict=True)]  # fmt: skip
        fig.add_trace(go.Bar(y=[str(k) for k in h["key"]], x=x, orientation="h", hovertext=text,
                             marker={"color": colour if panel == "b" else MUTED},
                             hoverinfo="text", showlegend=False), row=at[0], col=at[1])  # fmt: skip
        fig.update_yaxes(autorange="reversed", row=at[0], col=at[1])
    stacks = (("c", (1, 3), (("placed", PLACED), ("unmapped", UNMAPPED))),
              ("d", (2, 1), (("measured", PLACED), ("imputed", IMPUTED))))  # fmt: skip
    for panel, at, groups in stacks:
        edges = np.unique(floats(p[panel]["key"]))  # (d): the bins' left edges
        width = float(np.min(np.diff(edges))) if panel == "d" and len(edges) > 1 else 1.0
        for group, hue in groups:
            g = p[panel][p[panel]["group"] == group]
            if len(g):
                fig.add_trace(go.Bar(x=floats(g["key"]), y=floats(g["value"]), name=group,
                                     width=width, offset=0 if panel == "d" else None,
                                     marker={"color": hue}), row=at[0], col=at[1])  # fmt: skip
    for (r, c), xt, yt in (((1, 1), None, "checks"), ((1, 2), "length (bp)", None),
                           ((1, 3), "length (aa)", "distinct peptides"),
                           ((2, 1), "log2 abundance", "genes"),
                           ((2, 2), "missing (%)", None)):  # fmt: skip
        titles(fig, {"row": r, "col": c}, xt, yt)
    fig.update_xaxes(type="log", row=1, col=2)
    fig.update_xaxes(range=[0, 100], row=2, col=2)
    return finish(fig, title, 640, barmode="stack", legend={"orientation": "h", "y": -0.1})


# ---------------------------------------------------------------- V18


@interactive("V18", V18_COLUMNS)
def v18(data: pd.DataFrame, title: str) -> Any:
    """(a) the reliability points and y = x; (b) the TSV's p-value bins per family, the families
    ranked as the static plot draws them (at most 16)."""
    b = data[(data["panel"] == "b") & data["metric"].isna()]
    out = data[(data["panel"] == "b") & (data["metric"] == OUTSIDE_METRIC)]
    cells = zip(out["table"], out["family"], out["value"], strict=True)
    outside = {(str(t), str(f)): int(v) for t, f, v in cells}
    fams = list(dict.fromkeys(zip(b["table"].astype(str), b["family"].astype(str), strict=True)))
    hist = {f: b[(b["table"] == f[0]) & (b["family"] == f[1])] for f in fams}
    total = {f: float(np.nansum(floats(hist[f]["n"]))) for f in fams}
    order = {t: i for i, (t, _) in enumerate(FAMILIES)}
    keep = sorted(fams, key=lambda f: (-total[f], order.get(f[0], len(order)), f[1]))[:MAX_DRAWN]
    drawn = [f for f in fams if f in keep]
    n_rows = -(-len(drawn) // 4)
    names = ["a  detectability reliability"]
    for t, f in drawn:
        extra = f", {outside[t, f]} outside [0, 1]" if outside.get((t, f)) else ""
        names.append(f"{t}<br>{f}, n = {total[t, f]:g}{extra}")
    specs: list[list[dict[str, Any] | None]] = [[{"colspan": 2}, None, None, None]]
    fig = make_subplots(rows=1 + n_rows, cols=4, specs=specs + [[{}] * 4] * n_rows,
                        subplot_titles=names, vertical_spacing=0.3 / (1 + n_rows))  # fmt: skip
    diagonal = {"dash": "dash", "color": OKABE_ITO["grey"]}
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", showlegend=False, hoverinfo="skip",
                             line=diagonal), row=1, col=1)  # fmt: skip
    a = data[(data["panel"] == "a") & data["metric"].isna()]
    if len(a):
        text = [DOT.join((f"bin {fmt(k)}", f"mean p {fmt(x)}", f"observed {fmt(y)}",
                          f"n {fmt(n)}")) for k, x, y, n in
                zip(a["bin"], a["x"], a["y"], a["n"], strict=True)]  # fmt: skip
        ns = np.nan_to_num(floats(a["n"]))
        size = 4 + 12 * np.sqrt(ns / ns.max()) if ns.max() > 0 else 8
        fig.add_trace(go.Scatter(x=floats(a["x"]), y=floats(a["y"]), mode="lines+markers",
                                 hovertext=text, hoverinfo="text", showlegend=False,
                                 marker={"color": OKABE_ITO["blue"], "size": size}),
                      row=1, col=1)  # fmt: skip
    m = data[(data["panel"] == "a") & data["metric"].notna()]
    if len(m):
        stats = "<br>".join(f"{k} {fmt(v)}" for k, v in zip(m["metric"], m["value"], strict=True))
        note(fig, stats, "", 0.15, 0.85)
    titles(fig, {"row": 1, "col": 1}, "mean predicted p_detect", "observed detection rate",
           range=[0, 1])  # fmt: skip
    for i, f in enumerate(drawn):
        at = {"row": 2 + i // 4, "col": 1 + i % 4}
        fig.add_trace(go.Bar(x=floats(hist[f]["bin"]), y=floats(hist[f]["n"]), offset=0,
                             width=1 / N_BINS, marker={"color": OKABE_ITO["sky_blue"]},
                             showlegend=False), **at)  # fmt: skip
        titles(fig, at, "p-value", "count")
    if not drawn:
        note(fig, "no p-values in this run", x=0.75, y=0.85)
    return finish(fig, title, 340 + 220 * n_rows, bargap=0)
