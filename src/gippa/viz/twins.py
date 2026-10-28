"""Plotly twins of V2, V5, V6, V7 and V8 (M9-15; GIPPA_SPEC.md §10.1, §10.2), and the helpers
every twin shares (``twins_panels`` has V10, V11, V17, V18).

A twin reads only its plot's data TSV (1-based closed coordinates, NA cells as NaN) and the
optional run ``meta`` (``gippa.viz.interactive``); it recomputes nothing: histograms are the TSV's
pre-binned rows and V2's 20,000-element cap is already in its rows. Hover text holds the IDs and
the numbers as the TSV writes them; an NA CI reads "CI unestimated". Pure: no file I/O.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from gippa.config.model import AlphaConfig
from gippa.io.tsv import format_value
from gippa.viz import detection as det
from gippa.viz import style
from gippa.viz.interactive import interactive
from gippa.viz.landscape import COLUMNS as V02_ALL
from gippa.viz.mapping_fdr import COLUMNS as V06_COLUMNS
from gippa.viz.mapping_fdr import DECOY, TARGET
from gippa.viz.mapping_sankey import CLASSES, NODE_GREY, OUTCOME_COLOUR, OUTCOMES, ROUTES
from gippa.viz.mapping_sankey import COLUMNS as V05_COLUMNS

Meta = Mapping[str, Any] | None
DOT = " \N{MIDDLE DOT} "
BLACK, MUTED = "#000000", "#555555"
GL = 5_000  # above this many points a trace is WebGL
V02_COLUMNS = tuple(c for c in V02_ALL if c != "product")  # a TSV without it: product NA


def isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def fmt(v: object) -> str:
    """A cell as the TSV writes it (``%.6g``, ``NA``)."""
    v = v.item() if isinstance(v, np.generic) else v
    return "NA" if isna(v) else format_value(v)


def ci(lo: object, hi: object) -> str:
    return "CI unestimated" if isna(lo) or isna(hi) else f"CI {fmt(lo)} to {fmt(hi)}"


def floats(col: Iterable[object]) -> np.ndarray:
    return np.array([math.nan if isna(v) else float(str(v)) for v in col], dtype=float)


def flag_list(v: object) -> list[str]:
    return [] if isna(v) or not str(v) else str(v).split(";")


def alpha(meta: Meta, name: str) -> tuple[float, bool]:
    """The run's ``stats.alpha.<name>`` from ``meta`` (``True``), else the default (``False``)."""
    got = (meta or {}).get("alpha", {}).get(name)
    return (float(got), True) if got is not None else (float(getattr(AlphaConfig(), name)), False)


def contigs(data: pd.DataFrame, meta: Meta) -> list[str]:
    """The TSV's contigs in FASTA order (``meta``), then in TSV order."""
    reps, present = (meta or {}).get("replicons", {}), dict.fromkeys(map(str, data["contig"]))
    return [c for c in reps if c in present] + [c for c in present if c not in reps]


def rgba(hex_colour: str, a: float) -> str:
    r, g, b = (int(hex_colour[i : i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{a:g})"


def points(n: int, **kw: Any) -> Any:
    return (go.Scattergl if n > GL else go.Scatter)(hoverinfo="text", **kw)


def titles(fig: Any, at: Mapping[str, int], x: str | None, y: str | None, **kw: Any) -> None:
    fig.update_xaxes(title_text=x, **at)
    fig.update_yaxes(title_text=y, **at, **kw)


def finish(fig: Any, title: str, height: int, **layout: Any) -> Any:
    """The shared layout: the §10.1 template, the finding title, the height."""
    fig.update_layout(template=style.plotly_template(), title={"text": title}, height=height,
                      margin={"t": 70, "l": 60, "r": 20, "b": 50}, **layout)  # fmt: skip
    return fig


def note(fig: Any, text: str, ref: str = "paper", x: float = 0.5, y: float = 0.5) -> None:
    xref, yref = (ref, ref) if ref == "paper" else (f"x{ref} domain", f"y{ref} domain")
    fig.add_annotation(text=text, x=x, y=y, xref=xref, yref=yref, showarrow=False,
                       font={"color": MUTED})  # fmt: skip


# ---------------------------------------------------------------- V2


@interactive("V02", V02_COLUMNS)
def v02(data: pd.DataFrame, title: str, meta: Meta = None) -> Any:
    """Per replicon: measured points, the status rug as a marker row, the trend, island shapes."""
    order = contigs(data, meta)
    fig = make_subplots(rows=len(order), cols=1, subplot_titles=order,
                        vertical_spacing=min(0.1, 0.5 / len(order)))  # fmt: skip
    data = data.reindex(columns=V02_ALL)
    for i, contig in enumerate(order, start=1):
        at, first = {"row": i, "col": 1}, i == 1
        d = data[data["contig"].astype(str) == contig]
        kind = {k: d[d["kind"] == k] for k in ("point", "rug", "trend", "island", "oric", "ter")}
        top = float(np.nanmax(floats([*d["position_bp"], *d["end"]])))
        length = float((meta or {}).get("replicons", {}).get(contig, top))
        scale, unit = (1e6, "Mb") if length > 2_000_000 else (1e3, "kb")

        def xs(rows: pd.DataFrame, scale: float = scale) -> np.ndarray:
            return (floats(rows["position_bp"]) - 1) / scale

        def hover(rows: pd.DataFrame, last: Iterable[str]) -> list[str]:
            return [DOT.join((fmt(t), fmt(p), z)) for t, p, z in
                    zip(rows["locus_tag"], rows["product"], last, strict=True)]  # fmt: skip

        p, tr = kind["point"], kind["trend"]
        if len(p):
            fig.add_trace(points(len(p), x=xs(p), y=floats(p["value"]), mode="markers",
                                 name="detected (measured)", showlegend=first,
                                 hovertext=hover(p, (f"log2 {fmt(v)}" for v in p["value"])),
                                 marker={"color": style.DETECTION["detected"], "size": 4}),
                          **at)  # fmt: skip
        ys = floats([*p["value"], *tr["value"]])
        rug_y = float(np.nanmin(ys)) - 0.5 if np.isfinite(ys).any() else 0.0
        for status in sorted(set(kind["rug"]["category"].astype(str))):
            r = kind["rug"][kind["rug"]["category"].astype(str) == status]
            colour = style.STATUS_COLOUR.get(status, style.DETECTION["not_detected"])
            fig.add_trace(points(len(r), x=xs(r), y=np.full(len(r), rug_y), mode="markers",
                                 name=status.replace("_", " "), showlegend=first,
                                 hovertext=hover(r, [status] * len(r)),
                                 marker={"color": colour, "symbol": "line-ns-open", "size": 8}),
                          **at)  # fmt: skip
        if len(tr):
            fig.add_trace(points(len(tr), x=xs(tr), y=floats(tr["value"]), mode="lines",
                                 name="trend (LOESS)", showlegend=first,
                                 line={"color": BLACK, "width": 1.5}), **at)  # fmt: skip
        for _, r in kind["island"].iterrows():
            colour = style.ISLAND.get(str(r["category"]), style.ISLAND["high"])[0]
            x0, x1 = float(r["start"]) - 1, float(r["end"])
            cut = [(x0, min(x1, length)), (0.0, x1 - length)] if x1 > length else [(x0, x1)]
            for a, b in cut:  # an origin-crossing island is two spans (replicon length in meta)
                fig.add_vrect(x0=a / scale, x1=b / scale, fillcolor=colour, opacity=0.25,
                              line_width=0, layer="below", **at)  # fmt: skip
        for k, dash in (("oric", "solid"), ("ter", "dot")):
            for x in xs(kind[k]):
                fig.add_vline(x=x, line={"color": MUTED, "dash": dash, "width": 1}, **at)
        titles(fig, at, f"position ({unit})", "log2 abundance")
    return finish(fig, title, 120 + 230 * len(order), legend={"orientation": "h", "y": -0.15})


# ---------------------------------------------------------------- V5


@interactive("V05", V05_COLUMNS)
def v05(data: pd.DataFrame, title: str) -> Any:
    """The Sankey in the fixed M9-05 node order; one link per flow and stage, by outcome."""
    nodes = [*CLASSES, *ROUTES, *OUTCOMES]
    total = dict.fromkeys(nodes, 0)
    link: dict[str, list[Any]] = {"source": [], "target": [], "value": [], "color": []}
    for r in data.itertuples(index=False):
        flow, n = (str(r.input_class), str(r.route), str(r.outcome)), r.n
        if isna(n) or int(n) <= 0 or not all(f in total for f in flow):
            continue  # zero flows and the conflicts row are not links
        for f in flow:
            total[f] += int(n)
        for a, b in ((flow[0], flow[1]), (flow[1], flow[2])):
            for key, v in zip(link, (nodes.index(a), nodes.index(b), int(n),
                                     rgba(OUTCOME_COLOUR[flow[2]], 0.5)), strict=True):  # fmt: skip
                link[key].append(v)
    node = {"label": [f"{n} ({total[n]})" for n in nodes], "pad": 12,
            "color": [OUTCOME_COLOUR.get(n, NODE_GREY) for n in nodes]}  # fmt: skip
    return finish(go.Figure(go.Sankey(node=node, link=link, arrangement="fixed")), title, 420)


# ---------------------------------------------------------------- V6


@interactive("V06", V06_COLUMNS)
def v06(data: pd.DataFrame, title: str, meta: Meta = None) -> Any:
    """(a) targets/decoys on a log axis, (b) stratum FDR with the alpha line, (c) precision."""
    level, known = alpha(meta, "placement")
    sub = {p: data[data["panel"] == p] for p in "abc"}
    panels = ["a", "b", *(["c"] if len(sub["c"]) else [])]
    names = {"a": "(a) targets vs. decoys", "c": "(c) de novo precision",
             "b": "(b) stratum FDR" + (f" (dashed: {level:g})" if known else "")}  # fmt: skip
    fig = make_subplots(rows=1, cols=len(panels), subplot_titles=[names[p] for p in panels],
                        column_widths=[2.6, 1.5, 1.3][: len(panels)])  # fmt: skip
    for metric, colour in (("targets", TARGET), ("decoys", DECOY)):
        m = sub["a"][sub["a"]["metric"] == metric]
        x = [f"{fmt(c)} aa, {fmt(n)}" for c, n in zip(m["length_class"], m["n_diff"], strict=True)]
        text = [f"{s}{DOT}{metric} {fmt(v)}" for s, v in zip(x, m["value"], strict=True)]
        fig.add_trace(go.Bar(x=x, y=floats(m["value"]), name=metric, marker={"color": colour},
                             hovertext=text, hoverinfo="text"), row=1, col=1)  # fmt: skip
    b = sub["b"].pivot_table(index="n_diff", columns="metric", values="value", aggfunc="first",
                             dropna=False)  # fmt: skip
    b = b.reindex(columns=["fdr_stratum", "n_accepted"])
    for na in (False, True):  # an NA stratum: an open marker at the floor
        s = b[b["fdr_stratum"].isna() == na]
        text = [f"n_diff {fmt(m)}{DOT}FDR {fmt(f)}{' (no targets)' if na else ''}{DOT}{fmt(n)} "
                "accepted" for m, f, n in zip(s.index, s["fdr_stratum"], s["n_accepted"],
                                              strict=True)]  # fmt: skip
        y = np.zeros(len(s)) if na else floats(s["fdr_stratum"])
        fig.add_trace(go.Scatter(x=floats(s.index), y=y, mode="markers", hovertext=text,
                                 hoverinfo="text", showlegend=False,
                                 marker={"color": MUTED if na else TARGET, "size": 7,
                                         "symbol": "circle-open" if na else "circle"}),
                      row=1, col=2)  # fmt: skip
    fig.add_hline(y=level, line={"dash": "dash", "color": MUTED, "width": 1}, row=1, col=2)
    if "c" in panels:
        c = sub["c"]
        y, lo, hi = floats(c["value"]), floats(c["ci_low"]), floats(c["ci_high"])
        cells = zip(c["score_bin"], c["value"], c["ci_low"], c["ci_high"], strict=True)
        text = [f"bin {fmt(s)}{DOT}precision {fmt(v)} ({ci(a, z)})" for s, v, a, z in cells]
        sym = ["circle-open" if "c_assumed" in flag_list(f) else "circle" for f in c["flags"]]
        fig.add_trace(go.Scatter(x=floats(c["score_bin"]), y=y, mode="markers", hovertext=text,
                                 hoverinfo="text", showlegend=False,
                                 error_y={"array": hi - y, "arrayminus": y - lo, "color": MUTED},
                                 marker={"color": TARGET, "size": 7, "symbol": sym}),
                      row=1, col=3)  # fmt: skip
        titles(fig, {"row": 1, "col": 3}, "de novo score bin", "precision (CI)")
    titles(fig, {"row": 1, "col": 1}, "differences per length class (aa)",
           "placements (0 at the floor)", type="log")  # fmt: skip
    titles(fig, {"row": 1, "col": 2}, "differences (n_diff)", "stratum FDR", rangemode="tozero")
    return finish(fig, title, 380, barmode="group")


# ---------------------------------------------------------------- V7


@interactive("V07", det.V07_COLUMNS)
def v07(data: pd.DataFrame, title: str) -> Any:
    """The O - E forest in drawing order: one trace per row, region rows by region type."""
    fig, labels = go.Figure(), []
    for r in data.to_dict("records"):
        labels.append(f"{r['group']} (n = {fmt(r['n'])})")
        region = r["block"] == "region"
        colour = style.REGION.get(str(r["region_type"]), BLACK) if region else BLACK
        o, lo, hi, flags = (math.nan if isna(r["o_minus_e"]) else float(r["o_minus_e"]),
                            r["ci_low"], r["ci_high"], flag_list(r["flags"]))  # fmt: skip
        has_ci = not (isna(lo) or isna(hi) or lo == hi or "ci_unestimated" in flags)
        err = {"array": [float(hi) - o], "arrayminus": [o - float(lo)], "color": colour}
        text = DOT.join((str(r["group"]), f"O \N{MINUS SIGN} E {fmt(r['o_minus_e'])} (95 % "
                         f"{ci(lo, hi) if has_ci else 'CI unestimated'})",
                         f"q {fmt(r['q_value'])}"))  # fmt: skip
        sym = "circle-open" if "small_region" in flags else "circle"
        fig.add_trace(go.Scatter(x=[o], y=[labels[-1]], mode="markers", hovertext=[text],
                                 error_x=err if has_ci else None, hoverinfo="text",
                                 showlegend=False, marker={"color": colour, "size": 7,
                                                           "symbol": sym}))  # fmt: skip
    fig.add_vline(x=0.0, line={"color": "#777777", "width": 1})
    fig.update_yaxes(type="category", categoryorder="array", categoryarray=labels,
                     autorange="reversed")  # fmt: skip
    fig.update_xaxes(title_text=det.V07_XLABEL)
    return finish(fig, title, 140 + 22 * len(labels))


# ---------------------------------------------------------------- V8


@interactive("V08", det.V08_COLUMNS)
def v08(data: pd.DataFrame, title: str) -> Any:
    """Per COG category: a violin with points above 50 measured values, a strip otherwise; the
    detected fractions as bars when nothing is measured."""
    summ, genes = data[data["locus_tag"].isna()], data[data["locus_tag"].notna()]
    measured = bool(np.isfinite(floats(genes["abundance_log2"])).any())
    fig, letters = go.Figure(), [str(c) for c in summ["category"]]
    for x, (letter, r) in enumerate(zip(letters, summ.to_dict("records"), strict=True)):
        colour = det.NONE_GREY if letter == "none" else style.cog_colour(letter)
        count = f"{fmt(r['n_detected'])}/{fmt(r['n_genes'])}"
        fig.add_annotation(x=x, y=1.0, yref="paper", yanchor="bottom", showarrow=False, text=count)
        if not measured:
            frac = float(r["n_detected"]) / float(r["n_genes"]) if r["n_genes"] else 0.0
            fig.add_trace(go.Bar(x=[x], y=[frac], name=letter, marker={"color": colour},
                                 hovertext=[f"{letter}{DOT}{count} detected"],
                                 hoverinfo="text"))  # fmt: skip
            continue
        g = genes[genes["category"].astype(str) == letter]
        ok = np.isfinite(floats(g["abundance_log2"]))
        vals, tags = floats(g["abundance_log2"])[ok], g["locus_tag"].astype(str)[ok]
        text = [f"{t}{DOT}log2 {fmt(v)}" for t, v in zip(tags, vals, strict=True)]
        common: dict[str, Any] = {"y": vals, "name": letter, "hovertext": text, "hoverinfo": "text"}
        if len(vals) > det.MAX_STRIP:
            fig.add_trace(go.Violin(x=[x] * len(vals), points="all", fillcolor=rgba(colour, 0.6),
                                    line={"color": colour}, marker={"color": colour, "size": 3},
                                    box={"visible": False}, **common))  # fmt: skip
        elif len(vals):
            fig.add_trace(go.Scatter(x=[x + d for d in det._jitter(len(vals))], mode="markers",
                                     marker={"color": colour, "size": 5}, **common))  # fmt: skip
    fig.update_xaxes(title_text="COG category", tickvals=list(range(len(letters))),
                     ticktext=letters, range=[-0.6, len(letters) - 0.4])  # fmt: skip
    if measured:
        fig.update_yaxes(title_text=det.V08_YLABEL)
    else:
        fig.update_yaxes(title_text="detected fraction", range=[0.0, 1.0])
    return finish(fig, title, 380, showlegend=False)
