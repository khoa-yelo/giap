"""Plotly twins of V2, V6, V8, V10, V11, V17 and V18 (M9-15 behaviour 2-4, example 3, lead
decision 1: V2's data TSV gains ``product``; GIPPA_SPEC.md §10.2). Each twin's trace values are
checked against the data TSV columns it was built from, and its axis titles against the static
plot's. V5, V7, the template and the examples 1, 2, 4 are in ``test_twins.py``.

Frames are the data TSVs as the ``report`` step reads them (``helpers.twins``); coordinates are
the TSV's 1-based closed ones.

Interface assumed beyond the spec text (M9-15 test author; see also ``test_twins.py``):
1. Multi-panel twins are Plotly subplots; a trace's panel is found by the title of its x axis
   (and y axis where two panels share an x title), which is the static panel's axis label.
2. V2: one subplot per contig (its title, an annotation, names the contig); x is
   ``position (kb)`` or ``position (Mb)`` (or ``position (bp)``) with the positions divided
   accordingly; points, rug and trend are ``scatter``/``scattergl`` traces; the trend is a
   lines-only trace; islands are ``layout.shapes`` rectangles filled with the ``ISLAND`` colour;
   point hover holds the locus tag and the product (``NA`` when missing).
3. V6: panel (a) bars are traces named ``targets`` and ``decoys`` on a log y axis; panel (c)
   CIs are ``error_y`` (``array``/``arrayminus``); the alpha line is a ``layout.shapes`` line
   with ``y0 == y1`` on panel (b)'s y axis.
4. V8: one trace per COG category, named by its letter, in the TSV's category order; a
   ``violin`` (with ``points``) above 50 measured values, a ``scatter`` strip otherwise; without
   measured values, ``bar`` heights ``n_detected / n_genes`` on a "detected fraction" axis.
5. V10: envelopes are filled scatter traces (``fill`` ``toself``, or ``tonexty`` over the
   previous trace) coloured as the static bands; the observed line is the TU-block value.
6. V11, V17, V18: pre-binned rows are ``bar`` traces whose heights are the TSV's ``density`` /
   ``n`` / ``value`` cells; V17's warning list is a ``table`` trace; V17 (e) shows ``100 *
   value`` on the "missing (%)" axis; V18 (b) bars sit at each bin's left edge or centre.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import twins as tw
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")

SCATTER = ("scatter", "scattergl")


def _finite(xs: list[float]) -> list[float]:
    return [x for x in xs if not math.isnan(x)]


def _on_xy(s: dict[str, Any], xtitle: str, ytitle: str, *types: str) -> list[dict[str, Any]]:
    return [t for t in tw.on_axes(s, xtitle, *types) if tw.axis_title(s, t, "y") == ytitle]


def _marker(t: dict[str, Any], key: str, n: int) -> list[Any]:
    v = (t.get("marker") or {}).get(key)
    vals = tw.values(v)
    return vals * n if len(vals) == 1 and n > 1 else vals


# ---------------------------------------------------------------- V2


def _scale(s: dict[str, Any]) -> float:
    titles = {tw.title_of(v) for k, v in s["layout"].items() if k.startswith("xaxis")}
    units = {t for t in titles if t.startswith("position (")}
    assert len(units) == 1, titles
    return {"position (bp)": 1.0, "position (kb)": 1e3, "position (Mb)": 1e6}[units.pop()]


def _close(a: float, b: float, scale: float) -> bool:
    return abs(a * scale - b) <= 1.0 + 1e-6 * abs(b)  # 1-based TSV vs 0-based drawing: 1 bp


def test_v2_axes_and_trace_types() -> None:
    s = tw.spec("V02", tw.v02())
    assert {tw.ttype(t) for t in tw.traces(s)} <= set(SCATTER)
    assert _scale(s) == 1e3  # a 100-kb replicon is drawn in kb, as the static plot
    assert "log2 abundance" in tw.axis_titles(s)


def test_v2_points_equal_tsv_and_hover_has_locus_tag_and_product() -> None:
    data = tw.v02()
    s = tw.spec("V02", data)
    scale = _scale(s)
    pts = data[data["kind"] == "point"]
    drawn = tw.point_pairs(t for t in tw.traces(s) if "markers" in tw.mode(t))
    for _, r in pts.iterrows():
        assert any(_close(x, r["position_bp"], scale) and y == pytest.approx(r["value"])
                   for x, y in drawn), r.to_dict()  # fmt: skip
    hover = tw.all_hovers(s)
    g1 = [h for h in hover if "g1" in h and "DNA polymerase III subunit beta" in h]
    assert g1, hover[:10]
    assert any("g2" in h and "NA" in h for h in hover), "a missing product shows NA"


def test_v2_trend_is_a_line_through_the_trend_rows() -> None:
    data = tw.v02()
    s = tw.spec("V02", data)
    scale = _scale(s)
    trend = data[data["kind"] == "trend"]
    lines = [t for t in tw.traces(s) if tw.mode(t) == "lines"]
    assert len(lines) == 1
    xs, ys = tw.floats(lines[0]["x"]), tw.floats(lines[0]["y"])
    assert ys == pytest.approx(list(trend["value"]))
    assert all(_close(x, p, scale) for x, p in zip(xs, trend["position_bp"], strict=True))


def test_v2_rug_is_a_marker_row_in_status_colours() -> None:
    from gippa.viz.style import STATUS_COLOUR

    data = tw.v02()
    s = tw.spec("V02", data)
    scale = _scale(s)
    marks: list[tuple[float, float, str]] = []
    for t in tw.traces(s):
        if "markers" not in tw.mode(t):
            continue
        xs, ys = tw.floats(t.get("x")), tw.floats(t.get("y"))
        cols = _marker(t, "color", len(xs))
        marks += [(x, y, str(tw.to_hex(str(c)))) for x, y, c in zip(xs, ys, cols, strict=True)]
    rug = data[data["kind"] == "rug"]
    ys_rug = set()
    for _, r in rug.iterrows():
        want = STATUS_COLOUR[str(r["category"])].upper()
        hit = [m for m in marks if _close(m[0], r["position_bp"], scale) and m[2] == want]
        assert hit, (r.to_dict(), marks)
        ys_rug.add(hit[0][1])
    assert len(ys_rug) == 1, "the rug is one marker row"


def test_v2_islands_are_shapes_in_island_colours() -> None:
    from gippa.viz.style import ISLAND

    data = tw.v02()
    s = tw.spec("V02", data)
    scale = _scale(s)
    rects = [sh for sh in s["layout"].get("shapes", []) if sh.get("type") == "rect"]
    isl = data[data["kind"] == "island"]
    assert len(rects) == len(isl) == 2
    for _, r in isl.iterrows():
        colour = ISLAND[str(r["category"])][0].upper()
        hit = [sh for sh in rects if _close(float(sh["x0"]), r["start"] - 1, scale)
               and _close(float(sh["x1"]), r["end"], scale)]  # fmt: skip
        assert len(hit) == 1, (r.to_dict(), rects)
        assert tw.to_hex(str(hit[0]["fillcolor"])) == colour


def test_v2_one_subplot_per_contig() -> None:
    rows = tw.v02_rows("chr1") + [r for r in tw.v02_rows("p1") if r[1] in ("point", "rug")]
    s = tw.spec("V02", tw.frame(tw.V02_COLUMNS, rows))
    xaxes = {t.get("xaxis", "x") for t in tw.traces(s)}
    assert len(xaxes) == 2
    notes = " ".join(tw.annotations(s))
    assert "chr1" in notes and "p1" in notes


def test_v2_element_cap_twin_draws_the_tsvs_binned_rows_only() -> None:
    """Behaviour 4: above 20,000 points the TSV already holds the static plot's bin medians;
    the twin draws those rows and no more, in Mb on a 5-Mb replicon, within a size budget."""
    n = 20_000
    rows: list[tw.Row] = []
    for i in range(n):
        pos = 125 + 250 * i
        rows.append(("chr1", "point", None, None, pos, None, None, 2.0 + (i % 7) * 0.25,
                     "detected"))  # fmt: skip
        rows.append(("chr1", "trend", None, None, pos, None, None, 2.5, None))
    rows += [("chr1", "rug", f"r{i}", None, 1001 + 1000 * i, None, None, None, "not_detected")
             for i in range(4_000)]  # fmt: skip
    data = tw.frame(tw.V02_COLUMNS, rows)
    text = tw.figure_text("V02", data)
    assert text is not None
    s = json.loads(text)
    assert _scale(s) == 1e6
    drawn = sum(len(tw.values(t.get("x"))) for t in tw.traces(s))
    assert n <= drawn <= len(data)
    assert len(text.encode()) < 6_000_000, f"V2 twin JSON is {len(text.encode())} bytes"


def test_v2_data_tsv_gains_product_from_genes_tsv(tmp_path: Path) -> None:
    """Lead decision 1: the static V2 writes ``product`` (from ``genes.tsv``) into its data TSV;
    the twin built from that TSV shows it."""
    from helpers.viz_landscape import gene, landscape_out, plot_spec

    from gippa.config import resolve_config
    from gippa.io.run_view import RunView
    from gippa.viz import style

    style.apply()
    genes = []
    for i, (mid, log2) in enumerate(zip([10_000, 30_000, 50_000, 70_000, 90_000],
                                        [1.0, 2.0, 3.0, 2.0, 1.0], strict=True)):  # fmt: skip
        row = gene(f"g{i + 1}", mid, log2)
        row["product"] = f"protein number {i + 1}" if i != 2 else None
        genes.append(row)
    undet = gene("g9", 60_000, None, status="not_detected")
    undet["product"] = "transposase"
    out = landscape_out(tmp_path, genes=[*genes, undet])
    (inst,) = plot_spec("V02").fn(RunView.from_outdir(out), resolve_config())
    assert "product" in inst.data.columns
    pts = inst.data[inst.data["kind"] == "point"]
    got = {
        str(t): (None if pd.isna(p) else str(p))
        for t, p in zip(pts["locus_tag"], pts["product"], strict=True)
    }
    assert got == {"g1": "protein number 1", "g2": "protein number 2", "g3": None,
                   "g4": "protein number 4", "g5": "protein number 5"}  # fmt: skip
    s = tw.spec("V02", tw.instance_frame(inst.data, inst.coords))
    hover = tw.all_hovers(s)
    assert any("g4" in h and "protein number 4" in h for h in hover)
    assert any("g3" in h and "NA" in h for h in hover)


# ---------------------------------------------------------------- V6

A_X, B_X, C_X = "differences per length class (aa)", "differences (n_diff)", "de novo score bin"


def test_v6_three_subplots_with_static_axis_titles() -> None:
    s = tw.spec("V06", tw.v06())
    titles = tw.axis_titles(s)
    assert {A_X, B_X, C_X, "stratum FDR", "precision (CI)"} <= titles
    assert any(t.startswith("placements") for t in titles)
    assert len({t.get("xaxis", "x") for t in tw.traces(s)}) == 3


def test_v6_without_denovo_rows_has_two_subplots() -> None:
    s = tw.spec("V06", tw.v06(with_c=False))
    assert C_X not in tw.axis_titles(s)
    assert len({t.get("xaxis", "x") for t in tw.traces(s)}) == 2


def test_v6_panel_a_bars_equal_targets_and_decoys_on_log_axis() -> None:
    s = tw.spec("V06", tw.v06())
    bars = {str(t.get("name")): t for t in tw.on_axes(s, A_X, "bar")}
    assert {"targets", "decoys"} <= set(bars)
    for name, want in (("targets", [120.0, 30.0]), ("decoys", [2.0, 1.0])):
        ys = [y for y in _finite(tw.floats(bars[name]["y"])) if y > 0]
        assert ys == want, (name, tw.values(bars[name]["y"]))
    assert tw.axis(s, tw.axis_key(bars["targets"], "y")).get("type") == "log"


def test_v6_panel_b_points_alpha_shape_and_hover() -> None:
    s = tw.spec("V06", tw.v06())
    pts = tw.on_axes(s, B_X, *SCATTER)
    pairs = {(x, y) for x, y in tw.point_pairs(pts) if not math.isnan(y)}
    assert {(0.0, 0.004), (1.0, 0.03)} <= pairs
    hover = [h for t in pts for h in tw.hovers(t)]
    assert any("118" in h for h in hover), hover
    assert any("NA" in h for h in hover), "the NA stratum shows NA"
    yref = str(pts[0].get("yaxis", "y"))
    lines = [sh for sh in s["layout"].get("shapes", [])
             if sh.get("type") == "line" and float(sh["y0"]) == float(sh["y1"]) > 0]  # fmt: skip
    assert any(str(sh.get("yref", "y")) == yref for sh in lines), (yref, lines)


def test_v6_panel_c_precision_with_ci_error_bars() -> None:
    s = tw.spec("V06", tw.v06())
    (t,) = [t for t in tw.on_axes(s, C_X, *SCATTER) if t.get("error_y")]
    assert tw.floats(t["y"]) == pytest.approx([0.9, 0.95])
    assert tw.floats(t["x"]) == pytest.approx([1.0, 2.0])
    assert tw.floats(t["error_y"]["array"]) == pytest.approx([0.04, 0.03])
    assert tw.floats(t["error_y"]["arrayminus"]) == pytest.approx([0.05, 0.05])
    hover = [h for x in tw.on_axes(s, C_X) for h in tw.hovers(x)]
    assert any("CI 0.85 to 0.94" in h for h in hover), hover


# ---------------------------------------------------------------- V8


def _v8_traces(s: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(t["name"]): t for t in tw.traces(s) if t.get("name") is not None}


def test_v8_violin_above_50_strips_otherwise_values_equal_tsv() -> None:
    data = tw.v08()
    s = tw.spec("V08", data)
    ts = _v8_traces(s)
    genes = data[data["locus_tag"].notna()]
    for letter, kind in (("J", "violin"), ("C", "scatter"), ("none", "scatter")):
        t = ts[letter]
        assert tw.ttype(t) in ((kind,) if kind == "violin" else SCATTER), (letter, tw.ttype(t))
        want = sorted(genes.loc[genes["category"] == letter, "abundance_log2"].dropna())
        assert sorted(_finite(tw.floats(t["y"]))) == pytest.approx(want)
    assert ts["J"].get("points") in ("all", "outliers", "suspectedoutliers")
    assert tw.title_of(s["layout"]["xaxis"]) == "COG category"
    assert "log2 abundance" in tw.axis_titles(s)


def test_v8_category_order_and_colours_follow_the_tsv() -> None:
    from gippa.viz.style import cog_colour

    data = tw.v08()
    s = tw.spec("V08", data)
    names = [str(t["name"]) for t in tw.traces(s) if t.get("name") in ("J", "C", "none")]
    assert names == [str(c) for c in data.loc[data["locus_tag"].isna(), "category"]]
    ts = _v8_traces(s)
    for letter in ("J", "C"):
        t = ts[letter]
        col = (t.get("marker") or {}).get("color") or t.get("fillcolor") or t.get("line", {})
        col = col.get("color") if isinstance(col, dict) else col
        assert {tw.to_hex(str(c)) for c in tw.values(col)} == {cog_colour(letter).upper()}


def test_v8_strip_hover_names_the_locus_tag() -> None:
    hover = tw.hovers(_v8_traces(tw.spec("V08", tw.v08()))["C"])
    assert len(hover) == 5 and all(any(f"c{i}" in h for h in hover) for i in range(5))


def test_v8_without_measured_values_draws_detected_fractions() -> None:
    data = tw.v08(measured=False)
    s = tw.spec("V08", data)
    bars = tw.traces(s, "bar")
    assert bars and "detected fraction" in tw.axis_titles(s)
    got = sorted(y for t in bars for y in _finite(tw.floats(t["y"])))
    summ = data[data["locus_tag"].isna()]
    want = sorted(float(d) / float(n) for d, n in zip(summ["n_detected"], summ["n_genes"],
                                                      strict=True))  # fmt: skip
    assert got == pytest.approx(want)


# ---------------------------------------------------------------- V10


def _envelopes(s: dict[str, Any]) -> list[dict[str, Any]]:
    ts = tw.traces(s)
    out: list[dict[str, Any]] = []
    for i, t in enumerate(ts):
        fill = t.get("fill", "none")
        if fill == "toself":
            out.append(t)
        elif fill in ("tonexty", "tonextx") and i > 0:
            out += [ts[i - 1], t]
    return out


def test_v10_axis_titles_and_types() -> None:
    s = tw.spec("V10", tw.v10())
    assert {"lag (genes)", "Moran's I", "BB join statistic (Pearson residual)"} <= tw.axis_titles(s)
    assert {tw.ttype(t) for t in tw.traces(s)} <= set(SCATTER)


def test_v10_filled_envelopes_per_null_hold_null_lo_and_hi() -> None:
    data = tw.v10()
    s = tw.spec("V10", data)
    env = _envelopes(s)
    assert env, "no filled envelope traces"
    ys = {round(y, 9) for t in env for y in _finite(tw.floats(t.get("y")))}
    want = {round(float(v), 9) for c in ("null_lo", "null_hi") for v in data[c]}
    assert want <= ys
    fills = {tw.to_hex(str(t["fillcolor"])) for t in tw.traces(s) if t.get("fillcolor")}
    assert {"#CCCCCC", "#56B4E9"} <= fills


def test_v10_observed_line_is_the_tu_block_value_with_q_in_hover() -> None:
    data = tw.v10()
    s = tw.spec("V10", data)
    tu = data[(data["statistic"] == "morans_i") & (data["null_model"] == "tu_block")]
    want = [(float(lag), float(v)) for lag, v in zip(tu["lag"], tu["value"], strict=True)]
    obs = [t for t in tw.on_axes(s, "lag (genes)", *SCATTER) if t.get("fill", "none") == "none"
           and "lines" in tw.mode(t) and tw.point_pairs([t]) == pytest.approx(want)]  # fmt: skip
    assert obs, "no line through the TU-block Moran's I values"
    hover = [h for t in tw.traces(s) for h in tw.hovers(t)]
    assert any("q 0.001" in h for h in hover)


def test_v10_flagged_rows_draw_nothing() -> None:
    s = tw.spec("V10", tw.v10_not_estimable())
    drawn = [y for t in tw.traces(s) for y in _finite(tw.floats(t.get("y")))]
    assert drawn == []
    assert any("not estimable" in a for a in tw.annotations(s))


# ---------------------------------------------------------------- V11

DX = "|Δ| partner log2 abundance"
UP, DOWN = "upstream partner log2 abundance", "downstream partner log2 abundance"


def test_v11_axis_titles() -> None:
    s = tw.spec("V11", tw.v11())
    assert {DX, "density", "local fdr", UP, DOWN} <= tw.axis_titles(s)


def test_v11_density_bars_are_the_tsv_bins() -> None:
    s = tw.spec("V11", tw.v11())
    bars = _on_xy(s, DX, "density", "bar")
    by = {str(t.get("name", "")).split(" ")[0]: tw.floats(t["y"]) for t in bars}
    assert by.get("OP") == pytest.approx(tw.V11_OP_DENSITY)
    assert by.get("NO") == pytest.approx(tw.V11_NO_DENSITY)
    assert not {"CV", "SW", "DV"} & set(by), "classes without bins draw no bars"


def test_v11_lfdr_line_equals_panel_b_rows() -> None:
    data = tw.v11()
    s = tw.spec("V11", data)
    b = data[data["panel"] == "b"]
    want = sorted(zip(b["abs_delta_log2"].astype(float), b["lfdr"].astype(float), strict=True))
    lines = [t for t in _on_xy(s, DX, "local fdr", *SCATTER) if "lines" in tw.mode(t)]
    assert lines
    assert sorted(tw.point_pairs(lines[:1])) == pytest.approx(want)


def test_v11_partner_scatter_with_flag_symbols() -> None:
    data = tw.v11()
    s = tw.spec("V11", data)
    c = data[data["panel"] == "c"]
    ts = tw.on_axes(s, UP, *SCATTER)
    sym: dict[tuple[float, float], str] = {}
    for t in ts:
        xs, ys = tw.floats(t.get("x")), tw.floats(t.get("y"))
        for x, y, m in zip(xs, ys, _marker(t, "symbol", len(xs)) or ["circle"] * len(xs),
                           strict=True):  # fmt: skip
            sym[x, y] = str(m)
    want = {(float(u), float(d)) for u, d in zip(c["upstream_log2"], c["downstream_log2"],
                                                 strict=True)}  # fmt: skip
    assert want <= set(sym)
    plain = {sym[5.0, 5.25], sym[2.0, 3.5]}
    flagged = [sym[6.0, 5.5], sym[7.0, 4.0]]
    assert plain <= {"circle", "0"}
    assert all(f not in ("circle", "0") for f in flagged) and flagged[0] != flagged[1]
    hover = [h for t in ts for h in tw.hovers(t)]
    assert any("P0002" in h and "q 0.001" in h for h in hover), hover


# ---------------------------------------------------------------- V17


def test_v17_five_subplots_and_warning_table() -> None:
    data = tw.v17()
    s = tw.spec("V17", data)
    cart = [t for t in tw.traces(s) if tw.ttype(t) != "table"]
    assert len({t.get("xaxis", "x") for t in cart}) == 5
    (table,) = tw.traces(s, "table")
    cells = [str(v) for col in table["cells"]["values"] for v in tw.values(col)]
    msgs = data.loc[data["group"] == "message", "key"]
    assert all(any(m in c for c in cells) for m in msgs), cells
    assert {"checks", "length (bp)", "length (aa)", "distinct peptides", "log2 abundance",
            "genes", "missing (%)"} <= tw.axis_titles(s)  # fmt: skip


def test_v17_panel_a_level_counts() -> None:
    s = tw.spec("V17", tw.v17())
    (t,) = [t for t in tw.traces(s, "bar") if tw.axis_title(s, t, "y") == "checks"]
    assert [str(x) for x in tw.values(t["x"])] == ["ok", "info", "warn", "error"]
    assert tw.floats(t["y"]) == [7.0, 2.0, 3.0, 1.0]


def test_v17_panel_b_contig_lengths_on_log_axis() -> None:
    s = tw.spec("V17", tw.v17())
    (t,) = tw.on_axes(s, "length (bp)", "bar")
    assert t.get("orientation") == "h"
    assert dict(zip([str(y) for y in tw.values(t["y"])], tw.floats(t["x"]), strict=True)) == {
        "chr1": 4_600_000.0,
        "p1": 50_000.0,
    }
    assert tw.axis(s, tw.axis_key(t, "x")).get("type") == "log"


def test_v17_panels_c_d_e_values_equal_tsv() -> None:
    s = tw.spec("V17", tw.v17())
    c = {str(t.get("name")): dict(zip(tw.floats(t["x"]), tw.floats(t["y"]), strict=True))
         for t in tw.on_axes(s, "length (aa)", "bar")}  # fmt: skip
    assert c == {"placed": {8.0: 3.0, 9.0: 5.0}, "unmapped": {9.0: 2.0}}
    d = {str(t.get("name")): sorted(tw.floats(t["y"])) for t in tw.on_axes(s, "log2 abundance",
                                                                          "bar")}  # fmt: skip
    assert d == {"measured": [4.0, 6.0], "imputed": [1.0]}
    (e,) = tw.on_axes(s, "missing (%)", "bar")
    got = dict(zip([str(y) for y in tw.values(e["y"])], tw.floats(e["x"]), strict=True))
    assert got == pytest.approx({"abundance_raw": 10.0, "cai": 0.0, "mw_da": 25.0})


# ---------------------------------------------------------------- V18: example 3

RX, RY = "mean predicted p_detect", "observed detection rate"


def _hist_bars(s: dict[str, Any]) -> list[dict[str, Any]]:
    return tw.on_axes(s, "p-value", "bar")


def test_example3_v18_bar_heights_equal_the_prebinned_counts() -> None:
    s = tw.spec("V18", tw.v18())
    bars = _hist_bars(s)
    got = Counter(tuple(tw.floats(t["y"])) for t in bars)
    want = Counter(tuple(float(c) for c in counts) for counts in tw.V18_FAMILIES.values())
    assert got == want
    left = [k / 20 for k in range(20)]
    for t in bars:
        xs = tw.floats(t["x"])
        assert xs == pytest.approx(left) or xs == pytest.approx([x + 0.025 for x in left])
        assert "histfunc" not in t
    assert {"p-value", "count", RX, RY} <= tw.axis_titles(s)
    assert any("outside" in a for a in tw.annotations(s)), "n outside [0, 1] is reported"


def test_v18_reliability_points_and_diagonal() -> None:
    data = tw.v18()
    s = tw.spec("V18", data)
    a = data[(data["panel"] == "a") & data["metric"].isna()]
    want = {(float(x), float(y)) for x, y in zip(a["x"], a["y"], strict=True)}
    pts = [t for t in _on_xy(s, RX, RY, *SCATTER) if "markers" in tw.mode(t)]
    assert want <= set(tw.point_pairs(pts))
    diag = [t for t in _on_xy(s, RX, RY, *SCATTER) if tw.point_pairs([t]) == [(0, 0), (1, 1)]]
    shapes = [sh for sh in s["layout"].get("shapes", []) if sh.get("type") == "line"
              and (float(sh["x0"]), float(sh["y0"]), float(sh["x1"]), float(sh["y1"]))
              == (0.0, 0.0, 1.0, 1.0)]  # fmt: skip
    assert diag or shapes, "no y = x diagonal"


_counts = st.lists(st.integers(min_value=0, max_value=500), min_size=20, max_size=20)


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(_counts, min_size=1, max_size=4))
def test_property_v18_bars_equal_any_prebinned_counts(families: list[list[int]]) -> None:
    fams = {("pairs", f"test_{i}"): c for i, c in enumerate(families)}
    data = tw.v18(fams)
    text = tw.figure_text("V18", data)
    assert text is not None and text == tw.figure_text("V18", data)
    got = Counter(tuple(tw.floats(t["y"])) for t in _hist_bars(json.loads(text)))
    assert got == Counter(tuple(float(x) for x in c) for c in families)
