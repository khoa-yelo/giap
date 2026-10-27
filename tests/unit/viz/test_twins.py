"""Plotly twins of the static plots (M9-15; GIPPA_SPEC.md §10.1 colours and fonts, §10.2 column
"Interactive in report"; M9-14 behaviour 1: a twin is built from the plot's data TSV only).

This file: the registered-twin list, ``gippa.viz.style.plotly_template()``, examples 1-4, the
degenerate rows (header-only TSV -> no twin; NA hover fields -> "NA"), determinism, and the V5
(Sankey) and V7 (forest) twins. The other twins are in ``test_twins_panels.py``; the report-level
checks (``--no-interactive``, skipped plots, size budget)
are in ``tests/unit/report/test_twins_report.py``.

Frames are the data TSVs exactly as the ``report`` step reads them (``helpers.twins.frame``:
written with the TSV cell formatter, read with ``na_values=["NA"], keep_default_na=False``).

Interface assumed beyond the spec text (M9-15 test author):
1. ``gippa.viz.interactive.INTERACTIVE_PLOTS`` is the tuple of §10.2 plot IDs whose "Interactive
   in report" column is Plotly (``Vnn`` form, V16 included); after every ``gippa.viz`` module is
   imported, each of them except V14, V15, V16 has a twin (``has_twin``), and V01, V03, V12, V13
   have none.
2. ``figure_json(pid, data, title)`` returns ``None`` for a header-only data TSV (the empty state:
   no twin; ``gippa.steps.plots.twins`` already skips ``None``).
3. Hover text is fully formatted per point in the trace's ``hovertext`` (or ``text``); numbers
   appear as the TSV writes them (``%.6g``); a Unicode minus and thin spaces are allowed (the
   tests compare on ASCII ``-`` and plain spaces). An NA field shows ``NA``; an NA CI shows
   ``CI unestimated``.
4. Every twin's layout carries ``plotly_template()`` (its font and colorway are in the JSON).
5. V5: one ``sankey`` trace; ``node.label`` holds all 11 nodes in the fixed M9-05 order (class,
   route, outcome), each label starting with the node name; ``node.color`` is the M9-05 colour
   (outcomes in ``OUTCOME_COLOUR``, other nodes ``NODE_GREY``); links are per flow and per stage
   (no aggregation, as the static ribbons are), so every link of example 1 has value 1 and each
   link is coloured by its flow's outcome.
6. V7: horizontal error bars are ``error_x`` with ``array = ci_high - o_minus_e`` and
   ``arrayminus = o_minus_e - ci_low`` (symmetric arrays not used); each region row is its own
   trace (example 2: "2 region traces"), coloured by its region type; the y axis is categorical
   and lists the rows top to bottom in TSV order (``categoryarray`` with ``autorange =
   "reversed"``, or reversed ``categoryarray``); each category label starts with the row's group.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import twins as tw
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[3]
SPEC = REPO / "GIPPA_SPEC.md"

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")


# ---------------------------------------------------------------- the registered twins


def _spec_plotly_ids() -> set[str]:
    """§10.2 catalogue rows whose last column (Interactive in report) names Plotly."""
    text = SPEC.read_text(encoding="utf-8")
    body = text.split("### 10.2 Plot catalogue", 1)[1].split("### 10.3", 1)[0]
    out = set()
    for line in body.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        is_row = len(cells) >= 6 and cells[0].startswith("V") and cells[0][1:].isdigit()
        if is_row and "Plotly" in cells[-1]:
            out.add(f"V{int(cells[0][1:]):02d}")
    return out


def test_interactive_plots_list_is_derived_from_spec_10_2() -> None:
    interactive = tw.load()
    listed = tuple(interactive.INTERACTIVE_PLOTS)
    assert set(listed) == _spec_plotly_ids() == set(tw.PLOTLY_PLOTS)
    assert len(listed) == len(set(listed))


def test_every_plotly_plot_but_v14_v15_v16_has_a_twin() -> None:
    interactive = tw.load()
    missing = [p for p in interactive.INTERACTIVE_PLOTS if p not in tw.LATER
               and not interactive.has_twin(p)]  # fmt: skip
    assert missing == []
    assert set(tw.M915) <= set(interactive.INTERACTIVE_PLOTS)


@pytest.mark.parametrize("pid", tw.STATIC_ONLY)
def test_static_only_plots_have_no_twin(pid: str) -> None:
    interactive = tw.load()
    assert not interactive.has_twin(pid)
    assert pid not in interactive.INTERACTIVE_PLOTS


# ---------------------------------------------------------------- the template (§10.1)


def _template() -> Any:
    from gippa.viz.style import plotly_template

    return plotly_template()


def test_plotly_template_type_fonts_and_background() -> None:
    import plotly.graph_objects as go

    t = _template()
    assert isinstance(t, go.layout.Template)
    lay = t.layout
    family = [f.strip().strip("'\"") for f in str(lay.font.family).split(",")]
    assert family[0] == "DejaVu Sans" and family[-1] == "sans-serif", lay.font.family
    from gippa.viz.style import RC

    assert lay.font.size == RC["font.size"] == 7
    assert lay.title.font.size == RC["axes.titlesize"] == 8
    assert lay.xaxis.tickfont.size == lay.yaxis.tickfont.size == RC["xtick.labelsize"] == 6
    assert lay.legend.font.size == RC["legend.fontsize"] == 6
    assert tw.to_hex(str(lay.paper_bgcolor)) == "#FFFFFF"
    assert tw.to_hex(str(lay.plot_bgcolor)) == "#FFFFFF"


def test_plotly_template_colorway_is_the_static_category_cycle() -> None:
    from gippa.viz.style import OKABE_ITO

    cycle = [v.upper() for k, v in OKABE_ITO.items() if k not in ("black", "grey")]
    assert [tw.to_hex(str(c)) for c in _template().layout.colorway] == cycle


def test_template_colour_set_is_within_10_1() -> None:
    import plotly.graph_objects as go

    s = json.loads(go.Figure(layout={"template": _template()}).to_json())
    assert tw.disallowed(s) == set()


@pytest.mark.parametrize("pid", sorted(tw.EXAMPLES))
def test_every_twin_uses_the_template(pid: str) -> None:
    t = _template()
    lay = tw.spec(pid, tw.EXAMPLES[pid]())["layout"]
    tpl = lay["template"]["layout"]
    assert tpl["font"]["family"] == t.layout.font.family
    assert [tw.to_hex(c) for c in tpl["colorway"]] == [tw.to_hex(str(c)) for c in t.layout.colorway]


# ---------------------------------------------------------------- example 4: colours


@pytest.mark.parametrize("pid", sorted(tw.EXAMPLES))
def test_example4_twin_colours_within_10_1_and_viridis_cividis(pid: str) -> None:
    s = tw.spec(pid, tw.EXAMPLES[pid]())
    assert tw.disallowed(s) == set(), f"{pid} uses colours outside §10.1"
    assert any(True for _ in tw.colours(s["data"])), f"{pid}: traces set no colour"


def test_colour_checker_rejects_plotlys_default_palette() -> None:
    """The example 4 checker is not vacuous: the stock ``plotly`` template fails it."""
    import plotly.graph_objects as go

    s = json.loads(go.Figure(go.Scatter(x=[1], y=[1]), layout={"template": "plotly"}).to_json())
    assert tw.disallowed(s)


# ---------------------------------------------------------------- empty state, never recompute


@pytest.mark.parametrize("pid", sorted(tw.HEADERS))
def test_header_only_tsv_gives_no_twin(pid: str) -> None:
    assert tw.load().has_twin(pid), f"{pid} has no twin"
    empty = tw.frame(tw.HEADERS[pid], [])
    assert len(empty) == 0
    assert tw.figure_text(pid, empty) is None


@pytest.mark.parametrize("pid", sorted(tw.EXAMPLES))
def test_no_twin_bins_or_aggregates_in_plotly(pid: str) -> None:
    """Behaviour 3: histograms are the TSV's pre-binned rows, never Plotly's binning."""
    for t in tw.traces(tw.spec(pid, tw.EXAMPLES[pid]())):
        assert tw.ttype(t) not in ("histogram", "histogram2d", "histogram2dcontour"), pid
        assert not {"histfunc", "nbinsx", "nbinsy", "xbins", "ybins", "transforms"} & set(t)


@pytest.mark.parametrize("pid", sorted(tw.EXAMPLES))
def test_twin_json_has_data_and_layout_and_no_external_urls(pid: str) -> None:
    text = tw.figure_text(pid, tw.EXAMPLES[pid]())
    assert text is not None
    s = json.loads(text)
    assert {"data", "layout"} <= set(s) and s["data"]
    assert "http://" not in text and "https://" not in text


# ---------------------------------------------------------------- determinism


@pytest.mark.parametrize("pid", sorted(tw.EXAMPLES))
def test_figure_json_is_byte_identical_across_calls(pid: str) -> None:
    a = tw.figure_text(pid, tw.EXAMPLES[pid](), "V title")
    b = tw.figure_text(pid, tw.EXAMPLES[pid](), "V title")
    assert a is not None and a == b
    assert "uid" not in json.dumps(json.loads(a)["data"])


_HASHES = """
import hashlib, json, sys
from helpers import twins as tw
out = {p: hashlib.sha256(tw.figure_text(p, f(), "t").encode()).hexdigest()
       for p, f in sorted(tw.EXAMPLES.items())}
print(json.dumps(out, sort_keys=True))
"""


def test_figure_json_is_identical_across_processes_and_hash_seeds() -> None:
    env_base = {**os.environ, "PYTHONPATH": os.pathsep.join(
        [str(REPO / "src"), str(REPO / "tests" / "sim"), str(REPO / "tests"),
         os.environ.get("PYTHONPATH", "")])}  # fmt: skip
    runs = []
    for seed in ("0", "12345"):
        res = subprocess.run([sys.executable, "-c", _HASHES], env={**env_base,
                             "PYTHONHASHSEED": seed}, capture_output=True, text=True,
                             timeout=300, check=False)  # fmt: skip
        assert res.returncode == 0, res.stderr[-2000:]
        runs.append(json.loads(res.stdout.strip().splitlines()[-1]))
    assert runs[0] == runs[1]
    assert set(runs[0]) == set(tw.EXAMPLES)


# ---------------------------------------------------------------- V5: example 1 (Sankey)

NODES = ("peptide", "protein", "exact", "tolerant", "miniprot", "tblastn", "unmapped", "placed",
         "contaminant", "rejected by FDR", "unexplained")  # fmt: skip


def _sankey(data: pd.DataFrame) -> dict[str, Any]:
    (t,) = tw.traces(tw.spec("V05", data))
    assert tw.ttype(t) == "sankey"
    return t


def _node_names(t: dict[str, Any]) -> list[str]:
    labels = [str(x) for x in tw.values(t["node"]["label"])]
    names = []
    for lab in labels:
        hit = [n for n in sorted(NODES, key=len, reverse=True) if lab.startswith(n)]
        assert hit, f"node label {lab!r} names no node"
        names.append(hit[0])
    return names


def _links(t: dict[str, Any]) -> list[tuple[str, str, float]]:
    names = _node_names(t)
    link = t["link"]
    src, dst, val = (tw.values(link[k]) for k in ("source", "target", "value"))
    assert len(src) == len(dst) == len(val)
    return [(names[int(a)], names[int(b)], float(v)) for a, b, v in zip(src, dst, val, strict=True)]


def test_example1_sankey_nodes_2_5_4_in_fixed_order() -> None:
    t = _sankey(tw.v05())
    assert _node_names(t) == list(NODES)


def test_example1_sankey_five_flows_with_value_1() -> None:
    links = _links(_sankey(tw.v05()))
    assert all(v == 1.0 for _, _, v in links), links
    routes, outcomes = set(NODES[2:7]), set(NODES[7:])
    second = sorted((a, b) for a, b, _ in links if a in routes and b in outcomes)
    assert second == sorted((r, o) for (_, r, o) in tw.V05_EX1)
    first = sorted((a, b) for a, b, _ in links if a in ("peptide", "protein"))
    assert first == sorted((c, r) for (c, r, _) in tw.V05_EX1)
    assert len(links) == 2 * len(tw.V05_EX1)


def test_sankey_link_values_equal_the_tsv_n_column() -> None:
    counts = {("peptide", "exact", "placed"): 40, ("protein", "tblastn", "placed"): 3,
              ("protein", "unmapped", "unexplained"): 7}  # fmt: skip
    links = _links(_sankey(tw.v05(counts, conflicts=2)))
    got: dict[tuple[str, str], float] = {}
    for a, b, v in links:
        got[a, b] = got.get((a, b), 0.0) + v
    want: dict[tuple[str, str], float] = {}
    for (c, r, o), n in counts.items():
        want[c, r] = want.get((c, r), 0.0) + n
        want[r, o] = want.get((r, o), 0.0) + n
    assert got == want
    assert all(v > 0 for _, _, v in links), "zero flows and the conflicts row are not links"


def test_sankey_node_and_link_colours_are_m9_05s() -> None:
    from gippa.viz.mapping_sankey import NODE_GREY, OUTCOME_COLOUR

    t = _sankey(tw.v05())
    names = _node_names(t)
    node_cols = [tw.to_hex(str(c)) for c in tw.values(t["node"]["color"])]
    want = [OUTCOME_COLOUR.get(n, NODE_GREY).upper() for n in names]
    assert node_cols == want
    link_cols = {tw.to_hex(str(c)) for c in tw.values(t["link"]["color"])}
    assert link_cols <= {c.upper() for c in OUTCOME_COLOUR.values()}
    assert len(link_cols) == 4  # placed, contaminant, rejected by FDR, unexplained


# ---------------------------------------------------------------- V7: example 2 (forest)

XLABEL = "observed \N{MINUS SIGN} expected detected genes"


def _v7_points(s: dict[str, Any]) -> list[dict[str, Any]]:
    """Every V7 point: trace index, x, y label, error arrays and hover."""
    out = []
    for i, t in enumerate(tw.traces(s, "scatter", "scattergl")):
        xs, ys, hv = tw.values(t.get("x")), tw.values(t.get("y")), tw.hovers(t)
        ex = t.get("error_x") or {}
        plus, minus = tw.values(ex.get("array")), tw.values(ex.get("arrayminus"))
        for k, (x, y) in enumerate(zip(xs, ys, strict=True)):
            out.append({"trace": i, "x": x, "y": str(y), "hover": hv[k] if hv else "",
                        "plus": plus[k] if plus else None,
                        "minus": minus[k] if minus else None})  # fmt: skip
    return out


def _top_to_bottom(s: dict[str, Any]) -> list[str]:
    ya = s["layout"].get("yaxis", {})
    arr = [str(v) for v in tw.values(ya.get("categoryarray"))]
    assert arr, "the y axis lists the rows in a categoryarray"
    return arr if ya.get("autorange") == "reversed" else arr[::-1]


def _point_of(points: list[dict[str, Any]], group: str) -> dict[str, Any]:
    hit = [p for p in points if p["y"].startswith(group) and p["x"] is not None]
    assert len(hit) == 1, (group, hit)
    return hit[0]


def test_example2_two_region_traces_and_prophage_1_hover_q() -> None:
    s = tw.spec("V07", tw.v07())
    pts = _v7_points(s)
    region = {p["trace"] for p in pts if p["y"].startswith("prophage_")}
    assert len(region) == 2
    hover = _point_of(pts, "prophage_1")["hover"]
    assert "prophage_1" in hover and "q 0.003" in hover
    assert "-4.75" in hover and "CI -7.1 to -2.2" in hover


def test_v7_values_and_error_bars_equal_the_tsv() -> None:
    data = tw.v07()
    pts = _v7_points(tw.spec("V07", data))
    for r in data.to_dict("records"):
        p = _point_of(pts, str(r["group"]))
        assert p["x"] == pytest.approx(r["o_minus_e"])
        if not (pd.isna(r["ci_low"]) or pd.isna(r["ci_high"])):
            assert p["plus"] == pytest.approx(r["ci_high"] - r["o_minus_e"])
            assert p["minus"] == pytest.approx(r["o_minus_e"] - r["ci_low"])


def test_v7_rows_in_drawing_order_and_axis_title() -> None:
    data = tw.v07()
    s = tw.spec("V07", data)
    order = _top_to_bottom(s)
    groups = [str(g) for g in data["group"]]
    assert len(order) == len(groups)
    assert all(lab.startswith(g) for lab, g in zip(order, groups, strict=True)), order
    assert tw.title_of(s["layout"]["xaxis"]) == XLABEL


def test_v7_region_colour_by_type_and_black_basic_rows() -> None:
    from gippa.viz.style import REGION

    s = tw.spec("V07", tw.v07())
    ts = tw.traces(s, "scatter", "scattergl")
    for t in ts:
        ys = [str(y) for y in tw.values(t.get("y"))]
        cols = {tw.to_hex(str(c)) for c in tw.values((t.get("marker") or {}).get("color"))}
        if ys and all(y.startswith("prophage_") for y in ys):
            assert cols == {REGION["prophage"].upper()}
        elif ys and any(y.startswith(("chr1", "cds")) for y in ys):
            assert cols <= {"#000000"}


def test_v7_na_ci_shows_ci_unestimated_and_na_q_shows_na() -> None:
    rows = [*tw.V07_ROWS[:2], ("region", "plasmid_1", "plasmid", 10, 4.0, 4.5, -0.5, None,
                               None, None, "ci_unestimated")]  # fmt: skip
    pts = _v7_points(tw.spec("V07", tw.v07(rows)))
    assert "CI unestimated" in _point_of(pts, "cds")["hover"]
    hover = _point_of(pts, "plasmid_1")["hover"]
    assert "CI unestimated" in hover and "q NA" in hover


_num = st.floats(min_value=-50, max_value=50, allow_nan=False).map(lambda v: round(v, 2))


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.tuples(_num, _num, _num, st.sampled_from([0.001, 0.02, 0.5, None])),
                min_size=1, max_size=8))  # fmt: skip
def test_property_v7_twin_points_equal_any_tsv(cells: list[tuple[float, ...]]) -> None:
    rows = []
    for i, (a, b, c, q) in enumerate(cells):
        lo, mid, hi = sorted((a, b, c))
        rows.append(("region", f"prophage_{i + 1}", "prophage", 10, 3.0, 3.0 - mid, mid, lo,
                     hi, q, None))  # fmt: skip
    data = tw.v07(rows)
    text = tw.figure_text("V07", data)
    assert text is not None and text == tw.figure_text("V07", data)
    pts = _v7_points(json.loads(text))
    for r in data.to_dict("records"):
        p = _point_of(pts, str(r["group"]))
        assert p["x"] == pytest.approx(r["o_minus_e"])
        q = r["q_value"]
        assert (f"q {q:g}" if not (q is None or math.isnan(q)) else "q NA") in p["hover"]
