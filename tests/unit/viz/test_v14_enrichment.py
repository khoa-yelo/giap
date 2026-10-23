"""V14 ``enrichment``: which functions are high or low? (M9-19; GIAP_SPEC.md §10.2 V14, §7.7.3-
§7.7.5, ADR-0025/0026/0027/0029/0035, §9.3 ``enrichment_gsea.tsv`` / ``enrichment_detection.tsv``;
M9-17 ``plots/support/gsea_running``; lead decisions 1-3 of 2026-10-23).

Tests the plotted data (``Instance.data`` = ``plots/data/V14_enrichment.tsv``), the three panels,
the ADR-0025/0026 markings, the absence of any CI on NES (ADR-0027) and on ORA (decision 2), the
decision-1 state (nothing significant: the 20 smallest q with the alpha line), the title, the
degenerate cases, examples 1-5, the §14.4 checklist, determinism, registration (plot, Plotly twin,
``plots`` step). No ``--mpl`` baselines (Khoa approves those at M9-25, as for V08 / V12).

Interface assumed where the spec leaves it open (M9-19 test author):

- Importing every ``giap.viz`` module (``giap.steps.plots.load_plot_modules``) registers ``V14``
  (name ``enrichment``, ``module="function"``, ``layers=()``: the no-resources case is an empty
  state, not a skip); ``enrichment_gsea`` and ``enrichment_detection`` are among its ``needs`` /
  ``optional``; with both tables absent it skips ``input_absent:...``. One ``Instance``,
  ``suffix=None``, ``coords=()``, double width (7.2 in).
- Panel axes carry ``Axes.get_label()``: the six running-score axes ``"a"`` (one per support set,
  in the support file's order read row-major: top to bottom, then left to right; nothing else is
  labelled ``"a"``), the dot plot ``"b"``, the ORA plot ``"c"`` (any extra panel-(c) axes, e.g. a
  q axis, has a label starting with ``"c"``). A panel-level message ("no abundance input",
  "running scores not available") is a visible text of the figure.
- (a): the running sum is one ``Line2D`` with x = ``rank_position`` and y = ``running_es``; the
  leading edge of a positive set is one filled artist (``PolyCollection`` from ``fill_between`` or
  a ``Rectangle``) whose x extent ends at the peak rank (± 0.5); the panel title holds the term,
  the name, "(<collection>)", "NES" and "q"; a flagged set has a text "not calibrated (<flag>)" in
  its axes and a dashed outline (a dashed spine, or a dashed patch of the axes or figure).
- (b): x = ``nes``, y = -log10 ``q_value``, one marker per row with finite ``nes`` and ``q``; the
  collection colours are six distinct Okabe-Ito hues (not black or grey), fixed per collection
  whatever subset is present; flagged rows hollow; redundant rows achromatic and smaller than any
  non-redundant marker; the alpha line is an ``axhline`` (dashed) at -log10 alpha; no error bars
  (no ``LineCollection``, no ``ErrorbarContainer``); a row label is an ``ax.texts`` entry holding
  its term or name. A q = 0 row is drawn at a finite height and a figure text mentions "cap".
- (c): x = O/E = ``stat / (stat - effect)`` on a log axis of base 2, one marker per drawn row; the
  x label is exactly "observed / expected detected TUs"; "CI unestimated" is a figure text; no
  error bars; with nothing at q ≤ alpha the alpha line is a reference line (``axhline`` /
  ``axvline``) at alpha or -log10 alpha in an axes labelled "c...".
- Data columns, exactly: ``panel, collection, term, name, direction, rank_position, running_es,
  in_set, nes, q_value, o_over_e, coverage, flags``. Panel (b) has a row for every table row with
  finite ``nes`` and ``q_value``; panel (c) exactly the drawn rows.
"""

from __future__ import annotations

import logging
import math
import random
import tempfile
from collections.abc import Sequence
from pathlib import Path

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import twins as tw
from helpers import viz_enrichment as ve
from helpers.viz_artists import Marker, markers, texts
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import LineCollection, PathCollection, PolyCollection
from matplotlib.container import ErrorbarContainer
from matplotlib.patches import Rectangle

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

COLUMNS = [
    "panel", "collection", "term", "name", "direction", "rank_position", "running_es", "in_set",
    "nes", "q_value", "o_over_e", "coverage", "flags",
]  # fmt: skip
OVER, UNDER = "#0072B2", "#D55E00"
XLABEL_C = "observed / expected detected TUs"
FOOT_GSEA = "GSEA is conditional on detection, with no p_detect adjustment"
FOOT_COVERAGE = "coverage = share of eligible genes with ≥ 1 term"
NO_RESOURCES = "No annotation resources installed (giap setup --go --kegg --cog)"
EX3_TITLE = (
    "3 sets enriched by abundance at q ≤ 0.05 (conditional on detection); "
    "2 categories detected more or less than expected"
)
CATEGORICAL = {c.upper() for k, c in style.OKABE_ITO.items() if k not in ("black", "grey")}

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- calling the plot -------------------------------------------------------------------------


def _v14() -> base.PlotSpec:
    from giap.steps.plots import load_plot_modules

    load_plot_modules()
    reg = base.registered()
    assert "V14" in reg, f"V14 is not registered: {sorted(reg)}"
    return reg["V14"]


def _result(out: Path, cfg: GiapConfig) -> list[Instance] | Skip:
    return _v14().fn(RunView.from_outdir(out), cfg)


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _result(out, cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _main(tmp: Path, cfg: GiapConfig, **kw: object) -> Instance:
    return _render(ve.main_out(tmp, **kw), cfg)


def _isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _num(v: object) -> float:
    return math.nan if _isna(v) else float(str(v))


def _bool(v: object) -> bool:
    return str(v).strip().lower() in ("true", "1", "1.0")


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"].astype(str) == panel].reset_index(drop=True)


def _axes(inst: Instance, label: str) -> list[Axes]:
    """Axes labelled ``label``, row-major (top to bottom, then left to right)."""
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    return sorted(found, key=lambda ax: (-round(ax.get_position().y0, 3), ax.get_position().x0))


def _one(inst: Instance, label: str) -> Axes:
    found = _axes(inst, label)
    assert len(found) == 1, [ax.get_label() for ax in inst.figure.axes]
    return found[0]


def _c_axes(inst: Instance) -> list[Axes]:
    return [ax for ax in inst.figure.axes if ax.get_label().startswith("c")]


def _ax_texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.texts if t.get_visible() and t.get_text()]


def _all_text(inst: Instance) -> str:
    return " ".join(" ".join(texts(inst.figure)).split())


def _hex(m: Marker) -> set[str]:
    """The marker's colours without black / white edges."""
    return m.colours() - {"#000000", "#FFFFFF"}


def _sized(ax: Axes) -> list[tuple[float, float, float, Marker]]:
    """``(x, y, area pt^2, marker)`` per drawn point (scatter ``s`` or ``markersize^2``)."""
    out = []
    for line in ax.get_lines():
        if line.get_marker() in (None, "", "None", " ", "none", "|", "_") or not line.get_visible():
            continue
        size = float(line.get_markersize()) ** 2
        for m in markers_of_line(ax, line):
            out.append((m.x, m.y, size, m))
    for coll in ax.collections:
        if not isinstance(coll, PathCollection) or not coll.get_visible():
            continue
        offs = np.asarray(coll.get_offsets(), dtype=float).reshape(-1, 2)
        sizes = np.asarray(coll.get_sizes(), dtype=float)
        ms = [m for m in markers(ax) if any(abs(m.x - x) < 1e-12 and abs(m.y - y) < 1e-12
                                             for x, y in offs)]  # fmt: skip
        for i, (x, y) in enumerate(offs):
            if not (np.isfinite(x) and np.isfinite(y)):
                continue
            s = float(sizes[i if len(sizes) > 1 else 0]) if len(sizes) else 36.0
            m = next(mm for mm in ms if abs(mm.x - x) < 1e-12 and abs(mm.y - y) < 1e-12)
            out.append((float(x), float(y), s, m))
    return out


def markers_of_line(ax: Axes, line: object) -> list[Marker]:
    xs = np.asarray(line.get_xdata(), dtype=float)  # type: ignore[attr-defined]
    ys = np.asarray(line.get_ydata(), dtype=float)  # type: ignore[attr-defined]
    pts = {(float(x), float(y)) for x, y in zip(xs, ys, strict=True)}
    return [m for m in markers(ax) if (m.x, m.y) in pts]


def _at(ax: Axes, x: float, y: float, tol: float = 1e-6) -> list[Marker]:
    return [m for m in markers(ax) if abs(m.x - x) <= tol and abs(m.y - y) <= tol]


def _no_error_bars(ax: Axes) -> None:
    assert not any(isinstance(c, ErrorbarContainer) for c in ax.containers)
    assert not any(isinstance(c, LineCollection) for c in ax.collections), "error bars drawn"


def _dashed(artist: object) -> bool:
    ls = artist.get_linestyle()  # type: ignore[attr-defined]
    if isinstance(ls, str):
        return ls in ("--", "dashed", "-.", "dashdot", ":", "dotted")
    if isinstance(ls, tuple) and len(ls) == 2:
        return bool(ls[1])
    return False


def _outline_dashed(inst: Instance, ax: Axes) -> bool:
    if any(s.get_visible() and _dashed(s) for s in ax.spines.values()):
        return True
    patches = [*ax.patches, *inst.figure.patches]
    return any(p.get_visible() and p.get_linewidth() > 0 and _dashed(p) for p in patches)


def _ref_lines(ax: Axes) -> list[tuple[str, float]]:
    """``("h", y)`` / ``("v", x)`` for every axhline / axvline of ``ax``."""
    out = []
    for ln in ax.lines:
        tr = ln.get_transform()
        if tr == ax.get_yaxis_transform():
            out.append(("h", float(np.asarray(ln.get_ydata(), dtype=float)[0])))
        elif tr == ax.get_xaxis_transform():
            out.append(("v", float(np.asarray(ln.get_xdata(), dtype=float)[0])))
    return out


def _x_extents(ax: Axes) -> list[tuple[float, float]]:
    """x extents (data units) of the filled artists of ``ax``."""
    out = []
    for c in ax.collections:
        if isinstance(c, PolyCollection) and c.get_visible():
            for p in c.get_paths():
                v = np.asarray(p.vertices, dtype=float)
                v = v[np.isfinite(v).all(axis=1)]
                if len(v):
                    out.append((float(v[:, 0].min()), float(v[:, 0].max())))
    for p in ax.patches:
        if isinstance(p, Rectangle) and p.get_visible() and p.get_width() > 0:
            x0 = float(p.get_x())
            out.append((x0, x0 + float(p.get_width())))
    return out


def _sets() -> list[dict[str, object]]:
    return ve.main_sets()


# --- registration and skips -------------------------------------------------------------------


def test_v14_is_registered_as_enrichment() -> None:
    spec = _v14()
    assert spec.pid == "V14"
    assert spec.name == "enrichment"
    assert spec.module == "function"
    assert spec.layers == ()
    declared = {*spec.needs, *spec.optional}
    assert {"enrichment_gsea", "enrichment_detection"} <= declared


def test_function_module_off_is_skipped(tmp_path: Path, cfg: GiapConfig) -> None:
    out = ve.main_out(tmp_path, modules=["core"])
    assert _result(out, cfg) == Skip("module_off:function")


def test_no_enrichment_tables_is_input_absent(tmp_path: Path, cfg: GiapConfig) -> None:
    result = _result(ve.build_out(tmp_path, None, None), cfg)
    assert isinstance(result, Skip)
    assert result.reason.startswith("input_absent:enrichment_")


def test_plotly_twin_is_registered() -> None:
    assert tw.load().has_twin("V14")


# --- data TSV ---------------------------------------------------------------------------------


def test_columns_and_instance(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    assert list(inst.data.columns) == COLUMNS
    assert set(inst.data["panel"].astype(str)) == {"a", "b", "c"}


def test_panel_a_rows_copy_the_support(tmp_path: Path, cfg: GiapConfig) -> None:
    a = _rows(_main(tmp_path, cfg), "a")
    assert len(a) == 6 * ve.N_TUS
    for s in _sets():
        rows = a[(a["collection"] == s["collection"]) & (a["term"] == s["term"])]
        rows = rows.sort_values("rank_position", key=lambda c: c.map(_num))
        assert [int(_num(v)) for v in rows["rank_position"]] == list(range(1, ve.N_TUS + 1))
        got = [_num(v) for v in rows["running_es"]]
        assert got == pytest.approx([ve.g6(v) for v in s["running"]], rel=1e-5, abs=1e-9)  # type: ignore[attr-defined]
        members = set(s["members"])  # type: ignore[call-overload]
        assert [_bool(v) for v in rows["in_set"]] == [i in members for i in range(1, ve.N_TUS + 1)]


def test_panel_b_rows_are_every_finite_gsea_row(tmp_path: Path, cfg: GiapConfig) -> None:
    b = _rows(_main(tmp_path, cfg), "b")
    finite = [r for r in ve.main_gsea() if not _isna(r["nes"]) and not _isna(r["q_value"])]
    assert len(finite) == 8
    keys = {(str(c), str(t)) for c, t in zip(b["collection"], b["term"], strict=True)}
    assert keys == {(str(r["collection"]), str(r["term"])) for r in finite}
    by = b.set_index("term")
    for r in finite:
        assert _num(by.loc[r["term"], "nes"]) == pytest.approx(r["nes"], rel=1e-5)
        assert _num(by.loc[r["term"], "q_value"]) == pytest.approx(r["q_value"], rel=1e-5)
    assert "genomic_cluster" in str(by.loc["GO:0006412", "flags"])


def test_example_1_o_over_e(tmp_path: Path, cfg: GiapConfig) -> None:
    """ORA row stat 12, effect 4: O/E = 12 / 8 = 1.5."""
    c = _rows(_main(tmp_path, cfg), "c")
    j = c[(c["collection"] == "cog") & (c["term"] == "J")]
    assert len(j) == 1
    assert _num(j.iloc[0]["o_over_e"]) == pytest.approx(1.5)
    assert str(j.iloc[0]["direction"]) == "over"
    assert _num(j.iloc[0]["coverage"]) == pytest.approx(0.75)


def test_panel_c_rows_are_the_significant_ones(tmp_path: Path, cfg: GiapConfig) -> None:
    c = _rows(_main(tmp_path, cfg), "c")
    keys = sorted(zip(c["collection"], c["term"], c["direction"], strict=True))
    assert keys == [("cog", "J", "over"), ("kegg", "ko00010", "under")]
    k = c.set_index("term")
    assert _num(k.loc["ko00010", "o_over_e"]) == pytest.approx(0.25)
    assert "ci_unestimated" in str(k.loc["ko00010", "flags"])


def test_data_sorted_by_panel_collection_term_direction_rank(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    d = _main(tmp_path, cfg).data
    order = {c: i for i, c in enumerate(ve.COLLECTIONS)}

    def key(i: int) -> tuple[object, ...]:
        r = d.iloc[i]
        rank = _num(r["rank_position"])
        return (str(r["panel"]), order[str(r["collection"])], str(r["term"]),
                "" if _isna(r["direction"]) else str(r["direction"]),
                -1.0 if math.isnan(rank) else rank)  # fmt: skip

    keys = [key(i) for i in range(len(d))]
    assert keys == sorted(keys)


# --- panel (a) --------------------------------------------------------------------------------


def test_panel_a_six_small_multiples_in_support_order(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    axes = _axes(inst, "a")
    assert len(axes) == 6
    for ax, s in zip(axes, _sets(), strict=True):
        title = ax.get_title()
        for part in (str(s["term"]), str(s["name"]), f"({s['collection']})", "NES", "q"):
            assert part in title, (part, title)


def test_panel_a_running_sum_line(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    for ax, s in zip(_axes(inst, "a"), _sets(), strict=True):
        want = [ve.g6(v) for v in s["running"]]  # type: ignore[attr-defined]
        found = False
        for ln in ax.get_lines():
            x = np.asarray(ln.get_xdata(), dtype=float)
            y = np.asarray(ln.get_ydata(), dtype=float)
            ok = np.isfinite(x) & np.isfinite(y)
            x, y = x[ok], y[ok]
            pts = {round(a): float(b) for a, b in zip(x, y, strict=True)}
            if all(i in pts for i in range(1, ve.N_TUS + 1)):
                got = [pts[i] for i in range(1, ve.N_TUS + 1)]
                if np.allclose(got, want, rtol=1e-5, atol=1e-9):
                    found = True
        assert found, f"no running-sum line for {s['term']}"


def test_panel_a_leading_edge_shaded_up_to_the_peak(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    for ax, s in zip(_axes(inst, "a"), _sets(), strict=True):
        peak = ve.peak_rank(s["running"])  # type: ignore[arg-type]
        if float(s["es"]) <= 0 or peak < 3:  # type: ignore[arg-type]
            continue  # negative sets: the leading edge is after the peak; rank 1: no width
        ext = _x_extents(ax)
        assert any(abs(hi - peak) <= 0.5 and lo <= 1.5 for lo, hi in ext), (s["term"], peak, ext)


@pytest.mark.parametrize(("term", "flag"), [("GO:0006412", "genomic_cluster"),
                                            ("GO:0005840", "detection_confounded")])  # fmt: skip
def test_example_2_flagged_set_banner_and_dashed_outline(
    tmp_path: Path, cfg: GiapConfig, term: str, flag: str
) -> None:
    inst = _main(tmp_path, cfg)
    axes = dict(zip((str(s["term"]) for s in _sets()), _axes(inst, "a"), strict=True))
    assert f"not calibrated ({flag})" in " ".join(_ax_texts(axes[term]))
    assert _outline_dashed(inst, axes[term])
    plain = axes["ko00010"]
    assert not any("not calibrated" in t for t in _ax_texts(plain))
    assert not any(s.get_visible() and _dashed(s) for s in plain.spines.values())


def test_fewer_than_six_support_sets(tmp_path: Path, cfg: GiapConfig) -> None:
    sets = _sets()[:2]
    out = ve.build_out(tmp_path, ve.main_gsea(), ve.main_ora(), ve.running_rows(sets))
    inst = _render(out, cfg)
    with_lines = [ax for ax in _axes(inst, "a") if ax.get_visible() and any(
        len(ln.get_xdata()) >= ve.N_TUS for ln in ax.get_lines())]  # fmt: skip
    assert len(with_lines) == 2
    assert len(_rows(inst, "a")) == 2 * ve.N_TUS


def test_plain_tsv_support_file_is_read(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg, gz=False)
    assert len(_rows(inst, "a")) == 6 * ve.N_TUS


# --- examples 4 and 5: the support must match the table --------------------------------------


def _one_set_out(tmp: Path, peak: float, es: float) -> Path:
    r = [0.2, peak, 0.4, 0.1]
    rows = [["kegg", "ko00010", i, f"TU{i:04d}", 1.0 - 0.5 * i, i in (1, 2), v]
            for i, v in enumerate(r, start=1)]  # fmt: skip
    gsea = [ve.gsea_row("kegg", "ko00010", name="Glycolysis", es=es, nes=1.9, p=0.001, q=0.002)]
    return ve.build_out(tmp, gsea, ve.main_ora(), rows)


def test_example_4_support_peak_equal_to_es_is_accepted(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_one_set_out(tmp_path, 0.62, 0.62), cfg)
    assert len(_rows(inst, "a")) == 4


def test_example_5_support_peak_differs_from_es_raises(tmp_path: Path, cfg: GiapConfig) -> None:
    with pytest.raises(ValueError):
        _result(_one_set_out(tmp_path, 0.60, 0.62), cfg)


def test_mismatch_on_a_negative_set_raises(tmp_path: Path, cfg: GiapConfig) -> None:
    gsea = ve.main_gsea()
    j = next(r for r in gsea if r["term"] == "J")
    j["es"] = float(j["es"]) - 0.05  # type: ignore[arg-type]
    out = ve.build_out(tmp_path, gsea, ve.main_ora(), ve.running_rows(_sets()))
    with pytest.raises(ValueError):
        _result(out, cfg)


# --- panel (b) --------------------------------------------------------------------------------


def _b_point(r: dict[str, object]) -> tuple[float, float]:
    return float(r["nes"]), -math.log10(float(r["q_value"]))  # type: ignore[arg-type]


def test_panel_b_one_marker_per_finite_row(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "b")
    finite = [r for r in ve.main_gsea() if not _isna(r["nes"]) and not _isna(r["q_value"])]
    pts = sorted((round(m.x, 6), round(m.y, 6)) for m in markers(ax))
    assert pts == sorted((round(x, 6), round(y, 6)) for x, y in map(_b_point, finite))


def test_panel_b_colours_by_collection_fixed_and_distinct(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path / "all", cfg), "b")
    colour: dict[str, str] = {}
    for r in ve.main_gsea():
        if _isna(r["nes"]) or r["redundant"]:
            continue
        (m,) = _at(ax, *_b_point(r))
        hexes = _hex(m)
        assert len(hexes) == 1 and hexes <= CATEGORICAL, (r["term"], hexes)
        c = hexes.pop()
        assert colour.setdefault(str(r["collection"]), c) == c
    assert set(colour) == set(ve.COLLECTIONS) - {"go_mf"}  # go_mf's only row is redundant
    assert len(set(colour.values())) == len(colour)
    # the same collection keeps its colour when the others are absent (fixed, not by appearance)
    only = [ve.gsea_row("go_slim", "GO:0008150", nes=1.5, q=0.01)]
    ax2 = _one(_render(ve.build_out(tmp_path / "one", only, ve.main_ora()), cfg), "b")
    (m,) = _at(ax2, 1.5, -math.log10(0.01))
    assert _hex(m) == {colour["go_slim"]}


def test_panel_b_flagged_rows_hollow(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "b")
    for r in ve.main_gsea():
        if _isna(r["nes"]):
            continue
        (m,) = _at(ax, *_b_point(r))
        flagged = r["flags"] in ("genomic_cluster", "detection_confounded")
        assert m.hollow == flagged, r["term"]


def test_panel_b_redundant_small_and_grey(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "b")
    red = next(r for r in ve.main_gsea() if r["redundant"])
    x, y = _b_point(red)
    sized = _sized(ax)
    (rs,) = [s for s in sized if abs(s[0] - x) < 1e-6 and abs(s[1] - y) < 1e-6]
    others = [s for s in sized if not (abs(s[0] - x) < 1e-6 and abs(s[1] - y) < 1e-6)]
    assert all(rs[2] < o[2] for o in others), (rs[2], [o[2] for o in others])
    assert rs[3].colours() and all(mcolors.to_hex(c)[1:3] == mcolors.to_hex(c)[3:5]
                                   == mcolors.to_hex(c)[5:7] for c in rs[3].colours())  # fmt: skip


def test_panel_b_alpha_line(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "b")
    want = -math.log10(0.05)
    lines = [ln for ln in ax.lines if ln.get_transform() == ax.get_yaxis_transform()]
    assert any(abs(float(np.asarray(ln.get_ydata())[0]) - want) < 1e-9 and _dashed(ln)
               for ln in lines), _ref_lines(ax)  # fmt: skip


def test_panel_b_alpha_line_follows_the_config(tmp_path: Path) -> None:
    ax = _one(_main(tmp_path, resolve_config(overrides=["stats.alpha.enrichment=0.1"])), "b")
    assert ("h", pytest.approx(1.0)) in _ref_lines(ax)


def test_panel_b_no_error_bars(tmp_path: Path, cfg: GiapConfig) -> None:
    """ADR-0027: NES has no CI."""
    _no_error_bars(_one(_main(tmp_path, cfg), "b"))


def test_panel_b_labels_ten_most_significant_non_redundant(tmp_path: Path, cfg: GiapConfig) -> None:
    rows = [ve.gsea_row("kegg", f"ko9{i:04d}", name=f"pathway number {i:02d}", nes=1.0 + i / 10,
                        q=0.001 * (i + 1)) for i in range(14)]  # fmt: skip
    rows[0]["redundant"] = True  # the most significant row is redundant: never labelled
    inst = _render(ve.build_out(tmp_path, rows, ve.main_ora()), cfg)
    labels = " | ".join(_ax_texts(_one(inst, "b")))

    def labelled(i: int) -> bool:
        return f"ko9{i:04d}" in labels or f"pathway number {i:02d}" in labels

    assert not labelled(0)
    assert all(labelled(i) for i in range(1, 11)), labels
    assert not any(labelled(i) for i in range(11, 14)), labels


def test_q_zero_is_capped_and_noted(tmp_path: Path, cfg: GiapConfig) -> None:
    rows = [
        ve.gsea_row("kegg", "ko00010", nes=2.0, p=0.0, q=0.0),
        ve.gsea_row("cog", "J", nes=1.0, p=0.001, q=0.004),
    ]
    inst = _render(ve.build_out(tmp_path, rows, ve.main_ora()), cfg)
    ax = _one(inst, "b")
    ms = markers(ax)
    assert len(ms) == 2  # markers() drops non-finite points: the q = 0 row is drawn finite
    top = [m for m in ms if abs(m.x - 2.0) < 1e-9]
    assert len(top) == 1 and math.isfinite(top[0].y)
    assert top[0].y >= -math.log10(0.004)
    assert "cap" in _all_text(inst).lower()


# --- panel (c) --------------------------------------------------------------------------------


def test_panel_c_points_on_a_log2_axis(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "c")
    assert ax.get_xscale() == "log"
    assert ax.xaxis.get_transform().base == 2  # type: ignore[attr-defined]
    assert sorted(round(m.x, 9) for m in markers(ax)) == [0.25, 1.5]
    assert ax.get_xlabel() == XLABEL_C


def test_panel_c_coloured_by_direction(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _one(_main(tmp_path, cfg), "c")
    by_x = {round(m.x, 9): m for m in markers(ax)}
    assert OVER in by_x[1.5].colours()
    assert UNDER in by_x[0.25].colours()


def test_decision_2_ora_points_without_bars_labelled_ci_unestimated(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    inst = _main(tmp_path, cfg)
    for ax in _c_axes(inst):
        _no_error_bars(ax)
    assert "CI unestimated" in _all_text(inst)


def test_panel_c_coverage_printed(tmp_path: Path, cfg: GiapConfig) -> None:
    ora = ve.main_ora()
    ora[0]["coverage"], ora[1]["coverage"] = 0.62, 0.31
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), ora, ve.running_rows(_sets())), cfg)
    found = " ".join(t for ax in _c_axes(inst) for t in _ax_texts(ax))
    for cov in (0.62, 0.31):
        assert f"{cov:.2f}" in found or f"{round(cov * 100)} %" in found, found


def test_only_under_rows_significant(tmp_path: Path, cfg: GiapConfig) -> None:
    ora = [
        ve.ora_row("cog", "C", "under", 3, -3, 0.01),
        ve.ora_row("cog", "E", "under", 1, -4, 0.02),
        ve.ora_row("cog", "J", "over", 12, 4, 0.5),
    ]
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), ora, ve.running_rows(_sets())), cfg)
    ax = _one(inst, "c")
    ms = markers(ax)
    assert sorted(round(m.x, 9) for m in ms) == [0.2, 0.5]
    assert all(UNDER in m.colours() and OVER not in m.colours() for m in ms)


def _ora_none_significant(n: int = 25) -> list[dict[str, object]]:
    return [ve.ora_row(ve.COLLECTIONS[i % 6], f"T{i:03d}", "over" if i % 2 else "under",
                       5 + i % 4, 1 if i % 2 else -2, round(0.06 + 0.03 * ((i * 7) % n), 4))
            for i in range(n)]  # fmt: skip


def test_decision_1_nothing_significant_shows_the_20_smallest_q(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    ora = _ora_none_significant()
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), ora, ve.running_rows(_sets())), cfg)
    c = _rows(inst, "c")
    want = sorted(ora, key=lambda r: float(r["q_value"]))[:20]  # type: ignore[arg-type]
    assert sorted(c["term"].astype(str)) == sorted(str(r["term"]) for r in want)
    assert len(markers(_one(inst, "c"))) == 20


def test_decision_1_alpha_line_in_panel_c(tmp_path: Path, cfg: GiapConfig) -> None:
    ora = _ora_none_significant()
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), ora, ve.running_rows(_sets())), cfg)
    refs = [v for ax in _c_axes(inst) for _, v in _ref_lines(ax)]
    assert any(abs(v - 0.05) < 1e-9 or abs(v + math.log10(0.05)) < 1e-9 for v in refs), refs


# --- footnote and title -----------------------------------------------------------------------


def test_footnotes(tmp_path: Path, cfg: GiapConfig) -> None:
    text = _all_text(_main(tmp_path, cfg))
    assert FOOT_GSEA in text
    assert FOOT_COVERAGE in text


def test_example_3_title(tmp_path: Path, cfg: GiapConfig) -> None:
    """3 non-redundant GSEA rows and 1 redundant at q <= 0.05; 2 ORA rows."""
    assert _main(tmp_path, cfg).title == EX3_TITLE


def test_title_follows_alpha(tmp_path: Path) -> None:
    inst = _main(tmp_path, resolve_config(overrides=["stats.alpha.enrichment=0.25"]))
    assert inst.title == (
        "4 sets enriched by abundance at q ≤ 0.25 (conditional on detection); "
        "2 categories detected more or less than expected"
    )


def test_title_when_nothing_passes(tmp_path: Path, cfg: GiapConfig) -> None:
    gsea = [ve.gsea_row("kegg", "ko00020", nes=0.9, q=0.5)]
    inst = _render(ve.build_out(tmp_path, gsea, _ora_none_significant()), cfg)
    assert inst.title == "No set or category passes q ≤ 0.05"
    assert "no enrichment" not in inst.title.lower()


# --- degenerate cases -------------------------------------------------------------------------


def test_no_resources_is_an_empty_state(tmp_path: Path, cfg: GiapConfig) -> None:
    from giap.viz import checklist

    flags = ["layer_missing:go", "layer_missing:kegg", "layer_missing:cog"]
    inst = _render(ve.build_out(tmp_path, [], [], flags=flags), cfg)
    assert inst.title == NO_RESOURCES
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0
    assert checklist.is_empty_state(inst)


def test_no_abundance_panels_a_b_note_and_c_drawn(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(ve.build_out(tmp_path, [], ve.main_ora()), cfg)
    assert "no abundance input" in _all_text(inst)
    assert sorted(round(m.x, 9) for m in markers(_one(inst, "c"))) == [0.25, 1.5]
    assert set(inst.data["panel"].astype(str)) == {"c"}
    assert "2 categories detected more or less than expected" in inst.title


def test_support_absent_panel_a_note(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), ve.main_ora()), cfg)
    assert "running scores not available" in _all_text(inst)
    assert len(markers(_one(inst, "b"))) == 8
    assert len(markers(_one(inst, "c"))) == 2
    assert "a" not in set(inst.data["panel"].astype(str))
    assert inst.title == EX3_TITLE


def test_header_only_detection_with_gsea(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(ve.build_out(tmp_path, ve.main_gsea(), [], ve.running_rows(_sets())), cfg)
    assert len(markers(_one(inst, "b"))) == 8
    assert "c" not in set(inst.data["panel"].astype(str))


def test_double_width(tmp_path: Path, cfg: GiapConfig) -> None:
    assert _main(tmp_path, cfg).figure.get_figwidth() == pytest.approx(style.DOUBLE_IN)


# --- the §14.4 checklist ----------------------------------------------------------------------


def _check(inst: Instance, tmp: Path) -> list[str]:
    from giap.viz import checklist

    files = write_instance(inst, "V14", "enrichment", tmp, ["svg", "png"], 100)
    return [f"{r.item}: {r.detail}" for r in checklist.check_instance(inst, "V14", files)
            if r.passed is False]  # fmt: skip


@pytest.mark.parametrize("case", ["main", "nothing_significant", "no_abundance", "no_support",
                                  "no_resources"])  # fmt: skip
def test_checklist_gated_items_pass(tmp_path: Path, cfg: GiapConfig, case: str) -> None:
    run = tmp_path / "run"
    if case == "main":
        out = ve.main_out(run)
    elif case == "nothing_significant":
        out = ve.build_out(run, [ve.gsea_row("kegg", "ko00020", nes=0.9, q=0.5)],
                           _ora_none_significant())  # fmt: skip
    elif case == "no_abundance":
        out = ve.build_out(run, [], ve.main_ora())
    elif case == "no_support":
        out = ve.build_out(run, ve.main_gsea(), ve.main_ora())
    else:
        out = ve.build_out(run, [], [], flags=["layer_missing:go", "layer_missing:kegg",
                                               "layer_missing:cog"])  # fmt: skip
    assert _check(_render(out, cfg), tmp_path / "plots") == []


# --- determinism ------------------------------------------------------------------------------


def test_two_renders_write_identical_files(tmp_path: Path, cfg: GiapConfig) -> None:
    out = ve.main_out(tmp_path / "run")
    blobs = []
    for name in ("x", "y"):
        paths = write_instance(_render(out, cfg), "V14", "enrichment", tmp_path / name,
                               ["svg", "png"], 100)  # fmt: skip
        blobs.append({p.name: p.read_bytes() for p in paths})
    assert blobs[0] == blobs[1]
    tsv = blobs[0]["V14_enrichment.tsv"].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)


def test_input_row_order_does_not_change_the_data(tmp_path: Path, cfg: GiapConfig) -> None:
    a = _main(tmp_path / "a", cfg)
    gsea, ora = ve.main_gsea(), ve.main_ora()
    rng = random.Random(7)
    rng.shuffle(gsea)
    rng.shuffle(ora)
    b = _render(ve.build_out(tmp_path / "b", gsea, ora, ve.running_rows(_sets())), cfg)
    pd.testing.assert_frame_equal(a.data.reset_index(drop=True), b.data.reset_index(drop=True))
    assert a.title == b.title


# --- the Plotly twin --------------------------------------------------------------------------


def test_twin_builds_from_the_data_tsv(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _main(tmp_path, cfg)
    frame = tw.instance_frame(inst.data)
    s = tw.spec("V14", frame, inst.title)
    assert tw.traces(s)
    assert tw.disallowed(s) == set()
    for t in tw.traces(s):
        for key in ("error_x", "error_y"):
            assert not (t.get(key) or {}).get("visible", bool(t.get(key))), (key, t.get(key))
    assert tw.figure_text("V14", frame, inst.title) == tw.figure_text("V14", frame, inst.title)


def test_twin_header_only_is_none(tmp_path: Path, cfg: GiapConfig) -> None:
    flags = ["layer_missing:go", "layer_missing:kegg", "layer_missing:cog"]
    inst = _render(ve.build_out(tmp_path, [], [], flags=flags), cfg)
    assert tw.figure_text("V14", tw.instance_frame(inst.data), inst.title) is None


# --- the plots step ---------------------------------------------------------------------------


def test_plots_step_renders_v14_with_the_support_file(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    from helpers.report_fixtures import cfg as small_cfg
    from helpers.report_fixtures import write_preflight

    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    _v14()
    out = ve.main_out(tmp_path / "run")
    write_preflight(out)
    inputs = {
        "preflight": out / "preflight" / "preflight.json",
        "enrichment_gsea": out / "tables" / "enrichment_gsea.tsv",
        "enrichment_detection": out / "tables" / "enrichment_detection.tsv",
        "plots_support": out / "plots" / "support",
    }
    plots = tmp_path / "work" / "plots"
    caplog.set_level(logging.INFO, logger="giap.steps")
    reg.run_step("plots", small_cfg(), inputs, {"plots": plots})
    assert (plots / "V14_enrichment.svg").is_file()
    data = pd.read_csv(plots / "data" / "V14_enrichment.tsv", sep="\t", dtype=str,
                       keep_default_na=False)  # fmt: skip
    assert list(data.columns) == COLUMNS
    assert (data["panel"] == "a").sum() == 6 * ve.N_TUS
    rendered = [r for r in caplog.records if getattr(r, "event", None) == "plots.rendered"]
    assert any(dict(r.data).get("plot") == "V14" for r in rendered)  # type: ignore[attr-defined]


# --- properties -------------------------------------------------------------------------------

_ORA = st.lists(
    st.tuples(
        st.sampled_from(ve.COLLECTIONS),
        st.integers(1, 30),  # stat (k detected TUs, >= 1 so O/E > 0 on the log axis)
        st.integers(-10, 29),  # effect, clipped below so that E = stat - effect >= 1
        st.sampled_from([0.001, 0.01, 0.04, 0.05, 0.051, 0.2, 0.5, 0.9]),
        st.booleans(),  # over
    ),
    min_size=1,
    max_size=30,
)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_ORA, n_gsea_sig=st.integers(0, 3))
def test_property_panel_c_selection_o_over_e_and_title(
    rows: Sequence[tuple[str, int, int, float, bool]], n_gsea_sig: int
) -> None:
    style.apply()
    ora = []
    for i, (coll, stat, effect, q, over) in enumerate(rows):
        eff = min(effect, stat - 1)
        ora.append(ve.ora_row(coll, f"T{i:03d}", "over" if over else "under", stat, eff,
                              round(q + i * 1e-6, 7)))  # fmt: skip
    gsea = [ve.gsea_row("kegg", f"ko{i:05d}", nes=1.0 + i, q=0.01) for i in range(n_gsea_sig)]
    gsea.append(ve.gsea_row("cog", "K", nes=0.3, q=0.9))
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(ve.build_out(Path(tmp), gsea, ora), resolve_config())
    sig = [r for r in ora if float(r["q_value"]) <= 0.05]  # type: ignore[arg-type]
    want = sig or sorted(ora, key=lambda r: float(r["q_value"]))[:20]  # type: ignore[arg-type]
    c = _rows(inst, "c").set_index("term")
    assert sorted(c.index.astype(str)) == sorted(str(r["term"]) for r in want)
    for r in want:
        stat, eff = float(r["stat"]), float(r["effect"])  # type: ignore[arg-type]
        assert _num(c.loc[r["term"], "o_over_e"]) == pytest.approx(stat / (stat - eff), rel=1e-5)
    if not sig and not n_gsea_sig:
        assert inst.title == "No set or category passes q ≤ 0.05"
    else:
        assert inst.title.startswith(f"{n_gsea_sig} set")
        assert f"; {len(sig)} categor" in inst.title
