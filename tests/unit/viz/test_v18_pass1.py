"""V18 ``stat_diagnostics`` pass 1: (a) reliability, (b) p-value histograms (M9-07; GIAP_SPEC.md
§10.2 V18, §7.3.3, §8.1, §8.2, §9.2).

Tests the plotted data (``Instance.data`` = ``plots/data/V18_stat_diagnostics.tsv``), the axes
structure, the title, the family table, the bin edges and the degenerate cases. No ``--mpl``
baselines (M9-25 owns those).

Interface assumed where the spec leaves it open:

- The plot lives in ``giap.viz.diagnostics``; importing it registers ``V18`` (name
  ``stat_diagnostics``, ``needs=()`` because a missing ``detection_summary`` still draws panel (b),
  ``optional`` ⊇ {``detection_summary``, ``regions_stats``, ``autocorrelation``, ``pairs``,
  ``enrichment_detection``, ``enrichment_gsea``}, ``module="core"``). The tests call
  ``registered()["V18"].fn(view, cfg)``: a list of one ``Instance``, ``suffix=None``,
  ``coords=()``.
- ``RunView.PATHS`` gains ``enrichment_detection`` -> ``tables/enrichment_detection.tsv`` and
  ``enrichment_gsea`` -> ``tables/enrichment_gsea.tsv`` (M8-06a writes them; they are not yet keys).
- The module exposes the family table as its single source: ``FAMILIES``, a sequence of
  ``(table, key_columns)`` in the spec's order, and ``NOT_DRAWN``, the tables with a ``p_value``
  that are deliberately not histogrammed (⊇ {``islands``}).
- Panel (a) is the axes labelled ``"a"``. Each histogram of panel (b) is the axes labelled
  ``"b:<table> / <family>"``; with no family, one axes labelled ``"b"`` carries the text. No other
  axes label starts with ``"b"``. A histogram axes title holds the table name and ``"n = <n>"``.
- ``family`` in the data TSV is the key values joined by ``" / "`` (``o_minus_e``,
  ``morans_i / gene_label``, ``cog / over``); ``table`` is the table name.
- Data rows: panel (a) bin rows have ``bin`` = the bin number (1..10), ``x`` = mean_p,
  ``y`` = observed_rate, ``n`` = the bin's ``n`` metric, ``metric``/``value`` NA; overall rows have
  ``metric`` and ``value`` and ⊇ {brier, calibration_slope, calibration_intercept}. Panel (b): 20
  rows per family (zero counts included), ``bin`` = left edge k * 0.05, ``n`` = count, ``x``,
  ``y``, ``metric``, ``value`` NA. A family with no finite p has no rows.
- Bins: a p equal to a decimal edge k/20 (the permutation p-values of B = 19, 39, ...) lands in
  the bin with that left edge; 1 lands in the last bin.
- Marker area is proportional to ``n``.
"""

from __future__ import annotations

import math
import os
import random
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from helpers.viz_out_tree import ISLAND_ROW, SCHEMAS, make_out, schema_fields, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.collections import PathCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.plot_files import write_instance
from giap.io.run_view import PATHS, RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "m4"
VAL_OUT = Path(os.environ.get("GIAP_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = ["panel", "table", "family", "bin", "x", "y", "n", "metric", "value"]
TABLE_ORDER = [
    "regions_stats",
    "autocorrelation",
    "pairs",
    "enrichment_detection",
    "enrichment_gsea",
]
FAMILY_KEYS: dict[str, tuple[str, ...]] = {
    "regions_stats": ("stat",),
    "autocorrelation": ("statistic", "null_model"),
    "pairs": ("test",),
    "enrichment_detection": ("collection", "direction"),
    "enrichment_gsea": ("collection",),
}
OPTIONAL = {"detection_summary", *TABLE_ORDER}
CONSTANT_TEXT = "detectability constant: no calibration curve"
NO_P_TEXT = "no p-values in this run"
COLLECTIONS = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def _ds(group_type: str, group: object, metric: str, value: float | None) -> dict[str, object]:
    return {
        "section": "detectability",
        "group_type": group_type,
        "group": str(group),
        "metric": metric,
        "value": value,
    }


def _bins(points: Sequence[tuple[float, float, int]]) -> list[dict[str, object]]:
    """(mean_p, observed_rate, n) per bin 1.. in the producer's long format."""
    rows = []
    for i, (mp, rate, n) in enumerate(points, start=1):
        rows += [_ds("bin", i, "mean_p", mp), _ds("bin", i, "observed_rate", rate)]
        rows.append(_ds("bin", i, "n", n))
    return rows


def _overall(
    slope: float | None, intercept: float | None, brier: float | None, constant: bool = False
) -> list[dict[str, object]]:
    rows = [
        _ds("overall", "all", "brier", brier),
        _ds("overall", "all", "calibration_intercept", intercept),
        _ds("overall", "all", "calibration_slope", slope),
        _ds("overall", "all", "recalibration_separated", 0.0),
        _ds("overall", "all", "n_detected", 40.0),
        _ds("overall", "all", "n_eligible", 100.0),
    ]
    if constant:
        rows.append(_ds("overall", "all", "detectability_constant", 1.0))
    return rows


#: example 1
EX1_SUMMARY = [*_bins([(0.1, 0.12, 50), (0.6, 0.55, 50)]), *_overall(0.97, 0.02, 0.11)]
BASIC_ROW = {
    "section": "basic",
    "group_type": "replicon",
    "group": "chr1",
    "metric": "n_genes",
    "value": 10.0,
}


def _region(p: float | None, stat: str = "o_minus_e", subset: str = "all", i: int = 0) -> dict:
    return {
        "region_id": f"R{i:04d}",
        "region_type": "prophage",
        "subset": subset,
        "contig": "chr1",
        "start": 101 + 1000 * i,
        "end": 900 + 1000 * i,
        "n_genes": 5,
        "stat": stat,
        "p_value": p,
        "method": "x",
        "null_model": "circular_shift",
    }


def _autocorr(p: float | None, statistic: str, null_model: str, lag: int = 1) -> dict:
    return {
        "contig": "chr1",
        "statistic": statistic,
        "lag": lag,
        "null_model": null_model,
        "value": 0.1,
        "p_value": p,
    }


def _pair(p: float | None, test: str, i: int = 0) -> dict:
    return {
        "pair_id": f"P{i:04d}",
        "class": "OP",
        "gene_a": f"a{i}",
        "gene_b": f"b{i}",
        "strand": "+",
        "intergenic_bp": 20,
        "test": test,
        "p_value": p,
    }


def _ora(p: float | None, collection: str, direction: str, i: int = 0) -> dict:
    return {
        "collection": collection,
        "term": f"T{i:04d}",
        "direction": direction,
        "p_value": p,
        "method": "cpb_tu_ora",
        "null_model": "conditional_poisson_binomial_tu",
    }


def _gsea(p: float | None, collection: str, i: int = 0) -> dict:
    return {
        "collection": collection,
        "term": f"T{i:04d}",
        "p_value": p,
        "method": "fgsea_simple",
        "null_model": "tu_set_permutation",
    }


_SCHEMA = {
    "detection_summary": "detection_summary",
    "regions_stats": "regions_stats",
    "autocorrelation": "autocorrelation",
    "pairs": "pairs",
    "enrichment_detection": "enrichment_detection",
    "enrichment_gsea": "enrichment_gsea",
    "islands": "islands",
}


def _out(tmp: Path, **tables: Sequence[Mapping[str, object]] | None) -> Path:
    out = make_out(tmp)
    for name, rows in tables.items():
        if rows is not None:
            write_schema_table(out / "tables" / f"{name}.tsv", _SCHEMA[name], rows)
    return out


# --- calling the plot -------------------------------------------------------------------------


def _v18() -> base.PlotSpec:
    import giap.viz.diagnostics  # noqa: F401  (registers V18)

    return base.registered()["V18"]


def _render(out: Path, cfg: GiapConfig) -> Instance:
    result = _v18().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _rows(inst: Instance, panel: str) -> pd.DataFrame:
    return inst.data[inst.data["panel"] == panel].reset_index(drop=True)


def _num(value: object) -> float:
    return float(str(value))


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, f"axes {label!r}: {[a.get_label() for a in inst.figure.axes]}"
    return found[0]


def _hist_axes(inst: Instance) -> dict[str, Axes]:
    """``"<table> / <family>"`` -> histogram axes."""
    return {ax.get_label()[2:]: ax for ax in inst.figure.axes if ax.get_label().startswith("b:")}


def _texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _drawn(ax: Axes) -> list[Rectangle]:
    return [
        p
        for p in ax.patches
        if isinstance(p, Rectangle) and p.get_width() > 0 and p.get_height() > 0
    ]


def _bar_counts(ax: Axes) -> dict[int, float]:
    """bin index (0..19, from the bar centre) -> drawn height."""
    out: dict[int, float] = {}
    for r in _drawn(ax):
        k = math.floor((r.get_x() + r.get_width() / 2) / 0.05)
        out[k] = out.get(k, 0.0) + r.get_height()
    return out


def _family_counts(inst: Instance, table: str, family: str) -> dict[int, int]:
    """bin index -> count from the data rows of one family (all 20 bins)."""
    b = _rows(inst, "b")
    sub = b[(b["table"] == table) & (b["family"] == family)]
    return {round(_num(e) / 0.05): int(_num(n)) for e, n in zip(sub["bin"], sub["n"], strict=True)}


def _families(inst: Instance) -> list[tuple[str, str]]:
    b = _rows(inst, "b")
    return list(dict.fromkeys(zip(b["table"].astype(str), b["family"].astype(str), strict=True)))


def _metric(inst: Instance, metric: str) -> float:
    a = _rows(inst, "a")
    hit = a[a["metric"].astype(str) == metric]
    assert len(hit) == 1, f"metric {metric}: {len(hit)} rows\n{a}"
    return _num(hit.iloc[0]["value"])


def _is_grey(colour: Any) -> bool:
    r, g, b, _ = mcolors.to_rgba(colour)
    return abs(r - g) < 1e-6 and abs(g - b) < 1e-6 and 0.2 < r < 0.9


def _expected_bin(p: float) -> int:
    """The bin of a p away from the edges: exact arithmetic, 1 in the last bin."""
    return min(19, math.floor(Fraction(p) * 20))


# --- registration, family table ----------------------------------------------------------------


def test_v18_is_registered_as_stat_diagnostics() -> None:
    spec = _v18()
    assert spec.pid == "V18"
    assert spec.name == "stat_diagnostics"
    assert spec.needs == ()
    assert set(spec.optional) >= OPTIONAL
    assert spec.module == "core"


def test_run_view_knows_the_enrichment_tables() -> None:
    assert PATHS["enrichment_detection"] == "tables/enrichment_detection.tsv"
    assert PATHS["enrichment_gsea"] == "tables/enrichment_gsea.tsv"


def test_family_table_is_the_spec_table() -> None:
    import giap.viz.diagnostics as diag

    got = [(str(t), tuple(k)) for t, k in diag.FAMILIES]
    assert got == [(t, FAMILY_KEYS[t]) for t in TABLE_ORDER]
    assert "islands" in set(diag.NOT_DRAWN)
    assert not set(diag.NOT_DRAWN) & set(TABLE_ORDER)


def test_every_table_with_p_value_is_a_family_or_not_drawn() -> None:
    """A new test table cannot be forgotten silently (acceptance criterion)."""
    import giap.viz.diagnostics as diag

    with_p = {
        p.name.split(".")[0]
        for p in SCHEMAS.glob("*.schema.json")
        if any(f["name"] == "p_value" for f in schema_fields(p.name.split(".")[0]))
    }
    assert with_p >= {"regions_stats", "islands"}  # the glob found the schemas
    drawn = {str(t) for t, _ in diag.FAMILIES}
    forgotten = with_p - drawn - set(diag.NOT_DRAWN)
    assert not forgotten, f"tables with p_value in neither FAMILIES nor NOT_DRAWN: {forgotten}"
    for table, keys in diag.FAMILIES:
        cols = {f["name"] for f in schema_fields(str(table))}
        assert "p_value" in cols, table
        assert set(keys) <= cols, (table, keys)


# --- structure, title -------------------------------------------------------------------------


def _full_out(tmp: Path) -> Path:
    return _out(
        tmp,
        detection_summary=[BASIC_ROW, *EX1_SUMMARY],
        regions_stats=[_region(p, i=i) for i, p in enumerate([0.01, 0.04, 0.5, 1.0])]
        + [_region(0.3, "hl_shift", i=9)],
        autocorrelation=[
            _autocorr(0.2, "morans_i", "tu_block"),
            _autocorr(0.7, "morans_i", "gene_label"),
            _autocorr(0.05, "join_count_bb", "gene_label"),
        ],
        pairs=[_pair(0.02, "both_detected_t", 0), _pair(None, "none", 1)],
        enrichment_detection=[_ora(0.9, "kegg", "under", 0), _ora(0.1, "cog", "over", 1)],
        enrichment_gsea=[_gsea(0.33, "go_bp")],
    )


def test_one_instance_with_exact_columns(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_full_out(tmp_path), cfg)
    assert inst.suffix is None
    assert inst.coords == ()
    assert list(inst.data.columns) == COLUMNS
    assert set(inst.data["panel"]) == {"a", "b"}
    width = round(float(inst.figure.get_size_inches()[0]), 6)
    assert width in (style.SINGLE_IN, style.DOUBLE_IN)


def test_example1_title_and_bin_rows(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, detection_summary=EX1_SUMMARY), cfg)
    expected = "Detectability calibration slope 0.97 (intercept 0.02), Brier 0.11"
    assert inst.title == expected
    assert expected in [t.get_text() for t in inst.figure.findobj(Text)]
    a = _rows(inst, "a")
    bins = a[a["metric"].map(_isna)]
    assert len(bins) == 2
    got = [
        (int(_num(b)), _num(x), _num(y), int(_num(n)))
        for b, x, y, n in zip(bins["bin"], bins["x"], bins["y"], bins["n"], strict=True)
    ]
    assert got == [(1, 0.1, 0.12, 50), (2, 0.6, 0.55, 50)]
    assert all(_isna(v) for v in bins["value"])
    assert _metric(inst, "calibration_slope") == pytest.approx(0.97)
    assert _metric(inst, "calibration_intercept") == pytest.approx(0.02)
    assert _metric(inst, "brier") == pytest.approx(0.11)


def test_title_numbers_use_two_significant_digits(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.2, 0.3, 10), (0.7, 0.6, 10)]), *_overall(1.234, -0.0456, 0.1)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert inst.title == "Detectability calibration slope 1.2 (intercept -0.046), Brier 0.1"


def test_constant_title_counts_drawn_families(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.4, 0.4, 87)]), *_overall(None, None, 0.24, constant=True)]
    regions = [_region(0.2, i=0), _region(0.6, "hl_shift", i=1)]
    inst = _render(_out(tmp_path, detection_summary=summary, regions_stats=regions), cfg)
    assert inst.title == "Detectability constant; 2 test families shown"


def test_titles_never_say_calibrated_or_express(tmp_path: Path, cfg: GiapConfig) -> None:
    cases = [
        {"detection_summary": EX1_SUMMARY},
        {"detection_summary": _overall(None, None, 0.2, constant=True)},
        {"regions_stats": [_region(0.5)]},
    ]
    for i, tables in enumerate(cases):
        inst = _render(_out(tmp_path / str(i), **tables), cfg)
        assert "calibrated" not in inst.title.lower()
        assert "express" not in inst.title.lower()


def test_every_text_is_at_least_6_pt(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_full_out(tmp_path), cfg)
    sizes = [t.get_fontsize() for t in inst.figure.findobj(Text) if t.get_text()]
    assert sizes
    assert min(sizes) >= 6


# --- panel (a) --------------------------------------------------------------------------------


def test_reliability_axes_labels_and_limits(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _ax(_render(_out(tmp_path, detection_summary=EX1_SUMMARY), cfg), "a")
    assert ax.get_xlabel() == "mean predicted p_detect"
    assert ax.get_ylabel() == "observed detection rate"
    assert tuple(round(v, 6) for v in ax.get_xlim()) == (0.0, 1.0)
    assert tuple(round(v, 6) for v in ax.get_ylim()) == (0.0, 1.0)


def test_reliability_diagonal_is_grey_dashed(tmp_path: Path, cfg: GiapConfig) -> None:
    ax = _ax(_render(_out(tmp_path, detection_summary=EX1_SUMMARY), cfg), "a")
    diag = [
        ln
        for ln in ax.get_lines()
        if [round(float(v), 6) for v in ln.get_xdata()] == [0.0, 1.0]
        and [round(float(v), 6) for v in ln.get_ydata()] == [0.0, 1.0]
    ]
    assert len(diag) == 1
    assert diag[0].get_linestyle() == "--"
    assert _is_grey(diag[0].get_color())


def _joined_line(ax: Axes, xs: Sequence[float]) -> Line2D:
    want = [round(x, 6) for x in xs]
    found = [ln for ln in ax.get_lines() if [round(float(v), 6) for v in ln.get_xdata()] == want]
    assert len(found) == 1, [list(ln.get_xdata()) for ln in ax.get_lines()]
    return found[0]


def test_points_are_joined_in_numeric_bin_order(tmp_path: Path, cfg: GiapConfig) -> None:
    """Ten bins: bin 10 is last (numeric, not string order '1', '10', '2')."""
    pts = [(0.05 + 0.09 * i, 0.04 + 0.09 * i, 10 + i) for i in range(10)]
    rows = _bins(pts)
    random.Random(7).shuffle(rows)  # the producer's row order is not relied on
    inst = _render(_out(tmp_path, detection_summary=[*rows, *_overall(1.0, 0.0, 0.2)]), cfg)
    line = _joined_line(_ax(inst, "a"), [p[0] for p in pts])
    assert [round(float(v), 6) for v in line.get_ydata()] == [round(p[1], 6) for p in pts]
    a = _rows(inst, "a")
    bins = [int(_num(b)) for b, m in zip(a["bin"], a["metric"], strict=True) if _isna(m)]
    assert bins == list(range(1, 11))


def test_marker_area_is_proportional_to_n(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.2, 0.25, 10), (0.8, 0.7, 40)]), *_overall(0.9, 0.05, 0.15)]
    ax = _ax(_render(_out(tmp_path, detection_summary=summary), cfg), "a")
    colls = [c for c in ax.collections if isinstance(c, PathCollection) and len(c.get_offsets())]
    assert len(colls) == 1
    offsets = [(round(float(x), 6), round(float(y), 6)) for x, y in colls[0].get_offsets()]
    assert offsets == [(0.2, 0.25), (0.8, 0.7)]
    sizes = [float(s) for s in colls[0].get_sizes()]
    assert len(sizes) == 2
    assert sizes[1] / sizes[0] == pytest.approx(4.0)


def test_overall_numbers_are_in_the_corner(tmp_path: Path, cfg: GiapConfig) -> None:
    text = " ".join(_texts(_ax(_render(_out(tmp_path, detection_summary=EX1_SUMMARY), cfg), "a")))
    for word, value in (("slope", "0.97"), ("intercept", "0.02"), ("Brier", "0.11")):
        assert word in text and value in text, text


def test_example4_constant_model_draws_the_text(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.0, 0.0, 9)] * 10), *_overall(None, None, 0.0, constant=True)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    ax = _ax(inst, "a")
    assert CONSTANT_TEXT in _texts(ax)
    assert not [c for c in ax.collections if isinstance(c, PathCollection) and len(c.get_offsets())]
    assert inst.title == "Detectability constant; 0 test families shown"


def test_constant_value_zero_is_not_constant(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*EX1_SUMMARY, _ds("overall", "all", "detectability_constant", 0.0)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert CONSTANT_TEXT not in _texts(_ax(inst, "a"))
    assert inst.title.startswith("Detectability calibration slope 0.97")


def test_other_sections_are_ignored(tmp_path: Path, cfg: GiapConfig) -> None:
    other = {**BASIC_ROW, "group_type": "bin", "group": "1", "metric": "mean_p", "value": 0.5}
    inst = _render(_out(tmp_path, detection_summary=[other, *EX1_SUMMARY]), cfg)
    a = _rows(inst, "a")
    assert len(a[a["metric"].map(_isna)]) == 2


def test_detection_summary_absent_panel_a_not_available(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, regions_stats=[_region(0.01), _region(0.5, i=1)]), cfg)
    assert any("not available" in t for t in _texts(_ax(inst, "a")))
    assert len(_rows(inst, "a")) == 0
    assert list(_hist_axes(inst)) == ["regions_stats / o_minus_e"]
    assert _family_counts(inst, "regions_stats", "o_minus_e")[0] == 1


# --- panel (b) --------------------------------------------------------------------------------


def test_example2_regions_histogram(tmp_path: Path, cfg: GiapConfig) -> None:
    regions = [_region(p, i=i) for i, p in enumerate([0.01, 0.04, 0.5, 1.0])]
    inst = _render(_out(tmp_path, detection_summary=EX1_SUMMARY, regions_stats=regions), cfg)
    assert _families(inst) == [("regions_stats", "o_minus_e")]
    counts = _family_counts(inst, "regions_stats", "o_minus_e")
    assert counts == {k: {0: 2, 10: 1, 19: 1}.get(k, 0) for k in range(20)}
    b = _rows(inst, "b")
    assert [round(_num(e), 6) for e in b["bin"]] == [round(0.05 * k, 6) for k in range(20)]
    for col in ("x", "y", "metric", "value"):
        assert all(_isna(v) for v in b[col]), col
    ax = _hist_axes(inst)["regions_stats / o_minus_e"]
    assert _bar_counts(ax) == {0: 2.0, 10: 1.0, 19: 1.0}
    assert "regions_stats" in ax.get_title()
    assert "n = 4" in ax.get_title()
    assert tuple(round(v, 6) for v in ax.get_xlim()) == (0.0, 1.0)


@pytest.mark.parametrize("k", range(21))
def test_decimal_bin_edges_are_closed_on_the_left(tmp_path: Path, cfg: GiapConfig, k: int) -> None:
    """p = k/20 as written by a permutation test with B + 1 = 20 lands in bin k; 1 in bin 19."""
    p = float(f"{k / 20:.2f}")
    inst = _render(_out(tmp_path, regions_stats=[_region(p)]), cfg)
    counts = _family_counts(inst, "regions_stats", "o_minus_e")
    assert counts[min(k, 19)] == 1, (p, counts)
    assert sum(counts.values()) == 1


def test_one_p_is_a_one_count_histogram(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, pairs=[_pair(0.42, "one_detected_pmiss")]), cfg)
    ax = _hist_axes(inst)["pairs / one_detected_pmiss"]
    assert _bar_counts(ax) == {8: 1.0}
    assert "n = 1" in ax.get_title()


def test_p_zero_is_drawn_in_the_first_bin(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, regions_stats=[_region(0.0), _region(0.3, i=1)]), cfg)
    assert _family_counts(inst, "regions_stats", "o_minus_e")[0] == 1
    assert _bar_counts(_hist_axes(inst)["regions_stats / o_minus_e"])[0] == 1.0


def test_family_rows_filters(tmp_path: Path, cfg: GiapConfig) -> None:
    """regions: subset = all only; pairs: test none excluded even with a p; NA p dropped."""
    regions = [
        _region(0.1, i=0),
        _region(0.2, subset="integrase", i=1),
        _region(None, i=2),
        _region(None, "hl_shift", i=3),
    ]
    pairs = [_pair(0.3, "none", 0), _pair(0.6, "missed_operon_decoy", 1), _pair(None, "none", 2)]
    inst = _render(_out(tmp_path, regions_stats=regions, pairs=pairs), cfg)
    assert _families(inst) == [("regions_stats", "o_minus_e"), ("pairs", "missed_operon_decoy")]
    assert sum(_family_counts(inst, "regions_stats", "o_minus_e").values()) == 1
    assert sum(_family_counts(inst, "pairs", "missed_operon_decoy").values()) == 1
    assert set(_hist_axes(inst)) == {
        "regions_stats / o_minus_e",
        "pairs / missed_operon_decoy",
    }


def test_composite_family_keys_and_sort_order(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_full_out(tmp_path), cfg)
    assert _families(inst) == [
        ("regions_stats", "hl_shift"),
        ("regions_stats", "o_minus_e"),
        ("autocorrelation", "join_count_bb / gene_label"),
        ("autocorrelation", "morans_i / gene_label"),
        ("autocorrelation", "morans_i / tu_block"),
        ("pairs", "both_detected_t"),
        ("enrichment_detection", "cog / over"),
        ("enrichment_detection", "kegg / under"),
        ("enrichment_gsea", "go_bp"),
    ]
    assert len(_hist_axes(inst)) == 9
    data = inst.data
    assert list(data["panel"]) == sorted(data["panel"])
    b = _rows(inst, "b")
    keys = [
        (TABLE_ORDER.index(str(t)), str(f), _num(e))
        for t, f, e in zip(b["table"], b["family"], b["bin"], strict=True)
    ]
    assert keys == sorted(keys)
    assert len(b) == 9 * 20


def test_example3_islands_are_never_a_histogram(tmp_path: Path, cfg: GiapConfig) -> None:
    island = {**ISLAND_ROW, "stat": "scan_max", "p_value": 0.03}
    inst = _render(_out(tmp_path, detection_summary=EX1_SUMMARY, islands=[island]), cfg)
    assert "islands" not in set(inst.data["table"].dropna().astype(str))
    assert _hist_axes(inst) == {}
    assert NO_P_TEXT in _texts(_ax(inst, "b"))


def test_no_finite_p_draws_the_empty_text(tmp_path: Path, cfg: GiapConfig) -> None:
    acf = [_autocorr(None, "morans_i", "gene_label", lag) for lag in range(1, 5)]
    inst = _render(_out(tmp_path, detection_summary=EX1_SUMMARY, autocorrelation=acf), cfg)
    assert len(_rows(inst, "b")) == 0
    assert _hist_axes(inst) == {}
    ax = _ax(inst, "b")
    assert NO_P_TEXT in _texts(ax)
    assert _drawn(ax) == []
    assert len(_rows(inst, "a")) > 0


def test_enrichment_tables_are_read_when_present(tmp_path: Path, cfg: GiapConfig) -> None:
    ora = [_ora(0.01 * (i + 1), "go_bp", "over", i) for i in range(5)]
    gsea = [_gsea(0.5, "kegg", 0), _gsea(0.97, "kegg", 1)]
    inst = _render(_out(tmp_path, enrichment_detection=ora, enrichment_gsea=gsea), cfg)
    assert _families(inst) == [
        ("enrichment_detection", "go_bp / over"),
        ("enrichment_gsea", "kegg"),
    ]
    assert _family_counts(inst, "enrichment_detection", "go_bp / over")[0] == 4
    assert _family_counts(inst, "enrichment_gsea", "kegg") == {
        k: {10: 1, 19: 1}.get(k, 0) for k in range(20)
    }


def _twenty_families() -> tuple[dict[str, list[dict]], list[tuple[str, str]]]:
    """20 families with distinct sizes 1..20 (in the listing order below)."""
    fams: list[tuple[str, str, Any]] = []
    for c in COLLECTIONS:
        for d in ("over", "under"):
            fams.append(("enrichment_detection", f"{c} / {d}", (c, d)))
    for c in COLLECTIONS:
        fams.append(("enrichment_gsea", c, c))
    fams += [("regions_stats", "o_minus_e", "o_minus_e"), ("regions_stats", "hl_shift", "hl_shift")]
    tables: dict[str, list[dict]] = {
        "enrichment_detection": [],
        "enrichment_gsea": [],
        "regions_stats": [],
    }
    i = 0
    for size, (table, _, key) in enumerate(fams, start=1):
        for j in range(size):
            p = (j + 0.5) / (size + 1)
            if table == "enrichment_detection":
                tables[table].append(_ora(p, key[0], key[1], i))
            elif table == "enrichment_gsea":
                tables[table].append(_gsea(p, key, i))
            else:
                tables[table].append(_region(p, key, i=i))
            i += 1
    return tables, [(t, f) for t, f, _ in fams]


def test_example5_twenty_families_draw_the_16_largest(tmp_path: Path, cfg: GiapConfig) -> None:
    tables, fams = _twenty_families()
    inst = _render(_out(tmp_path, detection_summary=EX1_SUMMARY, **tables), cfg)
    assert set(_families(inst)) == set(fams)
    assert len(_rows(inst, "b")) == 20 * 20
    drawn = set(_hist_axes(inst))
    assert drawn == {f"{t} / {f}" for t, f in fams[4:]}
    for size, (t, f) in enumerate(fams, start=1):
        assert sum(_family_counts(inst, t, f).values()) == size


def test_grid_has_at_most_four_columns(tmp_path: Path, cfg: GiapConfig) -> None:
    tables, _ = _twenty_families()
    inst = _render(_out(tmp_path, **tables), cfg)
    x0 = {round(float(ax.get_position().x0), 3) for ax in _hist_axes(inst).values()}
    assert 1 <= len(x0) <= 4
    assert not [ax for ax in inst.figure.axes if ax.get_label() == "b"]


def test_constant_title_with_twenty_families_counts_16(tmp_path: Path, cfg: GiapConfig) -> None:
    tables, _ = _twenty_families()
    summary = _overall(None, None, 0.2, constant=True)
    inst = _render(_out(tmp_path, detection_summary=summary, **tables), cfg)
    assert inst.title == "Detectability constant; 16 test families shown"


# --- properties -------------------------------------------------------------------------------

_P_AWAY_FROM_EDGES = st.floats(min_value=0.0, max_value=1.0).filter(
    lambda p: p in (0.0, 1.0) or abs(p * 20 - round(p * 20)) > 1e-6
)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    ps=st.lists(st.one_of(st.none(), _P_AWAY_FROM_EDGES), min_size=1, max_size=40),
    seed=st.integers(0, 2**16),
)
def test_property_counts_match_exact_bins_and_ignore_row_order(
    ps: list[float | None], seed: int
) -> None:
    style.apply()
    cfg = resolve_config()
    rows = [_region(p, i=i) for i, p in enumerate(ps)]
    shuffled = list(rows)
    random.Random(seed).shuffle(shuffled)
    finite = [p for p in ps if p is not None]
    with tempfile.TemporaryDirectory() as tmp:
        first = _render(_out(Path(tmp) / "x", regions_stats=rows), cfg)
        second = _render(_out(Path(tmp) / "y", regions_stats=shuffled), cfg)
    pd.testing.assert_frame_equal(first.data, second.data)
    if not finite:
        assert len(_rows(first, "b")) == 0
        return
    counts = _family_counts(first, "regions_stats", "o_minus_e")
    expected = {k: 0 for k in range(20)}
    for p in finite:
        expected[_expected_bin(p)] += 1
    assert counts == expected
    assert sum(_bar_counts(_hist_axes(first)["regions_stats / o_minus_e"]).values()) == len(finite)


@settings(max_examples=6, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(sizes=st.permutations(list(range(1, 21))), n_fam=st.integers(1, 20))
def test_property_top_16_by_size(sizes: list[int], n_fam: int) -> None:
    style.apply()
    pool = [("enrichment_detection", c, d) for c in COLLECTIONS for d in ("over", "under")]
    pool += [("enrichment_gsea", c, None) for c in COLLECTIONS]
    pool += [
        ("autocorrelation", "morans_i", "gene_label"),
        ("autocorrelation", "morans_i", "tu_block"),
    ]
    chosen = list(zip(pool[:n_fam], sizes[:n_fam], strict=True))
    tables: dict[str, list[dict]] = {
        "enrichment_detection": [],
        "enrichment_gsea": [],
        "autocorrelation": [],
    }
    labels = {}
    i = 0
    for (table, a, b), size in chosen:
        label = f"{a} / {b}" if b is not None else a
        labels[f"{table} / {label}"] = size
        for j in range(size):
            p = (j + 1) / (size + 1)
            if table == "enrichment_detection":
                tables[table].append(_ora(p, a, b, i))
            elif table == "enrichment_gsea":
                tables[table].append(_gsea(p, a, i))
            else:
                tables[table].append(_autocorr(p, a, b, lag=j + 1))
            i += 1
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), **tables), resolve_config())
    top = sorted(labels, key=labels.__getitem__, reverse=True)[:16]
    assert set(_hist_axes(inst)) == set(top)
    assert {f"{t} / {f}" for t, f in _families(inst)} == set(labels)
    assert inst.title.startswith("Detectability")


# --- files, determinism, fixtures ---------------------------------------------------------------


def test_data_tsv_and_images_are_byte_identical(tmp_path: Path, cfg: GiapConfig) -> None:
    out = _full_out(tmp_path / "run")
    written = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(
            inst, "V18", "stat_diagnostics", tmp_path / name, ["svg", "png", "pdf"], 100
        )
        written.append({p.relative_to(tmp_path / name): p.read_bytes() for p in paths})
    assert set(written[0]) == {
        Path("V18_stat_diagnostics.svg"),
        Path("V18_stat_diagnostics.png"),
        Path("V18_stat_diagnostics.pdf"),
        Path("data/V18_stat_diagnostics.tsv"),
    }
    assert written[0] == written[1]
    tsv = written[0][Path("data/V18_stat_diagnostics.tsv")].decode("utf-8").splitlines()
    assert tsv[0] == "\t".join(COLUMNS)
    assert "b\tregions_stats\to_minus_e\t0\tNA\tNA\t2\tNA\tNA" in tsv
    assert "b\tregions_stats\to_minus_e\t0.95\tNA\tNA\t1\tNA\tNA" in tsv


def test_render_does_not_mutate_between_calls(tmp_path: Path, cfg: GiapConfig) -> None:
    view = RunView.from_outdir(_full_out(tmp_path))
    first = _v18().fn(view, cfg)
    second = _v18().fn(view, cfg)
    assert isinstance(first, list) and isinstance(second, list)
    pd.testing.assert_frame_equal(first[0].data, second[0].data)
    assert first[0].title == second[0].title


def test_renders_on_the_sim50k_golden_tables(tmp_path: Path, cfg: GiapConfig) -> None:
    out = make_out(tmp_path)
    for name in ("detection_summary", "regions_stats"):
        shutil.copyfile(GOLDEN / f"{name}.tsv", out / "tables" / f"{name}.tsv")
    inst = _render(out, cfg)
    assert CONSTANT_TEXT in _texts(_ax(inst, "a"))  # the golden fit is constant
    assert inst.title == "Detectability constant; 2 test families shown"
    assert _families(inst) == [("regions_stats", "hl_shift"), ("regions_stats", "o_minus_e")]
    for fam in ("hl_shift", "o_minus_e"):
        assert sum(_family_counts(inst, "regions_stats", fam).values()) == 2


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_renders_on_the_validation_run(tmp_path: Path, cfg: GiapConfig) -> None:
    """Example 4 and the 'every autocorrelation p NA' case on the real validation run."""
    inst = _render(VAL_OUT, cfg)
    assert CONSTANT_TEXT in _texts(_ax(inst, "a"))
    assert NO_P_TEXT in _texts(_ax(inst, "b"))
    assert len(_rows(inst, "b")) == 0
    assert inst.title == "Detectability constant; 0 test families shown"
    paths = write_instance(inst, "V18", "stat_diagnostics", tmp_path / "plots", ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)
