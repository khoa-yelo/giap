"""V18 ``stat_diagnostics``: review round of PR #203 (M9-07; GIAP_SPEC.md §10.2 V18, §8.2).

Pins what the review found unpinned (lead decisions 2-6 of
docs/specs/M9-07-v18-diagnostics-pass1.md) and the lead's rulings on the review:

- **B1, p outside [0, 1].** Never clamped, never silently dropped. Per family, values outside
  [0, 1] (including +-inf) are left out of the bars; one panel-b data row ``metric =
  n_outside_0_1``, ``value = k`` is written only when k > 0; the histogram's axes title gains
  ``", k outside [0, 1]"``. p == 1 is in bin 19; NaN is missing and dropped.
- **N3, lead decision 6 wins.** With no bin rows and ``detectability_constant`` not 1, panel (a)
  shows "detectability not available", not the constant text.

p-values are written to the TSVs as their ``repr`` strings so that the exact float reaches the plot
(the default TSV writer rounds to 6 significant digits).
"""

from __future__ import annotations

import math
import random
import tempfile
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path

import pandas as pd
import pytest
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip

COLLECTIONS = ("cog", "kegg", "go_bp", "go_mf", "go_cc", "go_slim")
CONSTANT_TEXT = "detectability constant: no calibration curve"
ABSENT_TEXT = "detectability not available"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GiapConfig:
    return resolve_config()


# --- builders (self-contained; mirror test_v18_pass1) -----------------------------------------


def _ds(
    group_type: str, group: object, metric: str, value: float | None, section: str = "detectability"
) -> dict[str, object]:
    return {
        "section": section,
        "group_type": group_type,
        "group": str(group),
        "metric": metric,
        "value": value,
    }


def _bins(points: Sequence[tuple[float, float, int]]) -> list[dict[str, object]]:
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


EX1_SUMMARY = [*_bins([(0.1, 0.12, 50), (0.6, 0.55, 50)]), *_overall(0.97, 0.02, 0.11)]
BASIC_ROW = {
    "section": "basic",
    "group_type": "replicon",
    "group": "chr1",
    "metric": "n_genes",
    "value": 10.0,
}


def _region(p: object, stat: str = "o_minus_e", subset: str = "all", i: int = 0) -> dict:
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


def _autocorr(p: object, statistic: str, null_model: str, lag: int = 1) -> dict:
    return {
        "contig": "chr1",
        "statistic": statistic,
        "lag": lag,
        "null_model": null_model,
        "value": 0.1,
        "p_value": p,
    }


def _ora(p: object, collection: str, direction: str, i: int = 0) -> dict:
    return {
        "collection": collection,
        "term": f"T{i:04d}",
        "direction": direction,
        "p_value": p,
        "method": "cpb_tu_ora",
        "null_model": "conditional_poisson_binomial_tu",
    }


def _gsea(p: object, collection: str, i: int = 0) -> dict:
    return {
        "collection": collection,
        "term": f"T{i:04d}",
        "p_value": p,
        "method": "fgsea_simple",
        "null_model": "tu_set_permutation",
    }


def _exact(p: float) -> str:
    """The cell that carries ``p`` exactly (``repr``; inf as the writer's ``INF``)."""
    if math.isinf(p):
        return "INF" if p > 0 else "-INF"
    return repr(float(p))


def _out(tmp: Path, **tables: Sequence[Mapping[str, object]] | None) -> Path:
    out = make_out(tmp)
    for name, rows in tables.items():
        if rows is not None:
            write_schema_table(out / "tables" / f"{name}.tsv", name, rows)
    return out


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


def _hist_labels(inst: Instance) -> list[str]:
    """Histogram axes, in the order they were added to the figure."""
    return [ax.get_label()[2:] for ax in inst.figure.axes if ax.get_label().startswith("b:")]


def _texts(ax: Axes) -> list[str]:
    return [t.get_text() for t in ax.findobj(Text) if t.get_text()]


def _bar_counts(ax: Axes) -> dict[int, float]:
    out: dict[int, float] = {}
    for r in ax.patches:
        if isinstance(r, Rectangle) and r.get_width() > 0 and r.get_height() > 0:
            k = math.floor((r.get_x() + r.get_width() / 2) / 0.05)
            out[k] = out.get(k, 0.0) + r.get_height()
    return out


def _bin_rows(inst: Instance, table: str, family: str) -> pd.DataFrame:
    b = _rows(inst, "b")
    sub = b[(b["table"] == table) & (b["family"] == family)]
    return sub[sub["metric"].map(_isna)]


def _family_counts(inst: Instance, table: str, family: str) -> dict[int, int]:
    sub = _bin_rows(inst, table, family)
    return {round(_num(e) / 0.05): int(_num(n)) for e, n in zip(sub["bin"], sub["n"], strict=True)}


def _outside_rows(inst: Instance) -> pd.DataFrame:
    b = _rows(inst, "b")
    return b[b["metric"].astype(str) == "n_outside_0_1"].reset_index(drop=True)


# --- B1: p outside [0, 1] ------------------------------------------------------------------------

B1_PS = [0.01, 1.0, 1.5, -0.2, math.inf, -math.inf]


def _b1_out(tmp: Path) -> Path:
    regions = [_region(_exact(p), i=i) for i, p in enumerate(B1_PS)]
    regions.append(_region(None, i=len(B1_PS)))  # NaN: missing, not outside
    regions.append(_region(0.3, "hl_shift", i=99))  # a clean family next to it
    return _out(tmp, regions_stats=regions)


def test_b1_outside_values_are_not_in_the_bars(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_b1_out(tmp_path), cfg)
    counts = _family_counts(inst, "regions_stats", "o_minus_e")
    assert counts == {k: {0: 1, 19: 1}.get(k, 0) for k in range(20)}, counts
    ax = _ax(inst, "b:regions_stats / o_minus_e")
    assert _bar_counts(ax) == {0: 1.0, 19: 1.0}
    assert "n = 2" in ax.get_title()


def test_b1_one_n_outside_row_per_affected_family(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_b1_out(tmp_path), cfg)
    out = _outside_rows(inst)
    assert len(out) == 1, out
    row = out.iloc[0]
    assert (row["table"], row["family"]) == ("regions_stats", "o_minus_e")
    assert _num(row["value"]) == 4
    for col in ("x", "y"):
        assert _isna(row[col]), col
    # the clean family has no such row and still its 20 bin rows
    assert len(_bin_rows(inst, "regions_stats", "hl_shift")) == 20
    assert len(_bin_rows(inst, "regions_stats", "o_minus_e")) == 20


def test_b1_histogram_title_counts_outside(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_b1_out(tmp_path), cfg)
    title = _ax(inst, "b:regions_stats / o_minus_e").get_title()
    assert title.endswith(", 4 outside [0, 1]"), title
    clean = _ax(inst, "b:regions_stats / hl_shift").get_title()
    assert "outside" not in clean, clean


def test_b1_no_outside_row_when_every_p_is_in_range(tmp_path: Path, cfg: GiapConfig) -> None:
    regions = [_region(p, i=i) for i, p in enumerate([0.0, 0.5, 1.0])]
    inst = _render(_out(tmp_path, regions_stats=regions), cfg)
    assert len(_outside_rows(inst)) == 0
    assert "outside" not in _ax(inst, "b:regions_stats / o_minus_e").get_title()
    assert _family_counts(inst, "regions_stats", "o_minus_e")[19] == 1


@pytest.mark.parametrize(
    "p", [-0.2, -1e-300, -5e-324, 1.0000000000000002, 1.5, 7.0, math.inf, -math.inf]
)
def test_b1_low_and_high_side_are_both_outside(tmp_path: Path, cfg: GiapConfig, p: float) -> None:
    """No clamp on either side: -0.2 is neither bin 0 nor (by negative indexing) bin 16."""
    regions = [_region(_exact(p), i=0), _region(0.52, i=1)]
    inst = _render(_out(tmp_path, regions_stats=regions), cfg)
    counts = _family_counts(inst, "regions_stats", "o_minus_e")
    assert counts == {k: int(k == 10) for k in range(20)}, (p, counts)
    out = _outside_rows(inst)
    assert len(out) == 1 and _num(out.iloc[0]["value"]) == 1


def test_b1_counts_are_per_family(tmp_path: Path, cfg: GiapConfig) -> None:
    acf = [
        _autocorr(_exact(1.3), "morans_i", "gene_label", 1),
        _autocorr(0.4, "morans_i", "gene_label", 2),
        _autocorr(_exact(-0.1), "morans_i", "tu_block", 1),
        _autocorr(_exact(math.inf), "morans_i", "tu_block", 2),
        _autocorr(0.6, "morans_i", "tu_block", 3),
    ]
    inst = _render(_out(tmp_path, autocorrelation=acf), cfg)
    got = {
        (str(t), str(f)): int(_num(v))
        for t, f, v in zip(
            *(_outside_rows(inst)[c] for c in ("table", "family", "value")), strict=True
        )
    }
    assert got == {
        ("autocorrelation", "morans_i / gene_label"): 1,
        ("autocorrelation", "morans_i / tu_block"): 2,
    }
    assert (
        _ax(inst, "b:autocorrelation / morans_i / tu_block")
        .get_title()
        .endswith(", 2 outside [0, 1]")
    )


_ANY_P = st.one_of(
    st.floats(allow_nan=True, allow_infinity=True, width=64),
    st.floats(min_value=0.0, max_value=1.0),
    st.sampled_from([0.0, 1.0, -0.0, 1.0000000000000002, -5e-324]),
)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(ps=st.lists(_ANY_P, min_size=1, max_size=30), seed=st.integers(0, 2**16))
def test_property_b1_bars_plus_outside_is_every_non_nan(ps: list[float], seed: int) -> None:
    """Every non-NaN p is in exactly one bar (iff in [0, 1]) or in the outside count."""
    style.apply()
    cfg = resolve_config()
    rows = [_region(None if math.isnan(p) else _exact(p), i=i) for i, p in enumerate(ps)]
    random.Random(seed).shuffle(rows)
    inside = [p for p in ps if not math.isnan(p) and 0.0 <= p <= 1.0]
    outside = [p for p in ps if not math.isnan(p) and not 0.0 <= p <= 1.0]
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), regions_stats=rows), cfg)
    out = _outside_rows(inst)
    if outside:
        assert len(out) == 1 and _num(out.iloc[0]["value"]) == len(outside)
    else:
        assert len(out) == 0
    if inside:
        expected = {k: 0 for k in range(20)}
        for p in inside:
            expected[min(19, math.floor(Fraction(repr(p)) * 20))] += 1
        assert _family_counts(inst, "regions_stats", "o_minus_e") == expected


# --- lead decision 5: top-16 tie-break -----------------------------------------------------------


def _seventeen_equal(size: int = 2) -> tuple[dict[str, list[dict]], list[str]]:
    """17 families of equal size, listed in spec table order then name."""
    tables: dict[str, list[dict]] = {
        "regions_stats": [],
        "autocorrelation": [],
        "enrichment_detection": [],
        "enrichment_gsea": [],
    }
    labels: list[str] = []
    i = 0
    for stat in ("hl_shift", "o_minus_e"):
        labels.append(f"regions_stats / {stat}")
        for j in range(size):
            tables["regions_stats"].append(_region((j + 1) / (size + 2), stat, i=i))
            i += 1
    labels.append("autocorrelation / morans_i / gene_label")
    for j in range(size):
        tables["autocorrelation"].append(_autocorr(0.3, "morans_i", "gene_label", lag=j + 1))
    ed = sorted(f"{c} / {d}" for c in COLLECTIONS for d in ("over", "under"))
    for fam in ed:
        labels.append(f"enrichment_detection / {fam}")
        c, d = fam.split(" / ")
        for _ in range(size):
            tables["enrichment_detection"].append(_ora(0.4, c, d, i))
            i += 1
    for c in ("cog", "kegg"):
        labels.append(f"enrichment_gsea / {c}")
        for _ in range(size):
            tables["enrichment_gsea"].append(_gsea(0.6, c, i))
            i += 1
    assert len(labels) == 17
    return tables, labels


def test_top16_ties_go_to_table_order_then_name(tmp_path: Path, cfg: GiapConfig) -> None:
    """17 equal families: the last by (table order, name) is dropped: enrichment_gsea / kegg."""
    tables, labels = _seventeen_equal()
    inst = _render(_out(tmp_path, **tables), cfg)
    assert set(_hist_labels(inst)) == set(labels) - {"enrichment_gsea / kegg"}
    b = _rows(inst, "b")
    assert len(set(zip(b["table"], b["family"], strict=True))) == 17


def test_top16_ties_name_order_within_a_table(tmp_path: Path, cfg: GiapConfig) -> None:
    """17 equal families inside one table: the largest name is dropped."""
    names = [f"c{k:02d}" for k in range(17)]
    random.Random(3).shuffle(names)
    gsea = [_gsea(0.5, c, i) for i, c in enumerate(names)]
    inst = _render(_out(tmp_path, enrichment_gsea=gsea), cfg)
    assert set(_hist_labels(inst)) == {f"enrichment_gsea / c{k:02d}" for k in range(16)}


# --- drawn order: spec table order, not rank order ------------------------------------------------


def test_families_are_drawn_in_spec_order_not_rank_order(tmp_path: Path, cfg: GiapConfig) -> None:
    tables = {
        "regions_stats": [_region(0.1, i=0)],
        "autocorrelation": [_autocorr(0.2, "morans_i", "gene_label", k) for k in range(1, 4)],
        "enrichment_gsea": [_gsea(0.5, "kegg", i) for i in range(6)]
        + [_gsea(0.5, "cog", 10 + i) for i in range(2)],
    }
    inst = _render(_out(tmp_path, **tables), cfg)
    want = [
        "regions_stats / o_minus_e",
        "autocorrelation / morans_i / gene_label",
        "enrichment_gsea / cog",
        "enrichment_gsea / kegg",
    ]
    assert _hist_labels(inst) == want
    pos = [_ax(inst, f"b:{lab}").get_position() for lab in want]
    x0 = [round(float(p.x0), 4) for p in pos]
    assert x0 == sorted(x0) and len(set(x0)) == 4  # one grid row, left to right in spec order


# --- titles (decision 2, 7) -----------------------------------------------------------------------


def test_na_slope_prints_na_not_nan(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.1, 0.12, 50), (0.6, 0.55, 50)]), *_overall(None, 0.02, 0.11)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert inst.title == "Detectability calibration slope NA (intercept 0.02), Brier 0.11"
    corner = " ".join(_texts(_ax(inst, "a")))
    assert "slope NA" in corner
    assert "nan" not in corner.lower()


def test_na_intercept_prints_na(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = [*_bins([(0.1, 0.12, 50), (0.6, 0.55, 50)]), *_overall(0.97, None, None)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert inst.title == "Detectability calibration slope 0.97 (intercept NA), Brier NA"


def test_absent_detection_summary_title(tmp_path: Path, cfg: GiapConfig) -> None:
    regions = [_region(0.2, i=0), _region(0.6, "hl_shift", i=1)]
    inst = _render(_out(tmp_path, regions_stats=regions), cfg)
    assert inst.title == "Detectability not available; 2 test families shown"


def test_no_inputs_title(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    assert inst.title == "Detectability not available; 0 test families shown"
    assert ABSENT_TEXT in _texts(_ax(inst, "a"))


# --- panel (a) rows (decision 3, 4) ---------------------------------------------------------------


def test_overall_rows_follow_bins_sorted_by_metric(tmp_path: Path, cfg: GiapConfig) -> None:
    overall = [*_overall(0.97, 0.02, 0.11), _ds("overall", "all", "auc", 0.7)]
    random.Random(11).shuffle(overall)
    bins = _bins([(0.1, 0.12, 50), (0.6, 0.55, 50)])
    inst = _render(_out(tmp_path, detection_summary=[*overall, *bins]), cfg)
    a = _rows(inst, "a")
    is_bin = [_isna(m) for m in a["metric"]]
    assert is_bin == [True, True] + [False] * (len(a) - 2)
    metrics = [str(m) for m in a["metric"] if not _isna(m)]
    want = sorted(str(r["metric"]) for r in overall)
    assert metrics == want
    assert all(_isna(f) for f in a["family"])
    assert set(a["table"]) == {"detection_summary"}


def test_constant_fit_keeps_bin_rows_in_tsv(tmp_path: Path, cfg: GiapConfig) -> None:
    pts = [(0.05 + 0.1 * i, 0.04 + 0.1 * i, 9) for i in range(10)]
    summary = [*_bins(pts), *_overall(None, None, 0.0, constant=True)]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert CONSTANT_TEXT in _texts(_ax(inst, "a"))
    a = _rows(inst, "a")
    bins = a[a["metric"].map(_isna)]
    got = [
        (int(_num(b)), round(_num(x), 6), round(_num(y), 6), int(_num(n)))
        for b, x, y, n in zip(bins["bin"], bins["x"], bins["y"], bins["n"], strict=True)
    ]
    assert got == [(i + 1, round(p[0], 6), round(p[1], 6), 9) for i, p in enumerate(pts)]
    assert _num(a[a["metric"].astype(str) == "detectability_constant"].iloc[0]["value"]) == 1


# --- decision 6 and N3: "not available" vs "constant" --------------------------------------------


def test_basic_section_only_is_not_available(tmp_path: Path, cfg: GiapConfig) -> None:
    inst = _render(_out(tmp_path, detection_summary=[BASIC_ROW]), cfg)
    texts = _texts(_ax(inst, "a"))
    assert ABSENT_TEXT in texts
    assert CONSTANT_TEXT not in texts
    assert inst.title == "Detectability not available; 0 test families shown"
    assert len(_rows(inst, "a")) == 0


def test_n3_no_bin_rows_without_constant_flag_is_not_available(
    tmp_path: Path, cfg: GiapConfig
) -> None:
    """Lead decision 6: the constant text needs detectability_constant == 1."""
    inst = _render(_out(tmp_path, detection_summary=_overall(0.9, 0.1, 0.2)), cfg)
    texts = _texts(_ax(inst, "a"))
    assert ABSENT_TEXT in texts
    assert CONSTANT_TEXT not in texts
    assert inst.title.startswith("Detectability not available;"), inst.title


def test_n3_no_bin_rows_with_constant_flag_is_constant(tmp_path: Path, cfg: GiapConfig) -> None:
    summary = _overall(None, None, 0.2, constant=True)
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    assert CONSTANT_TEXT in _texts(_ax(inst, "a"))
    assert inst.title == "Detectability constant; 0 test families shown"


# --- section filter -------------------------------------------------------------------------------


def test_section_filter_keeps_other_sections_out(tmp_path: Path, cfg: GiapConfig) -> None:
    other_bin = _ds("bin", 1, "mean_p", 0.5, section="basic")
    other_overall = _ds("overall", "all", "calibration_slope", 5.0, section="basic")
    summary = [other_bin, *EX1_SUMMARY, other_overall]
    inst = _render(_out(tmp_path, detection_summary=summary), cfg)
    a = _rows(inst, "a")
    bins = a[a["metric"].map(_isna)]
    first = bins[bins["bin"].map(_num) == 1]
    assert len(first) == 1
    assert _num(first.iloc[0]["x"]) == pytest.approx(0.1)
    slope = a[a["metric"].astype(str) == "calibration_slope"]
    assert len(slope) == 1 and _num(slope.iloc[0]["value"]) == pytest.approx(0.97)
    assert inst.title == "Detectability calibration slope 0.97 (intercept 0.02), Brier 0.11"


# --- exact bin rule -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("p", "k"),
    [
        ("0.44999999999999996", 8),
        ("0.8999999999999999", 17),
        ("0.30000000000000004", 6),
        ("0.45", 9),
        ("0.9", 18),
    ],
)
def test_exact_bin_rule_near_edges(tmp_path: Path, cfg: GiapConfig, p: str, k: int) -> None:
    """Naive ``floor(p * 20)`` puts 0.44999999999999996 in 9 and 0.8999999999999999 in 18."""
    import giap.viz.diagnostics as diag

    assert diag.p_bin(float(p)) == k
    inst = _render(_out(tmp_path, regions_stats=[_region(p)]), cfg)
    counts = _family_counts(inst, "regions_stats", "o_minus_e")
    assert counts == {j: int(j == k) for j in range(20)}, counts
    assert _bar_counts(_ax(inst, "b:regions_stats / o_minus_e")) == {k: 1.0}
