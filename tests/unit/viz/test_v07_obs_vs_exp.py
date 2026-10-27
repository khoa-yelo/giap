"""V7 ``obs_vs_exp``: the O - E forest (M9-06; GIPPA_SPEC.md §10.2 V7, §7.3.5(a)(b), ADR-0018).

Tests the plotted data (``Instance.data`` = ``plots/data/V07_obs_vs_exp.tsv``), the row labels
and their top-to-bottom order, the right-hand labels, marker fill and colour, the title and the
degenerate cases. No ``--mpl`` baselines (M9-25 owns those).

Interface assumed where the spec leaves it open:

- The plots live in ``gippa.viz.detection``; importing it registers ``V07`` (name
  ``obs_vs_exp``, ``needs=("detection_summary",)``, ``"regions_stats"`` in ``optional``,
  ``module="core"``, no ``layers``: ``layer_missing:cog`` drops the COG block, it does not skip).
  The tests call ``registered()["V07"].fn(view, cfg)``: a list of one ``Instance``,
  ``suffix=None``.
- Data columns, exactly: ``block, group, region_type, n, observed, expected, o_minus_e, ci_low,
  ci_high, q_value, flags``. Basic rows: ``region_type`` and ``q_value`` NA, ``n`` = the group's
  ``n_genes`` value, ``observed`` = ``n_detected``, ``flags`` = the ``o_minus_e`` row's flags.
  Region rows: ``group`` = ``region_id``, ``n`` = ``n_eligible``, ``o_minus_e`` = ``effect``.
  Within a basic block, rows keep the table order (the producer sorts by group).
- The forest is the axes whose x label is ``XLABEL`` (with U+2212 MINUS SIGN). Each row
  carries a text "``<group> (n = <n>)``" (tick label or text); rows are read top to bottom by
  screen position. Points are markers (``errorbar``/``plot`` markers or ``scatter``); caps,
  if any, use ``|``; nothing else in the axes is a marker.
- Title numbers: q as ``%.2g``; alpha as ``%g``; "tested" regions = region rows with a q value.
"""

from __future__ import annotations

import math
import os
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.viz_artists import markers, markers_at_x, texts, texts_top_to_bottom
from helpers.viz_out_tree import make_out, write_events, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes

from gippa.config import GippaConfig, resolve_config
from gippa.io.plot_files import write_instance
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

REPO = Path(__file__).resolve().parents[3]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k" / "m4"
VAL_OUT = Path(os.environ.get("GIPPA_ROOT", REPO / "local")) / "runs" / "val-mvp-2026-10-22" / "OUT"
COLUMNS = [
    "block",
    "group",
    "region_type",
    "n",
    "observed",
    "expected",
    "o_minus_e",
    "ci_low",
    "ci_high",
    "q_value",
    "flags",
]
XLABEL = "observed \u2212 expected detected genes"
BLACK = "#000000"
REGION_TYPES = ("prophage", "plasmid", "bgc", "crispr", "amr")


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


# --- builders -------------------------------------------------------------------------------


def basic(
    gtype: str,
    group: str,
    *,
    n_genes: int,
    n_detected: int,
    expected: float,
    ome: float | None = None,
    ci: tuple[float | None, float | None] = (None, None),
    flags: str | None = None,
) -> list[dict[str, object]]:
    """The seven ``basic`` rows of one group, as ``gippa.context.summary.basic_rows`` writes them."""
    value = n_detected - expected if ome is None else ome

    def row(metric: str, v: object, lo: object = None, hi: object = None, f: object = None) -> dict:
        return {
            "section": "basic",
            "group_type": gtype,
            "group": group,
            "metric": metric,
            "value": v,
            "ci_low": lo,
            "ci_high": hi,
            "n": n_genes,
            "flags": f,
        }

    return [
        row("n_genes", float(n_genes)),
        row("n_detected", float(n_detected)),
        row("pct_detected", 100.0 * n_detected / n_genes if n_genes else None),
        row("expected", expected),
        row("o_minus_e", value, ci[0], ci[1], flags),
        row("o_over_e", n_detected / expected if expected else None),
        row("n_tu", float(n_genes)),
    ]


def region(
    region_id: str,
    region_type: str,
    *,
    n_eligible: int = 25,
    observed: float | None = 10.0,
    expected: float | None = 10.0,
    effect: float | None = 0.0,
    ci: tuple[float | None, float | None] = (-1.0, 1.0),
    p: float | None = 0.5,
    q: float | None = 0.5,
    flags: str | None = None,
    label: str | None = None,
    subset: str = "all",
    stat: str = "o_minus_e",
) -> dict[str, object]:
    return {
        "region_id": region_id,
        "region_type": region_type,
        "subset": subset,
        "contig": "chr1",
        "start": 1001,
        "end": 9000,
        "n_genes": n_eligible,
        "n_eligible": n_eligible,
        "observed": observed,
        "expected": expected,
        "o_minus_e": effect,
        "o_over_e": None,
        "hl_shift_log2": None,
        "stat": stat,
        "effect": effect,
        "ci_low": ci[0],
        "ci_high": ci[1],
        "p_value": p,
        "q_value": q,
        "method": "replicon_contrast" if region_type == "plasmid" else "circular_shift_o_minus_e",
        "null_model": "none" if region_type == "plasmid" else "circular_shift",
        "n_perm": 0 if region_type == "plasmid" else 999,
        "label": label,
        "flags": flags,
    }


#: example 3 and two more basic groups (one per block)
BASIC = [
    *basic("cog_category", "J", n_genes=12, n_detected=3, expected=6.4, ci=(-5.0, -1.5)),
    *basic("feature_class", "cds", n_genes=40, n_detected=24, expected=19.0, ci=(-1.0, 9.0)),
    *basic("replicon", "chr1", n_genes=40, n_detected=20, expected=17.9, ome=2.1, ci=(-1.0, 5.2)),
]
#: examples 1 and 2
PROPHAGE_1 = region(
    "prophage_1",
    "prophage",
    n_eligible=25,
    observed=3.0,
    expected=7.75,
    effect=-4.75,
    ci=(-7.1, -2.2),
    p=0.0008,
    q=0.003,
)
PROPHAGE_2 = region("prophage_2", "prophage", effect=1.25, ci=(-0.5, 3.0), p=0.2, q=0.4)
PLASMID = region(
    "plasmid_1",
    "plasmid",
    n_eligible=6,
    effect=0.8,
    ci=(-0.5, 2.1),
    p=None,
    q=None,
    label="dosage_plus_regulation",
    flags="dosage_plus_regulation",
)
REGIONS = [
    PLASMID,
    PROPHAGE_2,
    PROPHAGE_1,
    # rows V7 must not draw
    region("prophage_1", "prophage", stat="hl_shift", effect=0.3, q=0.01),
    region("prophage_1", "prophage", subset="structural", effect=-2.0, q=None, p=None),
]


def _out(
    tmp: Path,
    *,
    summary: Sequence[Mapping[str, object]] | None = BASIC,
    regions: Sequence[Mapping[str, object]] | None = REGIONS,
    flags: Sequence[str] = (),
) -> Path:
    out = make_out(tmp)
    if summary is not None:
        write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", summary)
    if regions is not None:
        write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats", regions)
    if flags:
        records = [{"event": "layer_missing", "data": {"flag": f}} for f in flags]
        write_events(out / "logs" / "events.jsonl", records)
    return out


# --- calling the plot -------------------------------------------------------------------------


def _v07() -> base.PlotSpec:
    import gippa.viz.detection  # noqa: F401  (registers V07 and V08)

    return base.registered()["V07"]


def _render(out: Path, cfg: GippaConfig) -> Instance:
    result = _v07().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _forest(inst: Instance) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_xlabel() == XLABEL]
    assert len(found) == 1, [ax.get_xlabel() for ax in inst.figure.axes]
    return found[0]


def _isna(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and math.isnan(value))


def _num(value: object) -> float:
    return math.nan if _isna(value) else float(str(value))


def _row(inst: Instance, block: str, group: str) -> pd.Series:
    d = inst.data
    hit = d[(d["block"].astype(str) == block) & (d["group"].astype(str) == group)]
    assert len(hit) == 1, f"{block}/{group}:\n{d}"
    return hit.iloc[0]


def _keys(inst: Instance) -> list[tuple[str, str]]:
    return [(str(b), str(g)) for b, g in zip(inst.data["block"], inst.data["group"], strict=True)]


# --- registration and skips -------------------------------------------------------------------


def test_v07_is_registered_as_obs_vs_exp() -> None:
    spec = _v07()
    assert spec.pid == "V07"
    assert spec.name == "obs_vs_exp"
    assert spec.needs == ("detection_summary",)
    assert "regions_stats" in spec.optional
    assert spec.module == "core"
    assert spec.layers == ()


def test_no_detection_summary_is_skipped(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _out(tmp_path, summary=None)
    assert _v07().fn(RunView.from_outdir(out), cfg) == Skip("input_absent:detection_summary")


# --- data TSV: examples 1-3, blocks and order -------------------------------------------------


def test_columns_and_one_instance(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    assert inst.suffix is None
    assert list(inst.data.columns) == COLUMNS


def test_blocks_in_order_and_region_rows_filtered(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    assert _keys(inst) == [
        ("replicon", "chr1"),
        ("feature_class", "cds"),
        ("cog_category", "J"),
        # prophage before plasmid (context.region_types); within prophage by q
        ("region", "prophage_1"),
        ("region", "prophage_2"),
        ("region", "plasmid_1"),
    ]


def test_example3_basic_replicon_row(tmp_path: Path, cfg: GippaConfig) -> None:
    row = _row(_render(_out(tmp_path), cfg), "replicon", "chr1")
    assert _num(row["o_minus_e"]) == pytest.approx(2.1)
    assert _num(row["ci_low"]) == pytest.approx(-1.0)
    assert _num(row["ci_high"]) == pytest.approx(5.2)
    assert _num(row["n"]) == 40
    assert _num(row["observed"]) == 20
    assert _num(row["expected"]) == pytest.approx(17.9)
    assert _isna(row["region_type"])
    assert _isna(row["q_value"])


def test_example1_region_row(tmp_path: Path, cfg: GippaConfig) -> None:
    row = _row(_render(_out(tmp_path), cfg), "region", "prophage_1")
    assert str(row["region_type"]) == "prophage"
    assert _num(row["n"]) == 25
    assert _num(row["observed"]) == 3
    assert _num(row["expected"]) == pytest.approx(7.75)
    assert _num(row["o_minus_e"]) == pytest.approx(-4.75)
    assert (_num(row["ci_low"]), _num(row["ci_high"])) == pytest.approx((-7.1, -2.2))
    assert _num(row["q_value"]) == pytest.approx(0.003)


def test_region_flags_are_carried(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [region("prophage_1", "prophage", flags="small_region;few_shifts", q=0.2)]
    row = _row(_render(_out(tmp_path, regions=regions), cfg), "region", "prophage_1")
    assert "small_region" in str(row["flags"]).split(";")


def test_basic_ci_few_units_flag_is_carried(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = basic("replicon", "chr1", n_genes=8, n_detected=2, expected=3.0,
                    ci=(-2.0, 1.0), flags="ci_few_units")  # fmt: skip
    row = _row(_render(_out(tmp_path, summary=summary), cfg), "replicon", "chr1")
    assert "ci_few_units" in str(row["flags"]).split(";")


def test_region_order_type_then_q_na_last_then_id(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [
        region("bgc_2", "bgc", q=None, p=None),
        region("amr_1", "amr", q=0.01),
        region("bgc_1", "bgc", q=None, p=None),
        region("prophage_3", "prophage", q=0.2),
        region("prophage_1", "prophage", q=0.2),
        region("bgc_3", "bgc", q=0.9),
        region("prophage_2", "prophage", q=0.05),
    ]
    inst = _render(_out(tmp_path, regions=regions), cfg)
    assert [g for b, g in _keys(inst) if b == "region"] == [
        "prophage_2",
        "prophage_1",
        "prophage_3",
        "bgc_3",
        "bgc_1",
        "bgc_2",
        "amr_1",
    ]


def test_region_type_order_follows_context_region_types(tmp_path: Path) -> None:
    cfg = resolve_config(overrides=["context.region_types=[amr, bgc, crispr, plasmid, prophage]"])
    regions = [region("prophage_1", "prophage", q=0.01), region("amr_1", "amr", q=0.5)]
    inst = _render(_out(tmp_path, regions=regions), cfg)
    assert [g for b, g in _keys(inst) if b == "region"] == ["amr_1", "prophage_1"]


def test_layer_missing_cog_drops_the_cog_block(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = [
        *basic("cog_category", "none", n_genes=40, n_detected=20, expected=18.0, ci=(-1.0, 5.0)),
        *basic("replicon", "chr1", n_genes=40, n_detected=20, expected=18.0, ci=(-1.0, 5.0)),
    ]
    out = _out(tmp_path, summary=summary, flags=["layer_missing:cog"])
    inst = _render(out, cfg)
    assert "cog_category" not in set(inst.data["block"].astype(str))
    assert ("replicon", "chr1") in _keys(inst)


def test_regions_stats_absent_draws_basic_blocks_only(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, regions=None), cfg)
    assert [b for b, _ in _keys(inst)] == ["replicon", "feature_class", "cog_category"]


# --- drawing ----------------------------------------------------------------------------------


def test_row_labels_top_to_bottom(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    assert texts_top_to_bottom(inst.figure, "(n = ") == [
        "chr1 (n = 40)",
        "cds (n = 40)",
        "J (n = 12)",
        "prophage_1 (n = 25)",
        "prophage_2 (n = 25)",
        "plasmid_1 (n = 6)",
    ]


def test_vertical_line_at_zero(tmp_path: Path, cfg: GippaConfig) -> None:
    ax = _forest(_render(_out(tmp_path), cfg))
    zero = [
        line
        for line in ax.get_lines()
        if len(set(map(float, line.get_xdata()))) == 1 and float(line.get_xdata()[0]) == 0.0
    ]
    assert zero, "no vertical line at x = 0"


def test_one_point_per_row_at_its_estimate(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    ax = _forest(inst)
    xs = sorted(round(m.x, 9) for m in markers(ax))
    assert xs == sorted(round(_num(v), 9) for v in inst.data["o_minus_e"])


def test_ci_bars_span_the_interval(tmp_path: Path, cfg: GippaConfig) -> None:
    """Every drawn CI end (any line or collection vertex) is at a table CI bound."""
    inst = _render(_out(tmp_path), cfg)
    ax = _forest(inst)
    xs: set[float] = set()
    for line in ax.get_lines():
        xs.update(round(float(x), 9) for x in line.get_xdata())
    for coll in ax.collections:
        for seg in getattr(coll, "get_segments", list)():
            xs.update(round(float(x), 9) for x, _ in seg)
    for lo, hi in zip(inst.data["ci_low"], inst.data["ci_high"], strict=True):
        assert round(_num(lo), 9) in xs and round(_num(hi), 9) in xs, (lo, hi)


def test_example1_q_labels_and_example2_untested_label(tmp_path: Path, cfg: GippaConfig) -> None:
    found = texts(_render(_out(tmp_path), cfg).figure)
    assert "q = 0.003" in found
    assert "q = 0.4" in found
    assert "dosage_plus_regulation (untested)" in found


def test_q_label_uses_two_significant_digits(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [region("prophage_1", "prophage", q=0.012345, p=0.01)]
    assert "q = 0.012" in texts(_render(_out(tmp_path, regions=regions), cfg).figure)


def test_na_ci_is_labelled_unestimated(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = [
        *basic("feature_class", "cds", n_genes=5, n_detected=1, expected=2.0, flags="ci_few_units"),
        *basic("replicon", "chr1", n_genes=5, n_detected=1, expected=2.0, flags="ci_unestimated"),
    ]
    regions = [region("prophage_1", "prophage", q=0.5)]
    found = texts(_render(_out(tmp_path, summary=summary, regions=regions), cfg).figure)
    assert found.count("CI unestimated") == 2


def test_small_region_is_hollow_and_others_filled(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [
        region("prophage_1", "prophage", effect=-3.5, q=0.1, flags="small_region"),
        region("prophage_2", "prophage", effect=2.5, q=0.2),
    ]
    ax = _forest(_render(_out(tmp_path, regions=regions), cfg))
    small, big = markers_at_x(ax, -3.5), markers_at_x(ax, 2.5)
    assert small and all(m.hollow for m in small)
    assert big and not any(m.hollow for m in big)


def test_region_rows_coloured_by_type_basic_rows_black(tmp_path: Path, cfg: GippaConfig) -> None:
    ax = _forest(_render(_out(tmp_path), cfg))
    assert {style.REGION["prophage"]} <= markers_at_x(ax, -4.75)[0].colours()
    assert {style.REGION["plasmid"]} <= markers_at_x(ax, 0.8)[0].colours()
    assert markers_at_x(ax, 2.1)[0].colours() == {BLACK}


# --- titles -----------------------------------------------------------------------------------


def test_example1_title(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path), cfg)
    assert inst.title == "prophage_1 (prophage): 12 % detected vs 31 % expected (q = 0.003)"


def test_title_when_no_region_passes(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [
        region("prophage_1", "prophage", q=0.2),
        region("prophage_2", "prophage", q=0.4),
        PLASMID,
    ]
    inst = _render(_out(tmp_path, regions=regions), cfg)
    assert inst.title == "No region departs from detectability at q ≤ 0.05; 2 regions tested"


def test_title_alpha_from_stats_alpha_context(tmp_path: Path) -> None:
    regions = [
        region("prophage_1", "prophage", n_eligible=20, observed=2.0, expected=5.0, q=0.2),
        region("prophage_2", "prophage", q=0.4),
    ]
    out = _out(tmp_path, regions=regions)
    cfg = resolve_config(overrides=["stats.alpha.context=0.3"])
    assert _render(out, cfg).title == (
        "prophage_1 (prophage): 10 % detected vs 25 % expected (q = 0.2)"
    )
    cfg = resolve_config(overrides=["stats.alpha.context=0.1"])
    assert _render(out, cfg).title == (
        "No region departs from detectability at q ≤ 0.1; 2 regions tested"
    )


def test_title_picks_the_smallest_q_across_types(tmp_path: Path, cfg: GippaConfig) -> None:
    regions = [
        region("prophage_1", "prophage", n_eligible=10, observed=1.0, expected=4.0, q=0.04),
        region("bgc_1", "bgc", n_eligible=10, observed=9.0, expected=5.0, q=0.01),
    ]
    inst = _render(_out(tmp_path, regions=regions), cfg)
    assert inst.title == "bgc_1 (bgc): 90 % detected vs 50 % expected (q = 0.01)"


def test_basic_only_title_names_the_largest_ci_excluding_zero(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    # cds has the largest |O - E| (5.0) but its CI covers 0; J (-3.4) is the one that counts
    inst = _render(_out(tmp_path, regions=None), cfg)
    assert inst.title.startswith("J: ")
    assert "(95 % CI" in inst.title


def test_basic_only_title_when_every_ci_covers_zero(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = [
        *basic("feature_class", "cds", n_genes=40, n_detected=24, expected=19.0, ci=(-1.0, 9.0)),
        *basic("replicon", "chr1", n_genes=40, n_detected=20, expected=17.9, ci=(-1.0, 5.2)),
        # no CI: never named
        *basic("replicon", "plasmid_1", n_genes=4, n_detected=4, expected=0.5),
    ]
    inst = _render(_out(tmp_path, summary=summary, regions=None), cfg)
    assert inst.title == "Detection matches detectability in every group with a CI"


def test_validation_like_run_rows_unestimated_and_no_region_title(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """Every basic row has ``n_detected = 0`` and no CI; regions untested (spec degenerate row)."""
    summary = [
        *basic(
            "feature_class", "cds", n_genes=80, n_detected=0, expected=0.0, flags="ci_few_units"
        ),
        *basic("replicon", "chr1", n_genes=85, n_detected=0, expected=0.0, flags="ci_few_units"),
    ]
    regions = [
        region("plasmid_1", "plasmid", n_eligible=2, observed=None, expected=None, effect=None,
               ci=(None, None), p=None, q=None, label="dosage_plus_regulation",
               flags="small_region;dosage_plus_regulation"),
    ]  # fmt: skip
    inst = _render(_out(tmp_path, summary=summary, regions=regions), cfg)
    assert inst.title.startswith("No region departs from detectability at q ≤ 0.05")
    assert texts(inst.figure).count("CI unestimated") >= 2
    assert _num(_row(inst, "replicon", "chr1")["o_minus_e"]) == 0.0


# --- empty state and files --------------------------------------------------------------------


def test_no_rows_at_all_is_an_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, summary=[], regions=None), cfg)
    assert list(inst.data.columns) == COLUMNS
    assert len(inst.data) == 0
    assert all(not ax.axison for ax in inst.figure.axes)


def test_header_only_summary_with_regions_draws_regions(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path, summary=[]), cfg)
    assert [b for b, _ in _keys(inst)] == ["region"] * 3


def test_data_tsv_written_in_drawing_order(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(_out(tmp_path / "run"), cfg)
    paths = write_instance(inst, "V07", "obs_vs_exp", tmp_path / "plots", ["svg"], 100)
    tsv = tmp_path / "plots" / "data" / "V07_obs_vs_exp.tsv"
    assert tsv in paths
    lines = tsv.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "\t".join(COLUMNS)
    assert [ln.split("\t")[1] for ln in lines[1:]] == [
        "chr1", "cds", "J", "prophage_1", "prophage_2", "plasmid_1",
    ]  # fmt: skip


def test_render_is_deterministic(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _out(tmp_path / "run")
    blobs = []
    for name in ("x", "y"):
        inst = _render(out, cfg)
        paths = write_instance(inst, "V07", "obs_vs_exp", tmp_path / name, ["svg", "png"], 100)
        blobs.append({p.name: p.read_bytes() for p in paths})
    assert blobs[0] == blobs[1]


# --- properties -------------------------------------------------------------------------------

_REGION_ROWS = st.lists(
    st.tuples(
        st.sampled_from(REGION_TYPES),
        st.one_of(st.none(), st.sampled_from([0.001, 0.01, 0.04, 0.2, 0.5, 1.0])),
        st.floats(-10, 10, allow_nan=False).map(lambda v: round(v, 3)),
    ),
    min_size=1,
    max_size=12,
)


@settings(max_examples=15, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(rows=_REGION_ROWS)
def test_property_region_order_and_title(rows: list[tuple[str, float | None, float]]) -> None:
    style.apply()
    regions = [
        region(f"{t}_{i:02d}", t, n_eligible=20, observed=5.0, expected=4.0, effect=e,
               p=q, q=q)
        for i, (t, q, e) in enumerate(rows)
    ]  # fmt: skip
    cfg = resolve_config()
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render(_out(Path(tmp), regions=regions), cfg)

    def key(r: Mapping[str, Any]) -> tuple[int, int, float, str]:
        q = r["q_value"]
        return (REGION_TYPES.index(str(r["region_type"])), q is None, q or 0.0, str(r["region_id"]))

    want = [str(r["region_id"]) for r in sorted(regions, key=key)]
    got = [g for b, g in _keys(inst) if b == "region"]
    assert got == want
    tested = [r for r in regions if r["q_value"] is not None]
    passing = [r for r in tested if float(str(r["q_value"])) <= 0.05]
    if passing:
        # smallest q; ties go to the first in drawing order
        best = min(passing, key=lambda r: (float(str(r["q_value"])), *key(r)[:1], r["region_id"]))
        assert inst.title.startswith(f"{best['region_id']} ({best['region_type']}): 25 % detected")
    else:
        assert inst.title == (
            f"No region departs from detectability at q ≤ 0.05; {len(tested)} regions tested"
        )


# --- golden fixture and the validation run ----------------------------------------------------


def _copy_tables(src: Path, tmp: Path, names: Sequence[str]) -> Path:
    out = make_out(tmp)
    for name in names:
        (out / "tables" / name).write_bytes((src / name).read_bytes())
    return out


def test_renders_on_the_sim50k_golden_tables(tmp_path: Path, cfg: GippaConfig) -> None:
    out = _copy_tables(GOLDEN, tmp_path, ["detection_summary.tsv", "regions_stats.tsv"])
    inst = _render(out, cfg)
    assert _keys(inst) == [
        ("replicon", "chr1"),
        ("feature_class", "cds"),
        ("feature_class", "hypothetical"),
        ("cog_category", "none"),
        ("region", "prophage_1"),
        ("region", "amr_1"),
    ]
    assert inst.title == "No region departs from detectability at q ≤ 0.05; 2 regions tested"
    write_instance(inst, "V07", "obs_vs_exp", tmp_path / "plots", ["svg", "png"], 100)


@pytest.mark.slow
@pytest.mark.skipif(not VAL_OUT.is_dir(), reason=f"validation run absent: {VAL_OUT}")
def test_renders_on_the_validation_run(tmp_path: Path, cfg: GippaConfig) -> None:
    """2026-10-22 validation run: 49 ``basic`` summary rows (all O - E 0), regions all NA. Basic
    rows come first, then the two regions (lead fix: the summary is not header-only)."""
    inst = _render(VAL_OUT, cfg)
    keys = _keys(inst)
    assert keys[-2:] == [("region", "plasmid_1"), ("region", "crispr_1")]
    assert len(keys) > 2 and all(block != "region" for block, _ in keys[:-2])
    assert inst.title == "No region departs from detectability at q ≤ 0.05; 0 regions tested"
    assert "dosage_plus_regulation (untested)" in texts(inst.figure)
    paths = write_instance(inst, "V07", "obs_vs_exp", tmp_path / "plots", ["svg", "png"], 100)
    assert all(p.stat().st_size > 0 for p in paths)


@pytest.mark.parametrize("flags", ["", "ci_unestimated"])
def test_basic_only_title_never_headlines_an_unestimated_ci(
    tmp_path: Path, cfg: GippaConfig, flags: str
) -> None:
    """Review B1 (#202): a zero-width CI (one TU) or a flagged CI reads "CI unestimated" on the
    plot, so the title must not quote it."""
    ci = (-0.5, -0.5) if not flags else (-3.0, -1.0)
    summary = [
        *basic(
            "replicon", "chr2", n_genes=4, n_detected=0, expected=0.5, ci=ci, flags=flags or None
        ),
        *basic("feature_class", "cds", n_genes=40, n_detected=24, expected=19.0, ci=(-1.0, 9.0)),
    ]
    inst = _render(_out(tmp_path, summary=summary, regions=None), cfg)
    assert not inst.title.startswith("chr2")
    assert inst.title == "Detection matches detectability in every group with a CI"
