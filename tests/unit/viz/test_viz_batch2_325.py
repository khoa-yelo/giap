"""Viz batch 2 (PR #325, review S2): branches of the batch-2 changes no earlier test reached.

Specs: the 2026-10-24 amendments at the end of docs/specs/M9-05 (V06), M9-07 (V18), M9-20 (V15),
M9-18 (V12), M9-08 (V01) and M9-04 (V17).

- V06: with de novo rows but no placements the title adds "; de novo precision in N score bins";
  with no ``mapping_fdr`` rows at all the figure is the whole-figure empty state.
- V18: a run with no V18 input writes one explanatory row (``panel`` NA, ``metric``
  ``empty_state``, ``value`` = the three panel texts); a run with any input keeps zero rows for
  each empty panel.
- V15: a header-only ``determinants_model`` takes the "too few" empty-state branch.
- V12: the panel (b) site line prints the class error with ``base.sig``; title range ends that
  round alike print once.
- V01: the peptide-density ring is ``circos.DENSITY_COLOUR`` = ``style.NEUTRAL_GREY["envelope"]``
  and no #999999 is drawn.
- V17: the matplotlib (e) title is two lines; the plotly twin's subplot title is one line.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from helpers import viz_determinants as vd
from helpers.viz_mapping import EXAMPLE4_TOTALS, denovo_row, full_placement, mapping_out
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, checklist, style
from gippa.viz.base import Instance, Skip

SNAPSHOT = Path(__file__).resolve().parents[2] / "data" / "snapshot" / "sim_plots" / "OUT"
needs_snapshot = pytest.mark.skipif(not SNAPSHOT.is_dir(), reason="M9-25 snapshot fixture absent")
EN_DASH = "\u2013"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _spec(pid: str) -> base.PlotSpec:
    from gippa.steps.plots import load_plot_modules

    load_plot_modules()
    return base.registered()[pid]


def _render(pid: str, out: Path, cfg: GippaConfig) -> Instance:
    result = _spec(pid).fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    inst = result[0]
    assert isinstance(inst, Instance)
    return inst


def _snapshot(pid: str) -> Instance:
    return _render(pid, SNAPSHOT, resolve_config(SNAPSHOT / "resolved_config.yaml"))


def _isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


# --- V06 (M9-05 amendment) --------------------------------------------------------------------

V06_EMPTY = "No peptide placements or de novo precision in this run"
NO_PLACEMENTS = "Placement FDR unestimated: no peptide placements"


@pytest.mark.parametrize(
    ("bins", "tail"), [(1, "; de novo precision in 1 score bin"),
                       (3, "; de novo precision in 3 score bins")]
)  # fmt: skip
def test_v06_denovo_only_title_names_the_score_bins(
    tmp_path: Path, cfg: GippaConfig, bins: int, tail: str
) -> None:
    rows = [denovo_row(i, 0.8 + 0.05 * i, 0.7, 0.99) for i in range(1, bins + 1)]
    inst = _render("V06", mapping_out(tmp_path, fdr=rows), cfg)
    assert inst.title == NO_PLACEMENTS + tail
    assert set(inst.data["panel"]) == {"c"}
    assert len(inst.data) == bins
    assert not checklist.is_empty_state(inst)  # panel (c) is drawn
    assert "c" in [ax.get_label() for ax in inst.figure.axes]


def test_v06_no_rows_at_all_is_the_whole_figure_empty_state(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = _render("V06", mapping_out(tmp_path, fdr=[]), cfg)
    assert checklist.is_empty_state(inst)
    assert inst.title == V06_EMPTY
    assert len(inst.data) == 0
    assert "score bin" not in inst.title


def test_v06_shuffle_fallback_names_the_decoys_on_a_second_line(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """Batch 2: "(shuffle decoys)" goes on its own line (one line ran into the next title)."""
    rows = [*full_placement({}), *full_placement(EXAMPLE4_TOTALS, decoy_type="shuffle")]
    inst = _render("V06", mapping_out(tmp_path, fdr=rows), cfg)
    titles = {ax.get_label(): ax.get_title(loc="left") for ax in inst.figure.axes}
    assert titles["a"] == "(a) targets vs. decoys\n(shuffle decoys)"
    assert titles["b"].endswith("\n(shuffle decoys)")
    assert titles["b"].count("\n") == 1


# --- V18 (M9-07 amendment) --------------------------------------------------------------------

V18_NOTE = "detectability not available; no p-values in this run; null draws not available"


def _region(p: float, i: int = 0) -> dict[str, object]:
    return {"region_id": f"R{i:04d}", "region_type": "prophage", "subset": "all",
            "contig": "chr1", "start": 101 + 1000 * i, "end": 900 + 1000 * i, "n_genes": 5,
            "stat": "o_minus_e", "p_value": p, "method": "x",
            "null_model": "circular_shift"}  # fmt: skip


def _detection_summary() -> list[dict[str, object]]:
    rows = []
    for i, (mp, rate, n) in enumerate([(0.1, 0.12, 50), (0.6, 0.55, 50)], start=1):
        for metric, value in (("mean_p", mp), ("observed_rate", rate), ("n", n)):
            rows.append({"section": "detectability", "group_type": "bin", "group": str(i),
                         "metric": metric, "value": value})  # fmt: skip
    for metric, value in (("brier", 0.11), ("calibration_intercept", 0.02),
                          ("calibration_slope", 0.97)):  # fmt: skip
        rows.append({"section": "detectability", "group_type": "overall", "group": "all",
                     "metric": metric, "value": value})  # fmt: skip
    return rows


def test_v18_no_input_writes_one_explanatory_row(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render("V18", make_out(tmp_path), cfg)
    assert len(inst.data) == 1
    row = inst.data.iloc[0]
    assert _isna(row["panel"])
    assert row["metric"] == "empty_state"
    assert row["value"] == V18_NOTE
    for text in ("detectability not available", "no p-values in this run",
                 "null draws not available"):  # fmt: skip
        assert text in str(row["value"])
    others = [c for c in inst.data.columns if c not in ("metric", "value")]
    assert all(_isna(row[c]) for c in others), row.to_dict()


def test_v18_no_input_row_is_the_only_empty_state_row(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render("V18", make_out(tmp_path), cfg)
    assert list(inst.data["metric"]) == ["empty_state"]


def test_v18_p_values_only_keeps_panel_a_at_zero_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats",
                       [_region(0.01, 0), _region(0.5, 1)])  # fmt: skip
    inst = _render("V18", out, cfg)
    assert "empty_state" not in set(inst.data["metric"].dropna())
    assert int((inst.data["panel"] == "a").sum()) == 0
    assert inst.data["panel"].map(_isna).sum() == 0
    assert int((inst.data["panel"] == "b").sum()) > 0


def test_v18_detectability_only_keeps_panel_b_at_zero_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary",
                       _detection_summary())  # fmt: skip
    inst = _render("V18", out, cfg)
    assert "empty_state" not in set(inst.data["metric"].dropna())
    assert int((inst.data["panel"] == "b").sum()) == 0
    assert inst.data["panel"].map(_isna).sum() == 0
    assert int((inst.data["panel"] == "a").sum()) > 0


# --- V15 (M9-20 amendment) --------------------------------------------------------------------


def test_v15_header_only_model_is_the_too_few_empty_state(tmp_path: Path, cfg: GippaConfig) -> None:
    out = vd.build_out(
        tmp_path,
        model=[],
        attribution=vd.too_few_attribution(60),
        features=vd.feature_rows(60, tu=3, extra_non_model=3),
        support=False,
    )
    inst = _render("V15", out, cfg)
    assert checklist.is_empty_state(inst)
    assert inst.title == "Too few genes for the determinants model (60 rows, 20 TUs)"
    assert len(inst.data) == 0


def test_v15_header_only_model_with_full_inputs_is_still_too_few(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    inst = _render("V15", vd.build_out(tmp_path, model=[]), cfg)
    assert checklist.is_empty_state(inst)
    assert inst.title.startswith("Too few genes for the determinants model (")


# --- V12 (M9-18 amendment) --------------------------------------------------------------------


def _site(rate: object, **kw: object) -> dict[str, object]:
    row: dict[str, object] = {"gene": "g1", "aa_pos": 81, "contig": "chr1", "codon_start": 100,
                              "strand": "+", "ref_aa": "Y", "alt_aa": "D",
                              "alt_type": "substitution", "confidence": "high", "n_alt": 8,
                              "n_ref": 0, "class_error_rate": rate,
                              "class_error_n": 1805}  # fmt: skip
    row.update(kw)
    return row


@pytest.mark.parametrize(
    ("rate", "shown"),
    [(0.872576, "0.87"), (1.0, "1"), (0.933333, "0.93"), (0.0234567, "0.023"), (0.0, "0"),
     ("0.872576", "0.87"), ("unestimated", "unestimated")],
)  # fmt: skip
def test_v12_site_line_rounds_the_class_error(rate: object, shown: str) -> None:
    from gippa.viz.variants import _header

    text = _header(_site(rate))
    if shown == "unestimated":
        assert text.endswith("class error unestimated"), text
    else:
        star = "*" if shown == "1" else ""
        assert text.endswith(f"class error {shown}{star} (pooled, class of 1805 paired reads)"), (
            text
        )


@settings(max_examples=200, deadline=None)
@given(st.floats(min_value=1e-3, max_value=1.0, allow_nan=False))
def test_v12_site_line_class_error_is_base_sig(rate: float) -> None:
    from gippa.viz.variants import _header

    star = "*" if rate == 1 else ""
    want = f"class error {base.sig(rate)}{star} (pooled, class of 1805 paired reads)"
    assert _header(_site(rate)).endswith(want)


@needs_snapshot
def test_v12_snapshot_panel_b_site_lines_are_rounded() -> None:
    inst = _snapshot("V12")
    lines = [
        t.get_text()
        for ax in inst.figure.axes
        for t in ax.texts
        if "class error" in t.get_text() and not t.get_text().startswith("*")
    ]  # site lines, not the rate-1 footnote; fmt: skip
    assert lines
    tails = {
        line.rsplit("class error ", 1)[1].split(" (pooled", 1)[0].rstrip("*") for line in lines
    }
    assert all("(pooled, class of " in line for line in lines if "unestimated" not in line)
    assert "0.87" in tails
    for tail in tails:
        assert tail == "unestimated" or tail == base.sig(float(tail)), tail
    assert not any("0.872576" in line for line in lines)


def _alt(rates: list[object]) -> pd.DataFrame:
    n = len(rates)
    return pd.DataFrame({"call": ["ALT"] * n, "confidence": ["high"] * n,
                         "class_error_rate": rates})  # fmt: skip


@pytest.mark.parametrize(
    ("rates", "rng"),
    [([0.872576, 1.0], f"0.87{EN_DASH}1*"), ([1.0, 1.0], "1*"), ([0.996, 1.0], "1*"),
     ([0.5], "0.5"), ([0.8721, 0.8749, 0.87], "0.87")],
)  # fmt: skip
def test_v12_title_equal_ends_print_once(rates: list[object], rng: str) -> None:
    from gippa.viz.variants import title

    assert title(_alt(rates)).endswith(
        f"(class error {rng}, pooled per class over ALT and mixed sites)"
    )


@settings(max_examples=200, deadline=None)
@given(st.lists(st.floats(min_value=1e-3, max_value=1.0, allow_nan=False), min_size=1,
                max_size=6))  # fmt: skip
def test_v12_title_range_property(rates: list[float]) -> None:
    from gippa.viz.variants import title

    lo, hi = base.sig(min(rates)), base.sig(max(rates))
    want = lo if lo == hi else f"{lo}{EN_DASH}{hi}"
    want += "*" if max(rates) == 1 else ""
    got = title(_alt(list(rates)))
    assert got.endswith(f"(class error {want}, pooled per class over ALT and mixed sites)"), got
    assert (EN_DASH in got) == (lo != hi)


# --- V01 (M9-08 amendment) --------------------------------------------------------------------


def _colours(fig: Any) -> set[str]:
    out: set[str] = set()
    for artist in fig.findobj():
        for name in ("get_facecolor", "get_edgecolor", "get_color"):
            getter = getattr(artist, name, None)
            if getter is None:
                continue
            try:
                arr = np.atleast_2d(mcolors.to_rgba_array(getter()))
            except (TypeError, ValueError):
                continue
            out |= {mcolors.to_hex(r) for r in arr if r[3] > 0}
    return out


def test_v01_density_colour_is_the_envelope_grey() -> None:
    from gippa.viz import circos

    assert style.NEUTRAL_GREY["envelope"] == circos.DENSITY_COLOUR
    assert mcolors.to_hex(circos.DENSITY_COLOUR) == "#cccccc"
    assert mcolors.to_hex(circos.DENSITY_COLOUR) != mcolors.to_hex(style.OKABE_ITO["grey"])


@needs_snapshot
def test_v01_snapshot_density_ring_is_envelope_grey_and_no_okabe_grey() -> None:
    inst = _snapshot("V01")
    ring = [c for ax in inst.figure.axes for c in ax.collections
            if c.get_gid() == "peptide_density"]  # fmt: skip
    assert len(ring) == 1
    faces = {mcolors.to_hex(f) for f in ring[0].get_facecolor()}
    assert faces == {mcolors.to_hex(style.NEUTRAL_GREY["envelope"])}
    assert "#999999" not in _colours(inst.figure)


# --- V17 (M9-04 amendment) --------------------------------------------------------------------

E_TWO_LINES = "(e) missingness\n(non-partial genes)"
E_ONE_LINE = "(e) missingness (non-partial genes)"


@needs_snapshot
def test_v17_matplotlib_panel_e_title_is_two_lines() -> None:
    inst = _snapshot("V17")
    (ax,) = [a for a in inst.figure.axes if a.get_label() == "e"]
    assert ax.get_title(loc="left") == E_TWO_LINES
    assert E_ONE_LINE not in [a.get_title(loc="left") for a in inst.figure.axes]


@needs_snapshot
def test_v17_plotly_twin_panel_e_title_is_one_line() -> None:
    from helpers import twins as tw

    inst = _snapshot("V17")
    spec = tw.spec("V17", inst.data, inst.title)
    notes = [str(a.get("text", "")) for a in spec.get("layout", {}).get("annotations", [])]
    assert E_ONE_LINE in notes, notes
    assert not any("missingness\n" in n or "missingness<br>" in n for n in notes)
    json.dumps(spec)  # the twin stays plain JSON
