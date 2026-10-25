"""V11 title: what is flagged versus what is shown (viz batch 2, 2026-10-24).

``operon_flags.tsv`` rows are drawn in the data TSV and panel (c) whatever their q (M9-11
decision 4). When some rows miss ``stats.alpha.operon_flags`` the title names both counts
("0 of 1 candidate flags at q ≤ 0.1") instead of "0 pairs flagged", which read as if no flag row
existed; with every row at or below alpha (or none) the M9-11 wording is unchanged.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import style
from giap.viz.operon_pairs import operon_pairs, title

MEDS = {"OP": 0.4, "NO": 1.3}
SNAPSHOT = Path(__file__).resolve().parents[2] / "data" / "snapshot" / "sim_plots" / "OUT"


def test_rows_above_alpha_are_counted_as_candidates() -> None:
    got = title(MEDS, 0, 0.1, 1)
    assert got == "Median |Δ| OP 0.4 vs NO 1.3 log2; 0 of 1 candidate flags at q ≤ 0.1"
    assert "pairs flagged" not in got


def test_some_rows_pass_alpha() -> None:
    assert title(MEDS, 2, 0.05, 5).endswith("; 2 of 5 candidate flags at q ≤ 0.05")


def test_every_row_passing_keeps_the_m9_11_wording() -> None:
    assert title(MEDS, 3, 0.1, 3).endswith("; 3 pairs flagged at q ≤ 0.1")
    assert title(MEDS, 0, 0.1, 0).endswith("; 0 pairs flagged at q ≤ 0.1")
    assert title(MEDS, 1, 0.1).endswith("; 1 pairs flagged at q ≤ 0.1")  # n not given


def test_too_few_form_uses_the_same_clause() -> None:
    meds = {"OP": 0.4, "NO": math.nan}
    got = title(meds, 0, 0.1, 2)
    assert got == "0 of 2 candidate flags at q ≤ 0.1; too few measured pairs for densities"


@pytest.mark.skipif(not SNAPSHOT.is_dir(), reason="M9-25 snapshot fixture absent")
def test_snapshot_fixture_terminator_conflict_row_at_q_0_89() -> None:
    style.apply()
    view = RunView.from_outdir(SNAPSHOT)
    flags = view.table("operon_flags")
    assert flags is not None and len(flags) == 1  # the fixture's one row, q 0.89
    result = operon_pairs(view, resolve_config(SNAPSHOT / "resolved_config.yaml"))
    assert isinstance(result, list)
    assert result[0].title.endswith("; 0 of 1 candidate flags at q ≤ 0.1")
