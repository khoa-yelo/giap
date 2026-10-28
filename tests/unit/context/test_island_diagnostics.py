"""Per-replicon island diagnostics: l-hat, l' and the lag-1 Moran's I of the detrended TU scores
(ADR-0021 Amendment 1, decision 2a; M5-06 decision 48; docs/methods/islands.md "validity
domain").

Amendment 1 dropped the ``uncalibrated`` flag and put **these three numbers per replicon** in
its place: "no per-replicon statistic certifies the domain, so a flag would give false
reassurance ... Instead: docs/methods/islands.md states the domain above, and the §9.4 run
report carries l-hat, l' (``block_len``) and the lag-1 Moran's I per replicon so a reader can
see where their data sit." M5-09 is the step that has to persist them, because ``IslandScan``
dies with the step.

**The thing that must not happen.** ``IslandScan.l_hat`` is ``0`` both when the replicon was not
scanned and when ``context.islands.block_len`` fixes the block length. Printed as a number,
l-hat = 0 reads "no autocorrelation" — the false reassurance Amendment 1 refuses. So the persisted
table must carry **NA, never 0**, and the renderer must say "not estimated"; the two reasons must
stay distinguishable (fixed block: l' > 0 and Î₁ finite; not scanned: both missing and Î₁ NaN),
and with a fixed block the ADR-0021 calibration does not apply at all.

Interface assumed (**new in M5-09; no lead decision names it yet — see the test-author report**)::

    gippa.context.diagnostics.DIAGNOSTICS_COLUMNS
        ("contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags")
    gippa.context.diagnostics.BLOCK_RULES == ("auto", "fixed", "not_scanned")
    gippa.context.diagnostics.island_diagnostics(scan, cfg) -> pd.DataFrame  # one row per
        replicon, C1 order, nullable integers for l_hat / block_len
    gippa.context.diagnostics.render_l_hat(l_hat, block_rule) -> str
    gippa.context.diagnostics.calibration_applies(block_rule) -> bool

The step persists the frame as ``work/context/island_diagnostics.parquet`` (checked in
``tests/integration/test_context_stats_m5.py``).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as isd
from region_tests_data import Contig

SEED = 5


def _contigs() -> list[Contig]:
    """One scannable chromosome and one replicon too short to scan (decisions 40-41)."""
    rng = np.random.default_rng(SEED)
    return [
        isd.sim_contig(rng, "chr1", 80, rho=0.5),
        isd.sim_contig(rng, "p1", 3, detected_frac=0.6),
    ]


def _scan(cfg: Any) -> Any:
    """``IslandScan`` of the two replicons (M5-06, merged)."""
    from gippa.context.islands import scan_replicons

    genes, tus, replicons = isd.frames(_contigs())
    series = isd.build_series(genes, tus, replicons, cfg)
    return scan_replicons(series, cfg, np.random.default_rng(SEED))


def _diagnostics(cfg: Any) -> pd.DataFrame:
    from gippa.context.diagnostics import island_diagnostics

    return island_diagnostics(_scan(cfg), cfg)


def _row(df: pd.DataFrame, contig: str) -> pd.Series:
    rows = df[df["contig"] == contig]
    assert len(rows) == 1, (contig, list(df["contig"]))
    return rows.iloc[0]


# --------------------------------------------------------------------------- the frame


def test_columns_and_one_row_per_replicon_in_c1_order() -> None:
    from gippa.context.diagnostics import DIAGNOSTICS_COLUMNS

    df = _diagnostics(isd.config(bootstrap=39, permutations=19))
    assert tuple(df.columns) == DIAGNOSTICS_COLUMNS
    assert list(df["contig"]) == ["chr1", "p1"]


def test_auto_block_length_reports_an_estimated_l_hat() -> None:
    """A scanned replicon with ``block_len = "auto"``: l-hat ≥ 1, l' = min(20, 2l-hat + 1, T)
    (ADR-0021 decision 1) and a finite Î₁."""
    cfg = isd.config(bootstrap=39, permutations=19, block_len="auto")
    scan = _scan(cfg)
    from gippa.context.diagnostics import island_diagnostics

    row = _row(island_diagnostics(scan, cfg), "chr1")
    assert row["scanned"] is True or bool(row["scanned"])
    assert row["block_rule"] == "auto"
    assert int(row["l_hat"]) == scan.l_hat[0] >= 1
    assert int(row["block_len"]) == scan.block_len[0] == min(20, 2 * scan.l_hat[0] + 1, 80)
    assert math.isfinite(float(row["morans_i_lag1"]))
    assert not set(row["flags"])


def test_fixed_block_length_reports_l_hat_as_not_estimated() -> None:
    """``context.islands.block_len = 4``: ``IslandScan.l_hat`` is 0 because l-hat was never
    estimated. The table must say NA (never 0), keep l' = 4 and keep Î₁ finite."""
    cfg = isd.config(bootstrap=39, permutations=19, block_len=4)
    scan = _scan(cfg)
    assert scan.l_hat[0] == 0 and scan.block_len[0] == 4  # the trap in the producer
    from gippa.context.diagnostics import island_diagnostics

    row = _row(island_diagnostics(scan, cfg), "chr1")
    assert row["block_rule"] == "fixed"
    assert pd.isna(row["l_hat"])
    assert int(row["block_len"]) == 4
    assert math.isfinite(float(row["morans_i_lag1"]))


def test_unscanned_replicon_reports_everything_missing() -> None:
    """A replicon that is not scanned (decisions 40-41) has no l-hat, no l' and no Î₁, and keeps the
    ``few_values`` run flag — distinguishable from the fixed-block row above."""
    cfg = isd.config(bootstrap=39, permutations=19, block_len=4)
    scan = _scan(cfg)
    assert scan.scanned == (True, False)
    from gippa.context.diagnostics import island_diagnostics

    row = _row(island_diagnostics(scan, cfg), "p1")
    assert not bool(row["scanned"])
    assert row["block_rule"] == "not_scanned"
    assert pd.isna(row["l_hat"]) and pd.isna(row["block_len"])
    assert pd.isna(row["morans_i_lag1"])
    assert set(row["flags"]) == {"few_values"}


@pytest.mark.parametrize("block_len", ["auto", 4])
def test_zero_never_reaches_the_table(block_len: int | str) -> None:
    """l-hat = 0 and l' = 0 are "not estimated" sentinels inside ``IslandScan``, never values."""
    cfg = isd.config(bootstrap=39, permutations=19, block_len=block_len)
    df = _diagnostics(cfg)
    for column in ("l_hat", "block_len"):
        values = [v for v in df[column] if not pd.isna(v)]
        assert all(int(v) >= 1 for v in values), (column, list(df[column]))


def test_the_three_cases_are_distinguishable() -> None:
    """Fixed block and not scanned must not collapse onto the same row."""
    cfg = isd.config(bootstrap=39, permutations=19, block_len=4)
    df = _diagnostics(cfg)
    fixed, unscanned = _row(df, "chr1"), _row(df, "p1")
    assert fixed["block_rule"] != unscanned["block_rule"]
    assert pd.isna(fixed["l_hat"]) and pd.isna(unscanned["l_hat"])
    assert not pd.isna(fixed["block_len"]) and pd.isna(unscanned["block_len"])
    assert math.isfinite(float(fixed["morans_i_lag1"])) and pd.isna(unscanned["morans_i_lag1"])


def test_block_rules_are_the_three_documented_values() -> None:
    from gippa.context.diagnostics import BLOCK_RULES

    assert set(BLOCK_RULES) == {"auto", "fixed", "not_scanned"}


# --------------------------------------------------------------------------- the rendering


@pytest.mark.parametrize("rule", ["fixed", "not_scanned"])
@pytest.mark.parametrize("value", [0, None, float("nan")])
def test_l_hat_without_an_estimate_renders_as_not_estimated(rule: str, value: Any) -> None:
    """ADR-0021 Amendment 1: printed as a number, l-hat = 0 reads "no autocorrelation"."""
    from gippa.context.diagnostics import render_l_hat

    text = render_l_hat(value, rule)
    assert isinstance(text, str)
    assert "not estimated" in text.lower()
    assert "0" not in text


def test_the_two_reasons_render_differently() -> None:
    from gippa.context.diagnostics import render_l_hat

    assert render_l_hat(None, "fixed") != render_l_hat(None, "not_scanned")


def test_an_estimated_l_hat_renders_as_its_number() -> None:
    from gippa.context.diagnostics import render_l_hat

    assert render_l_hat(3, "auto") == "3"
    assert render_l_hat(1, "auto") == "1"


def test_a_fixed_block_voids_the_calibration() -> None:
    """ADR-0021 decision 1 fixes l' = min(20, 2l-hat + 1, T); the FWER evidence of Amendment 1 is
    conditional on that rule, so with an overridden ``block_len`` nothing is calibrated."""
    from gippa.context.diagnostics import calibration_applies

    assert calibration_applies("auto") is True
    assert calibration_applies("fixed") is False
    assert calibration_applies("not_scanned") is False
