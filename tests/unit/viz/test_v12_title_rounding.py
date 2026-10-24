"""V12 title: the class-error range is printed with ``giap.viz.base.sig`` (2 significant digits),
not ``%g`` (viz batch 2, 2026-10-24: the snapshot fixture printed "0.872576" as the low end).
A range whose ends round to the same text prints once.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from giap.config import resolve_config
from giap.io.run_view import RunView
from giap.viz import style
from giap.viz.variants import EN_DASH, title, variants

SNAPSHOT = Path(__file__).resolve().parents[2] / "data" / "snapshot" / "sim_plots" / "OUT"


def _frame(rates: list[object], conf: list[str] | None = None) -> pd.DataFrame:
    conf = conf or ["high"] * len(rates)
    return pd.DataFrame(
        {"call": ["ALT"] * len(rates), "confidence": conf, "class_error_rate": rates}
    )


def test_range_is_rounded_to_two_significant_digits() -> None:
    got = title(_frame([0.872576, 1.0, 0.933333], ["high", "candidate", "high"]))
    assert got == f"2 high-confidence and 1 candidate ALT sites (class error 0.87{EN_DASH}1)"


def test_small_rates_keep_their_digits() -> None:
    got = title(_frame([0.004, 0.0234567]))
    assert got.endswith(f"(class error 0.004{EN_DASH}0.023)")


def test_ends_rounding_alike_print_once() -> None:
    assert title(_frame([0.8721, 0.8749])).endswith("(class error 0.87)")


def test_unestimated_rates_are_left_out() -> None:
    assert title(_frame(["unestimated", "unestimated"])).endswith("(class error unestimated)")
    assert title(_frame(["unestimated", 0.123456])).endswith("(class error 0.12)")


@pytest.mark.skipif(not SNAPSHOT.is_dir(), reason="M9-25 snapshot fixture absent")
def test_snapshot_fixture_title_has_no_unrounded_rate() -> None:
    style.apply()
    view = RunView.from_outdir(SNAPSHOT)
    result = variants(view, resolve_config(SNAPSHOT / "resolved_config.yaml"))
    assert isinstance(result, list)
    assert result[0].title.endswith(f"(class error 0.87{EN_DASH}1)")
    assert "0.872576" not in result[0].title
