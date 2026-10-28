"""The K-12 values recorded in ``calibration.yaml`` (ADR-0022), pinned so that a changed file
fails a fast test, not only the M4 golden (review of #292). Other tests read these values from
the same file the steps read, so they cannot catch a wrong value."""

from __future__ import annotations

from gippa.config import GippaConfig
from gippa.data.calibration import load as load_calibration


def test_recorded_k12_values() -> None:
    cal = load_calibration()
    assert cal["operon_distance"]["d0"] == 102.68667744326773
    assert cal["operon_distance"]["s"] == 47.16037265004493
    assert cal["uniop_pair_threshold"]["value"] == 0.5  # Khoa, 2026-10-24
    assert cal["terminatornet_yes"]["value"] == 0.9999999999961204
    assert all(
        cal[k]["calibrated"]
        for k in ("operon_distance", "uniop_pair_threshold", "terminatornet_yes")
    )


def test_defaults_mirror_the_recorded_thresholds() -> None:
    ann = GippaConfig().annotation
    assert ann.uniop.pair_threshold == 0.5
    assert ann.terminators.terminatornet_yes == 0.9999999999961204
