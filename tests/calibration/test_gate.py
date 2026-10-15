"""Unit tests of ``_gate.calibrate_levels`` (M4-05 decision 7; unmarked: run every commit).

The helper must only route a multi-level criterion through ``_criteria.calibrate``: one simulation
per replicate however many levels read it, the level's column per check, and the ADR-0015 verdicts
of ``calibrate`` unchanged.
"""

from __future__ import annotations

from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import NOT_FAILED, calibrate_levels

LEVELS = (0.01, 0.05, 0.10)
R = 40


def test_one_simulation_per_replicate_across_levels(tmp_path: Path) -> None:
    calls: list[float] = []

    def fn(rng: np.random.Generator) -> float:
        u = float(rng.uniform())
        calls.append(u)
        return 0.5 + 0.5 * u  # never rejects: every level passes in block 1

    gate = calibrate_levels(
        "regions", "level", LEVELS, fn, R, stream="t.gate.once", log_dir=tmp_path
    )
    assert len(calls) == R
    assert [r["status"] for r in gate.records] == ["pass"] * 3
    assert [r["level"] for r in gate.records] == list(LEVELS)
    assert len(gate.values) == R
    gate.assert_ok()


def test_values_equal_unmemoised_calibrate(tmp_path: Path) -> None:
    def fn(rng: np.random.Generator) -> float:
        return float(rng.uniform())

    gate = calibrate_levels("regions", "level", LEVELS, fn, R, stream="t.gate.eq", log_dir=tmp_path)
    for a, rec in zip(LEVELS, gate.records, strict=True):
        ref = C.calibrate("regions", "level", a, fn, R, stream="t.gate.eq", log_dir=tmp_path / "x")
        assert rec["values"] == ref["values"] and rec["status"] == ref["status"]


def test_vector_output_gives_each_level_its_column(tmp_path: Path) -> None:
    def fn(rng: np.random.Generator) -> np.ndarray:
        rng.uniform()
        return np.array([0.0, 0.25, 0.5])

    gate = calibrate_levels("gene_fdr", "fdr", LEVELS, fn, R, stream="t.gate.col", log_dir=tmp_path)
    assert [r["values"][0] for r in gate.records] == [0.0, 0.25, 0.5]
    assert [r["status"] for r in gate.records] == ["pass", "fail", "fail"]
    assert gate.mean() == pytest.approx([0.0, 0.25, 0.5])
    with pytest.raises(AssertionError):
        gate.assert_ok()


def test_confirmed_fail_runs_block_two_once(tmp_path: Path) -> None:
    calls = 0

    def fn(rng: np.random.Generator) -> float:
        nonlocal calls
        calls += 1
        rng.uniform()
        return 0.0  # rejects at every level in both blocks

    gate = calibrate_levels(
        "regions", "level", LEVELS, fn, R, stream="t.gate.fail", log_dir=tmp_path
    )
    assert calls == 2 * R  # block 2 simulated once, shared by the three levels
    assert all(r["status"] == "fail" and r["block"] == 2 for r in gate.records)
    assert "fail" not in NOT_FAILED
    with pytest.raises(AssertionError):
        gate.assert_ok()


def test_wrong_vector_length_refused(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="2 values for 3 levels"):
        calibrate_levels(
            "regions", "level", LEVELS, lambda g: np.array([0.5, 0.5]), R, stream="t.gate.len",
            log_dir=tmp_path,
        )  # fmt: skip


def test_block_one_seconds_kept_when_a_level_goes_to_block_two(tmp_path: Path) -> None:
    """PR #69 decision 2: the simulation time of block 1 survives a level that needs block 2."""
    import time

    def fn(rng: np.random.Generator) -> float:
        time.sleep(0.002)
        rng.uniform()
        return 0.0  # every level fails block 1 and runs block 2

    gate = calibrate_levels(
        "regions", "level", LEVELS, fn, R, stream="t.gate.time", log_dir=tmp_path
    )
    assert all(r["block"] == 2 for r in gate.records)
    assert gate.block_seconds[1] >= R * 0.002
    assert gate.block_seconds[2] >= R * 0.002


def test_block_seconds_single_block(tmp_path: Path) -> None:
    gate = calibrate_levels(
        "regions", "level", LEVELS, lambda g: 0.5 + 0.5 * g.uniform(), R, stream="t.gate.t1",
        log_dir=tmp_path,
    )  # fmt: skip
    assert set(gate.block_seconds) == {1} and gate.block_seconds[1] > 0.0
