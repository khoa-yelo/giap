"""Nightly task 2's trap rows carry the ``calibration_nightly`` label (issue M11-03b, lead decision
Q5).

A killed or early-failing task 2 used to write ``task 2: …`` (label ``task 2``), so readiness R3
kept reading an older passing ``calibration_nightly`` row. Now its EXIT/TERM trap rows read
``calibration_nightly: …`` with their non-zero exit; the MODE=record skip row stays ``task 2: …``
(exit 0, not a calibration result). ``nightly.sbatch`` runs in the sandbox of
``tests/helpers/slurm_sandbox.py`` exactly as tests/unit/test_m2_21_nightly_task_reporting.py runs
it (fake ``python -m pytest``; ``NIGHTLY_POST=dry``).
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.release_readiness import VERSION, Fixture, load_module, nightly_row, passing_nightly
from helpers.slurm_sandbox import REPO, Sandbox, assert_honours_giap_repo, make_sandbox

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
JOB = "4242"


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_giap_repo(NIGHTLY.read_text(), "nightly.sbatch")
    return make_sandbox(tmp_path / "sb")


def _task2(sb: Sandbox, **env: str) -> list[list[str]]:
    base = {"SLURM_ARRAY_TASK_ID": "2", "SLURM_ARRAY_JOB_ID": JOB, "SLURM_JOB_ID": "5002",
            "SLURM_ARRAY_TASK_COUNT": "4", "SLURM_CPUS_PER_TASK": "2",
            "SLURM_JOB_NAME": "giap-nightly", "NIGHTLY_POST": "dry"}  # fmt: skip
    sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")], env={**base, **env})
    rows = [r for r in sb.tsv_rows() if len(r) >= 5 and r[1] == f"{JOB}_2"]
    assert len(rows) == 1, rows
    return rows


def _label(summary: str) -> str:
    return summary.split(":", 1)[0].strip()


def test_killed_task2_row_is_a_failing_calibration_nightly_row(sb: Sandbox) -> None:
    (row,) = _task2(sb, FAKE_PYTEST="term")
    assert row[3] == "143", row
    assert _label(row[4]) == "calibration_nightly", row
    assert "TERM" in row[4], row


@pytest.mark.parametrize("env", [{"NIGHTLY_REF": "origin/no-such-branch"}, {"MODE": "bogus"}],
                         ids=["checkout", "bad-mode"])  # fmt: skip
def test_task2_that_ends_early_writes_a_failing_calibration_nightly_row(
    sb: Sandbox, env: dict[str, str]
) -> None:
    (row,) = _task2(sb, **env)
    assert row[3] != "0", row
    assert _label(row[4]) == "calibration_nightly", row


def test_record_mode_skip_row_keeps_the_task_label(sb: Sandbox) -> None:
    (row,) = _task2(sb, MODE="record")
    assert row[3] == "0" and _label(row[4]) == "task 2", row


def test_green_task2_row_is_unchanged(sb: Sandbox) -> None:
    (row,) = _task2(sb)
    assert row[3] == "0" and _label(row[4]) == "calibration_nightly", row


def test_a_killed_night_turns_r3_red(sb: Sandbox, tmp_path: Path) -> None:
    """The killed row, newer than a passing row of the same code, decides R3."""
    (row,) = _task2(sb, FAKE_PYTEST="term")
    (tmp_path / "fx").mkdir()
    fx = Fixture.build(tmp_path / "fx")
    short = fx.short(fx.r)
    fx.set_nightly(passing_nightly(short, date="2026-11-02")
                   + nightly_row("2026-11-03", f"{JOB}_2", short, int(row[3]), row[4]))  # fmt: skip
    rd = load_module()
    rd.code_equal.cache_clear()
    r3 = rd.r3(rd.Ctx(repo=fx.repo, version=VERSION, logs=fx.logs))
    assert r3["status"] == "fail", r3
