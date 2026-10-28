"""``scripts/slurm/weekly.sbatch`` run in the sandbox (M2-21 lead decisions 22, 24 and 27).

As for the nightly (``test_m2_21_nightly_task_reporting.py``): fake ``sbatch``/``sacct``/``gh``/
``gippa``/``hostname`` and a fake ``python -m pytest`` first on ``PATH``, a sandbox ``GIPPA_REPO``
with a local bare ``origin``, a temporary ``GIPPA_ROOT`` with no genomes or recordings (so the run
itself fails: no live Bakta JSON). Nothing is submitted.

* Decision 27: every exit writes exactly one ``weekly.tsv`` row for ``<job>_<task>`` and prints
  ``WEEKLY_EXIT=<rc>``: a normal (here failing) run, SIGTERM during the capability gate (the task
  stops: no ``gippa run`` afterwards; exit 143), a checkout failure and a bad ``MODE``.
* Decision 22: ``WEEKLY_PROFILE`` with a ``farmshare`` (Slurm executor) profile without ``_local``
  exits 2 with a message naming R14 or ``farmshare_local``, runs no ``gippa``, writes one row.
* Decision 24: works with ``USER`` unset (scrontab).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.slurm_sandbox import REPO, Sandbox, assert_honours_gippa_repo, make_sandbox

WEEKLY = REPO / "scripts" / "slurm" / "weekly.sbatch"
JOB = "4300"


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_gippa_repo(WEEKLY.read_text(), "weekly.sbatch")
    return make_sandbox(tmp_path)


def _weekly(sb: Sandbox, task: int = 0, unset: tuple[str, ...] = (), **env: str):  # type: ignore[no-untyped-def]
    base = {
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": str(5100 + task),
        "SLURM_ARRAY_TASK_COUNT": "2",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "gippa-weekly",
    }
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "weekly.sbatch")],
                  env={**base, **env}, unset=unset)  # fmt: skip


def _rows(sb: Sandbox, task: int = 0) -> list[list[str]]:
    return [r for r in sb.tsv_rows("weekly.tsv") if len(r) >= 4 and r[1] == f"{JOB}_{task}"]


def _exits(stdout: str) -> list[str]:
    return re.findall(r"WEEKLY_EXIT=(\S+)", stdout)


def _gippa_runs(sb: Sandbox) -> list[str]:
    return [c for c in sb.calls_log() if c.startswith("gippa ")]


def test_normal_run_writes_one_row(sb: Sandbox) -> None:
    res = _weekly(sb)
    rows = _rows(sb)
    assert len(rows) == 1, f"{rows}\n{res.stdout}\n{res.stderr}"
    assert rows[0][3] == str(res.returncode), (rows, res.returncode)
    assert _exits(res.stdout) == [rows[0][3]], res.stdout


def test_sigterm_during_the_gate_writes_one_row_and_stops(sb: Sandbox) -> None:
    res = _weekly(sb, FAKE_PYTEST="term")
    assert not _gippa_runs(sb), f"the task went on after SIGTERM: {_gippa_runs(sb)}"
    rows = _rows(sb)
    assert len(rows) == 1, f"expected one row, got {rows}\n{res.stdout}\n{res.stderr}"
    assert rows[0][3] == "143", rows
    assert _exits(res.stdout) == ["143"], res.stdout
    assert res.returncode != 0


def test_checkout_failure_writes_one_row(sb: Sandbox) -> None:
    res = _weekly(sb, WEEKLY_REF="origin/no-such-branch")
    assert res.returncode != 0
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] != "0", f"{rows}\n{res.stdout}\n{res.stderr}"
    assert not _gippa_runs(sb)


def test_bad_mode_writes_one_row(sb: Sandbox) -> None:
    res = _weekly(sb, MODE="bogus")
    assert res.returncode == 2, res.stderr
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] != "0", f"{rows}\n{res.stdout}\n{res.stderr}"


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "apptainer,farmshare", "farmshare"])
def test_slurm_executor_profile_is_refused(sb: Sandbox, profile: str) -> None:
    res = _weekly(sb, WEEKLY_PROFILE=profile)
    out = res.stdout + res.stderr
    assert res.returncode == 2, out
    assert re.search(r"R14|farmshare_local", out), f"no clear refusal message:\n{out}"
    assert not _gippa_runs(sb), _gippa_runs(sb)
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == "2", rows


def test_runs_without_user(sb: Sandbox) -> None:
    res = _weekly(sb, unset=("USER", "LOGNAME"))
    assert "unbound variable" not in res.stderr, res.stderr
    assert _gippa_runs(sb), f"the run stopped before gippa run:\n{res.stderr}"
    assert len(_rows(sb)) == 1
