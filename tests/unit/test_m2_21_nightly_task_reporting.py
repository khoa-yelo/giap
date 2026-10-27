"""``scripts/slurm/nightly.sbatch`` run for real in a sandbox (M2-21 lead decisions 14 and 15).

The script runs under bash with fake ``sbatch``/``sacct``/``gh``/``hostname`` and a fake
``python -m pytest`` first on ``PATH`` (``tests/helpers/slurm_sandbox.py``); its ``GIPPA_REPO`` is a
sandbox git repository with a local bare ``origin``, its ``GIPPA_ROOT`` a temporary directory.
``NIGHTLY_POST=dry`` keeps ``post_summary.py`` offline. Nothing is submitted.

* Decision 14: every task reports. A task killed by Slurm (``SIGTERM`` to the job's processes at
  the time limit; simulated by the fake pytest sending ``SIGTERM`` to its process group) and a
  task that exits early (checkout failure) still append their ``nightly.tsv`` row with a non-zero
  exit and print ``NIGHTLY_EXIT=<non-zero>``. A failed benchmark resubmission (``sbatch`` exit 1)
  writes a failure row instead of exiting 0 silently.
* Decision 26: every failure path writes exactly one row for its task.
* Decision 15: the benchmark resubmission passes ``NIGHTLY_REF=<commit>`` (the commit this task
  checked out: a hex prefix of at least 7 characters of the sandbox ``HEAD``, never
  ``origin/main``), either in ``--export`` or in the environment of an ``--export=ALL`` call.
  After a successful resubmission the task writes no row: the pinned job reports the benchmarks.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.slurm_sandbox import (
    REPO,
    Sandbox,
    assert_honours_gippa_repo,
    export_pairs,
    make_sandbox,  # fmt: skip
    option,
)

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
JOB = "4242"


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_gippa_repo(NIGHTLY.read_text(), "nightly.sbatch")
    return make_sandbox(tmp_path)


def _nightly(sb: Sandbox, task: int, **env: str):  # type: ignore[no-untyped-def]
    base = {
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": str(5000 + task),
        "SLURM_ARRAY_TASK_COUNT": "4",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "gippa-nightly",
        "NIGHTLY_POST": "dry",
    }
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
                  env={**base, **env})  # fmt: skip


def _rows(sb: Sandbox, task: int) -> list[list[str]]:
    return [r for r in sb.tsv_rows() if len(r) >= 4 and r[1] == f"{JOB}_{task}"]


def _exits(stdout: str) -> list[str]:
    return re.findall(r"NIGHTLY_EXIT=(\S+)", stdout)


def test_sandbox_harness_a_green_calibration_task_reports_exit_0(sb: Sandbox) -> None:
    res = _nightly(sb, 2)
    rows = _rows(sb, 2)
    assert len(rows) == 1, (
        f"expected one nightly.tsv row for {JOB}_2: {rows}\n{res.stdout}\n{res.stderr}"
    )
    assert rows[0][3] == "0", rows
    assert "calibration_nightly" in rows[0][-1], rows
    assert res.returncode == 0, res.stderr
    assert _exits(res.stdout)[-1:] == ["0"], res.stdout


def test_task_killed_by_sigterm_still_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: the Slurm time limit sends SIGTERM; a trap reports before the task dies."""
    res = _nightly(sb, 2, FAKE_PYTEST="term")
    rows = _rows(sb, 2)
    assert len(rows) == 1, (
        f"a task killed by SIGTERM wrote {len(rows)} nightly.tsv rows, not 1: {rows}\n"
        f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
    )
    assert rows[0][3] != "0", f"killed task reported as passed: {rows}"
    exits = _exits(res.stdout)
    assert exits and exits[-1] != "0", f"no NIGHTLY_EXIT=<non-zero> line:\n{res.stdout}"
    assert res.returncode != 0


def test_task_that_cannot_check_out_its_ref_still_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: "a failure trap writes its NIGHTLY_EXIT row on any exit"."""
    res = _nightly(sb, 0, NIGHTLY_REF="origin/no-such-branch")
    assert res.returncode != 0
    rows = _rows(sb, 0)
    assert len(rows) == 1, (
        f"an early exit (checkout failure) wrote {len(rows)} rows, not 1: {rows}\n"
        f"{res.stdout}\n{res.stderr}"
    )
    assert rows[0][3] != "0", rows
    exits = _exits(res.stdout)
    assert exits and exits[-1] != "0", res.stdout


def _bench_call(sb: Sandbox) -> dict:  # type: ignore[type-arg]
    calls = sb.sbatch_calls()
    assert len(calls) == 1, f"expected one benchmark resubmission, got {calls}"
    return calls[0]


def test_benchmark_resubmission_pins_the_commit(sb: Sandbox) -> None:
    """Decision 15: the resubmitted benchmarks test the commit the other tasks tested."""
    res = _nightly(sb, 1, FAKE_HOSTNAME="barley-01")
    call = _bench_call(sb)
    argv = call["argv"]
    assert "6148" in (option(argv, "constraint", "C") or ""), argv
    assert option(argv, "array", "a") == "1", argv
    exported = export_pairs(argv)
    export = option(argv, "export") or "ALL"
    ref = exported.get("NIGHTLY_REF")
    if ref is None and export.split(",")[0] == "ALL":
        ref = call["env"].get("NIGHTLY_REF")
    head = sb.head
    assert ref is not None, f"NIGHTLY_REF not passed to the resubmission: {argv} {call['env']}"
    assert re.fullmatch(r"[0-9a-f]{7,40}", ref) and head.startswith(ref), (
        f"NIGHTLY_REF={ref!r} is not the checked-out commit {head}"
    )
    assert exported.get("NIGHTLY_ARRAY_JOB", call["env"].get("NIGHTLY_ARRAY_JOB")) == JOB
    # decision 26: without NIGHTLY_EXPECT the pinned job would post with --expect 1, early
    assert exported.get("NIGHTLY_EXPECT", call["env"].get("NIGHTLY_EXPECT")) == "4", argv
    assert res.returncode == 0, res.stderr
    assert not _rows(sb, 1), f"the pinned job reports the benchmarks, not this task: {_rows(sb, 1)}"


def test_failed_benchmark_resubmission_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: sbatch failing on the resubmission path is reported, not silent."""
    res = _nightly(sb, 1, FAKE_HOSTNAME="barley-01", FAKE_SBATCH_FAIL="all")
    assert sb.sbatch_calls(), "the benchmark task did not try to resubmit itself"
    rows = _rows(sb, 1)
    assert len(rows) == 1, (
        f"failed resubmission wrote {len(rows)} rows, not 1: {rows}\n{res.stdout}\n{res.stderr}"
    )
    assert rows[0][3] != "0", rows
    exits = _exits(res.stdout)
    assert exits and exits[-1] != "0", res.stdout


# --- second review of PR #107: decisions 22-24, 26 ----------------------------------------------


def test_bad_mode_writes_exactly_one_row(sb: Sandbox) -> None:
    res = _nightly(sb, 2, MODE="bogus")
    assert res.returncode == 2, res.stderr
    rows = _rows(sb, 2)
    assert len(rows) == 1 and rows[0][3] != "0", rows


def test_sigterm_stops_the_task_and_the_trap_reports_once(sb: Sandbox) -> None:
    """Decision 23: SIGTERM during task 0's pytest ends the task at once (nf-test never starts)
    and the TERM trap writes the one row, exit 143, naming the signal. A trap body of ``true``
    or ``''`` lets the task carry on into nf-test."""
    res = _nightly(sb, 0, FAKE_PYTEST="term")
    nf = [c for c in sb.calls_log() if c.startswith("nf-test")]
    assert not nf, f"the task went on after SIGTERM: {nf}"
    rows = _rows(sb, 0)
    assert len(rows) == 1, f"expected one row, got {rows}\n{res.stdout}\n{res.stderr}"
    assert rows[0][3] == "143", rows
    assert "TERM" in rows[0][-1], f"the row does not say the task was killed by SIGTERM: {rows}"
    assert res.returncode != 0


def test_sigterm_during_report_keeps_the_real_row(sb: Sandbox) -> None:
    """Decision 23: the trap writes a row only when none was written yet. SIGTERM lands while
    report() runs post_summary.py, after this task's (green) line was written."""
    res = _nightly(sb, 2, FAKE_TERM_IN_POST="1")
    assert (sb.state / "term_in_post.done").exists(), "post_summary.py was never reached"
    rows = _rows(sb, 2)
    assert len(rows) == 1, f"SIGTERM in report() added a row: {rows}\n{res.stdout}\n{res.stderr}"
    assert rows[0][3] == "0" and "calibration_nightly" in rows[0][-1], rows


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "apptainer,farmshare", "farmshare"])
def test_slurm_executor_profile_is_refused(sb: Sandbox, profile: str) -> None:
    """Decision 22: a ``farmshare`` (Slurm executor) profile without ``_local`` exits 2 with a
    message naming R14 or ``farmshare_local``; nothing runs and the task writes one row."""
    res = _nightly(sb, 3, NIGHTLY_PROFILE=profile)
    out = res.stdout + res.stderr
    assert res.returncode == 2, out
    assert re.search(r"R14|farmshare_local", out), f"no clear refusal message:\n{out}"
    assert not [c for c in sb.calls_log() if c.startswith("gippa ")], sb.calls_log()
    rows = _rows(sb, 3)
    assert len(rows) == 1 and rows[0][3] == "2", rows


SID = "6100"


def _squeue_calls(sb: Sandbox) -> list[str]:
    return [c for c in sb.calls_log() if c.startswith("squeue")]


def test_resubmission_looks_only_at_pending_summary_jobs(sb: Sandbox) -> None:
    """Decision 26 (``-t PENDING``): a running summary job can no longer wait for anything."""
    sb.set_squeue("gippa-nightly-summary", [f"{SID} nightly:{JOB}"])
    res = _nightly(sb, 1, FAKE_HOSTNAME="barley-01")
    assert res.returncode == 0, res.stderr
    calls = _squeue_calls(sb)
    assert calls, "no squeue query"
    pend = r"(-t\s*|--states[= ])(PENDING|PD)\b"
    assert all(re.search(pend, c) for c in calls), f"squeue without a PENDING filter: {calls}"


def test_resubmission_works_without_user(sb: Sandbox) -> None:
    """Decision 24: under scrontab USER may be unset; ``set -u`` must not abort task 1."""
    import subprocess

    sb.set_squeue("gippa-nightly-summary", [f"{SID} nightly:{JOB}"])
    res = sb.run(
        ["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
        env={"SLURM_ARRAY_TASK_ID": "1", "SLURM_ARRAY_JOB_ID": JOB, "SLURM_JOB_ID": "5001",
             "SLURM_ARRAY_TASK_COUNT": "4", "SLURM_CPUS_PER_TASK": "2",
             "FAKE_HOSTNAME": "barley-01", "NIGHTLY_POST": "dry"},
        unset=("USER", "LOGNAME"),
    )  # fmt: skip
    me = subprocess.run(["id", "-un"], capture_output=True, text=True, check=True).stdout.strip()
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    assert not _rows(sb, 1), f"spurious task-1 row: {_rows(sb, 1)}"
    calls = _squeue_calls(sb)
    pat = rf"-u\s*{re.escape(me)}\b|--user[= ]{re.escape(me)}\b"
    ok = bool(calls) and all(re.search(pat, c) for c in calls)
    assert ok, f"squeue not called with -u {me}: {calls}\n{res.stderr}"
    assert [c for c in sb.calls_log() if c.startswith("scontrol") and "update" in c]


def test_green_task_without_user(sb: Sandbox) -> None:
    res = sb.run(
        ["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
        env={"SLURM_ARRAY_TASK_ID": "2", "SLURM_ARRAY_JOB_ID": JOB, "SLURM_JOB_ID": "5002",
             "SLURM_ARRAY_TASK_COUNT": "4", "SLURM_CPUS_PER_TASK": "2", "NIGHTLY_POST": "dry"},
        unset=("USER", "LOGNAME"),
    )  # fmt: skip
    rows = _rows(sb, 2)
    assert res.returncode == 0 and len(rows) == 1 and rows[0][3] == "0", (rows, res.stderr)
