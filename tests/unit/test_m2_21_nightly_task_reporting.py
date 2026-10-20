"""``scripts/slurm/nightly.sbatch`` run for real in a sandbox (M2-21 lead decisions 14 and 15).

The script runs under bash with fake ``sbatch``/``sacct``/``gh``/``hostname`` and a fake
``python -m pytest`` first on ``PATH`` (``tests/helpers/slurm_sandbox.py``); its ``GIAP_REPO`` is a
sandbox git repository with a local bare ``origin``, its ``GIAP_ROOT`` a temporary directory.
``NIGHTLY_POST=dry`` keeps ``post_summary.py`` offline. Nothing is submitted.

* Decision 14: every task reports. A task killed by Slurm (``SIGTERM`` to the job's processes at
  the time limit; simulated by the fake pytest sending ``SIGTERM`` to its process group) and a
  task that exits early (checkout failure) still append their ``nightly.tsv`` row with a non-zero
  exit and print ``NIGHTLY_EXIT=<non-zero>``. A failed benchmark resubmission (``sbatch`` exit 1)
  writes a failure row instead of exiting 0 silently.
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
    assert_honours_giap_repo,
    export_pairs,
    make_sandbox,  # fmt: skip
    option,
)

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
JOB = "4242"


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_giap_repo(NIGHTLY.read_text(), "nightly.sbatch")
    return make_sandbox(tmp_path)


def _nightly(sb: Sandbox, task: int, **env: str):  # type: ignore[no-untyped-def]
    base = {
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": str(5000 + task),
        "SLURM_ARRAY_TASK_COUNT": "4",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "giap-nightly",
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
    assert rows, f"no nightly.tsv row for {JOB}_2\n{res.stdout}\n{res.stderr}"
    assert rows[-1][3] == "0", rows
    assert "calibration_nightly" in rows[-1][-1], rows
    assert res.returncode == 0, res.stderr
    assert _exits(res.stdout)[-1:] == ["0"], res.stdout


def test_task_killed_by_sigterm_still_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: the Slurm time limit sends SIGTERM; a trap reports before the task dies."""
    res = _nightly(sb, 2, FAKE_PYTEST="term")
    rows = _rows(sb, 2)
    assert rows, (
        f"a task killed by SIGTERM wrote no nightly.tsv row (failure trap missing)\n"
        f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
    )
    assert rows[-1][3] != "0", f"killed task reported as passed: {rows}"
    exits = _exits(res.stdout)
    assert exits and exits[-1] != "0", f"no NIGHTLY_EXIT=<non-zero> line:\n{res.stdout}"
    assert res.returncode != 0


def test_task_that_cannot_check_out_its_ref_still_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: "a failure trap writes its NIGHTLY_EXIT row on any exit"."""
    res = _nightly(sb, 0, NIGHTLY_REF="origin/no-such-branch")
    assert res.returncode != 0
    rows = _rows(sb, 0)
    assert rows, (
        f"an early exit (checkout failure) wrote no nightly.tsv row\n{res.stdout}\n{res.stderr}"
    )
    assert rows[-1][3] != "0", rows
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
    assert res.returncode == 0, res.stderr
    assert not _rows(sb, 1), f"the pinned job reports the benchmarks, not this task: {_rows(sb, 1)}"


def test_failed_benchmark_resubmission_writes_a_failure_row(sb: Sandbox) -> None:
    """Decision 14: sbatch failing on the resubmission path is reported, not silent."""
    res = _nightly(sb, 1, FAKE_HOSTNAME="barley-01", FAKE_SBATCH_FAIL="all")
    assert sb.sbatch_calls(), "the benchmark task did not try to resubmit itself"
    rows = _rows(sb, 1)
    assert rows, f"failed resubmission wrote no nightly.tsv row\n{res.stdout}\n{res.stderr}"
    assert rows[-1][3] != "0", rows
    exits = _exits(res.stdout)
    assert exits and exits[-1] != "0", res.stdout
