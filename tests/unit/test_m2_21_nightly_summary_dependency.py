"""The summary job also waits for a resubmitted benchmark task (M2-21 lead decision 19).

``nightly.sbatch`` task 1 and ``nightly_launch.sh`` run in the sandbox of
``tests/helpers/slurm_sandbox.py`` (fake ``sbatch``/``squeue``/``scontrol``/``gh``/``hostname``
first on ``PATH``). Nothing is submitted, updated or posted.

Conventions settled by the test author:

* Task 1 on a non-wheat node resubmits itself (fake job ID ``NEW``, distinct from the array ID).
  The fake ``squeue`` answers a query by name ``giap-nightly-summary`` with ``<id> <comment>``
  lines: one summary job of another night (``nightly:9999``) listed first, then the one of this
  night (``nightly:<array>``). The task then calls
  ``scontrol update JobId=<this night's id> Dependency=afterany:<ids>`` where the ``afterany``
  list contains both the array ID and ``NEW``; it never touches the other night's job, writes no
  ``<array>_1`` row and no ``<array>_note`` row, and posts nothing.
* With no matching summary job, the task writes one ``<array>_note`` row (a warning mentioning
  the resubmitted job ID), no ``<array>_1`` row, calls no ``scontrol update``, and posts nothing
  (posting is enabled here: ``NIGHTLY_POST`` unset, so any post would reach the fake ``gh``).
* The launcher's summary job carries ``--comment=nightly:<array id>`` (``--comment V`` accepted),
  the key the task uses to find it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers.slurm_sandbox import (
    REPO,
    Sandbox,
    assert_honours_giap_repo,
    make_sandbox,  # fmt: skip
    option,
)

SLURM = REPO / "scripts" / "slurm"
NIGHTLY = SLURM / "nightly.sbatch"
LAUNCH = SLURM / "nightly_launch.sh"
JOB = "4242"
NEW = 7001  # the fake sbatch's ID for the resubmission
SID = "6100"  # this night's pending summary job
OTHER = "6099"  # another night's summary job


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_giap_repo(NIGHTLY.read_text(), "nightly.sbatch")
    s = make_sandbox(tmp_path)
    s.set_next_job_id(NEW)
    return s


def _task1(sb: Sandbox, **env: str):  # type: ignore[no-untyped-def]
    base = {
        "SLURM_ARRAY_TASK_ID": "1",
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": "5001",
        "SLURM_ARRAY_TASK_COUNT": "4",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "giap-nightly",
        "FAKE_HOSTNAME": "barley-01",
    }
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
                  env={**base, **env})  # fmt: skip


def _rows(sb: Sandbox, suffix: str) -> list[list[str]]:
    return [r for r in sb.tsv_rows() if len(r) >= 4 and r[1] == f"{JOB}_{suffix}"]


def _updates(sb: Sandbox) -> list[str]:
    return [c for c in sb.calls_log() if c.startswith("scontrol") and " update" in c]


def _summary_queries(sb: Sandbox) -> list[str]:
    return [c for c in sb.calls_log()
            if c.startswith("squeue") and "giap-nightly-summary" in c]  # fmt: skip


def _resubmitted(sb: Sandbox, res) -> None:  # type: ignore[no-untyped-def]
    calls = [c for c in sb.sbatch_calls() if not c["failed"]]
    assert len(calls) == 1 and calls[0]["id"] == NEW, (
        f"expected one benchmark resubmission\n{calls}\n{res.stdout}\n{res.stderr}"
    )


def test_resubmission_extends_the_pending_summary_dependency(sb: Sandbox) -> None:
    sb.set_squeue("giap-nightly-summary", [f"{OTHER} nightly:9999", f"{SID} nightly:{JOB}"])
    res = _task1(sb, NIGHTLY_POST="dry")
    _resubmitted(sb, res)
    assert _summary_queries(sb), f"no squeue query for giap-nightly-summary: {sb.calls_log()}"
    ups = _updates(sb)
    assert len(ups) == 1, f"expected one scontrol update, got {ups}\n{res.stdout}\n{res.stderr}"
    up = ups[0]
    m_id = re.search(r"\bJob[Ii][Dd]=(\S+)", up)
    assert m_id and m_id.group(1) == SID, f"updated {up!r}, not this night's summary job {SID}"
    m_dep = re.search(r"\bDependency=(\S+)", up)
    assert m_dep, f"no Dependency= in {up!r}"
    dep = m_dep.group(1)
    assert dep.startswith("afterany:"), dep
    ids = dep.split(",")[0].split(":")[1:]
    assert JOB in ids and str(NEW) in ids, (
        f"Dependency={dep}: must be afterany on both the array {JOB} and the resubmission {NEW}"
    )
    assert OTHER not in up
    assert res.returncode == 0, res.stderr
    assert not _rows(sb, "1"), f"the pinned job reports the benchmarks: {_rows(sb, '1')}"
    assert not _rows(sb, "note"), f"summary job found, yet a warning row: {_rows(sb, 'note')}"
    assert not sb.gh_calls()


def test_no_summary_job_writes_a_warning_row_and_posts_nothing(sb: Sandbox) -> None:
    sb.set_squeue("giap-nightly-summary", [f"{OTHER} nightly:9999"])
    res = _task1(sb)  # posting enabled: any post would reach the fake gh
    _resubmitted(sb, res)
    assert _summary_queries(sb), f"no squeue query for giap-nightly-summary: {sb.calls_log()}"
    assert not _updates(sb), f"no summary job of {JOB}, yet scontrol update: {_updates(sb)}"
    notes = _rows(sb, "note")
    assert len(notes) == 1, f"expected one {JOB}_note row: {sb.tsv_rows()}\n{res.stdout}"
    assert str(NEW) in notes[0][-1], f"the warning does not name the resubmission: {notes}"
    assert not _rows(sb, "1"), f"no task-1 row expected: {_rows(sb, '1')}"
    assert not sb.gh_calls(), f"task 1 posted by itself: {sb.gh_calls()}"
    assert "NIGHTLY_EXIT" not in res.stdout, res.stdout
    assert res.returncode == 0, res.stderr


def test_launcher_summary_job_carries_the_nightly_comment(tmp_path: Path) -> None:
    assert_honours_giap_repo(LAUNCH.read_text(), "nightly_launch.sh")
    sb = make_sandbox(tmp_path)
    res = sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly_launch.sh")])
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    calls = sb.sbatch_calls()
    arrays = [c for c in calls if option(c["argv"], "array", "a") == "0-3"]
    assert len(arrays) == 1, calls
    array_id = str(arrays[0]["id"])
    summaries = [c for c in calls if "afterany" in (option(c["argv"], "dependency", "d") or "")]
    assert len(summaries) == 1, calls
    assert option(summaries[0]["argv"], "comment") == f"nightly:{array_id}", summaries[0]["argv"]
