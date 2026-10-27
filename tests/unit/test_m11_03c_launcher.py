# ruff: noqa: F811  (pytest fixture sb imported from test_m11_03b_launcher)
"""``calibration_full_launch.sh`` submits one job per suite (issue M11-03c PR 4, criterion 6).

Spec: docs/specs/M11-03c-calibration-full-sharding.md, "calibration_full_launch.sh (changed)",
examples 6-9. The launcher validates the **worktree's** plan (``calibration_lanes.py check``; exit
2, nothing submitted), refuses when ``jobs + 16 > 20`` or ``cpus + sum(lanes x cpus_per_lane) + 1 >
128``, refuses when ``now + max(time) + 20 min`` reaches the next nightly start, submits the 15
``calibration_suite.sbatch`` jobs in plan order (``gippa-cal-<suite>``, ``--cpus-per-task`` = lanes x
cpus_per_lane, ``--mem`` = lanes x mem_per_lane_mb, ``--time``), then the aggregator
(``afterany:<j0>:...:<j14>``, ``--job <j0>``, ``--expect 15``). A failed ``sbatch`` cancels exactly
the jobs already submitted. Fixtures: ``test_m11_03b_launcher`` (fake sbatch IDs 4242, 4243, ...).
"""

from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path
from typing import Any

import pytest
from helpers.slurm_sandbox import REPO, Sandbox, option, write_exe
from unit.test_m11_03b_launcher import (  # noqa: F401  (fixture sb)
    _export,
    _job_name,
    _launch,
    _out,
    _queue,
    at,
    sb,
)

PLAN = REPO / "scripts" / "slurm" / "calibration_plan.tsv"
IDS = [str(4242 + k) for k in range(15)]
ROWS = "4257"


def _plan(path: Path = PLAN) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _lanes(r: dict[str, str]) -> int:
    return int(r["shards"]) + int(r["groups"])


def _suite_jobs(sb: Sandbox) -> list[dict[str, Any]]:
    return [c for c in sb.sbatch_calls() if not c["failed"]
            and any(a.endswith("calibration_suite.sbatch") for a in c["argv"])]  # fmt: skip


def _rows_job(sb: Sandbox) -> dict[str, Any]:
    calls = [c for c in sb.sbatch_calls() if _job_name(c) == "gippa-calibration-rows"]
    assert len(calls) == 1, sb.sbatch_calls()
    return calls[0]


def _cancelled(sb: Sandbox) -> list[str]:
    return [i for c in sb.calls_log() if c.startswith("scancel")
            for i in re.findall(r"\b\d+\b", c)]  # fmt: skip


def _push(sb: Sandbox, edit: dict[str, dict[str, str]], keep_working_copy: bool) -> None:
    """Commit a plan edit to origin/main; optionally restore the sandbox checkout's own copy."""
    path = sb.repo / "scripts" / "slurm" / "calibration_plan.tsv"
    rows = _plan(path)
    for r in rows:
        r.update(edit.get(r["suite"], {}))
    head = "\t".join(rows[0]) + "\n"
    path.write_text(head + "".join("\t".join(r.values()) + "\n" for r in rows))
    for args in (("commit", "-q", "-am", "plan"), ("push", "-q", "origin", "main")):
        subprocess.run(["git", *args], cwd=sb.repo, env=sb.env, check=True, capture_output=True)
    if keep_working_copy:
        subprocess.run(["git", "checkout", "HEAD~1", "--", str(path)], cwd=sb.repo, env=sb.env,
                       check=True, capture_output=True)  # fmt: skip


def test_example6_one_job_per_suite_then_the_aggregator(sb: Sandbox) -> None:
    _queue(sb, ["1 gippa-nightly-launch"])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    plan, jobs = _plan(), _suite_jobs(sb)
    assert [_job_name(j) for j in jobs] == [f"gippa-cal-{r['suite']}" for r in plan], jobs
    assert not [c for c in sb.sbatch_calls() if option(c["argv"], "array", "a")]
    assert not [c for c in sb.sbatch_calls() if "calibration_full.sbatch" in " ".join(c["argv"])]
    wt = f"{sb.root}/ci/calibration-{sb.head[:12]}"
    log = _export(jobs[0], "GIPPA_CALIBRATION_LOG_DIR") or ""
    assert log.startswith(f"{sb.root}/logs/calibration/full-{sb.head[:12]}-"), log
    for r, j in zip(plan, jobs, strict=True):
        a = j["argv"]
        assert "--parsable" in a and f"{wt}/scripts/slurm/calibration_suite.sbatch" in a, a
        assert option(a, "cpus-per-task", "c") == str(_lanes(r) * int(r["cpus_per_lane"])), a
        mem = option(a, "mem") or ""
        assert re.fullmatch(rf"{_lanes(r) * int(r['mem_per_lane_mb'])}M?", mem), a
        assert option(a, "time", "t") == r["time"], a
        assert _export(j, "GIPPA_SUITE") == r["suite"] and _export(j, "GIPPA_WORKTREE") == wt, j
        assert _export(j, "GIPPA_CALIBRATION_LOG_DIR") == log, j
    first = jobs[0]["argv"]
    assert option(first, "cpus-per-task", "c") == "20" and option(first, "time", "t") == "06:00:00"
    rows = _rows_job(sb)
    assert option(rows["argv"], "dependency", "d") == "afterany:" + ":".join(IDS), rows["argv"]
    assert option(rows["argv"], "cpus-per-task", "c") == "1", rows["argv"]
    cmd = " ".join(rows["argv"]) + " " + (rows.get("script") or "")
    for pat in (rf"--job[= ]{IDS[0]}\b", r"--expect[= ]15\b", rf"--commit[= ]{sb.head}\b"):
        assert re.search(pat, cmd), (pat, cmd)
    assert all((_job_name(c) or "").startswith("gippa-") for c in sb.sbatch_calls())
    want = rf"jobs={','.join(IDS)} rows={ROWS} log={re.escape(log)} cpus=98\b"
    assert re.search(want, res.stdout), res.stdout
    seen = (sb.state / "squeue_seen").read_text().split()
    assert seen and seen[0] == "0", seen  # squeue before the first sbatch


@pytest.mark.parametrize(("lines", "ok"),
                         [(["29 gippa-x"], True),        # 29 + 98 + 1 = 128
                          (["30 gippa-x"], False),       # example 7: 129 > 128
                          (["15 gippa-a", "15 gippa-b", "64 other"], False),
                          (["1 gippa-j"] * 4, True),     # 4 + 16 = 20 jobs
                          (["0 gippa-j"] * 5, False)])  # 21 > 20 jobs  # fmt: skip
def test_r14_cpu_fit_and_job_count(sb: Sandbox, lines: list[str], ok: bool) -> None:
    _queue(sb, lines)
    res = _launch(sb)
    if ok:
        assert res.returncode == 0 and len(_suite_jobs(sb)) == 15, _out(res)
        return
    out = _out(res)
    assert res.returncode == 2 and sb.sbatch_calls() == [] and "R14" in out, out
    if sum(int(ln.split()[0]) for ln in lines if ln.split()[1].startswith("gippa-")) == 30:
        assert all(re.search(rf"\b{n}\b", out) for n in (30, 99, 128)), out


@pytest.mark.parametrize(("hour", "minute", "nightly", "ok"),
                         [(19, 9, "01:30", True), (19, 10, "01:30", False),
                          (22, 0, "01:30", False), (2, 0, "01:30", True),
                          (10, 0, "16:20", False), (10, 0, "16:21", True)])  # fmt: skip
def test_example8_window_is_max_time_plus_20_min(
    sb: Sandbox, hour: int, minute: int, nightly: str, ok: bool
) -> None:
    _queue(sb, ["8 gippa-other"])  # CPUs used no longer stretch the window: no waves
    res = _launch(sb, GIPPA_NOW=at(hour, minute), GIPPA_NIGHTLY_AT=nightly)
    if ok:
        assert res.returncode == 0 and len(_suite_jobs(sb)) == 15, _out(res)
    else:
        assert res.returncode == 2 and sb.sbatch_calls() == [] and nightly in _out(res), _out(res)


def test_example9_failed_suite_submission_cancels_exactly_the_submitted(sb: Sandbox) -> None:
    _queue(sb, [])
    res = _launch(sb, FAKE_SBATCH_FAIL="8")  # the 9th sbatch: ora
    assert res.returncode == 1, _out(res)
    assert sorted(_cancelled(sb)) == IDS[:8], sb.calls_log()
    assert len(sb.sbatch_calls()) == 9 and _job_name(sb.sbatch_calls()[8]) == "gippa-cal-ora"


def test_failed_aggregator_cancels_all_15(sb: Sandbox) -> None:
    _queue(sb, [])
    res = _launch(sb, FAKE_SBATCH_FAIL="15")
    assert res.returncode == 1, _out(res)
    assert sorted(_cancelled(sb)) == IDS, sb.calls_log()


def test_bad_worktree_plan_refuses_before_any_sbatch(sb: Sandbox) -> None:
    _push(sb, {"placement": {"cpus_per_lane": "3"}}, keep_working_copy=True)
    _queue(sb, [])
    res = _launch(sb)
    assert res.returncode == 2 and sb.sbatch_calls() == [], _out(res)


def test_jobs_are_sized_from_the_worktree_plan(sb: Sandbox) -> None:
    _push(sb, {"compare": {"cpus_per_lane": "4", "time": "03:00:00"}}, keep_working_copy=True)
    _queue(sb, [])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    last = _suite_jobs(sb)[-1]["argv"]
    assert option(last, "cpus-per-task", "c") == "4" and option(last, "time", "t") == "03:00:00"
    assert re.search(r"\bcpus=100\b", res.stdout), res.stdout


def test_a_plan_without_a_final_newline_keeps_its_last_suite(sb: Sandbox) -> None:
    """#399 review S1: ``check`` accepts it, so the bash read must not drop the last row."""
    path = sb.repo / "scripts" / "slurm" / "calibration_plan.tsv"
    path.write_text(path.read_text().rstrip("\n"))
    for args in (("commit", "-q", "-am", "plan"), ("push", "-q", "origin", "main")):
        subprocess.run(["git", *args], cwd=sb.repo, env=sb.env, check=True, capture_output=True)
    _queue(sb, [])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    assert len(_suite_jobs(sb)) == 15 and re.search(r"\bcpus=98\b", res.stdout), res.stdout


def test_a_failed_scancel_is_reported(sb: Sandbox) -> None:
    """#399 review S2: orphaned suite jobs would hold R14 resources unnoticed."""
    fail = '#!/bin/bash\necho "scancel $*" >> "$FAKE_SLURM_DIR/calls.log"\nexit 1\n'
    write_exe(sb.bin / "scancel", fail)
    _queue(sb, [])
    res = _launch(sb, FAKE_SBATCH_FAIL="15")
    assert res.returncode == 1, _out(res)
    assert "scancel FAILED" in _out(res) and "4242" in _out(res), _out(res)


@pytest.mark.parametrize("args", [("origin/main", "--fast"), ("--shard",)])
def test_harness_arguments_are_still_refused(sb: Sandbox, args: tuple[str, ...]) -> None:
    _queue(sb, [])
    res = _launch(sb, *args)
    assert res.returncode == 2 and sb.sbatch_calls() == [], _out(res)
