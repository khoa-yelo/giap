"""``scripts/slurm/calibration_full_launch.sh`` in the Slurm sandbox (issue M11-03b, acceptance
criterion 6, example 10).

``tests/helpers/slurm_sandbox.py`` provides the sandbox checkout (bare local ``origin``), fake
``sbatch`` (IDs 4242, 4243, …) and the logging ``scancel``. Its ``squeue`` fake answers only
``-n <name>`` queries, so this module replaces it with one that prints ``state/r14.out`` (lines
``<cpus> <name>``, as ``-o '%C %j'``), exits ``FAKE_SQUEUE_RC`` when set and records how many
``sbatch`` calls had happened when it ran. The sandbox commit also carries
``tests/calibration/tolerances.yaml`` (15 ``replicates_full`` keys).

R14 (M11-03c PR 4, superseding the M11-03b array and throttle): a launch queues 16 jobs (15
per-suite jobs and the aggregator), so refuse when ``jobs + 16 > 20``; refuse when
``cpus + Σ + 1 > 128`` (Σ = 98 from ``calibration_plan.tsv``); refuse when the worst-case end
``now + max(time) + 20 min`` (6 h 20 m) reaches the next nightly start (01:30). The clock is
``GIPPA_NOW``; ``_launch`` pins it to 10:00 local unless a test passes its own. The per-suite
details are pinned in ``test_m11_03c_launcher.py``.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pytest
from helpers.slurm_sandbox import (
    REPO,
    Sandbox,
    assert_honours_gippa_repo,
    export_pairs,
    make_sandbox,
    option,
    write_exe,
)

LAUNCH = REPO / "scripts" / "slurm" / "calibration_full_launch.sh"
FIRST, ROWS = "4242", "4257"
SUITE_IDS = [str(4242 + i) for i in range(15)]

R14_SQUEUE = r"""#!/bin/bash
printf '%s\n' "squeue $*" >> "$FAKE_SLURM_DIR/calls.log"
n=0
[ -f "$FAKE_SLURM_DIR/sbatch.jsonl" ] && n=$(wc -l < "$FAKE_SLURM_DIR/sbatch.jsonl")
echo "$n" >> "$FAKE_SLURM_DIR/squeue_seen"
if [ -n "${FAKE_SQUEUE_RC:-}" ]; then echo "squeue: error: fake" >&2; exit "$FAKE_SQUEUE_RC"; fi
[ -f "$FAKE_SLURM_DIR/r14.out" ] && cat "$FAKE_SLURM_DIR/r14.out"
exit 0
"""


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert LAUNCH.is_file(), f"{LAUNCH} does not exist (M11-03b)"
    assert_honours_gippa_repo(LAUNCH.read_text(), "calibration_full_launch.sh")
    box = make_sandbox(tmp_path)
    tol = box.repo / "tests" / "calibration" / "tolerances.yaml"
    tol.parent.mkdir(parents=True)
    shutil.copyfile(REPO / "tests" / "calibration" / "tolerances.yaml", tol)
    for args in (("add", "-A"), ("commit", "-q", "-m", "tolerances"),
                 ("push", "-q", "origin", "main")):  # fmt: skip
        subprocess.run(["git", *args], cwd=box.repo, env=box.env, check=True, capture_output=True)
    write_exe(box.bin / "squeue", R14_SQUEUE)
    return box


def _queue(sb: Sandbox, lines: list[str]) -> None:
    (sb.state / "r14.out").write_text("".join(f"{ln}\n" for ln in lines))


def at(hour: int, minute: int = 0) -> str:
    """Epoch seconds of ``hour:minute`` local time on a fixed day (2026-10-26)."""
    return str(int(time.mktime((2026, 10, 26, hour, minute, 0, 0, 0, -1))))


def _launch(sb: Sandbox, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
    env.setdefault("GIPPA_NOW", at(10))
    script = sb.repo / "scripts" / "slurm" / "calibration_full_launch.sh"
    return sb.run(["bash", str(script), *args], env=env)


def _suites(sb: Sandbox) -> list[dict[str, Any]]:
    calls = [c for c in sb.sbatch_calls() if any(a.endswith("calibration_suite.sbatch")
                                                for a in c["argv"])]  # fmt: skip
    assert len(calls) == 15, sb.sbatch_calls()
    return calls


def _rows_job(sb: Sandbox) -> dict[str, Any]:
    calls = [c for c in sb.sbatch_calls() if "afterany" in (option(c["argv"], "dependency", "d")
                                                            or "")]  # fmt: skip
    assert len(calls) == 1, sb.sbatch_calls()
    return calls[0]


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return res.stdout + res.stderr


def _job_name(call: dict[str, Any]) -> str | None:
    name = option(call["argv"], "job-name", "J")
    if name is None and call.get("script"):
        m = re.findall(r"^#SBATCH\s+(?:--job-name[= ]|-J\s*)(\S+)", call["script"], re.M)
        name = m[-1] if m else None
    return name


def _export(call: dict[str, Any], key: str) -> str | None:
    return export_pairs(call["argv"]).get(key, call["env"].get(key))


# ------------------------------------------------------------------ the happy path


def test_empty_queue_submits_the_suite_jobs_and_aggregator(sb: Sandbox) -> None:
    _queue(sb, [])
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    sha = sb.head
    arr, rows = _suites(sb)[0], _rows_job(sb)
    assert option(arr["argv"], "array", "a") is None, arr["argv"]
    assert "--parsable" in arr["argv"]
    assert (option(arr["argv"], "export") or "").split(",")[0] == "ALL", arr["argv"]
    wt = _export(arr, "GIPPA_WORKTREE")
    assert wt == f"{sb.root}/ci/calibration-{sha[:12]}", wt
    head = subprocess.run(["git", "-C", str(wt), "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True, env=sb.env).stdout.strip()  # fmt: skip
    assert head == sha
    log = _export(arr, "GIPPA_CALIBRATION_LOG_DIR") or ""
    assert re.fullmatch(rf"{re.escape(str(sb.root))}/logs/calibration/full-{sha[:12]}-"
                        r"\d{8}-\d{6}", log), log  # fmt: skip
    assert option(rows["argv"], "cpus-per-task", "c") == "1", rows["argv"]
    want = "afterany:" + ":".join(SUITE_IDS)
    assert option(rows["argv"], "dependency", "d") == want, rows["argv"]
    assert _job_name(rows) == "gippa-calibration-rows"
    assert (_job_name(arr) or "").startswith("gippa-"), arr
    cmd = " ".join(rows["argv"]) + " " + (rows.get("script") or "")
    assert "calibration_rows.py" in cmd, cmd
    for pat in (rf"--job[= ]{FIRST}\b", rf"--commit[= ]{sha}\b", r"--expect[= ]15\b",
                rf"--log-dir[= ]{re.escape(log)}(\s|$|['\"])",
                rf"--tsv[= ]{re.escape(str(sb.root))}/logs/nightly\.tsv"):  # fmt: skip
        assert re.search(pat, cmd), (pat, cmd)
    jobs = ",".join(SUITE_IDS)
    assert re.search(rf"jobs={jobs} rows={ROWS} log={re.escape(log)} cpus=98\b",
                     res.stdout), res.stdout  # fmt: skip


def test_squeue_runs_before_the_first_sbatch(sb: Sandbox) -> None:
    _queue(sb, [])
    assert _launch(sb).returncode == 0
    seen = (sb.state / "squeue_seen").read_text().split()
    assert seen and seen[0] == "0", seen
    q = next(c for c in sb.calls_log() if c.startswith("squeue"))
    assert re.search(r"(^|\s)-r(\s|$)|--array\b", q), q
    assert "RUNNING" in q and "PENDING" in q, q
    assert re.search(r"-u\s*\S+|--user[= ]\S+", q), q


@pytest.mark.parametrize(
    "lines",
    [["8 gippa-a", "4 gippa-b", "64 not-gippa"],  # gippa cpus 12: 12 + 98 + 1 <= 128; jobs 2 + 16
     ["1 gippa-j"] * 4,                         # jobs 4 + 16 = 20 (boundary)
     ["29 gippa-big"]],                         # 29 + 98 + 1 = 128 (boundary)
)  # fmt: skip
def test_launch_fits_the_r14_budget(sb: Sandbox, lines: list[str]) -> None:
    _queue(sb, lines)
    res = _launch(sb)
    assert res.returncode == 0, _out(res)
    assert len(_suites(sb)) == 15
    assert "cpus=98" in res.stdout


# ------------------------------------------------------------------ refusals


@pytest.mark.parametrize(
    "lines",
    [["16 gippa-n"] * 4 + ["3 gippa-x", "3 gippa-y"],   # example 10: 70 CPUs / 6 jobs
     ["112 gippa-big"],                               # 112 + 98 + 1 = 211 > 128
     ["1 gippa-j"] * 5],                              # 5 + 16 = 21 > 20 jobs
)  # fmt: skip
def test_over_budget_refuses_without_sbatch(sb: Sandbox, lines: list[str]) -> None:
    _queue(sb, lines)
    res = _launch(sb)
    assert res.returncode == 2, _out(res)
    assert sb.sbatch_calls() == []
    assert "R14" in _out(res)


def test_example10_message_names_the_numbers(sb: Sandbox) -> None:
    _queue(sb, ["16 gippa-n"] * 4 + ["3 gippa-x", "3 gippa-y"])
    res = _launch(sb)
    out = _out(res)
    assert res.returncode == 2, out
    # 6 + 16 > 20: the job count is checked before the CPU fit (test_m11_03c_launcher pins that)
    assert all(re.search(rf"\b{n}\b", out) for n in (6, 16, 20)), out


def test_squeue_failure_refuses(sb: Sandbox) -> None:
    res = _launch(sb, FAKE_SQUEUE_RC="1")
    assert res.returncode == 2, _out(res)
    assert sb.sbatch_calls() == []


def test_non_integer_cpu_count_refuses(sb: Sandbox) -> None:
    _queue(sb, ["N/A gippa-odd"])
    res = _launch(sb)
    assert res.returncode == 2, _out(res)
    assert sb.sbatch_calls() == []


@pytest.mark.parametrize("args", [("origin/main", "--fast"), ("origin/main", "--dry-run"),
                                  ("origin/main", "--shard", "0/2"), ("--fast",)])  # fmt: skip
def test_extra_arguments_are_refused(sb: Sandbox, args: tuple[str, ...]) -> None:
    _queue(sb, [])
    res = _launch(sb, *args)
    assert res.returncode == 2, _out(res)
    assert sb.sbatch_calls() == []


def test_explicit_ref_is_resolved_to_its_sha(sb: Sandbox) -> None:
    _queue(sb, [])
    sha = sb.head
    res = _launch(sb, sha[:10])
    assert res.returncode == 0, _out(res)
    assert re.search(rf"--commit[= ]{sha}\b", " ".join(_rows_job(sb)["argv"])
                     + " " + (_rows_job(sb).get("script") or ""))  # fmt: skip


def test_failed_aggregator_submission_cancels_exactly_the_suite_jobs(sb: Sandbox) -> None:
    _queue(sb, [])
    res = _launch(sb, FAKE_SBATCH_FAIL="15")  # the 16th sbatch call (the aggregator) fails
    assert res.returncode == 1, _out(res)
    cancels = [c for c in sb.calls_log() if c.startswith("scancel")]
    assert len(cancels) == 1, cancels
    assert re.findall(r"\b\d+\b", cancels[0]) == SUITE_IDS, cancels


def test_failed_array_submission_submits_nothing_else(sb: Sandbox) -> None:
    _queue(sb, [])
    res = _launch(sb, FAKE_SBATCH_FAIL="0")
    assert res.returncode != 0, _out(res)
    assert len(sb.sbatch_calls()) == 1
    assert not [c for c in sb.calls_log() if c.startswith("scancel")]
