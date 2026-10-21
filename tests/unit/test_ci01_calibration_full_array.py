"""``scripts/slurm/calibration_full.sbatch`` and its launcher (CI-01 decisions 1-4).

The scripts are parsed (``bash -n``), read as text and run for real in the sandbox of
``tests/helpers/slurm_sandbox.py`` (fake ``sbatch``/``sacct``/``gh`` and a fake ``python -m
pytest`` first on ``PATH``, sandbox ``GIAP_REPO``/``GIAP_ROOT``). Nothing is submitted, nothing is
posted and the real ``calibration_full`` tier never runs here: the fake pytest answers
``--collect-only`` from a list of test ids and prints a ``--durations`` table for the ids it is
given.

What the tests pin down:

* decision 1: one shard per array task, assigned by test id through ``scripts/ci/shard_tests.py``;
  the shards of one array together run every collected test exactly once;
* decision 2: the shard count is a variable (``CALFULL_SHARDS``) and the run records per-test
  durations; R14: array tasks x ``--cpus-per-task`` + the 1-CPU summary job <= 128;
* decision 3: one ``calibration_full.tsv`` row per shard in the shape of ``nightly.tsv`` (also
  after SIGTERM or an early exit), ``scripts/ci/post_summary.py`` posting once every shard has a
  line, and a launcher that adds an ``afterany`` summary job; ``NIGHTLY_POST=dry`` posts nothing;
* decision 4: no scrontab entry anywhere for this job.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest
from helpers.slurm_sandbox import (
    REPO,
    Sandbox,
    assert_honours_giap_repo,
    make_sandbox,  # fmt: skip
    option,
    write_exe,
)

SLURM = REPO / "scripts" / "slurm"
SBATCH = SLURM / "calibration_full.sbatch"
LAUNCH = SLURM / "calibration_full_launch.sh"
JOB = "4242"
SHARDS = 4
IDS = [
    f"tests/calibration/{d}test_mod{i % 5}.py::test_case{i}[full]"
    for i, d in enumerate(["", "attacks/"] * 9)
]

# A fake `python`: `-m pytest --collect-only` prints $FAKE_IDS, `-m pytest <ids>` prints a
# --durations table for the ids it was given, everything else is the real interpreter.
_PYTEST = r"""#!/bin/bash
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "pytest" ]; then
  printf '%s\n' "$*" >> "$FAKE_SLURM_DIR/pytest.argv"
  case " $* " in
    *" --collect-only "*)
      cat "$FAKE_IDS"
      echo
      echo "$(wc -l < "$FAKE_IDS")/9019 tests collected (8858 deselected) in 1.00s"
      exit "${FAKE_COLLECT_RC:-0}" ;;
  esac
  n=0
  for a in "$@"; do
    case "$a" in
      *::*) n=$((n + 1)); printf '%s\n' "$a" >> "$FAKE_SLURM_DIR/ran.txt"
            printf '%s\n' "1.50s call     $a" ;;
    esac
  done
  echo "=========================== slowest durations ==========================="
  case "${FAKE_PYTEST:-ok}" in
    term) kill -TERM 0; sleep 5; echo "$n failed in 5.00s"; exit 1 ;;
    fail) echo "1 failed, $((n - 1)) passed in 3.00s"; exit 1 ;;
    *) echo "$n passed in 3.00s"; exit 0 ;;
  esac
fi
exec {python} "$@"
"""


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    for path in (SBATCH, LAUNCH):
        assert path.is_file(), f"{path.relative_to(REPO)} does not exist (CI-01)"
        assert_honours_giap_repo(path.read_text(), path.name)
    box = make_sandbox(tmp_path)
    ids = tmp_path / "ids.txt"
    ids.write_text("\n".join(IDS) + "\n")
    for name in ("python", "python3"):
        write_exe(box.bin / name, _PYTEST.replace("{python}", sys.executable))
    box.env["FAKE_IDS"] = str(ids)
    return box


def _shard(sb: Sandbox, task: int, **env: str) -> subprocess.CompletedProcess[str]:
    base = {
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": str(5000 + task),
        "SLURM_ARRAY_TASK_COUNT": str(SHARDS),
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "giap-calfull",
    }
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "calibration_full.sbatch")],
                  env={**base, **env})  # fmt: skip


def _rows(sb: Sandbox, task: int | None = None) -> list[list[str]]:
    rows = [r for r in sb.tsv_rows("calibration_full.tsv") if len(r) >= 4]
    return [r for r in rows if task is None or r[1] == f"{JOB}_{task}"]


def _ran(sb: Sandbox) -> list[str]:
    path = sb.state / "ran.txt"
    return path.read_text().splitlines() if path.is_file() else []


# --- the script itself ---------------------------------------------------------------------------


def _directive(path: Path, key: str) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}[= ](\S+)", path.read_text(), re.M)
    return m.group(1) if m else None


@pytest.mark.parametrize("path", [SBATCH, LAUNCH])
def test_the_scripts_parse(path: Path) -> None:
    res = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr


def test_directives_name_the_job_giap_and_declare_an_array() -> None:
    assert (_directive(SBATCH, "job-name") or "").startswith("giap-"), "R14: giap-* job names"
    array = _directive(SBATCH, "array")
    assert array and re.fullmatch(r"0-\d+", array), f"--array={array}: an array, one shard per task"
    out = _directive(SBATCH, "output") or ""
    assert "%A" in out and "%a" in out, f"--output={out} does not separate the array tasks"


def test_worst_case_cpus_within_r14() -> None:
    """Decision 2 / R14: array tasks x cpus-per-task + the 1-CPU summary job <= 128."""
    array = _directive(SBATCH, "array") or "0-0"
    tasks = int(array.split("-")[1]) + 1
    cpus = int(_directive(SBATCH, "cpus-per-task") or 1)
    assert tasks * cpus + 1 <= 128, f"{tasks} tasks x {cpus} CPUs + 1 > 128 (R14)"


def test_the_shard_count_is_a_documented_variable() -> None:
    text = SBATCH.read_text()
    assert re.search(r"CALFULL_SHARDS:-", text), "no ${CALFULL_SHARDS:-…} default"
    array = int((_directive(SBATCH, "array") or "0-0").split("-")[1]) + 1
    header = "\n".join(ln for ln in text.splitlines() if ln.startswith("#") and "SBATCH" not in ln)
    assert "shard_tests.py" in text, "decision 1: the shards come from scripts/ci/shard_tests.py"
    assert re.search(r"\b" + str(array) + r"\b", header), (
        f"the header does not explain the default of {array} shards (decision 2)"
    )
    assert "durations" in text, "decision 2: the job records per-test durations"


def test_no_scrontab_schedules_the_full_tier() -> None:
    """Decision 4: scheduling needs Khoa; the job is submitted by hand."""
    for path in SLURM.glob("*.scrontab"):
        assert "calibration_full" not in path.read_text(), f"{path.name} schedules the full tier"


def test_the_harness_array_kept_its_own_script() -> None:
    """CI-01 renamed the M4-05 harness array to ``calibration_suite.sbatch`` (the name
    ``calibration_full.sbatch`` now belongs to the tier array): no stale path may be left."""
    assert (SLURM / "calibration_suite.sbatch").is_file(), "the M4-05 harness array is gone"
    suite = (SLURM / "calibration_suite.sbatch").read_text()
    assert "_criteria.py" in suite, "calibration_suite.sbatch is not the harness entry point"
    assert "--task-id" not in SBATCH.read_text(), "the tier array is not the harness"
    for name in ("M5-index.md", "M5-07-island-calibration-power.md"):
        text = (REPO / "docs" / "specs" / name).read_text()
        rows = [ln for ln in text.splitlines() if "calibration_suite.sbatch" in ln]
        assert rows, f"{name} does not point the harness suites at calibration_suite.sbatch"


def test_the_backlog_no_longer_claims_the_tier_runs() -> None:
    """Decision 5: the entry assumed a nightly full array that never existed."""
    text = (REPO / "docs" / "backlog.md").read_text()
    stale = [ln for ln in text.splitlines() if re.search(r"nightly full array", ln)]
    assert not stale, f"docs/backlog.md still says the full tier runs nightly: {stale}"
    assert "calibration_full.sbatch" in text, "the entry does not point at the job that runs it"


# --- the array task, run in the sandbox ---------------------------------------------------------


def test_a_green_shard_writes_one_row_and_its_durations(sb: Sandbox) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=str(SHARDS))
    rows = _rows(sb, 0)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    assert len(rows) == 1, f"expected one calibration_full.tsv row: {rows}\n{res.stdout}"
    date, job_task, commit, code, summary = rows[0][:5]
    assert (date, job_task, code) == (rows[0][0], f"{JOB}_0", "0")
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", date) and re.fullmatch(r"[0-9a-f]{7,40}", commit)
    assert "shard 0" in summary and "passed" in summary, summary
    durations = sb.tsv_rows("calibration_full_durations.tsv")
    assert durations, "decision 2: no per-test durations were recorded"
    ran = _ran(sb)
    assert {r[3] for r in durations} == set(ran), durations[:3]
    assert all(r[1] == f"{JOB}_0" and float(r[4]) > 0 for r in durations), durations[:3]
    assert not sb.gh_calls(), "NIGHTLY_POST defaults to dry: nothing may be posted"


def test_the_shards_of_one_array_run_every_test_exactly_once(sb: Sandbox) -> None:
    """Decision 1: sharding by test id, with no test lost and none run twice."""
    for task in range(SHARDS):
        res = _shard(sb, task, CALFULL_SHARDS=str(SHARDS))
        assert res.returncode == 0, f"shard {task}: {res.stdout}\n{res.stderr}"
    ran = _ran(sb)
    assert sorted(ran) == sorted(IDS), f"the shards did not cover the selection: {sorted(ran)}"
    assert len(_rows(sb)) == SHARDS
    sizes = [len(_ran_of(sb, t)) for t in range(SHARDS)]
    assert all(sizes), f"an empty shard with {len(IDS)} tests and {SHARDS} shards: {sizes}"


def _ran_of(sb: Sandbox, task: int) -> list[str]:
    """The ids a shard was given, read back from the pytest argv log."""
    calls = (sb.state / "pytest.argv").read_text().splitlines()
    runs = [c for c in calls if "--collect-only" not in c]
    return [tok for tok in runs[task].split() if "::" in tok]


def test_a_shard_is_the_same_set_on_every_run(sb: Sandbox) -> None:
    first = [_shard(sb, 1, CALFULL_SHARDS=str(SHARDS)) for _ in range(2)]
    assert all(r.returncode == 0 for r in first)
    ids = _ran(sb)
    half = len(ids) // 2
    assert ids[:half] == ids[half:], "the same shard ran different tests the second time"


def test_a_failing_shard_reports_a_non_zero_exit(sb: Sandbox) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=str(SHARDS), FAKE_PYTEST="fail")
    rows = _rows(sb, 0)
    assert res.returncode != 0
    assert len(rows) == 1 and rows[0][3] != "0", rows
    assert "failed" in rows[0][4], rows


def test_a_shard_killed_by_sigterm_still_writes_one_row(sb: Sandbox) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=str(SHARDS), FAKE_PYTEST="term")
    rows = _rows(sb, 0)
    assert len(rows) == 1, f"expected one row, got {rows}\n{res.stdout}\n{res.stderr}"
    assert rows[0][3] == "143" and "TERM" in rows[0][4], rows
    assert res.returncode != 0


def test_a_shard_that_cannot_check_out_its_ref_still_writes_one_row(sb: Sandbox) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=str(SHARDS), CALFULL_REF="origin/no-such-branch")
    rows = _rows(sb, 0)
    assert res.returncode != 0
    assert len(rows) == 1 and rows[0][3] != "0", f"{rows}\n{res.stdout}\n{res.stderr}"


def test_a_failed_collection_is_reported_not_run(sb: Sandbox) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=str(SHARDS), FAKE_COLLECT_RC="4")
    rows = _rows(sb, 0)
    assert res.returncode == 4, f"{res.stdout}\n{res.stderr}"
    assert len(rows) == 1 and rows[0][3] == "4" and "collection" in rows[0][4], rows
    assert not _ran(sb), "tests ran although collection failed"


@pytest.mark.parametrize("shards", ["0", "two", "-1"])  # empty means "unset": the array's size
def test_a_nonsense_shard_count_is_refused(sb: Sandbox, shards: str) -> None:
    res = _shard(sb, 0, CALFULL_SHARDS=shards)
    rows = _rows(sb, 0)
    assert res.returncode == 2, f"CALFULL_SHARDS={shards!r}: {res.stdout}\n{res.stderr}"
    assert len(rows) == 1 and rows[0][3] == "2", rows


def test_a_task_outside_the_shard_range_is_refused(sb: Sandbox) -> None:
    res = _shard(sb, 5, CALFULL_SHARDS=str(SHARDS))
    assert res.returncode == 2, f"{res.stdout}\n{res.stderr}"
    rows = _rows(sb, 5)
    assert len(rows) == 1 and rows[0][3] == "2", rows
    assert not _ran(sb)


def test_the_last_shard_posts_the_summary_once(sb: Sandbox) -> None:
    """Decision 3: post_summary.py waits for every shard, then posts one comment."""
    for task in range(SHARDS):
        res = _shard(sb, task, CALFULL_SHARDS=str(SHARDS), NIGHTLY_POST="yes")
        assert res.returncode == 0, res.stderr
        posted = [c for c in sb.gh_calls() if c["argv"][:2] == ["issue", "comment"]]
        assert len(posted) == (1 if task == SHARDS - 1 else 0), (
            f"after shard {task} of {SHARDS}: {len(posted)} comments"
        )
    (comment,) = [c for c in sb.gh_calls() if c["argv"][:2] == ["issue", "comment"]]
    body = comment["body"] or ""
    assert JOB in body and "green" in body.lower(), body
    created = [c for c in sb.gh_calls() if c["argv"][:2] == ["issue", "create"]]
    titles = [c["argv"][c["argv"].index("--title") + 1] for c in created]
    assert titles and "Nightly (FarmShare)" not in titles, (
        f"decision 4 leaves the nightly issue to Khoa; posted to {titles}"
    )


# --- the launcher --------------------------------------------------------------------------------


def _launch(sb: Sandbox, **env: str) -> subprocess.CompletedProcess[str]:
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "calibration_full_launch.sh")],
                  env=env)  # fmt: skip


def _calls(sb: Sandbox) -> tuple[dict[str, Any], dict[str, Any]]:
    calls = sb.sbatch_calls()
    arrays = [c for c in calls if any(a.endswith("calibration_full.sbatch") for a in c["argv"])]
    summaries = [c for c in calls if "afterany" in (option(c["argv"], "dependency", "d") or "")]
    assert len(arrays) == 1 and len(summaries) == 1, calls
    return arrays[0], summaries[0]


def test_the_launcher_submits_the_array_and_an_afterany_summary(sb: Sandbox) -> None:
    res = _launch(sb, CALFULL_SHARDS="6")
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    array, summary = _calls(sb)
    assert option(array["argv"], "array", "a") == "0-5", array["argv"]
    dep = option(summary["argv"], "dependency", "d") or ""
    assert re.search(rf"afterany:(?:[0-9_:]*:)?{array['id']}\b", dep), dep
    assert str(array["id"]) in summary["argv"], "the summary job is not told which array to report"
    name = option(summary["argv"], "job-name", "J") or ""
    assert name.startswith("giap-"), f"summary job name {name!r} (R14: giap-*)"
    assert (option(summary["argv"], "cpus-per-task", "c") or "1") == "1"
    assert "post_summary.py" in (summary["script"] or ""), "the summary job must reuse post_summary"
    assert option(summary["argv"], "expect") == "6", "the summary must wait for all six shards"
    tsv = option(summary["argv"], "tsv") or ""
    assert tsv.endswith("calibration_full.tsv"), f"the summary reads {tsv!r}"
    assert "--dry-run" in summary["argv"], "NIGHTLY_POST defaults to dry"


def test_the_launcher_submits_no_summary_for_a_refused_array(sb: Sandbox) -> None:
    res = _launch(sb, FAKE_SBATCH_FAIL="0")
    assert res.returncode != 0, "a refused array submission was hidden"
    later = [c for c in sb.sbatch_calls() if not c["failed"]]
    assert not later, f"summary submitted for an array that does not exist: {later}"


def test_the_launcher_refuses_a_nonsense_shard_count(sb: Sandbox) -> None:
    res = _launch(sb, CALFULL_SHARDS="-3")
    assert res.returncode == 2 and not sb.sbatch_calls(), f"{res.stdout}\n{res.stderr}"


def test_the_summary_job_posts_a_killed_shard_from_sacct(sb: Sandbox) -> None:
    """Decision 3: the appended options make nightly_summary.sh report this array's TSV."""
    commit = sb.head[:7]
    tsv = sb.root / "logs" / "calibration_full.tsv"
    tsv.write_text("".join(
        f"2026-10-20\t{JOB}_{t}\t{commit}\t0\tshard {t}/4: 3 passed\n" for t in (0, 1, 3)
    ))  # fmt: skip
    sb.set_sacct([
        {"JobID": f"{JOB}_{t}", "JobName": "giap-calfull",
         "State": "COMPLETED" if t != 2 else "TIMEOUT"}
        for t in range(SHARDS)
    ])  # fmt: skip
    res = sb.run(
        ["bash", str(sb.repo / "scripts" / "slurm" / "nightly_summary.sh"), JOB,
         "--tsv", str(tsv), "--expect", str(SHARDS), "--title", "Calibration full tier"],
        env={},
    )  # fmt: skip
    comments = [c for c in sb.gh_calls() if c["argv"][:2] == ["issue", "comment"]]
    assert len(comments) == 1, f"{comments}\n{res.stdout}\n{res.stderr}"
    body = comments[0]["body"] or ""
    row = re.search(r"^\|\s*2\s*\|.*$", body, re.M)
    assert row and "TIMEOUT" in row.group(0), f"the killed shard is missing:\n{body}"
