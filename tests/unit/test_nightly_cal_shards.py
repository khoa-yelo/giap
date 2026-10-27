"""Nightly task 2 in shards (``NIGHTLY_CAL_SHARDS``): one ``calibration_nightly`` row per night.

The unsharded task 2 (the whole fast suite with coverage, ADR-0052) took 5 h 43 on 2026-10-24
and was killed at the 6 h limit on 2026-10-27 (job 1900002_2, no failure names). Pinned here, in
the M2-21 sandbox (``tests/helpers/slurm_sandbox.py``) with fake ``pytest``/``coverage``:

1. Shard 1 (the array task) submits shards 2..N (``--array=2``, the same commit, its own worktree
   ``<W>-s<k>``), extends the summary job and writes no row; each shard runs its
   ``scripts/ci/shard_tests.py`` files with coverage data in ``<job>_2.shards/coverage.s<k>``; the
   last shard writes exactly one row: counts summed, the largest exit, the combined coverage gate.
   The row stays readable by ``scripts/release/readiness.py`` (label, exit, ``N failed``).
2. A failing shard, a refused shard submission or coverage below ``coverage_min`` fails the row.
3. A shard killed by SIGTERM records exit 143 with the FAILED lines pytest -v printed so far.
4. ``nightly_launch.sh`` asks for 2 shards; R14 worst case stays within 128 CPUs.
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from helpers.gnu_userland import require_gnu_userland
from helpers.slurm_sandbox import (
    _PYTHON,
    REPO,
    Sandbox,
    assert_honours_giap_repo,
    export_pairs,
    make_sandbox,
    option,
    write_exe,
)

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
JOB = "4242"
FILES = [f"tests/test_{c}.py" for c in "abcde"]

# Before the sandbox wrapper's own pytest fake: `python -m coverage` logs and prints a TOTAL line
# (FAKE_COV=low: under the gate, exit 2); `python -m pytest` logs its argv, writes $COVERAGE_FILE
# and prints FAKE_SHARD_PYTEST's summary (ok | fail | term: a FAILED line, then SIGTERM).
_HOOK = r"""
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "coverage" ]; then
  [ "${FAKE_COV:-}" = real ] && exec {python} "$@"
  printf 'coverage %s\n' "${*:3}" >> "$FAKE_SLURM_DIR/calls.log"
  if [ "${3:-}" = report ]; then
    echo "TOTAL 1000 50 95%"; [ "${FAKE_COV:-ok}" = low ] && exit 2
  fi
  exit 0
fi
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "pytest" ]; then
  printf 'pytest %s\n' "${*:3}" >> "$FAKE_SLURM_DIR/calls.log"
  if [ "${FAKE_COV:-}" = real ]; then   # real data: src/giap/mod.py run under coverage, branch on
    {python} -c 'import coverage, runpy; c = coverage.Coverage(branch=True); c.start()
runpy.run_path("src/giap/mod.py"); c.stop(); c.save()'
  elif [ -n "${COVERAGE_FILE:-}" ]; then echo data > "$COVERAGE_FILE"; fi
  case "${FAKE_SHARD_PYTEST:-ok}" in
    fail) echo "2 failed, 10 passed, 1 warning in 7.25s"; exit 1 ;;
    term) echo "[gw0] [ 50%] FAILED tests/test_a.py::test_x"; sleep 1
          kill -TERM 0; sleep 5; exit 1 ;;
    *) echo "10 passed, 1 skipped, 2 warnings in 2.50s"; exit 0 ;;
  esac
fi
"""


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    require_gnu_userland()
    assert_honours_giap_repo(NIGHTLY.read_text(), "nightly.sbatch")
    box = make_sandbox(tmp_path)
    head, _, rest = _PYTHON.replace("{python}", sys.executable).partition("\n")
    for name in ("python", "python3"):
        write_exe(box.bin / name, f"{head}\n{_HOOK.replace('{python}', sys.executable)}{rest}")
    env = {**os.environ, **box.env}
    srcs = {f: "def test_x():\n    pass\n" for f in FILES}  # shard_tests.py needs test files
    srcs["src/giap/mod.py"] = "def f(x):\n    return x + 1 if x else 0\n\n\nf(1)\nf(0)\n"
    for f, text in srcs.items():
        (box.repo / f).parent.mkdir(parents=True, exist_ok=True)
        (box.repo / f).write_text(text)
    for cmd in (
        ["add", "-A"],
        ["commit", "-q", "-m", "tests"],
        ["push", "-q", "origin", "HEAD:main"],
    ):
        subprocess.run(["git", *cmd], cwd=box.repo, env=env, check=True, capture_output=True)
    box.set_squeue("giap-nightly-summary", [f"9100 nightly:{JOB}"])
    return box


def _env(**extra: str) -> dict[str, str]:
    return {"SLURM_ARRAY_TASK_ID": "2", "SLURM_ARRAY_JOB_ID": JOB, "SLURM_JOB_ID": "5242",
            "SLURM_ARRAY_TASK_COUNT": "4", "SLURM_CPUS_PER_TASK": "2", "NIGHTLY_POST": "dry",
            "NIGHTLY_SLOW_NODE_GUARD": "0", "NIGHTLY_CAL_SHARDS": "2", **extra}  # fmt: skip


def _run(sb: Sandbox, **extra: str) -> subprocess.CompletedProcess[str]:
    return sb.run(
        ["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")], env=_env(**extra)
    )


def _shard_calls(sb: Sandbox) -> list[list[str]]:
    calls = [c for c in sb.sbatch_calls() if not c["failed"]]
    return [c["argv"] for c in calls if c["argv"] and c["argv"][-1].endswith("nightly.sbatch")]


def _sibling_env(argv: list[str]) -> dict[str, str]:
    """What Slurm hands the submitted shard: its --export pairs (ALL keeps the rest)."""
    return {k: v for k, v in export_pairs(argv).items() if k != "ALL"}


def _rows(sb: Sandbox) -> list[list[str]]:
    return [r for r in sb.tsv_rows() if len(r) >= 5 and r[1] == f"{JOB}_2"]


def _pytest_files(sb: Sandbox) -> list[list[str]]:
    return [[a for a in ln.split()[1:] if a.endswith(".py")]
            for ln in sb.calls_log() if ln.startswith("pytest ")]  # fmt: skip


def _diag(res: subprocess.CompletedProcess[str], sb: Sandbox) -> str:
    return f"\nstdout:\n{res.stdout}\nstderr:\n{res.stderr}\nrows: {sb.tsv_rows()}"


def _load(path: Path, name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def _readiness_view(row: list[str]) -> tuple[str, bool]:
    rd = _load(REPO / "scripts" / "release" / "readiness.py", "giap_readiness_shards")
    rows, bad = rd.parse_log("\t".join(row) + "\n")
    assert bad == 0 and len(rows) == 1, row
    return rows[0].label, rows[0].ok()


def _two_shards(sb: Sandbox, **shard2: str) -> subprocess.CompletedProcess[str]:
    res = _run(sb)
    (sub,) = _shard_calls(sb)
    assert not _rows(sb), "shard 1 wrote the row before shard 2 finished" + _diag(res, sb)
    return _run(sb, **_sibling_env(sub), **shard2)


# ---------------------------------------------------------------- 1. one merged row


def test_shard_one_submits_the_rest_on_the_same_commit(sb: Sandbox) -> None:
    sb.set_next_job_id(7001)
    res = _run(sb)
    diag = _diag(res, sb)
    subs = _shard_calls(sb)
    assert len(subs) == 1, f"expected one shard submission: {subs}" + diag
    assert option(subs[0], "array") == "2", subs[0]
    ex = _sibling_env(subs[0])
    assert (ex.get("NIGHTLY_CAL_SHARD"), ex.get("NIGHTLY_CAL_SHARDS")) == ("2", "2"), ex
    assert ex.get("NIGHTLY_ARRAY_JOB") == JOB and ex.get("NIGHTLY_EXPECT") == "4", ex
    assert sb.head.startswith(ex.get("NIGHTLY_REF", "-")), f"commit not pinned: {ex}"
    assert ex.get("NIGHTLY_WORKTREE", "").endswith("main-t2-s2"), ex
    assert any(re.search(rf"scontrol .*afterany:{JOB}:7001\b", ln) for ln in sb.calls_log()), diag
    assert res.returncode == 0 and not _rows(sb), diag


def test_two_shards_write_one_merged_row(sb: Sandbox) -> None:
    res = _two_shards(sb)
    diag = _diag(res, sb)
    rows = _rows(sb)
    assert len(rows) == 1, "expected exactly one calibration_nightly row" + diag
    assert rows[0][3] == "0", rows[0]
    assert rows[0][4].startswith(
        "calibration_nightly: 20 passed, 2 skipped, 4 warnings in 2.50s"
    ), rows[0]
    assert "coverage 95%" in rows[0][4], rows[0]
    assert _readiness_view(rows[0]) == ("calibration_nightly", True)
    shards = _pytest_files(sb)
    assert len(shards) == 2 and sorted(shards[0] + shards[1]) == FILES, shards
    comb = [ln for ln in sb.calls_log() if ln.startswith("coverage combine")]
    assert len(comb) == 1 and "coverage.s1" in comb[0] and "coverage.s2" in comb[0], comb
    rep = [ln for ln in sb.calls_log() if ln.startswith("coverage report")]
    assert rep and "--fail-under=90" in rep[0], rep


def test_real_coverage_keeps_the_shard_data_and_reports_the_figure(sb: Sandbox) -> None:
    """``coverage report`` combines ``<data-file>.*``: the shard files must survive (#420 B1)."""
    res = _two_shards(sb, FAKE_COV="real")
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == "0", _diag(res, sb)
    assert "coverage 100% (min 90%)" in rows[0][4], rows[0]
    d = sb.root / "ci" / "runs" / f"{JOB}_2.shards"
    for name in ("coverage.s1", "coverage.s2", "combined.coverage", "coverage_report.txt"):
        assert (d / name).is_file(), f"{name} is gone: {sorted(p.name for p in d.iterdir())}"
    assert "TOTAL" in (d / "coverage_report.txt").read_text()


# ---------------------------------------------------------------- 2. failures fail the row


@pytest.mark.parametrize(
    ("extra", "rc"),
    [({"FAKE_SHARD_PYTEST": "fail"}, "1"), ({"FAKE_COV": "low"}, "2")],
    ids=["failing-shard", "coverage-below-min"],
)
def test_a_failing_shard_or_low_coverage_fails_the_row(
    sb: Sandbox, extra: dict[str, str], rc: str
) -> None:
    res = _two_shards(sb, **extra)
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == rc, _diag(res, sb)
    assert _readiness_view(rows[0]) == ("calibration_nightly", False)


def test_failed_counts_are_summed(sb: Sandbox) -> None:
    res = _two_shards(sb, FAKE_SHARD_PYTEST="fail")
    (row,) = _rows(sb)
    assert row[4].startswith("calibration_nightly: 2 failed, 20 passed, 1 skipped"), _diag(res, sb)


def test_refused_shard_submission_still_writes_a_failing_row(sb: Sandbox) -> None:
    res = _run(sb, FAKE_SBATCH_FAIL="all")
    rows = _rows(sb)
    assert len(rows) == 1, "no row after a refused shard submission" + _diag(res, sb)
    assert rows[0][3] == "2" and "refused" in rows[0][4], rows[0]
    assert _readiness_view(rows[0]) == ("calibration_nightly", False)


def test_too_many_shards_are_refused(sb: Sandbox) -> None:
    res = _run(sb, NIGHTLY_CAL_SHARDS="5")
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == "2" and "R14" in rows[0][4], _diag(res, sb)
    assert not _shard_calls(sb) and not _pytest_files(sb), "5 shards were not refused"


def _summary(sb: Sandbox) -> subprocess.CompletedProcess[str]:
    script = sb.repo / "scripts" / "slurm" / "nightly_summary.sh"
    return sb.run(["bash", str(script), JOB, "--dry-run"], env={"NIGHTLY_EXPECT": "4"})


def test_summary_job_fails_the_row_of_a_shard_that_never_recorded(sb: Sandbox) -> None:
    """A shard SIGKILLed or cancelled leaves no record: no older pass may stand in (#420 B2)."""
    sb.set_next_job_id(7001)
    res = _run(sb)
    assert not _rows(sb), _diag(res, sb)
    sb.set_sacct([{"JobID": "7001", "JobName": "giap-nightly", "State": "CANCELLED"}])
    for _ in range(2):  # the second run must not add a row
        res = _summary(sb)
    rows = _rows(sb)
    assert len(rows) == 1, "expected one failing row from the summary job" + _diag(res, sb)
    assert rows[0][3] == "1" and sb.head.startswith(rows[0][2]), rows[0]
    assert "shard 2/2 left no record (sacct CANCELLED)" in rows[0][4], rows[0]
    assert _readiness_view(rows[0]) == ("calibration_nightly", False)


def test_summary_job_adds_nothing_when_the_shards_wrote_the_row(sb: Sandbox) -> None:
    res = _two_shards(sb)
    res = _summary(sb)
    assert len(_rows(sb)) == 1, _diag(res, sb)


# ---------------------------------------------------------------- 3. a killed shard


def test_killed_shard_names_its_failures_so_far(sb: Sandbox) -> None:
    res = _two_shards(sb, FAKE_SHARD_PYTEST="term")
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == "143", _diag(res, sb)
    assert re.search(r"SIGTERM.*1 failed so far", rows[0][4]), rows[0]
    so_far = next((sb.root / "ci" / "runs").glob(f"{JOB}_2-s2/failed_so_far.txt"))
    assert so_far.read_text().split() == ["FAILED", "tests/test_a.py::test_x"]


def test_unsharded_kill_also_names_failures_so_far(sb: Sandbox) -> None:
    res = _run(sb, NIGHTLY_CAL_SHARDS="1", FAKE_SHARD_PYTEST="term")
    assert not _shard_calls(sb), "an unsharded task 2 submitted shards"
    rows = _rows(sb)
    assert len(rows) == 1 and rows[0][3] == "143", _diag(res, sb)
    assert "1 failed so far" in rows[0][4], rows[0]


# ---------------------------------------------------------------- 4. launcher and R14


def test_launcher_asks_for_two_shards_within_r14(sb: Sandbox) -> None:
    res = sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly_launch.sh")])
    arrays = [c for c in sb.sbatch_calls() if option(c["argv"], "array", "a") == "0-3"]
    assert len(arrays) == 1, _diag(res, sb)
    shards = int(arrays[0]["env"].get("NIGHTLY_CAL_SHARDS", "1"))
    assert shards == 2, arrays[0]["env"]
    text = NIGHTLY.read_text()
    cpus = int(re.search(r"^#SBATCH --cpus-per-task=(\d+)", text, re.M)[1])  # type: ignore[index]
    arr = re.search(r"^#SBATCH --array=(\d+)-(\d+)", text, re.M)
    assert arr, "no #SBATCH --array range"
    tasks = int(arr[2]) - int(arr[1]) + 1
    assert tasks * cpus + (shards - 1) * cpus + 1 <= 128


# ---------------------------------------------------------------- the fold itself

NS = _load(REPO / "scripts" / "ci" / "nightly_shards.py", "giap_nightly_shards")


def test_fold_sums_counts_keeps_the_longest_time_and_notes() -> None:
    lab = "calibration_nightly: "
    rc, text = NS.fold([
        {"shard": 2, "exit": 143, "summary": lab + "task 2 killed by SIGTERM"},
        {"shard": 1, "exit": 1, "summary": lab + "1 failed, 5 passed, 1 warning in 3700.50s"},
        {"shard": 3, "exit": 0, "summary": lab + "7 passed, 2 warnings, 1 error in 10.00s"},
    ], 3)  # fmt: skip
    assert rc == 143
    assert text.startswith("1 failed, 12 passed, 3 warnings, 1 errors in 3700.50s (1:01:40"), text
    assert text.endswith("; shard 2/3: task 2 killed by SIGTERM"), text


def test_merge_waits_for_every_shard(tmp_path: Path) -> None:
    NS.record(tmp_path, 1, 0, "calibration_nightly: 3 passed in 1.00s")
    assert NS.merge(tmp_path, 2) is None
    NS.record(tmp_path, 2, 0, "calibration_nightly: 4 passed in 2.00s")
    assert NS.merge(tmp_path, 2) == (0, "7 passed in 2.00s (0:00:02, longest of 2 shards)")
