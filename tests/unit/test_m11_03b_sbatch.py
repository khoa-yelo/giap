"""``calibration_full.sbatch`` writes a JUnit file per pytest run and a task record
(issue M11-03b, acceptance criteria 4 and 5).

Static checks of the script text, then runs under bash with a stub ``uv`` (the pattern of
tests/unit/test_calibration_full_array.py): the stub logs its call and ``PYTEST_ADDOPTS``, writes a
one-testcase JUnit file wherever ``--junitxml=`` points (as pytest would), and exits by role. The
checkout (``GIAP_WORKTREE``) is a git repository in ``tmp_path`` holding a copy of
``tests/calibration/tolerances.yaml`` (the ``replicates_full`` lookup of the suite name).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest

REPO = Path(__file__).resolve().parents[2]
SBATCH = REPO / "scripts" / "slurm" / "calibration_full.sbatch"

STUB_UV = r"""#!/bin/bash
# stub uv: log the call, write a JUnit file at --junitxml=, exit by role
echo "ADDOPTS=${PYTEST_ADDOPTS-<unset>} ARGS=$*" >> "$STUB_LOG"
j=$(printf '%s\n' "${PYTEST_ADDOPTS-}" | grep -oE -- '--junitxml=[^ ]+' | tail -n 1)
j=${j#--junitxml=}
if [ -n "$j" ]; then
  mkdir -p "$(dirname "$j")"
  {
    printf '%s' '<?xml version="1.0" encoding="utf-8"?><testsuites>'
    printf '%s' '<testsuite name="pytest" tests="1"><testcase name="test_a" time="0.1"'
    printf '%s' ' classname="tests.calibration.attacks.test_attack_x" /></testsuite></testsuites>'
  } > "$j"
fi
if [ -n "${STUB_TERM:-}" ]; then kill -TERM "$PPID"; sleep 1; fi
args=" $* "
case "$args" in *" --merge "*|*" --merge="*) exit "${STUB_RC_MERGE:-0}" ;; esac
if [[ "$args" =~ --shard[=\ ]([0-9]+)/([0-9]+) ]]; then
  v="STUB_RC_SHARD_${BASH_REMATCH[1]}"; exit "${!v:-0}"
fi
exit "${STUB_RC_PLAIN:-0}"
"""

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture", "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture", "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
}  # fmt: skip


def _git(repo: Path, *args: str) -> str:
    out = subprocess.run(["git", "-C", str(repo), *args], env={**os.environ, **_GIT_ENV},
                         check=True, capture_output=True, text=True)  # fmt: skip
    return out.stdout.strip()


def make_checkout(tmp: Path) -> tuple[Path, str]:
    repo = tmp / "repo"
    (repo / "scripts").mkdir(parents=True)
    (repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    (repo / "tests" / "calibration").mkdir(parents=True)
    shutil.copyfile(REPO / "tests" / "calibration" / "tolerances.yaml",
                    repo / "tests" / "calibration" / "tolerances.yaml")  # fmt: skip
    (repo / "tests" / "sim").mkdir()
    (repo / "tests" / "sim" / "keep.py").write_text("# sim\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "fixture")
    return repo, _git(repo, "rev-parse", "HEAD")


def run_task(tmp: Path, task: int, *args: str, env: dict[str, str] | None = None
             ) -> tuple[subprocess.CompletedProcess[str], list[str], Path]:  # fmt: skip
    repo = tmp / "repo"
    bindir = tmp / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "uv").write_text(STUB_UV)
    (bindir / "uv").chmod(0o755)
    calls_log, log = tmp / "calls.txt", tmp / "log"
    calls_log.unlink(missing_ok=True)
    full = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp),
        "GIAP_REPO": str(repo),
        "GIAP_WORKTREE": str(repo),
        "GIAP_ROOT": str(tmp),
        "GIAP_CALIBRATION_LOG_DIR": str(log),
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": "1900000",
        "STUB_LOG": str(calls_log),
        **_GIT_ENV,
        **(env or {}),
    }
    if "UV_PROJECT_ENVIRONMENT" in os.environ:
        full.setdefault("UV_PROJECT_ENVIRONMENT", os.environ["UV_PROJECT_ENVIRONMENT"])
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run([bash, str(SBATCH), *args], env=full, capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip
    calls = calls_log.read_text().splitlines() if calls_log.exists() else []
    return proc, calls, log


def records(log: Path) -> list[Path]:
    return sorted(log.glob("task-*.json"))


def read_record(log: Path, task: int) -> dict[str, Any]:
    recs = records(log)
    assert [p.name for p in recs] == [f"task-{task}.json"], recs
    doc: dict[str, Any] = json.loads(recs[0].read_text())
    return doc


def junit_of(call: str) -> str | None:
    m = re.findall(r"--junitxml=(\S+)", call)
    return m[-1] if m else None


@pytest.fixture
def checkout(tmp_path: Path) -> tuple[Path, str]:
    return make_checkout(tmp_path)


# ------------------------------------------------------------------ static (AC4)


def test_array_and_cpus_unchanged() -> None:
    text = SBATCH.read_text()
    assert re.findall(r"^#SBATCH --array=(\S+)\s*$", text, re.M) == ["0-14"]
    assert re.findall(r"^#SBATCH --cpus-per-task=(\S+)\s*$", text, re.M) == ["4"]
    assert re.findall(r"^#SBATCH --time=(\S+)\s*$", text, re.M) == ["04:00:00"]


def test_script_sets_junitxml_under_the_log_dir() -> None:
    text = SBATCH.read_text()
    lines = [ln for ln in text.splitlines() if "PYTEST_ADDOPTS" in ln and not ln.lstrip()
             .startswith("#")]  # fmt: skip
    assert lines, "no PYTEST_ADDOPTS assignment"
    assert any("--junitxml=" in ln and "GIAP_CALIBRATION_LOG_DIR" in ln for ln in lines), lines
    assert "${PYTEST_ADDOPTS:+" in text, "an existing PYTEST_ADDOPTS must be kept"


def test_script_traps_exit_and_term_for_the_record() -> None:
    code = "\n".join(ln for ln in SBATCH.read_text().splitlines()
                     if not ln.lstrip().startswith("#"))  # fmt: skip
    traps = re.findall(r"^\s*trap\s+(.+)$", code, re.M)
    assert any(re.search(r"\bEXIT\b", t) for t in traps), traps
    assert any(re.search(r"\bTERM\b", t) for t in traps), traps
    assert "task-" in code and ".json" in code


# ------------------------------------------------------------------ stub runs (AC5)


def test_unsharded_task_writes_its_junit_and_record(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    _repo, head = checkout
    proc, calls, log = run_task(tmp_path, 3)
    assert proc.returncode == 0, proc.stderr
    assert len(calls) == 1 and junit_of(calls[0]) == f"{log}/task-3.run.junit.xml", calls
    rec = read_record(log, 3)
    assert {k: rec[k] for k in ("task", "suite", "job", "commit", "dirty", "exit", "args",
                                "sharded", "parts")} == {
        "task": 3, "suite": "regions", "job": "1900000", "commit": head, "dirty": False,
        "exit": 0, "args": [], "sharded": False, "parts": ["task-3.run.junit.xml"],
    }  # fmt: skip
    assert re.fullmatch(r"\d{4}-\d\d-\d\d", rec["date"]), rec
    assert {p.name for p in log.iterdir()} == {"task-3.json", "task-3.run.junit.xml"}


def test_task_exit_1_is_recorded(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    proc, _calls, log = run_task(tmp_path, 3, env={"STUB_RC_PLAIN": "1"})
    assert proc.returncode == 1, proc.stderr
    assert read_record(log, 3)["exit"] == 1


def test_sigterm_records_exit_143_once(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    proc, _calls, log = run_task(tmp_path, 3, env={"STUB_TERM": "1"})
    assert proc.returncode == 143, (proc.returncode, proc.stderr)
    rec = read_record(log, 3)
    assert rec["exit"] == 143, rec
    assert not [p for p in log.iterdir() if p.name not in ("task-3.json",
                                                           "task-3.run.junit.xml")]  # fmt: skip


def test_task13_gives_each_shard_and_the_merge_its_own_junit(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    proc, calls, log = run_task(tmp_path, 13, env={"STUB_RC_MERGE": "1"})
    assert proc.returncode == 1, proc.stderr
    assert len(calls) == 3, calls
    got = {junit_of(c) for c in calls}
    want = {f"{log}/task-13.shard0.junit.xml", f"{log}/task-13.shard1.junit.xml",
            f"{log}/task-13.merge.junit.xml"}  # fmt: skip
    assert got == want, calls
    merge = [c for c in calls if "--merge" in c.split()]
    assert merge and junit_of(merge[0]) == f"{log}/task-13.merge.junit.xml", calls
    rec = read_record(log, 13)
    assert rec["sharded"] is True and rec["suite"] == "variants" and rec["exit"] == 1, rec
    assert sorted(rec["parts"]) == sorted(Path(p).name for p in want), rec


def test_harness_arguments_are_recorded(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    proc, _calls, log = run_task(tmp_path, 3, "--fast")
    assert proc.returncode == 0, proc.stderr
    assert read_record(log, 3)["args"] == ["--fast"]


def test_existing_pytest_addopts_is_kept(tmp_path: Path, checkout) -> None:  # type: ignore[no-untyped-def]
    _proc, calls, log = run_task(tmp_path, 3, env={"PYTEST_ADDOPTS": "-x"})
    assert calls and calls[0].startswith("ADDOPTS=-x --junitxml="), calls
    assert junit_of(calls[0]) == f"{log}/task-3.run.junit.xml"


@pytest.mark.parametrize(("rel", "dirty"), [("tests/sim/new.py", True),
                                            ("tests/helpers/new.py", True),
                                            ("tests/calibration/new.py", True),
                                            ("notes.txt", False)])  # fmt: skip
def test_dirty_flag_covers_the_calibration_inputs(  # type: ignore[no-untyped-def]
    tmp_path: Path, checkout, rel: str, dirty: bool
) -> None:
    repo, _head = checkout
    (repo / rel).parent.mkdir(parents=True, exist_ok=True)
    (repo / rel).write_text("# change\n")
    proc, _calls, log = run_task(tmp_path, 3)
    assert proc.returncode == 0, proc.stderr
    assert read_record(log, 3)["dirty"] is dirty
