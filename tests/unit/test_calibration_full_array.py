"""The nightly calibration array covers every ``replicates_full`` suite (docs/specs/
M7-13-calibration-gate.md § Array extension; ADR-0042; backlog line 614; GIPPA_SPEC.md §8.4 "one
suite per task"): ``--array=0-(n-1)`` and the header comment "<n> suites" with n the number of
``replicates_full`` keys (15 today: task 13 ``variants``, 14 ``compare``)."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "calibration"))
import _criteria as C  # noqa: E402

SBATCH = REPO / "scripts" / "slurm" / "calibration_full.sbatch"


def test_array_spans_every_suite() -> None:
    text = SBATCH.read_text()
    suites = C.array_suites()
    assert len(suites) == 15 and suites[13:] == ("variants", "compare")
    arrays = re.findall(r"^#SBATCH --array=(\d+)-(\d+)\s*$", text, flags=re.M)
    assert arrays == [("0", str(len(suites) - 1))], arrays
    assert f"{len(suites)} suites" in text


# --------------------------------------------------------------------------- task-13 wrapper
# Lead decision 7 (statistician B1/B2, code review B1/B2/N4): task 13 runs shards 0/n .. (n-1)/n
# (``shard_replicates`` is 0-based), then ``--merge``; the task fails if any shard or the merge
# is non-zero; ``--time`` >= 04:00:00; ``--shard=K/N`` / ``--merge`` given by hand are not split.
# The wrapper runs under bash with a stub ``uv`` on PATH (no pytest, no GIPPA run).

STUB_UV = r"""#!/bin/bash
# stub uv: log the call, exit by role (STUB_RC_SHARD_<k>, STUB_RC_MERGE, STUB_RC_PLAIN)
echo "$*" >> "$STUB_LOG"
args=" $* "
case "$args" in *" --merge "*|*" --merge="*) exit "${STUB_RC_MERGE:-0}" ;; esac
if [[ "$args" =~ --shard[=\ ]([0-9]+)/([0-9]+) ]]; then
  v="STUB_RC_SHARD_${BASH_REMATCH[1]}"; exit "${!v:-0}"
fi
exit "${STUB_RC_PLAIN:-0}"
"""


def run_wrapper(
    tmp: Path, task: int, *args: str, env: dict[str, str] | None = None
) -> tuple[int, list[str]]:
    repo = tmp / "repo"
    (repo / "scripts").mkdir(parents=True, exist_ok=True)
    (repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    bindir = tmp / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "uv").write_text(STUB_UV)
    (bindir / "uv").chmod(0o755)
    log = tmp / "calls.txt"
    log.unlink(missing_ok=True)
    full = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp),
        "GIPPA_REPO": str(repo),
        "GIPPA_WORKTREE": str(repo),
        "GIPPA_ROOT": str(tmp),
        "GIPPA_CALIBRATION_LOG_DIR": str(tmp / "log"),
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": "1",
        "STUB_LOG": str(log),
        **(env or {}),
    }
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run([bash, str(SBATCH), *args], env=full, capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip
    calls = log.read_text().splitlines() if log.exists() else []
    return proc.returncode, calls


def shard_of(call: str) -> str | None:
    m = re.search(r"--shard[= ](\S+)", call)
    return m.group(1) if m else None


def test_array_time_fits_shards_and_block_two() -> None:
    (t,) = re.findall(r"^#SBATCH --time=(\S+)\s*$", SBATCH.read_text(), flags=re.M)
    parts = [int(x) for x in re.split(r"[-:]", t)]
    days, (h, m, s) = (parts[0], parts[1:]) if len(parts) == 4 else (0, parts)
    assert days * 86400 + h * 3600 + m * 60 + s >= 4 * 3600, t


def test_task13_runs_shards_from_zero_then_merges(tmp_path: Path) -> None:
    rc, calls = run_wrapper(tmp_path, 13)
    assert rc == 0, calls
    assert len(calls) == 3 and "--merge" in calls[-1].split(), calls
    assert sorted(shard_of(c) for c in calls[:-1]) == ["0/2", "1/2"], calls
    assert shard_of(calls[-1]) is None
    assert all("--task-id 13" in c and "_criteria.py" in c for c in calls), calls


def test_task13_shard_indices_match_shard_replicates(tmp_path: Path) -> None:
    """Every shard the wrapper starts is a valid ``shard_replicates`` block, and together they
    cover range(R) exactly (job 1872876 ran 50:100 as "shard 1/2"; 0:50 was never run)."""
    _rc, calls = run_wrapper(tmp_path, 13)
    R = C.replicates("variants", True)
    got: list[int] = []
    for c in calls[:-1]:
        k, n = (int(x) for x in str(shard_of(c)).split("/"))
        got.extend(C.shard_replicates(R, n, k))
    assert sorted(got) == list(range(R)), calls


def test_task13_passes_extra_arguments_to_every_call(tmp_path: Path) -> None:
    rc, calls = run_wrapper(tmp_path, 13, "--fast")
    assert rc == 0 and len(calls) == 3 and all("--fast" in c.split() for c in calls), calls


@pytest.mark.parametrize(
    ("env", "what"),
    [({"STUB_RC_SHARD_0": "1"}, "first shard pytest fails"),
     ({"STUB_RC_SHARD_1": "2"}, "second shard pytest usage error"),
     ({"STUB_RC_SHARD_0": "1", "STUB_RC_SHARD_1": "2"}, "both shards"),
     ({"STUB_RC_MERGE": "1"}, "merge: confirmed fail"),
     ({"STUB_RC_MERGE": "2"}, "merge error")],
)  # fmt: skip
def test_task13_fails_when_a_shard_or_the_merge_fails(
    tmp_path: Path, env: dict[str, str], what: str
) -> None:
    rc, calls = run_wrapper(tmp_path, 13, env=env)
    assert rc != 0, (what, calls)


def test_other_tasks_run_once_and_keep_the_pytest_exit(tmp_path: Path) -> None:
    rc, calls = run_wrapper(tmp_path, 3)
    assert rc == 0 and len(calls) == 1 and shard_of(calls[0]) is None, calls
    rc, calls = run_wrapper(tmp_path, 3, env={"STUB_RC_PLAIN": "1"})
    assert rc == 1 and len(calls) == 1, calls


@pytest.mark.parametrize(
    "args",
    [("--shard", "1/4"), ("--shard=1/4",), ("--merge",), ("--dry-run",),
     ("--shard=0/2", "--block", "2")],
)  # fmt: skip
def test_task13_hand_given_shard_or_merge_is_not_split(
    tmp_path: Path, args: tuple[str, ...]
) -> None:
    """N4: a ``--shard``/``--merge``/``--dry-run`` given by hand (either form) runs once, as
    given; the wrapper never appends its own ``--shard k/2`` (argparse keeps the last one)."""
    rc, calls = run_wrapper(tmp_path, 13, *args)
    assert rc == 0 and len(calls) == 1, calls
    assert calls[0].split()[-len(args) :] == list(args), calls
    assert len(re.findall(r"--shard", calls[0])) == sum(a.startswith("--shard") for a in args)
