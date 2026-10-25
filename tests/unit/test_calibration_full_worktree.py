"""``calibration_full.sbatch`` with ``GIAP_WORKTREE`` tests the worktree's code, not main's.

Slurm 1885838 ran M6-13's branch tests from a worktree but imported giap from the shared venv's
editable install (the main checkout's src) and failed with ``TypeError: discordance() got an
unexpected keyword argument 'sigma_model'``. Pinned here:

- with ``GIAP_WORKTREE`` set the job exports
  ``PYTHONPATH=$GIAP_WORKTREE/src:$GIAP_WORKTREE/tests/sim:$GIAP_WORKTREE/tests`` (as
  test_branch.sbatch does) before the harness runs; without it PYTHONPATH is left alone;
- fail closed: under ``GIAP_WORKTREE`` the harness (``_criteria.py``, every mode including
  ``--dry-run`` and ``--merge``) prints ``giap.__file__`` and exits 2 unless it lies under
  ``$GIAP_WORKTREE/src``.

The sbatch runs under bash with a stub ``uv`` (the Slurm script test kit of
test_calibration_full_array.py); the guard runs the real harness in a subprocess.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SBATCH = REPO / "scripts" / "slurm" / "calibration_full.sbatch"
HARNESS = REPO / "tests" / "calibration" / "_criteria.py"

STUB_UV = r"""#!/bin/bash
# stub uv: log PYTHONPATH (or <unset>) and the call
echo "PYTHONPATH=${PYTHONPATH-<unset>} ARGS=$*" >> "$STUB_LOG"
exit 0
"""


def run_sbatch(tmp: Path, worktree: Path | None, task: int = 3) -> list[str]:
    repo = tmp / "repo"
    (repo / "scripts").mkdir(parents=True, exist_ok=True)
    (repo / "scripts" / "farmshare_env.sh").write_text("# stub\n")
    bindir = tmp / "bin"
    bindir.mkdir(exist_ok=True)
    (bindir / "uv").write_text(STUB_UV)
    (bindir / "uv").chmod(0o755)
    log = tmp / "calls.txt"
    log.unlink(missing_ok=True)
    env = {
        "PATH": f"{bindir}{os.pathsep}{os.environ['PATH']}",
        "HOME": str(tmp),
        "GIAP_REPO": str(repo),
        "GIAP_ROOT": str(tmp),
        "GIAP_CALIBRATION_LOG_DIR": str(tmp / "log"),
        "SLURM_ARRAY_TASK_ID": str(task),
        "SLURM_ARRAY_JOB_ID": "1",
        "STUB_LOG": str(log),
    }
    if worktree is not None:
        env["GIAP_WORKTREE"] = str(worktree)
    bash = shutil.which("bash") or "/bin/bash"
    proc = subprocess.run([bash, str(SBATCH)], env=env, capture_output=True, text=True,
                          timeout=60, check=False)  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    return log.read_text().splitlines()


def test_worktree_job_exports_the_worktree_pythonpath(tmp_path: Path) -> None:
    wt = tmp_path / "wt"
    wt.mkdir()
    calls = run_sbatch(tmp_path, wt)
    want = f"PYTHONPATH={wt}/src:{wt}/tests/sim:{wt}/tests "
    assert calls and all(c.startswith(want) for c in calls), calls


def test_task13_shards_and_merge_get_the_worktree_pythonpath(tmp_path: Path) -> None:
    wt = tmp_path / "wt"
    wt.mkdir()
    calls = run_sbatch(tmp_path, wt, task=13)
    assert len(calls) == 3 and all(c.startswith(f"PYTHONPATH={wt}/src:") for c in calls), calls


def test_main_checkout_job_leaves_pythonpath_alone(tmp_path: Path) -> None:
    calls = run_sbatch(tmp_path, None)
    assert calls and all(c.startswith("PYTHONPATH=<unset> ") for c in calls), calls


def _harness(worktree: Path, pythonpath: str, *args: str) -> subprocess.CompletedProcess[str]:
    keep = {k: v for k, v in os.environ.items() if not k.startswith("GIAP_CALIBRATION_")}
    env = {**keep, "GIAP_WORKTREE": str(worktree), "PYTHONPATH": pythonpath,
           "GIAP_CALIBRATION_LOG_DIR": str(worktree / "log")}  # fmt: skip
    return subprocess.run([sys.executable, str(HARNESS), "--suite", "discordance", *args],
                          env=env, capture_output=True, text=True, timeout=120, check=False,
                          cwd=REPO)  # fmt: skip


def test_harness_refuses_giap_from_outside_the_worktree(tmp_path: Path) -> None:
    """``GIAP_WORKTREE`` names another checkout but giap comes from this one: exit 2."""
    other = tmp_path / "other_worktree"
    (other / "src").mkdir(parents=True)
    proc = _harness(other, f"{REPO}/src:{REPO}/tests/sim", "--dry-run")
    assert proc.returncode == 2, (proc.stdout, proc.stderr)
    assert "giap.__file__" in proc.stderr and str(other) in proc.stderr
    assert "pytest" not in proc.stdout  # nothing ran


@pytest.mark.parametrize("args", [("--dry-run",), ("--merge",)])
def test_harness_refusal_covers_merge_too(tmp_path: Path, args: tuple[str, ...]) -> None:
    other = tmp_path / "other_worktree"
    (other / "src").mkdir(parents=True)
    proc = _harness(other, f"{REPO}/src:{REPO}/tests/sim", *args)
    assert proc.returncode == 2 and "giap.__file__" in proc.stderr, (proc.stdout, proc.stderr)


def test_harness_accepts_giap_from_the_worktree(tmp_path: Path) -> None:
    proc = _harness(REPO, f"{REPO}/src:{REPO}/tests/sim", "--dry-run")
    assert proc.returncode == 0, proc.stderr
    assert f"giap.__file__={(REPO / 'src' / 'giap' / '__init__.py').resolve()}" in proc.stderr
