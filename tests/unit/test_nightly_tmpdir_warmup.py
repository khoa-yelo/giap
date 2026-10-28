"""``scripts/slurm/nightly.sbatch`` integration task: node-local TMPDIR and the Nextflow warm-up
(#323 review S1, N3; M9-13 decision 5; M2-21 decision 23).

Run for real in the M2-21 sandbox (``tests/helpers/slurm_sandbox.py``), with this module's own
``python`` wrapper in the sandbox ``bin/`` that records, at each ``python -m pytest`` call, the
``TMPDIR`` it sees and whether that directory exists. Pinned:

* TMPDIR is a fresh directory under ``/tmp`` (node-local), not under the task's ``$RUN`` on NFS
  (``$GIPPA_ROOT/ci/runs/...``), and it is removed when the task ends (also after SIGTERM);
* ``nextflow -version`` runs before the integration pytest; nf-test still starts only after the
  pytest run ends (decision 23).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from helpers.slurm_sandbox import (
    _PYTHON,
    REPO,
    Sandbox,
    assert_honours_gippa_repo,
    make_sandbox,
    write_exe,
)

NIGHTLY = REPO / "scripts" / "slurm" / "nightly.sbatch"
JOB = "4242"

_RECORD = r"""
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "pytest" ]; then
  if [ -d "${TMPDIR:-/nonexistent}" ]; then e=1; else e=0; fi
  printf 'pytest TMPDIR=%s exists=%s\n' "${TMPDIR:-}" "$e" >> "$FAKE_SLURM_DIR/calls.log"
fi
"""


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    assert_honours_gippa_repo(NIGHTLY.read_text(), "nightly.sbatch")
    box = make_sandbox(tmp_path)
    head, _, rest = _PYTHON.replace("{python}", sys.executable).partition("\n")
    for name in ("python", "python3"):
        write_exe(box.bin / name, f"{head}\n{_RECORD}{rest}")
    return box


def _integration(sb: Sandbox, **env: str):  # type: ignore[no-untyped-def]
    base = {
        "SLURM_ARRAY_TASK_ID": "0",
        "SLURM_ARRAY_JOB_ID": JOB,
        "SLURM_JOB_ID": "5000",
        "SLURM_ARRAY_TASK_COUNT": "4",
        "SLURM_CPUS_PER_TASK": "2",
        "SLURM_JOB_NAME": "gippa-nightly",
        "NIGHTLY_POST": "dry",
    }
    return sb.run(["bash", str(sb.repo / "scripts" / "slurm" / "nightly.sbatch")],
                  env={**base, **env})  # fmt: skip


def _pytest_tmpdirs(sb: Sandbox) -> list[tuple[str, str]]:
    return [
        (m[1], m[2])
        for ln in sb.calls_log()
        if (m := re.fullmatch(r"pytest TMPDIR=(.*) exists=([01])", ln))
    ]


def _index(lines: list[str], pattern: str) -> list[int]:
    return [i for i, ln in enumerate(lines) if re.match(pattern, ln)]


@pytest.mark.parametrize("pytest_mode", ["ok", "term"])
def test_tmpdir_is_node_local_and_removed_at_exit(sb: Sandbox, pytest_mode: str) -> None:
    res = _integration(sb, FAKE_PYTEST=pytest_mode)
    seen = _pytest_tmpdirs(sb)
    assert seen, f"the integration pytest never ran\n{res.stdout}\n{res.stderr}"
    run_root = str(sb.root / "ci" / "runs")
    for tmpdir, exists in seen:
        assert tmpdir, "TMPDIR is unset for pytest"
        assert exists == "1", f"TMPDIR {tmpdir} does not exist while pytest runs"
        p = Path(tmpdir)
        assert p.is_absolute() and p.parts[:2] == ("/", "tmp"), f"TMPDIR {tmpdir} is not under /tmp"
        assert not tmpdir.startswith(run_root), f"TMPDIR {tmpdir} is under $RUN (NFS)"
        assert not tmpdir.startswith(str(sb.root)), f"TMPDIR {tmpdir} is under $GIPPA_ROOT"
        assert p != Path("/tmp"), "TMPDIR is /tmp itself, not a task directory"
        assert not p.exists(), f"TMPDIR {tmpdir} is left behind after the task ({pytest_mode})"
    runs = sorted((sb.root / "ci" / "runs").glob(f"{JOB}_0"))
    assert runs, res.stdout
    assert not (runs[0] / "tmp").exists(), "the task still creates $RUN/tmp"


def test_nextflow_warm_up_runs_before_pytest_and_nf_test_after(sb: Sandbox) -> None:
    res = _integration(sb)
    log = sb.calls_log()
    warm = _index(log, r"nextflow\s+-version\b")
    tests = _index(log, r"pytest TMPDIR=")
    nf = _index(log, r"nf-test\b")
    assert tests, f"the integration pytest never ran\n{log}\n{res.stdout}\n{res.stderr}"
    assert warm, f"no `nextflow -version` warm-up\n{log}"
    assert min(warm) < min(tests), f"the warm-up runs after pytest starts\n{log}"
    assert nf, f"nf-test never ran\n{log}"
    assert min(nf) > max(tests), f"nf-test starts before the pytest run ends (decision 23)\n{log}"


def test_warm_up_is_not_nf_test(sb: Sandbox) -> None:
    """Decision 23: nothing named nf-test (``nf-test version`` included) runs before pytest."""
    _integration(sb)
    log = sb.calls_log()
    tests = _index(log, r"pytest TMPDIR=")
    assert tests
    assert not [ln for ln in log[: min(tests)] if ln.startswith("nf-test")], log
