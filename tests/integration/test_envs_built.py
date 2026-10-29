"""Built tool envs report their pinned version and have lock files (M2-01 lead decision 6).

Runs on FarmShare after ``sbatch scripts/slurm/build_envs.sbatch``; skipped without micromamba
(``needs_conda``). A tool whose env is not under ``$MAMBA_ROOT_PREFIX/envs`` is skipped
(``all`` leaves container-run and opt-in tools unbuilt, decision 4).
"""

from __future__ import annotations

import hashlib
import os
import shlex
import subprocess
from pathlib import Path

import pytest

from gippa.annotate.tools import TOOLS

pytestmark = pytest.mark.needs_conda

REPO = Path(__file__).resolve().parents[2]
_GIPPA_REPO = Path(os.environ.get("GIPPA_REPO", REPO))
_ROOT = Path(os.environ.get("GIPPA_ROOT", _GIPPA_REPO / "local"))
PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", _ROOT / "micromamba"))
LOCKS = REPO / "workflow" / "envs" / "locks"
CONDA_TOOLS = sorted(n for n, t in TOOLS.items() if not t.in_process)


def _env_dir(tool: str) -> Path:
    env = PREFIX / "envs" / TOOLS[tool].env_name
    if not env.is_dir():
        pytest.skip(f"{TOOLS[tool].env_name} not built under {PREFIX / 'envs'}")
    return env


@pytest.mark.parametrize("tool", CONDA_TOOLS)
def test_version_command_prints_pinned_version(tool: str) -> None:
    _env_dir(tool)
    spec = TOOLS[tool]
    proc = subprocess.run(
        f"micromamba run -n {spec.env_name} {shlex.join(spec.version_cmd)}",
        shell=True,
        capture_output=True,
        text=True,
        timeout=600,
    )
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out
    assert spec.version.removeprefix("git:") in out, out


@pytest.mark.parametrize("tool", CONDA_TOOLS)
def test_lock_files_present(tool: str) -> None:
    env = _env_dir(tool)
    explicit = LOCKS / f"{tool}.explicit.txt"
    pip = LOCKS / f"{tool}.pip.txt"
    assert explicit.is_file(), explicit
    assert "@EXPLICIT" in explicit.read_text()
    # fgsea: frozen explicit lock only, no pip layer (M8-17 spec §1); an env without bin/pip
    # (miniprot) gets no pip lock from build_envs.sh lock mode
    if tool != "fgsea" and (env / "bin" / "pip").exists():
        assert pip.is_file(), pip


@pytest.mark.parametrize("tool", CONDA_TOOLS)
def test_build_marker_matches_env_file(tool: str) -> None:
    env = _env_dir(tool)
    marker = env / ".gippa-built"
    assert marker.is_file(), marker
    sha = hashlib.sha256((REPO / TOOLS[tool].env_file).read_bytes()).hexdigest()
    assert sha in marker.read_text()
