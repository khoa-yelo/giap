"""#354 review S1-S2 (M11-09b): ``build_envs.sbatch`` exits with the env-build status, and the
three Slurm ``grep -P`` readers find a pin in both the bare and the legacy prefixed form.

S1 runs ``build_envs.sbatch`` with bash in a throw-away ``GIPPA_REPO`` whose stub
``scripts/build_envs.sh`` exits 3. S2 runs each script's own ``grep -oP`` line, copied verbatim
from the script, through bash on a one-line config, so a mandatory ``docker://`` group in any
reader fails here.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SLURM = REPO / "scripts" / "slurm"
DIGEST = "b" * 64
BAKTA = f"oschwengers/bakta@sha256:{DIGEST}"

FAKE_APPTAINER = r"""#!/bin/bash
if [ "$1" = pull ]; then printf 'fake sif\n' > "$2"; fi
exit 0
"""


def _build_envs(tmp: Path, config: str, env_status: int) -> subprocess.CompletedProcess[str]:
    repo, root = tmp / "repo", tmp / "root"
    (repo / "scripts" / "slurm").mkdir(parents=True)
    (repo / "workflow" / "conf").mkdir(parents=True)
    shutil.copy(REPO / "scripts" / "farmshare_env.sh", repo / "scripts")
    shutil.copy(SLURM / "build_envs.sbatch", repo / "scripts" / "slurm")
    (repo / "scripts" / "build_envs.sh").write_text(f"#!/bin/bash\nexit {env_status}\n")
    (repo / "workflow" / "conf" / "containers.config").write_text(config)
    (root / "bin").mkdir(parents=True)
    fake = root / "bin" / "apptainer"
    fake.write_text(FAKE_APPTAINER)
    fake.chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG"}}
    env.update(GIPPA_REPO=str(repo), GIPPA_ROOT=str(root), HOME=str(root / "home"))
    return subprocess.run(
        ["bash", str(repo / "scripts" / "slurm" / "build_envs.sbatch")],
        capture_output=True, text=True, cwd=repo, env=env, timeout=120, check=False,
    )  # fmt: skip


def _config(prefix: str) -> str:
    return f"process {{\n  withName: 'BAKTA' {{ container = '{prefix}{BAKTA}' }}\n}}\n"


# ------------------------------------------------------------------ S1


@pytest.mark.parametrize("config", ["process {\n}\n", _config(""), _config("docker://")])
@pytest.mark.parametrize("status", [0, 3])
def test_build_envs_exits_with_the_env_build_status(
    tmp_path: Path, config: str, status: int
) -> None:
    proc = _build_envs(tmp_path, config, status)
    assert proc.returncode == status, proc.stdout + proc.stderr


# ------------------------------------------------------------------ S2


def _grep_line(script: str, pattern: str) -> str:
    """The quoted ``grep -oP "…"`` pattern of ``script`` whose text contains ``pattern``."""
    text = (SLURM / script).read_text(encoding="utf-8")
    found = [m.group(1) for m in re.finditer(r'grep -oP ("[^"]*")', text) if pattern in m.group(1)]
    assert len(found) == 1, (script, found)
    return found[0]


def _run_grep(tmp: Path, quoted: str, config: str, arg: str = "BAKTA") -> str:
    conf = tmp / "containers.config"
    conf.write_text(config)
    script = f'f() {{ grep -oP {quoted} "$CONF"; }}; f {arg}'
    proc = subprocess.run(
        ["bash", "-c", script], capture_output=True, text=True, check=False,
        env={**os.environ, "CONF": str(conf)},
    )  # fmt: skip
    assert proc.returncode == 0, (quoted, config, proc.stderr)
    return proc.stdout.strip()


@pytest.mark.parametrize("script", ["fetch_bakta_db.sbatch", "record_bakta_full.sbatch"])
@pytest.mark.parametrize("prefix", ["", "docker://"])
def test_bakta_digest_readers_accept_both_forms(tmp_path: Path, script: str, prefix: str) -> None:
    quoted = _grep_line(script, "withName: 'BAKTA'")
    assert _run_grep(tmp_path, quoted, _config(prefix)) == DIGEST


@pytest.mark.parametrize("prefix", ["", "docker://"])
def test_record_tools_image_reader_returns_the_bare_ref(tmp_path: Path, prefix: str) -> None:
    quoted = _grep_line("record_tools.sbatch", "container = '(docker://)?")
    assert _run_grep(tmp_path, quoted, _config(prefix)) == BAKTA
