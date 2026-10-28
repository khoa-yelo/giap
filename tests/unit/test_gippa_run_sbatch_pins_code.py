"""``gippa_run.sbatch`` runs the code of ``$GIPPA_REPO``, not only its working directory.

The shared venv's editable ``gippa`` points at the main checkout, so a head job submitted with
``GIPPA_REPO=<worktree>`` used to import the main checkout's code and record its commit (pilot A2,
2026-10-25). The script now puts ``$GIPPA_REPO/src`` first on ``PYTHONPATH``.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RUN = REPO / "scripts" / "slurm" / "gippa_run.sbatch"


def _run(tmp_path: Path, pythonpath: str | None) -> list[str]:
    repo = tmp_path / "wt"
    (repo / "scripts").mkdir(parents=True)
    (repo / "scripts" / "farmshare_env.sh").write_text(": # stub\n")
    bin_ = tmp_path / "bin"
    bin_.mkdir()
    uv = bin_ / "uv"
    uv.write_text(
        '#!/bin/bash\necho "cwd=$PWD"\necho "pythonpath=${PYTHONPATH:-}"\necho "args=$*"\n'
    )
    uv.chmod(0o755)
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    env.update(PATH=f"{bin_}:{env['PATH']}", GIPPA_REPO=str(repo))
    if pythonpath is not None:
        env["PYTHONPATH"] = pythonpath
    proc = subprocess.run(
        ["bash", str(RUN), "--outdir", "x"], env=env, capture_output=True, text=True, check=False
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout.splitlines()


def test_gippa_repo_src_comes_first_on_pythonpath(tmp_path: Path) -> None:
    out = _run(tmp_path, None)
    assert f"cwd={tmp_path / 'wt'}" in out
    assert f"pythonpath={tmp_path / 'wt' / 'src'}" in out
    assert "args=run --no-sync gippa run --profile farmshare,apptainer --outdir x" in out


def test_an_existing_pythonpath_is_kept_after_it(tmp_path: Path) -> None:
    out = _run(tmp_path, "/opt/extra")
    assert f"pythonpath={tmp_path / 'wt' / 'src'}:/opt/extra" in out
