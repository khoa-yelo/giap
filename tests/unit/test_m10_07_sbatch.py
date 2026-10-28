"""``scripts/slurm/entrapment.sbatch`` run under plain bash (M10-07 criterion 12; rev-1007 B2, S5).

Binding: docs/specs/M10-07-free-entrapment.md, "Lead decisions after the code review (rev-1007)":

* B2: the default output path (``$GIPPA_ROOT/pilot/metrics``, with ``GIPPA_ROOT`` defaulting to
  ``$GIPPA_REPO/local``) must be git-ignored (``git check-ignore``), and the default run must not
  exit 2;
* S5: ``OUT_A1`` defaults to the M10-04 A1 directory (M10-05's input table:
  ``$GIPPA_ROOT/runs/pilot/pass1/out_a1``); an explicitly set ``OUT_A1``, ``DIANN_CMD`` or
  ``PROVENANCE`` that does not exist exits 2 before anything runs; the sbatch never skips the
  cross-checks silently (a skipped cross-check is logged).

The script sources ``scripts/farmshare_env.sh``, which sets ``UV_PROJECT_ENVIRONMENT`` to
``$GIPPA_ROOT/venv``. The tests point ``GIPPA_ROOT`` at a scratch directory and put a stub
``venv/bin/python`` there: it records its argv (one token per line) and exits 0, or (B2) execs
this interpreter so the real ``entrapment.py`` runs. ``#SBATCH`` lines are bash comments, so
nothing is submitted to Slurm.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from typing import Any

import pytest
from helpers.entrapment_m10_07 import (
    CONTAMINANTS,
    EXAMPLE2_ROWS,
    K12,
    NISSLE,
    REPO,
    coding_text,
    write_a1,
    write_faa,
    write_genome,
    write_runs_manifest,
)

SBATCH = REPO / "scripts" / "slurm" / "entrapment.sbatch"
OPTIONAL_VARS = ("OUT_A1", "DIANN_CMD", "PROVENANCE")
INPUT_VARS = ("PEP_A1", "K12_FASTA", "CONTAMINANTS", "GENOME", "STRAIN_DB", "RUNS_MANIFEST")


def _inputs(root: Path) -> dict[str, str]:
    """The six required inputs (the spec's example-2 fixture), outside the checkout."""
    return {
        "PEP_A1": str(write_a1(root / "a1.tsv.gz", list(EXAMPLE2_ROWS))),
        "K12_FASTA": str(write_faa(root / "UP000000625.fasta", K12)),
        "CONTAMINANTS": str(write_faa(root / "contaminants.fasta", CONTAMINANTS, "Cont_P")),
        "GENOME": str(
            write_genome(root / "nissle.fna", {"nissle": coding_text(*NISSLE)}, ["nissle"])
        ),
        "STRAIN_DB": str(write_faa(root / "strain_db.faa", ["MSAGLLKVEPTLDEKAAGGYWR"], "N_")),
        "RUNS_MANIFEST": str(write_runs_manifest(root / "pxd078459_files.tsv")),
    }


def _stub_python(gippa_root: Path, *, real: bool) -> Path:
    """``$GIPPA_ROOT/venv/bin/python``: records argv to ``$GIPPA_ROOT/argv.txt``; ``real`` execs
    this interpreter afterwards."""
    py = gippa_root / "venv" / "bin" / "python"
    py.parent.mkdir(parents=True, exist_ok=True)
    tail = f'exec "{sys.executable}" "$@"' if real else "exit 0"
    py.write_text(
        f'#!/bin/bash\nprintf "%s\\n" "$@" > "{gippa_root / "argv.txt"}"\n{tail}\n',
        encoding="utf-8",
    )
    py.chmod(0o755)
    return gippa_root / "argv.txt"


def _run(gippa_root: Path, **env_over: str) -> subprocess.CompletedProcess[str]:
    env = {
        k: v
        for k, v in os.environ.items()
        if k not in (*OPTIONAL_VARS, *INPUT_VARS, "METRICS_DIR", "GIPPA_ROOT", "GIPPA_REPO")
    }
    env |= {"GIPPA_REPO": str(REPO), "GIPPA_ROOT": str(gippa_root), **env_over}
    return subprocess.run(
        ["bash", str(SBATCH)], capture_output=True, text=True, env=env, check=False, timeout=300
    )


def _argv(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines() if path.is_file() else []


def _flag(argv: list[str], flag: str) -> str | None:
    return argv[argv.index(flag) + 1] if flag in argv else None


def _show(res: Any) -> str:
    return f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


@pytest.fixture
def root(tmp_path: Path) -> Path:
    return tmp_path / "gippa_root"


def test_s5_out_a1_defaults_to_the_m10_04_a1_directory(root: Path, tmp_path: Path) -> None:
    out_a1 = root / "runs" / "pilot" / "pass1" / "out_a1"
    out_a1.mkdir(parents=True)
    argv_file = _stub_python(root, real=False)
    res = _run(root, **_inputs(tmp_path / "in"))
    assert res.returncode == 0, _show(res)
    argv = _argv(argv_file)
    assert argv, f"entrapment.py was never started\n{_show(res)}"
    assert _flag(argv, "--out-a1") == str(out_a1), argv
    assert _flag(argv, "--out") == str(root / "pilot" / "metrics"), argv


@pytest.mark.parametrize("var", OPTIONAL_VARS)
def test_s5_explicit_missing_optional_path_exits_2_before_running(
    root: Path, tmp_path: Path, var: str
) -> None:
    argv_file = _stub_python(root, real=False)
    res = _run(root, **_inputs(tmp_path / "in"), **{var: str(tmp_path / "no-such-path")})
    assert res.returncode == 2, _show(res)
    assert not _argv(argv_file), "entrapment.py ran although an explicit input is missing"
    assert "no-such-path" in res.stdout + res.stderr


def test_s5_a_skipped_cross_check_is_logged(root: Path, tmp_path: Path) -> None:
    """No ``OUT_A1`` given and the default A1 directory absent: the run goes on without
    ``--out-a1`` but says that cross-checks 1 and 2 were skipped."""
    argv_file = _stub_python(root, real=False)
    res = _run(root, **_inputs(tmp_path / "in"))
    assert res.returncode == 0, _show(res)
    assert "--out-a1" not in _argv(argv_file)
    log = (res.stdout + res.stderr).lower()
    assert "cross-check" in log and "skip" in log, _show(res)


def test_b2_default_metrics_dir_inside_the_checkout_runs(tmp_path: Path) -> None:
    """``GIPPA_ROOT`` under the checkout's git-ignored ``local/`` (farmshare_env.sh's default
    layout) and ``METRICS_DIR`` left at its default: the real entrapment.py exits 0."""
    base = REPO / "local"
    made_base = not base.exists()
    root = base / f"_m10_07_sbatch_{uuid.uuid4().hex[:8]}"
    try:
        _stub_python(root, real=True)
        res = _run(root, **_inputs(tmp_path / "in"))
        assert res.returncode == 0, _show(res)
        assert (root / "pilot" / "metrics" / "entrapment.json").is_file(), _show(res)
    finally:
        shutil.rmtree(root, ignore_errors=True)
        if made_base:
            shutil.rmtree(base, ignore_errors=True)


def test_b2_the_ignored_check_is_git_check_ignore() -> None:
    """Lead decision B2 names the mechanism: ``git check-ignore`` (sbatch or entrapment.py)."""
    texts = [
        SBATCH.read_text(encoding="utf-8"),
        (REPO / "validation" / "pilot" / "entrapment.py").read_text(encoding="utf-8"),
    ]
    assert any("check-ignore" in t for t in texts)
