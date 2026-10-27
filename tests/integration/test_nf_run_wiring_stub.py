"""``gippa run --genome G`` reaches ANNOTATE in stub mode (M2-07 follow-up, review item 1 of PR #65).

``gippa run`` is invoked through the CLI; a ``nextflow`` wrapper first on ``PATH`` forwards the
runner's exact command line to the real Nextflow, swapping the profile for ``test,stub`` and adding
``-stub-run`` (no tool, container or database is touched). Before the follow-up, params.json had
no ``inputs.replicons`` and main.nf ran the one-task smoke branch instead.
"""

from __future__ import annotations

import csv
import os
import shutil
import stat
from pathlib import Path

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "workflow" / "tests" / "data" / "annotate"

WRAPPER = """#!/bin/sh
set -e
n=$#; i=0; prev=""
for a in "$@"; do
  i=$((i + 1))
  if [ "$prev" = "-profile" ]; then a="test,stub"; fi
  set -- "$@" "$a"
  prev="$a"
done
shift "$n"
exec "{real}" "$@" -stub-run
"""


@pytest.fixture(autouse=True)
def _fake_gippa_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Every pinned image and prebuilt env present (run-wiring decision 4a); never the real root."""
    root = fake_gippa_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    return root


def _db(path: Path) -> Path:
    path.mkdir()
    write_setup_manifest(path)  # real `gippa setup` format (decision 4c)
    return path


@pytest.fixture
def wrapped_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    real = shutil.which("nextflow")
    assert real is not None
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(WRAPPER.replace("{real}", real))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)


def test_gippa_run_reaches_annotate_in_stub_mode(
    wrapped_nf: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.cli import cli

    monkeypatch.setenv("GIPPA_DB_DIR", str(_db(tmp_path / "db")))
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(">prot1\nM" + "ACDEFGHIKLMNPQRSTVWY" * 8 + "\n")
    out = tmp_path / "OUT"
    argv = ["--genome", str(DATA / "genome.fa"), "--proteins", str(proteins), "--outdir", str(out)]
    res = CliRunner().invoke(cli, ["run", *argv, "--modules", "core"])
    log = res.output + (res.stderr if res.stderr_bytes is not None else "")
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, log
    with (out / "logs" / "trace.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    names = [r["name"] for r in rows]
    procs = {n.split(" (")[0].split(":")[-1] for n in names}
    assert "BAKTA" in procs, f"ANNOTATE did not run: {sorted(names)}"
    assert any("parse_bakta" in n for n in names), sorted(names)
    assert not any("(smoke)" in n for n in names), "main.nf took the smoke branch"
    assert all(r["status"] == "COMPLETED" for r in rows)
