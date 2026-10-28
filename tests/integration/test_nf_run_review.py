"""Code review of PR #78, lead decisions 3a and 3b, against real Nextflow (stub mode, no tools).

docs/specs/M2-07-nextflow-tool-processes.md, "Lead decisions for the run-wiring follow-up", item 3:

a. Two ``gippa run ... --resume`` stub runs into the same OUT: the second run's trace has cached
   tasks (preflight must not rewrite unchanged ``work/preflight/*`` files, whose mtime is part of
   Nextflow's task hash).
b. With ``params.db_dir`` set, the evaluated apptainer profile binds it read-only
   (``-B <db_dir>:<db_dir>:ro`` in ``apptainer.runOptions`` or a process ``containerOptions``);
   without it there is no ``null`` bind. A probe script is run with ``-c workflow/nextflow.config``
   and ``-params-file`` (``nextflow config`` cannot take params); no task or container starts.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import stat
import subprocess
from pathlib import Path

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "workflow" / "tests" / "data" / "annotate"
NF_CONFIG = REPO / "workflow" / "nextflow.config"

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

PROBE_NF = """
def collect(m, opts) {
    m.each { k, v ->
        if (v instanceof Map) collect(v, opts)
        else if (k == 'containerOptions') opts << (v instanceof Closure ? v.call() : v).toString()
    }
    return opts
}

workflow {
    def c = workflow.session.config
    def opts = collect(c.process ?: [:], [(c.apptainer?.runOptions ?: '').toString()])
    new File(System.getenv('PROBE_OUT')).text = groovy.json.JsonOutput.toJson(opts)
}
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


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_a_second_resume_run_caches_tasks(
    wrapped_nf: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gippa.cli import cli

    monkeypatch.setenv("GIPPA_DB_DIR", str(_db(tmp_path / "db")))
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(">prot1\nM" + "ACDEFGHIKLMNPQRSTVWY" * 8 + "\n")
    out = tmp_path / "OUT"
    argv = ["run", "--genome", str(DATA / "genome.fa"), "--proteins", str(proteins)]
    argv += ["--outdir", str(out), "--modules", "core", "--resume"]
    for attempt in (1, 2):
        res = CliRunner().invoke(cli, argv)
        log = res.output + (res.stderr if res.stderr_bytes is not None else "")
        if res.exception is not None and not isinstance(res.exception, SystemExit):
            raise res.exception
        assert res.exit_code == ExitCode.SUCCESS, f"run {attempt}: {log}"
        if attempt == 1:
            first = _trace(out)
            assert first and all(r["status"] == "COMPLETED" for r in first)
    rows = _trace(out)
    cached = [r["name"] for r in rows if r["status"] == "CACHED"]
    assert len(cached) > 0, f"--resume cached nothing: {[(r['name'], r['status']) for r in rows]}"
    bakta = [r for r in rows if r["name"].split(" (")[0].split(":")[-1] == "BAKTA"]
    assert bakta and all(r["status"] == "CACHED" for r in bakta), bakta


def _probe(tmp_path: Path, params: dict[str, object]) -> list[str]:
    nf = shutil.which("nextflow")
    assert nf is not None
    run = tmp_path / "probe"
    run.mkdir()
    (run / "probe.nf").write_text(PROBE_NF)
    (run / "params.json").write_text(json.dumps(params))
    out = run / "opts.json"
    env = {**os.environ, "PROBE_OUT": str(out), "NXF_ANSI_LOG": "false"}
    cmd = [nf, "run", "probe.nf", "-c", str(NF_CONFIG), "-params-file", "params.json"]
    cmd += ["-profile", "apptainer", "-work-dir", str(run / "work")]
    proc = subprocess.run(cmd, cwd=run, env=env, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    opts = json.loads(out.read_text())
    assert isinstance(opts, list)
    return [str(o) for o in opts]


def test_b_apptainer_profile_binds_db_dir_read_only(tmp_path: Path) -> None:
    db = tmp_path / "gippadb"
    db.mkdir()
    opts = _probe(tmp_path, {"db_dir": str(db)})
    bind = f"{db}:{db}:ro"
    assert any(bind in o for o in opts), f"no `-B {bind}` in apptainer options: {opts}"


def test_b_no_db_dir_means_no_bind(tmp_path: Path) -> None:
    opts = _probe(tmp_path, {})
    assert not any("null" in o for o in opts), opts
    assert not any(":ro" in o and "ps" not in o for o in opts), opts
