"""Step events reach ``OUT/logs/events.jsonl`` under Nextflow (M3-17 lead decision 15).

docs/specs/M3-17-mapping-steps.md, decision 15 after the code review of PR #96 (finding 1: the
``search_db`` event stayed in ``<task>/logs/events.jsonl``): the per-task ``logs/events.jsonl``
written by ``gippa _step`` is an optional GIPPA_STEP output; the workflow collects them and
``gippa run`` appends them to ``OUT/logs/events.jsonl`` after Nextflow ends, in task completion
order as listed in the trace, each event keeping its step name. The test reads the run directory
only, never a task directory.

``gippa run`` is invoked through the CLI with a ``nextflow`` wrapper that swaps the profile for
``test,stub`` and adds ``-stub-run`` (as tests/integration/test_nf_run_wiring_stub.py). Stub
tasks never execute ``gippa _step``, so the placement stub writes its ``search_db`` event (one
§6.7 JSON line, like the loci stub's §9.3 headers): this checks the collect-and-append path end to
end. Inputs: workflow/tests/data/annotate/genome.fa, one long protein and a peptide table; the
database directory is a tmp copy of the real-format ``gippa setup`` manifest.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest

from gippa.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "workflow" / "tests" / "data" / "annotate"
LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
PY_STEPS = ("sixframe_index", "map_exact", "map_tolerant", "placement", "protein_leftovers", "loci")

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
def _env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    root = fake_gippa_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    real = shutil.which("nextflow")
    assert real is not None
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(WRAPPER.replace("{real}", real))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    monkeypatch.setenv("GIPPA_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)


@pytest.fixture(scope="function")
def run_dir(tmp_path: Path) -> Path:
    from gippa.cli import cli

    proteins = tmp_path / "in" / "proteins.faa"
    proteins.parent.mkdir()
    proteins.write_text(f">P1\n{LONG_PROTEIN}\n")
    peptides = tmp_path / "in" / "peptides.tsv"
    peptides.write_text(
        "peptide\tprotein_ids\trun\tintensity\nPEPTIDEK\tP1\tr1\t100\nAAGLIEVQAPILSR\tP1\tr1\t50\n"
    )
    out = tmp_path / "OUT"
    argv = ["run", "--genome", str(DATA / "genome.fa"), "--proteins", str(proteins)]
    argv += ["--peptides", str(peptides), "--outdir", str(out), "--modules", "core"]
    res = CliRunner().invoke(cli, argv)
    log = res.output + (res.stderr if res.stderr_bytes is not None else "")
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, log
    return out


def _events(out: Path) -> list[dict[str, Any]]:
    path = out / "logs" / "events.jsonl"
    assert path.is_file(), f"no {path}; logs/: {sorted(p.name for p in (out / 'logs').glob('*'))}"
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _is_search_db(ev: dict[str, Any]) -> bool:
    data = ev.get("data") if isinstance(ev.get("data"), dict) else {}
    return "search_db" in str(ev.get("event", "")) or "search_db" in data


def _trace_steps(out: Path) -> list[str]:
    """Py step names (GIPPA_STEP alias tags) in trace order (task completion)."""
    with (out / "logs" / "trace.tsv").open() as fh:
        names = [r["name"] for r in csv.DictReader(fh, delimiter="\t")]
    tags = [n.rsplit("(", 1)[-1].rstrip(") ") for n in names if "(" in n]
    return [t for t in tags if t in PY_STEPS]


def test_run_dir_events_jsonl_has_the_placement_search_db_event(run_dir: Path) -> None:
    evs = _events(run_dir)
    for ev in evs:
        assert set(ev) >= {"ts", "step", "level", "event", "data"}, ev
    placed = [ev for ev in evs if _is_search_db(ev)]
    assert placed, f"no search_db event in OUT/logs/events.jsonl: {evs[-10:]}"
    assert any(ev["step"] == "placement" for ev in placed), f"step name not kept: {placed}"


def test_run_dir_step_events_follow_trace_order(run_dir: Path) -> None:
    order = _trace_steps(run_dir)
    assert "placement" in order, order
    seen = [ev["step"] for ev in _events(run_dir) if ev.get("step") in PY_STEPS]
    assert "placement" in seen, seen
    firsts = list(dict.fromkeys(seen))
    idx = [order.index(s) for s in firsts]
    assert idx == sorted(idx), f"events in step order {firsts}, trace order {order}"
