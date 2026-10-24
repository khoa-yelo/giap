"""Resume and failure paths of VARIANTS and NOVEL under Nextflow (docs/specs/M7-15-integration-
hardening.md item 3, examples 2 and 3; GIAP_SPEC.md §3.8 "Resume", §6.3 exit 4, §11.2).

``giap run`` through the CLI on the 50-kb M7 fixture's reference genome and de novo peptides,
in stub mode: the M2-21 ``nextflow`` wrapper (tests/integration/test_m2_21_resume_failure.py)
swaps the profile for ``test,stub`` and adds ``-stub-run``; ``M7_15_EXTRA_CONFIG`` appends an
injected fault with ``-c``. No ``--config`` is given, so no layer is precomputed and NOVEL runs
its tool tasks (PROMOTECH_NOVEL, OSTIR_NOVEL; ``params.stub.novel_windows_chunks`` defaults to 1).

* example 2: a completed run, then ``--resume`` with no change: every VARIANTS and NOVEL task is
  ``CACHED`` (and there is at least one of each). Observed before M7-15 (Slurm 1869582): the
  unchanged resume re-executes 8 of 39 tasks -- GENETABLE:GENE_TABLE, DETECTABILITY,
  CONTEXT_STATS, SENSITIVITY, OPERONS_STEP, NOVEL:NOVEL_ORFS, ... -- so some GENETABLE input
  changes between two identical ``giap run`` calls;
* item 3: ``--resume`` with ``--set variants.high_min_peptides=3`` re-runs exactly the
  ``variants`` step task; every other task is ``CACHED``. (Today the resolved config is read by
  path, not hashed, so a config change re-runs nothing: the step's config must enter its task
  hash, e.g. as the slice of the config the step reads.)
* example 3: OSTIR_NOVEL forced to exit 1 (``beforeScript = 'exit 1'``) -> ``giap run`` exits 4,
  its message names ``ostir_novel`` (any case) and an existing ``.command.log`` of the failed
  task (the M2-21 example 6 contract of ``giap.runner.nextflow``), and the trace's only FAILED
  task is OSTIR_NOVEL.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import stat
import sys
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m7_chain as m7

pytestmark = [pytest.mark.needs_nextflow, pytest.mark.slow]

SEED = 13
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
set -- "$@" -stub-run
if [ -n "${M7_15_EXTRA_CONFIG:-}" ]; then set -- "$@" -c "$M7_15_EXTRA_CONFIG"; fi
exec "{real}" "$@"
"""
OSTIR_FAIL = "process { withName: '.*OSTIR_NOVEL' { beforeScript = 'exit 1' } }\n"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    return m7.simulate(tmp_path_factory.mktemp("m7_15_nf") / "sim", SEED)


@pytest.fixture
def stub_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from helpers.fake_env import fake_giap_root, write_setup_manifest

    real = shutil.which("nextflow")
    if real is None:
        pytest.skip("nextflow not on PATH")
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "giap_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(WRAPPER.replace("{real}", real))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    monkeypatch.delenv("M7_15_EXTRA_CONFIG", raising=False)


def _argv(sim: Any, out: Path, *extra: str) -> list[str]:
    return ["run", "--genome", str(sim.reference_fasta), "--peptides", str(m7.denovo(sim)),
            "--outdir", str(out), "--modules", "core", "--seed", str(SEED), "--resume",
            *extra]  # fmt: skip


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _path(name: str) -> list[str]:
    """Workflow path of a trace name: ``NOVEL:OSTIR_NOVEL (novel)`` -> [NOVEL, OSTIR_NOVEL]."""
    return name.split(" (")[0].split(":")


def _tag(name: str) -> str:
    m = re.search(r"\((.*)\)\s*$", name)
    return m.group(1) if m else ""


def _in(sub: str, row: dict[str, str]) -> bool:
    return sub in _path(row["name"])[:-1]


def _first_run(sim: Any, out: Path) -> list[dict[str, str]]:
    code, log = m7.invoke(_argv(sim, out))
    assert code == 0, log[-4000:]
    rows = _trace(out)
    for sub in ("VARIANTS", "NOVEL"):
        assert any(_in(sub, r) for r in rows), f"no {sub} task in the stub run"
    assert any(_path(r["name"])[-1] == "OSTIR_NOVEL" for r in rows), "OSTIR_NOVEL did not run"
    return rows


def test_example_2_resume_without_change_caches_variants_and_novel(
    sim: Any, stub_env: None, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    _first_run(sim, out)
    code, log = m7.invoke(_argv(sim, out))
    assert code == 0, log[-4000:]
    rows = [r for r in _trace(out) if _in("VARIANTS", r) or _in("NOVEL", r)]
    assert rows
    ran = sorted(r["name"] for r in rows if r["status"] != "CACHED")
    assert not ran, f"VARIANTS / NOVEL tasks re-executed on -resume: {ran}"


def test_item_3_changed_high_min_peptides_reruns_only_variants(
    sim: Any, stub_env: None, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    _first_run(sim, out)
    code, log = m7.invoke(_argv(sim, out, "--set", "variants.high_min_peptides=3"))
    assert code == 0, log[-4000:]
    rows = _trace(out)
    ran = [r for r in rows if r["status"] != "CACHED"]
    names = sorted(r["name"] for r in ran)
    assert ran, "nothing re-ran after variants.high_min_peptides changed"
    assert all(_in("VARIANTS", r) and _tag(r["name"]) == "variants" for r in ran), names
    assert all(r["status"] == "COMPLETED" for r in ran), names


def test_example_3_failed_ostir_novel_exits_4_naming_it_and_its_log(
    sim: Any, stub_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.exitcodes import ExitCode

    fault = tmp_path / "fault.config"
    fault.write_text(OSTIR_FAIL)
    monkeypatch.setenv("M7_15_EXTRA_CONFIG", str(fault))
    out = tmp_path / "OUT"
    code, log = m7.invoke(_argv(sim, out))
    assert code == ExitCode.STEP_FAILED, log[-4000:]
    assert re.search(r"\bostir_novel\b", log, re.I), log[-4000:]
    m = re.search(r"(\S+/\.command\.log)\b", log)
    assert m, f"no .command.log path in the message: {log[-4000:]}"
    assert Path(m.group(1)).is_file(), m.group(1)
    failed = [r for r in _trace(out) if r["status"] == "FAILED"]
    assert failed and all(_path(r["name"])[-1] == "OSTIR_NOVEL" for r in failed), failed
