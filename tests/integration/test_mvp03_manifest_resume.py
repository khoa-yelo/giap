"""MVP-03a examples 11 and 12: the manifest after the M2-21 Promotech fault and after ``--resume``
(docs/specs/MVP-03-manifest-finalize.md "Examples" 11-12, edge case "-resume", D-B step 1, lead
decision 5; GIAP_SPEC.md §3.8 "Resume", §6.3 exit 4, §6.7).

Stub variants of tests/integration/test_m2_21_resume_failure.py (its wrapper and fault are copied:
test modules are not importable). A ``nextflow`` wrapper first on ``PATH`` swaps the profile for
``test,stub``, adds ``-stub-run`` and appends ``-c $M2_21_EXTRA_CONFIG`` when set. The fault is
``beforeScript = 'exit 137'`` on ``PROMOTECH`` (base.config retries 137 twice; every retry is killed
too). ``$GIAP_ROOT`` and ``$GIAP_DB_DIR`` are fakes from ``helpers.fake_env``.

* 11: exit 4; manifest written and valid; ``status failed``, ``exit_code 4``, ``failed_step
  "promotech"``, ``failed_log`` ends in ``.command.log``; ``steps[promotech].status ==
  "failed"``; no ``finalize`` step entry.
* 12: the same arguments plus ``--resume`` without the fault → exit 0; ``status completed``,
  ``resume true``; every step whose tasks all completed in run 1 is ``cached``;
  ``stable_view`` equals that of a fresh run (new OUT, no fault, no resume); ``logs/events.jsonl``
  holds no duplicate line and one preflight ``environment.resources`` event, and run 1's file is
  ``logs/events.previous.jsonl``.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import stat
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.manifest_data import validate

from giap.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

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
if [ -n "${M2_21_EXTRA_CONFIG:-}" ]; then set -- "$@" -c "$M2_21_EXTRA_CONFIG"; fi
exec "{real}" "$@"
"""
KILL = "process { withName: 'PROMOTECH' { beforeScript = 'exit 137' } }\n"
LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim

    return giapsim.simulate(
        giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim"
    )


@pytest.fixture
def stub_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from helpers.fake_env import fake_giap_root, write_setup_manifest

    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "giap_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)
    real = shutil.which("nextflow")
    assert real is not None
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(WRAPPER.replace("{real}", real))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    monkeypatch.delenv("M2_21_EXTRA_CONFIG", raising=False)


def _fault(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, on: bool) -> None:
    if not on:
        monkeypatch.delenv("M2_21_EXTRA_CONFIG", raising=False)
        return
    cfg = tmp_path / "fault.config"
    cfg.write_text(KILL)
    monkeypatch.setenv("M2_21_EXTRA_CONFIG", str(cfg))


def _argv(sim: Any, tmp_path: Path, out: Path, *extra: str) -> list[str]:
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(f">P1\n{LONG_PROTEIN}\n")
    argv = ["run", "--genome", str(sim.genome_fasta), "--proteins", str(proteins)]
    return [*argv, "--outdir", str(out), "--modules", "core", *extra]


def _invoke(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}; OUT has {sorted(p.name for p in out.glob('*'))}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _step(name: str) -> str:
    from giap.runner.nextflow import _step_name

    return _step_name(name)


def _lines(path: Path) -> list[str]:
    return [ln for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def _example_11(sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    out = tmp_path / "OUT"
    _fault(tmp_path, monkeypatch, True)
    code, log = _invoke(_argv(sim, tmp_path, out))
    assert code == ExitCode.STEP_FAILED, log
    return out


def test_example_11_failed_run_manifest(
    sim: Any, stub_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = _example_11(sim, tmp_path, monkeypatch)
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["failed_step"]) == ("failed", 4, "promotech")
    assert str(m["failed_log"]).endswith(".command.log")
    by = {s["name"]: s for s in m["steps"]}
    assert by["promotech"]["status"] == "failed"
    assert "finalize" not in by
    assert m["executor"] == "nextflow"


def test_example_12_resume_after_the_fault(
    sim: Any, stub_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.manifest import stable_view

    out = _example_11(sim, tmp_path, monkeypatch)
    first = _trace(out)
    run1_events = _lines(out / "logs" / "events.jsonl")
    tasks: dict[str, set[str]] = {}
    for r in first:
        tasks.setdefault(_step(r["name"]), set()).add(r["status"])
    done = sorted(s for s, st in tasks.items() if st == {"COMPLETED"} and s != "promotech")
    assert done, f"no step finished before the fault: {tasks}"

    _fault(tmp_path, monkeypatch, False)
    code, log = _invoke(_argv(sim, tmp_path, out, "--resume"))
    assert code == ExitCode.SUCCESS, log
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["resume"]) == ("completed", 0, True)
    by = {s["name"]: s for s in m["steps"]}
    not_cached = [s for s in done if by.get(s, {}).get("status") != "cached"]
    assert not not_cached, f"finished in run 1 but not cached: {not_cached}"

    # events: rotated at the start of run 2, every line once
    cur = _lines(out / "logs" / "events.jsonl")
    dup = [ln for ln, n in Counter(cur).items() if n > 1]
    assert not dup, f"duplicate events.jsonl lines after -resume: {dup[:3]}"
    res = [json.loads(ln) for ln in cur]
    pre = [e for e in res if e.get("step") == "preflight" and e.get("event") ==
           "environment.resources"]  # fmt: skip
    assert len(pre) == 1, pre
    assert _lines(out / "logs" / "events.previous.jsonl") == run1_events

    # the resumed manifest's stable view equals a fresh run's
    fresh = tmp_path / "FRESH" / "OUT"
    code, log = _invoke(_argv(sim, tmp_path, fresh))
    assert code == ExitCode.SUCCESS, log
    mf = _manifest(fresh)
    va, vb = stable_view(m), stable_view(mf)
    differ = sorted(k for k in set(va) | set(vb) if va.get(k) != vb.get(k))
    assert not differ, f"resumed and fresh manifests differ on {differ}"
