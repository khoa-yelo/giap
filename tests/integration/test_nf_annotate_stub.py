"""ANNOTATE in stub mode through ``main.nf`` and the runner (M2-07 edge cases; GIPPA_SPEC.md §6.3).

- every ``withName`` pin of an annotation tool matches a process (no "no process matching config
  selector" warning for it);
- a failing optional-module process fails the run with exit 4 naming the step, unless the module
  is off (no silent skip). The failure is injected with a ``beforeScript = 'exit 3'`` override for
  GENOMAD / ANTISMASH; nothing real runs, local executor only.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
from pathlib import Path

import pytest

from gippa.exitcodes import ExitCode
from gippa.runner.nextflow import launch_command, run_nextflow

pytestmark = pytest.mark.needs_nextflow

REPO = Path(__file__).resolve().parents[2]
MAIN = REPO / "workflow" / "main.nf"
DATA = REPO / "workflow" / "tests" / "data" / "annotate"
PINNED = ["BAKTA", "TRANSTERMHP", "GENOMAD", "ANTISMASH"]


def _params(tmp: Path, out: Path, modules: list[str]) -> Path:
    data = {
        "outdir": str(out),
        "modules": modules,
        "inputs": {"genome": str(DATA / "genome.fa"), "replicons": str(DATA / "replicons.tsv")},
        "annotation": json.loads((DATA / "annotation.json").read_text()),
    }
    path = tmp / "params.json"
    path.write_text(json.dumps(data))
    return path


def _env() -> dict[str, str]:
    drop = {"SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"}
    env = {k: v for k, v in os.environ.items() if k not in drop}
    env["NXF_ANSI_LOG"] = "false"
    return env


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _proc(name: str) -> str:
    return name.split(" (")[0].split(":")[-1]


def test_stub_run_has_every_pinned_process_and_no_selector_warning(tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    params = _params(tmp_path, out, ["core", "mobile", "bgc", "function"])
    cmd = [*launch_command(MAIN, params, "test,stub", out, resume=False), "-stub-run"]
    out.mkdir()
    res = subprocess.run(cmd, cwd=out, capture_output=True, text=True, env=_env(), timeout=1200)
    log = res.stdout + res.stderr
    assert res.returncode == 0, log
    for name in PINNED:
        assert f"config selector: {name}\n" not in log + "\n", f"withName {name} matches no process"
    ran = {_proc(r["name"]) for r in _trace(out)}
    assert set(PINNED) <= ran, sorted(ran)
    assert all(r["status"] == "COMPLETED" for r in _trace(out))


def _failing(tmp: Path, proc: str) -> Path:
    cfg = tmp / "fail.config"
    cfg.write_text(f"process {{ withName: '{proc}' {{ beforeScript = 'exit 3' }} }}\n")
    return cfg


@pytest.mark.parametrize(("proc", "module"), [("GENOMAD", "mobile"), ("ANTISMASH", "bgc")])
def test_failing_optional_module_fails_the_run_with_exit_4(
    tmp_path: Path,
    proc: str,
    module: str,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    out = tmp_path / "OUT"
    out.mkdir()
    params = _params(tmp_path, out, ["core", module])
    cmd = launch_command(MAIN, params, "test,stub", out, resume=False)
    cmd += ["-stub-run", "-c", str(_failing(tmp_path, proc))]
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    code = run_nextflow(cmd, out)
    err = capsys.readouterr().err
    assert code == ExitCode.STEP_FAILED, err
    assert f"step {proc.lower()} failed" in err, err


@pytest.mark.parametrize("proc", ["GENOMAD", "ANTISMASH"])
def test_optional_module_off_is_not_run(tmp_path: Path, proc: str) -> None:
    out = tmp_path / "OUT"
    out.mkdir()
    params = _params(tmp_path, out, ["core"])
    cmd = launch_command(MAIN, params, "test,stub", out, resume=False)
    cmd += ["-stub-run", "-c", str(_failing(tmp_path, proc))]
    assert run_nextflow(cmd, out) == ExitCode.SUCCESS
    ran = {_proc(r["name"]) for r in _trace(out)}
    assert proc not in ran
    assert "BAKTA" in ran, f"ANNOTATE did not run: {sorted(ran)}"
