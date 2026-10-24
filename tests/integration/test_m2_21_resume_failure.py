"""Resume and failure paths of ``giap run`` on sim100k (M2-21 examples 5 and 6; GIAP_SPEC.md §11.2
"resume and failure paths", §3.8 "Resume", §6.3 exit 4).

* 5: ``giap run`` on sim100k, then kill one Promotech task and ``-resume`` → completes; finished
  tasks cached (trace shows ``CACHED``).
* 6: a tool exits 1 (simulated with a broken env var) → ``giap run`` exits 4 naming the step and
  its ``.command.log``.

``giap run`` is invoked through the CLI. A ``nextflow`` wrapper first on ``PATH`` forwards the
runner's exact command line to the real Nextflow and appends ``-c $M2_21_EXTRA_CONFIG`` when set
(the injected fault); with ``M2_21_STUB=1`` it also swaps the profile for ``test,stub`` and adds
``-stub-run`` (no tool, container or database is touched, as
tests/integration/test_nf_run_wiring_stub.py).

Faults (settled by the test author):

* "kill one Promotech task": ``beforeScript = 'exit 137'`` on ``PROMOTECH`` (the status of a
  SIGKILLed task; base.config retries 137 twice, every retry is killed too, the run stops). The
  resumed run has no fault.
* "a tool exits 1 (broken env var)": ``PROMOTECH_HOME`` pointed at a missing directory. In stub
  mode the stub never reads it, so a ``beforeScript`` checks it and exits 1; the real variant
  only sets the env var and lets Promotech fail.

Stub variants run in the integration layer (``needs_nextflow``); real-tool variants also need
``needs_conda`` and ``needs_promotech_models`` and are ``slow`` (nightly, array task 0).
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from giap.exitcodes import ExitCode

pytestmark = pytest.mark.needs_nextflow

WRAPPER = """#!/bin/sh
set -e
n=$#; i=0; prev=""
for a in "$@"; do
  i=$((i + 1))
  if [ "$prev" = "-profile" ] && [ "${M2_21_STUB:-0}" = 1 ]; then a="test,stub"; fi
  set -- "$@" "$a"
  prev="$a"
done
shift "$n"
if [ "${M2_21_STUB:-0}" = 1 ]; then set -- "$@" -stub-run; fi
if [ -n "${M2_21_EXTRA_CONFIG:-}" ]; then set -- "$@" -c "$M2_21_EXTRA_CONFIG"; fi
exec "{real}" "$@"
"""

KILL = "process { withName: 'PROMOTECH' { beforeScript = 'exit 137' } }\n"
BROKEN_ENV_STUB = (
    "env { PROMOTECH_HOME = '/nonexistent/m2-21-broken' }\n"
    "process { withName: 'PROMOTECH' { beforeScript = "
    '\'[ -d "$PROMOTECH_HOME" ] || { echo "PROMOTECH_HOME=$PROMOTECH_HOME missing" >&2; '
    "exit 1; }' "
    "} }\n"
)
BROKEN_ENV_REAL = "env { PROMOTECH_HOME = '/nonexistent/m2-21-broken' }\n"
LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim

    return giapsim.simulate(
        giapsim.SIM100K, giapsim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim"
    )


def _wrap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, stub: bool) -> None:
    real = shutil.which("nextflow")
    assert real is not None
    bindir = tmp_path / "bin"
    bindir.mkdir(exist_ok=True)
    exe = bindir / "nextflow"
    exe.write_text(WRAPPER.replace("{real}", real))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    monkeypatch.setenv("M2_21_STUB", "1" if stub else "0")
    monkeypatch.delenv("M2_21_EXTRA_CONFIG", raising=False)


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
    _wrap(tmp_path, monkeypatch, stub=True)


@pytest.fixture
def real_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Real ``$GIAP_ROOT`` / ``$GIAP_DB_DIR`` (prebuilt envs, RF-HOT model); farmshare_local."""
    _wrap(tmp_path, monkeypatch, stub=False)


def _fault(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, text: str | None) -> None:
    if text is None:
        monkeypatch.delenv("M2_21_EXTRA_CONFIG", raising=False)
        return
    cfg = tmp_path / "fault.config"
    cfg.write_text(text)
    monkeypatch.setenv("M2_21_EXTRA_CONFIG", str(cfg))


def _argv(sim: Any, tmp_path: Path, out: Path, real: bool) -> list[str]:
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(f">P1\n{LONG_PROTEIN}\n")
    argv = ["run", "--genome", str(sim.genome_fasta), "--proteins", str(proteins)]
    argv += ["--outdir", str(out), "--modules", "core", "--resume"]
    if real:  # every layer precomputed except promoters, so Promotech is the only heavy tool
        argv += ["--config", str(sim.config_yaml), "--set", "annotation.precomputed.promoters=null"]
        argv += ["--profile", "farmshare_local"]
    return argv


def _invoke(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _trace(out: Path) -> list[dict[str, str]]:
    with (out / "logs" / "trace.tsv").open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _proc(name: str) -> str:
    return name.split(" (")[0].split(":")[-1]


def _example_5(sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, real: bool) -> None:
    out = tmp_path / "OUT"
    argv = _argv(sim, tmp_path, out, real)
    _fault(tmp_path, monkeypatch, KILL)
    code, log = _invoke(argv)
    assert code == ExitCode.STEP_FAILED, f"killed Promotech must fail the run: {log}"
    first = _trace(out)
    assert any(_proc(r["name"]) == "PROMOTECH" and r["status"] == "FAILED" for r in first), first
    done = {r["name"] for r in first if r["status"] == "COMPLETED"}
    assert done, "no task finished before Promotech was killed"

    _fault(tmp_path, monkeypatch, None)
    code, log = _invoke(argv)
    assert code == ExitCode.SUCCESS, f"resumed run: {log}"
    second = _trace(out)
    status = {r["name"]: r["status"] for r in second}
    not_cached = sorted(n for n in done if status.get(n) != "CACHED")
    assert not not_cached, f"finished tasks re-ran instead of CACHED: {not_cached}"
    promotech = [r for r in second if _proc(r["name"]) == "PROMOTECH"]
    assert promotech and all(r["status"] == "COMPLETED" for r in promotech), promotech
    assert all(r["status"] in ("COMPLETED", "CACHED") for r in second), second


def _example_6(sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, real: bool) -> None:
    out = tmp_path / "OUT"
    _fault(tmp_path, monkeypatch, BROKEN_ENV_REAL if real else BROKEN_ENV_STUB)
    code, log = _invoke(_argv(sim, tmp_path, out, real))
    assert code == ExitCode.STEP_FAILED, log
    assert re.search(r"\bstep promotech failed\b", log, re.I), log
    m = re.search(r"(\S+/\.command\.log)\b", log)
    assert m, f"no .command.log path in the message: {log}"
    cmd_log = Path(m.group(1))
    assert cmd_log.is_file(), cmd_log
    task = cmd_log.parent
    assert (task / ".exitcode").read_text().strip() != "0"
    failed = [r for r in _trace(out) if r["status"] == "FAILED"]
    assert failed and all(_proc(r["name"]) == "PROMOTECH" for r in failed), failed
    if real:
        text = cmd_log.read_text(errors="replace") + (task / ".command.err").read_text(
            errors="replace"
        )
        assert "/nonexistent/m2-21-broken" in text, text[-2000:]


def test_example_5_killed_promotech_then_resume_stub(
    sim: Any, stub_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _example_5(sim, tmp_path, monkeypatch, real=False)


def test_example_6_broken_env_var_exits_4_naming_step_and_log_stub(
    sim: Any, stub_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _example_6(sim, tmp_path, monkeypatch, real=False)


@pytest.mark.slow
@pytest.mark.needs_conda
@pytest.mark.needs_promotech_models
def test_example_5_killed_promotech_then_resume_real(
    sim: Any, real_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _example_5(sim, tmp_path, monkeypatch, real=True)


@pytest.mark.slow
@pytest.mark.needs_conda
@pytest.mark.needs_promotech_models
def test_example_6_broken_env_var_exits_4_naming_step_and_log_real(
    sim: Any, real_env: None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _example_6(sim, tmp_path, monkeypatch, real=True)
