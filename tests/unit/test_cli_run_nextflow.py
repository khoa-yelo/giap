"""``giap run`` launches the Nextflow runner (M2-07 owns the wiring; M2-06 lead decision 7).

GIAP_SPEC.md §6.3 "Launch": ``giap run`` validates inputs and config in Python, writes
``OUT/work/params.json`` and ``OUT/resolved_config.yaml``, then runs
``nextflow run <pkg>/workflow/main.nf -params-file … -profile … -work-dir OUT/work -with-trace …
[-resume]``; a failed step is exit 4 naming the step (§6.3 "Error strategy"). A missing tool is
exit 2 (§5.1).

A fake ``nextflow`` executable on ``PATH`` records its argv and working directory and exits with
``$FAKE_NF_RC``; with ``$FAKE_NF_TRACE`` it copies that file to the ``-with-trace`` path first.
Nothing real is launched.
"""

from __future__ import annotations

import json
import os
import random
import re
import stat
from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from giap.exitcodes import ExitCode
from giap.runner.nextflow import launch_command

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
pwd > "$FAKE_NF_ARGS.cwd"
if [ -n "$FAKE_NF_TRACE" ]; then
  prev=""
  for a in "$@"; do
    if [ "$prev" = "-with-trace" ]; then mkdir -p "$(dirname "$a")"; cp "$FAKE_NF_TRACE" "$a"; fi
    prev="$a"
  done
fi
exit "${FAKE_NF_RC:-0}"
"""

TRACE_HEADER = "task_id\thash\tnative_id\tname\ttag\tstatus\texit\tattempt\tworkdir\n"


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    rng = random.Random(3)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    path = tmp_path / "in" / "genome.fa"
    path.parent.mkdir()
    lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    path.write_text(">chr1 [topology=circular]\n" + "\n".join(lines) + "\n")
    return path


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    monkeypatch.delenv("FAKE_NF_RC", raising=False)
    monkeypatch.delenv("FAKE_NF_TRACE", raising=False)
    return args


def _run(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    result = CliRunner().invoke(cli, ["run", *argv])
    out = result.output + (result.stderr if result.stderr_bytes is not None else "")
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return result.exit_code, out


def _argv(path: Path) -> list[str]:
    assert path.is_file(), "giap run did not launch nextflow"
    return path.read_text().splitlines()


def test_run_launches_nextflow_with_the_spec_command(
    genome: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out)])
    assert code == ExitCode.SUCCESS, text
    args = _argv(fake_nf)
    workflow = Path(args[1])
    assert args[0] == "run"
    assert workflow.name == "main.nf" and workflow.parent.name == "workflow" and workflow.is_file()
    assert workflow.is_absolute()
    params = out.absolute() / "work" / "params.json"
    want = launch_command(workflow, params, "farmshare,apptainer", out.absolute(), resume=False)
    assert ["nextflow", *args] == want
    assert Path((fake_nf.parent / (fake_nf.name + ".cwd")).read_text().strip()).resolve() == (
        out.resolve()
    ), "nextflow must run with OUT as its working directory"


def test_run_writes_params_and_resolved_config(genome: Path, fake_nf: Path, tmp_path: Path) -> None:
    out = tmp_path / "OUT"
    code, text = _run(
        [
            "--genome",
            str(genome),
            "--outdir",
            str(out),
            "--threads",
            "4",
            "--modules",
            "core,mobile",
        ]
    )
    assert code == ExitCode.SUCCESS, text
    params = json.loads((out / "work" / "params.json").read_text())
    assert params["threads"] == 4
    assert params["modules"] == ["core", "mobile"]
    assert params["inputs"], "params.json has no inputs"
    for key, value in params["inputs"].items():
        assert Path(value).is_absolute() and Path(value).exists(), f"input {key}: {value}"
    resolved = yaml.safe_load((out / "resolved_config.yaml").read_text())
    assert resolved["threads"] == 4
    assert Path(params["config"]) == (out / "resolved_config.yaml").absolute()
    assert params["annotation"] == resolved["annotation"]


def test_run_memory_reaches_params_in_a_nextflow_form(
    genome: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--outdir", str(out), "--memory", "32.GB"])
    assert code == ExitCode.SUCCESS, text
    mem = json.loads((out / "work" / "params.json").read_text())["memory"]
    m = re.fullmatch(r"(\d+(?:\.\d+)?)\s*\.?\s*([KMGT]B)", str(mem))
    assert m, f"params.memory {mem!r} is not a Nextflow memory string"
    assert float(m.group(1)) * {"KB": 2**-20, "MB": 2**-10, "GB": 1, "TB": 2**10}[m.group(2)] == 32


@pytest.mark.parametrize("profile", ["farmshare_local", "farmshare,conda"])
def test_run_profile_and_resume(genome: Path, fake_nf: Path, tmp_path: Path, profile: str) -> None:
    out = tmp_path / "OUT"
    code, text = _run(
        ["--genome", str(genome), "--outdir", str(out), "--profile", profile, "--resume"]
    )
    assert code == ExitCode.SUCCESS, text
    args = _argv(fake_nf)
    assert args[args.index("-profile") + 1] == profile
    assert args[-1] == "-resume"
    assert args.count("-resume") == 1


def test_run_failed_step_exits_4_naming_the_step(
    genome: Path, fake_nf: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    trace = tmp_path / "trace_src.tsv"
    wd = tmp_path / "work" / "ab" / "cdef"
    trace.write_text(
        TRACE_HEADER
        + "1\taa/111111\t-\tANNOTATE:GIAP_STEP (parse_bakta)\tparse_bakta\tCOMPLETED\t0\t1\t"
        + f"{tmp_path}\n"
        + f"2\tab/cdef00\t-\tANNOTATE:TERMINATORNET (chr1)\tchr1\tFAILED\t1\t1\t{wd}\n"
    )
    monkeypatch.setenv("FAKE_NF_TRACE", str(trace))
    monkeypatch.setenv("FAKE_NF_RC", "1")
    code, text = _run(["--genome", str(genome), "--outdir", str(out)])
    assert code == ExitCode.STEP_FAILED, text
    assert "terminatornet" in text
    assert f"{wd}/.command.log" in text


def test_run_failure_without_trace_still_exits_4(
    genome: Path, fake_nf: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_NF_RC", "1")
    code, text = _run(["--genome", str(genome), "--outdir", str(tmp_path / "OUT")])
    assert code == ExitCode.STEP_FAILED, text


def test_run_without_nextflow_on_path_is_an_environment_error(
    genome: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "emptybin"
    empty.mkdir()
    py = Path(os.__file__).parents[2] / "bin"  # keep a python on PATH, drop nextflow
    monkeypatch.setenv("PATH", f"{empty}{os.pathsep}{py}")
    code, text = _run(["--genome", str(genome), "--outdir", str(tmp_path / "OUT")])
    assert code == ExitCode.ENVIRONMENT, text
    assert "nextflow" in text.lower()


def test_run_is_no_longer_a_stub(genome: Path, fake_nf: Path, tmp_path: Path) -> None:
    code, text = _run(["--genome", str(genome), "--outdir", str(tmp_path / "OUT")])
    assert "not implemented" not in text.lower()
    assert code != ExitCode.INTERNAL
