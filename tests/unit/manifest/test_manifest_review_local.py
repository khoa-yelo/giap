"""MVP-03a review round (PR #196), local executor and the builder: blocking findings B1 and B2 and
the surviving mutants M17 and M19 (docs/specs/MVP-03-manifest-finalize.md D-B steps 1 and 5,
D-C ``steps``, edge cases "Exit 5" and "-resume", example 12).

Interface pinned for the fix:

* B1 (lead decision): the launch rotation also moves ``OUT/logs/trace.tsv`` to
  ``OUT/logs/trace.previous.tsv`` under ``--executor local`` (same rotation as Nextflow).
* B2: an unexpected exception after the provenance core (``local.publish`` raising ``OSError``,
  ``local.run_local`` raising anything that is not ``StepError``/``WiringError``) exits 5
  (``ExitCode.INTERNAL``) and the tail still writes ``run_manifest.json`` with ``status failed``
  and ``exit_code 5``; the exception text is reported (stderr or a file under ``OUT/logs``). If
  the tail raises too, no manifest is written and the exit stays 5.
* M17: a step whose last-attempt tasks are a mix of CACHED and COMPLETED is ``completed``;
  ``cached`` only when every task of the step was CACHED (example 12: "the steps that finished
  in run 1 have status cached"; a step with any task run in this session did not finish in run 1).
* M19: ``command_line`` starts with ``["giap", "run"]`` (the raw giap argv, not ``sys.argv``).

The toy graph is the one of ``test_manifest_cli_local.py`` (cut down; copied, not imported).
"""

from __future__ import annotations

import importlib
import json
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.manifest_data import build, spec, trace_row, validate
from helpers.mapping_steps import write_test_db

from giap.exitcodes import ExitCode


def _toy(reg: ModuleType, *, bad_event_line: bool = False) -> None:
    def x(ctx: Any) -> None:
        out = Path(ctx.outputs["b.txt"])
        out.write_text("b\n")
        if bad_event_line:  # makes the tail's read_events raise
            logs = out.parents[1] / "logs"
            logs.mkdir(parents=True, exist_ok=True)
            with (logs / "events.jsonl").open("a", encoding="utf-8") as fh:
                fh.write("this is not json\n")

    def y(ctx: Any) -> None:
        ctx.outputs["c.txt"].write_text(ctx.inputs["b.txt"].read_text() + "y\n")

    reg.step("x", inputs=[], outputs=["b.txt"])(x)
    reg.step("y", inputs=["b.txt"], outputs=["c.txt"])(y)
    reg.external("term_tool", inputs=[], outputs=["term_tool.tsv"], precomputed="terminators")
    importlib.reload(importlib.import_module("giap.steps.manifest"))


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    genome = tmp_path / "in" / "g.fa"
    genome.parent.mkdir()
    genome.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    pep = tmp_path / "in" / "p.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    layer = tmp_path / "in" / "terminators.tsv"
    layer.write_text("contig\tstart\tend\n")
    monkeypatch.setenv("GIAP_DB_DIR", str(write_test_db(tmp_path / "test_db")))
    return {"genome": genome, "peptides": pep, "layer": layer}


def _argv(env: dict[str, Path], out: Path) -> list[str]:
    return ["run", "--genome", str(env["genome"]), "--outdir", str(out), "--executor", "local",
            "--modules", "core", "--seed", "5", "--peptides", str(env["peptides"]),
            "--set", f"annotation.precomputed.terminators={env['layer']}"]  # fmt: skip


def _invoke_raw(argv: list[str]) -> Any:
    from giap.cli import cli

    return CliRunner().invoke(cli, argv)


def _invoke(argv: list[str]) -> tuple[int, str]:
    res = _invoke_raw(argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}; OUT has {sorted(p.name for p in out.glob('*'))}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _reported(res: Any, out: Path) -> str:
    logs = out / "logs"
    files = sorted(p for p in logs.rglob("*") if p.is_file()) if logs.is_dir() else []
    return res.output + "".join(p.read_text(encoding="utf-8", errors="replace") for p in files)


@pytest.fixture
def toy(registry: ModuleType) -> Iterator[ModuleType]:
    _toy(registry)
    yield registry


# ----------------------------------------------------------------------------- B1 rotation


def test_b1_local_rotates_a_leftover_trace(
    toy: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    (out / "logs").mkdir(parents=True)
    stale = "task_id\tname\tstatus\n1\tANNOTATE:PROMOTECH (chunk_1)\tFAILED\n"
    (out / "logs" / "trace.tsv").write_text(stale)
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    assert not (out / "logs" / "trace.tsv").exists(), "the earlier run's trace was left in place"
    prev = out / "logs" / "trace.previous.tsv"
    assert prev.is_file(), sorted(p.name for p in (out / "logs").iterdir())
    assert prev.read_text(encoding="utf-8") == stale
    m = _manifest(out)
    assert [s["name"] for s in m["steps"]] == ["x", "y"]
    assert m["failed_step"] is None


# ----------------------------------------------------------------------------- B2 exit 5


def _check_exit5_manifest(res: Any, out: Path, needle: str) -> None:
    assert isinstance(res.exception, SystemExit), (
        f"the exception escaped giap run (main() would exit 5 with no manifest): {res.exception!r}"
    )
    assert res.exit_code == ExitCode.INTERNAL, res.output
    assert (out / "work" / "provenance.json").is_file()
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["executor"]) == ("failed", 5, "local")
    assert m["failed_step"] is None
    assert needle in _reported(res, out), "the exception is reported nowhere"


def test_b2_publish_oserror_exits_5_with_a_failed_manifest(
    toy: ModuleType, env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.runner import local

    def bad(*a: Any, **k: Any) -> None:
        raise OSError("disk full on purpose")

    monkeypatch.setattr(local, "publish", bad)
    out = tmp_path / "OUT"
    res = _invoke_raw(_argv(env, out))
    _check_exit5_manifest(res, out, "disk full on purpose")
    m = _manifest(out)
    assert [s["name"] for s in m["steps"]] == ["x", "y"]  # the steps known so far


def test_b2_unexpected_error_inside_run_local_exits_5_with_a_failed_manifest(
    toy: ModuleType, env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.runner import local

    def bad(*a: Any, **k: Any) -> None:
        raise local.LayerMissing("annotation.precomputed.on_purpose")

    monkeypatch.setattr(local, "run_local", bad)
    out = tmp_path / "OUT"
    res = _invoke_raw(_argv(env, out))
    _check_exit5_manifest(res, out, "on_purpose")


def test_b2_tail_failure_after_exit_5_writes_nothing_and_stays_5(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.runner import local

    _toy(registry, bad_event_line=True)

    def bad(*a: Any, **k: Any) -> None:
        raise OSError("disk full on purpose")

    monkeypatch.setattr(local, "publish", bad)
    out = tmp_path / "OUT"
    res = _invoke_raw(_argv(env, out))
    assert isinstance(res.exception, SystemExit), repr(res.exception)
    assert res.exit_code == ExitCode.INTERNAL, res.output
    assert not (out / "run_manifest.json").exists()


# ----------------------------------------------------------------------------- M19


def test_m19_local_command_line_is_the_giap_argv(
    toy: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    cl = _manifest(out)["command_line"]
    assert cl[:2] == ["giap", "run"], cl
    assert cl[2:] == _argv(env, out)[1:], cl


# ----------------------------------------------------------------------------- M17


def _mixed(statuses: list[str]) -> dict[str, Any]:
    trace = [trace_row(f"ANNOTATE:PROMOTECH (chunk_{i})", st, realtime="1s")
             for i, st in enumerate(statuses, 1)]  # fmt: skip
    m = build(executor="nextflow", trace=trace, local_steps=(), resume=True,
              registered=[spec("promotech", kind="ext")])  # fmt: skip
    (s,) = m["steps"]
    return s


@pytest.mark.parametrize(
    "statuses",
    [["CACHED", "COMPLETED"], ["COMPLETED", "CACHED"], ["CACHED", "CACHED", "COMPLETED"]],
)
def test_m17_mixed_cached_and_completed_step_is_completed(statuses: list[str]) -> None:
    s = _mixed(statuses)
    assert s["status"] == "completed", s
    assert s["n_tasks"] == len(statuses)
    assert s["wall_time_s"] == pytest.approx(float(len(statuses)))  # cached rows keep metrics


def test_m17_all_cached_step_is_cached() -> None:
    assert _mixed(["CACHED", "CACHED"])["status"] == "cached"


def test_m17_mixed_with_a_failed_task_is_failed() -> None:
    assert _mixed(["CACHED", "COMPLETED", "FAILED"])["status"] == "failed"
