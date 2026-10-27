"""MVP-03a review round (PR #196), Nextflow executor: blocking findings B1 and B2 and the surviving
mutants M7 and M19 (docs/specs/MVP-03-manifest-finalize.md D-B steps 1 and 5, edge cases
"trace.tsv absent" and "Exit 5"; Goal "every run that gets past preflight writes
run_manifest.json").

Interface pinned for the fix:

* B1 (lead decision): the launch rotation (D-B step 1, next to ``events.jsonl`` →
  ``events.previous.jsonl``) also moves ``OUT/logs/trace.tsv`` to ``OUT/logs/trace.previous.tsv``
  under both executors, before anything is launched. A run whose Nextflow writes no trace then
  has no ``logs/trace.tsv``: ``steps == []``, ``failed_step null`` and no step events copied by
  ``append_step_events`` from the earlier run.
* B2: an unexpected exception from ``gippa.cli._run_nextflow`` (after the provenance core) exits
  5 (``ExitCode.INTERNAL``), and the tail still writes ``run_manifest.json`` with ``status
  failed`` and ``exit_code 5``; the exception text is reported (stderr or a file under
  ``OUT/logs``).
* M7: ``gippa.runner.nextflow._step_name`` maps an alias named ``*_STEP`` to its tag.
* M19: ``command_line`` is the raw ``gippa`` argv (``["gippa", "run", ...]``), never
  ``sys.argv`` (pytest's argv under CliRunner).

The fake ``nextflow`` harness is that of ``test_manifest_cli_nextflow.py`` (copied, not imported:
test modules are not packages here).
"""

from __future__ import annotations

import json
import os
import random
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import fake_gippa_root, write_setup_manifest
from helpers.manifest_data import validate

from gippa.exitcodes import ExitCode

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
if [ -n "$FAKE_NF_TRACE" ]; then
  prev=""
  for a in "$@"; do
    if [ "$prev" = "-with-trace" ]; then mkdir -p "$(dirname "$a")"; cp "$FAKE_NF_TRACE" "$a"; fi
    prev="$a"
  done
fi
if [ -n "$FAKE_NF_STEP_EVENTS" ]; then
  mkdir -p logs/steps
  cp "$FAKE_NF_STEP_EVENTS" logs/steps/parse_bakta.events.jsonl
fi
exit "${FAKE_NF_RC:-0}"
"""
COLS = ["task_id", "hash", "native_id", "name", "tag", "status", "exit", "attempt", "realtime",
        "%cpu", "peak_rss", "workdir", "container"]  # fmt: skip
STEP_EVENT = {"ts": "2026-10-22T00:00:00+00:00", "step": "parse_bakta", "level": "WARNING",
              "event": "parse_bakta.run_flag", "data": {"flag": "anti_sd_default"}}  # fmt: skip


@pytest.fixture
def nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIPPA_ROOT", str(fake_gippa_root(tmp_path / "gippa_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    for var in ("SLURM_CPUS_ON_NODE", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE"):
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(tmp_path / "nf_args.txt"))
    for var in ("FAKE_NF_RC", "FAKE_NF_TRACE", "FAKE_NF_STEP_EVENTS"):
        monkeypatch.delenv(var, raising=False)
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome.write_text(">chr1 [topology=circular]\n" + "\n".join(seq[i:i + 60]
                      for i in range(0, 6000, 60)) + "\n")  # fmt: skip
    aa = "ACDEFGHIKLMNPQRSTVWY"
    proteins = tmp_path / "in" / "proteins.faa"
    proteins.write_text("".join(f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n"
                                for i in range(3)))  # fmt: skip
    events = tmp_path / "step_events.jsonl"
    events.write_text(json.dumps(STEP_EVENT, sort_keys=True) + "\n")
    return {"genome": genome, "proteins": proteins, "events": events,
            "args": tmp_path / "nf_args.txt"}  # fmt: skip


def _task(tmp_path: Path, name: str, versions: str = "") -> Path:
    d = tmp_path / "tasks" / name
    d.mkdir(parents=True)
    (d / "versions.yml").write_text(versions)
    (d / ".command.log").write_text("log\n")
    return d


def _trace(tmp_path: Path, rows: list[dict[str, str]], name: str = "trace_src.tsv") -> Path:
    path = tmp_path / name
    lines = ["\t".join(COLS)]
    for i, r in enumerate(rows, 1):
        full = {c: "-" for c in COLS} | {"task_id": str(i), "attempt": "1", "exit": "0"} | r
        lines.append("\t".join(full[c] for c in COLS))
    path.write_text("\n".join(lines) + "\n")
    return path


def _ok_trace(tmp_path: Path) -> Path:
    pb = _task(tmp_path, "pb", '"ANNOTATE:GIPPA_STEP":\n    gippa: 0.1.0\n')
    return _trace(tmp_path, [
        {"name": "ANNOTATE:GIPPA_STEP (parse_bakta)", "status": "COMPLETED", "realtime": "5.5s",
         "%cpu": "65.8%", "peak_rss": "269.9 MB", "workdir": str(pb)},
    ], "ok_trace.tsv")  # fmt: skip


def _failed_trace(tmp_path: Path) -> Path:
    pt = _task(tmp_path, "pt")
    return _trace(tmp_path, [
        {"name": "ANNOTATE:PROMOTECH (chunk_1)", "status": "FAILED", "exit": "1",
         "workdir": str(pt)},
    ], "failed_trace.tsv")  # fmt: skip


def _argv(nf: dict[str, Path], out: Path, *extra: str) -> list[str]:
    return ["run", "--genome", str(nf["genome"]), "--proteins", str(nf["proteins"]),
            "--outdir", str(out), *extra]  # fmt: skip


def _run(nf: dict[str, Path], out: Path, *extra: str) -> tuple[int, str]:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, _argv(nf, out, *extra))
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _manifest(out: Path) -> dict[str, Any]:
    path = out / "run_manifest.json"
    assert path.is_file(), f"no {path}; OUT has {sorted(p.name for p in out.glob('*'))}"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _events(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


# ----------------------------------------------------------------------------- B1 stale trace


def _second_run_without_trace(
    nf: dict[str, Path], out: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> tuple[int, str]:
    """Run 2 into the same OUT: Nextflow dies before the first task (no trace, no events)."""
    monkeypatch.delenv("FAKE_NF_TRACE", raising=False)
    monkeypatch.delenv("FAKE_NF_STEP_EVENTS", raising=False)
    monkeypatch.setenv("FAKE_NF_RC", "1")
    return _run(nf, out, *extra)


def test_b1_stale_failed_trace_is_not_this_runs_failed_step(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    monkeypatch.setenv("FAKE_NF_TRACE", str(_failed_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_RC", "1")
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    assert _manifest(out)["failed_step"] == "promotech"  # run 1 as before
    code, log = _second_run_without_trace(nf, out, monkeypatch)
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    validate(m)
    assert m["status"] == "failed"
    assert m["failed_step"] is None, f"stale failed_step from run 1: {m['failed_step']}"
    assert m["failed_log"] is None
    assert m["steps"] == [], f"stale steps from run 1: {m['steps']}"
    assert "step promotech failed" not in log, log


def test_b1_stale_successful_trace_gives_no_steps_events_or_seeds(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_STEP_EVENTS", str(nf["events"]))
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    assert "anti_sd_default" in _manifest(out)["flags"]  # run 1 as before
    code, log = _second_run_without_trace(nf, out, monkeypatch)
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    validate(m)
    assert m["steps"] == [], f"stale steps from run 1: {m['steps']}"
    assert m["failed_step"] is None
    assert "parse_bakta" not in m["seeds"], m["seeds"]
    assert "anti_sd_default" not in m["flags"], m["flags"]
    assert "parse_bakta:parse_bakta.run_flag" not in m["warnings"], m["warnings"]
    cur = _events(out / "logs" / "events.jsonl")
    stale = [e for e in cur if e.get("event") == "parse_bakta.run_flag"]
    assert not stale, f"run 1's step events copied into run 2's events.jsonl: {stale}"


def test_b1_trace_is_rotated_to_trace_previous(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    out = tmp_path / "OUT"
    src = _failed_trace(tmp_path)
    monkeypatch.setenv("FAKE_NF_TRACE", str(src))
    monkeypatch.setenv("FAKE_NF_RC", "1")
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    run1 = (out / "logs" / "trace.tsv").read_text(encoding="utf-8")
    code, log = _second_run_without_trace(nf, out, monkeypatch)
    assert code == ExitCode.STEP_FAILED, log
    assert not (out / "logs" / "trace.tsv").exists(), "run 1's trace left as this run's trace"
    prev = out / "logs" / "trace.previous.tsv"
    assert prev.is_file(), sorted(p.name for p in (out / "logs").iterdir())
    assert prev.read_text(encoding="utf-8") == run1


def test_b1_resume_without_a_new_trace_gives_no_stale_steps(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Under ``-resume`` Nextflow rewrites the trace from its cache, so the rotation loses
    nothing; a resumed launch that dies before the session still reports no steps."""
    out = tmp_path / "OUT"
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    code, log = _second_run_without_trace(nf, out, monkeypatch, "--resume")
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    assert m["resume"] is True
    assert m["steps"] == []
    assert m["failed_step"] is None


def test_b1_a_new_trace_replaces_the_rotated_one(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The rotation never hides this run's own trace."""
    out = tmp_path / "OUT"
    monkeypatch.setenv("FAKE_NF_TRACE", str(_failed_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_RC", "1")
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_RC", "0")
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    m = _manifest(out)
    assert [s["name"] for s in m["steps"]] == ["parse_bakta"]
    assert m["failed_step"] is None
    assert (out / "logs" / "trace.previous.tsv").is_file()


# ----------------------------------------------------------------------------- B2 exit 5


def test_b2_exception_in_run_nextflow_exits_5_with_a_failed_manifest(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gippa.cli as gcli

    def boom(*a: Any, **k: Any) -> int:
        raise RuntimeError("nextflow wrapper exploded on purpose")

    monkeypatch.setattr(gcli, "_run_nextflow", boom)
    out = tmp_path / "OUT"
    res = CliRunner().invoke(gcli.cli, _argv(nf, out))
    assert isinstance(res.exception, SystemExit), (
        f"the exception escaped gippa run (main() would exit 5 with no manifest): {res.exception!r}"
    )
    assert res.exit_code == ExitCode.INTERNAL, res.output
    assert (out / "work" / "provenance.json").is_file()
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["executor"]) == ("failed", 5, "nextflow")
    assert m["failed_step"] is None
    text = res.output + "".join(
        p.read_text(encoding="utf-8", errors="replace")
        for p in (out / "logs").rglob("*") if p.is_file()
    )  # fmt: skip
    assert "exploded on purpose" in text, "the exception is reported nowhere"


# ----------------------------------------------------------------------------- M7 / M19


@pytest.mark.parametrize(
    ("name", "step"),
    [
        ("OPERONS:OPERONS_STEP (operons)", "operons"),
        ("ANNOTATE:GIPPA_STEP (parse_bakta)", "parse_bakta"),
        ("ANNOTATE:PROMOTECH (chunk_1)", "promotech"),
    ],
)
def test_m7_step_name_of_a_star_step_alias_is_its_tag(name: str, step: str) -> None:
    from gippa.runner.nextflow import _step_name

    assert _step_name(name) == step


def test_m7_star_step_alias_groups_under_its_tag_in_the_manifest() -> None:
    from helpers.manifest_data import build, spec, trace_row

    trace = [trace_row("OPERONS:OPERONS_STEP (operons)", realtime="1s")]
    m = build(executor="nextflow", trace=trace, local_steps=(), registered=[spec("operons")])
    assert [s["name"] for s in m["steps"]] == ["operons"]


def test_m19_command_line_is_the_gippa_argv(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    cl = _manifest(out)["command_line"]
    assert cl[:2] == ["gippa", "run"], cl
    assert cl[2:] == _argv(nf, out)[1:], cl
