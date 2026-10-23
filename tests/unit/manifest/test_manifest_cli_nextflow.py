"""MVP-03a ``giap run`` (Nextflow executor) run tail: the manifest after success and after a failed
step, the events rotation and the trace → steps/tools gathering
(docs/specs/MVP-03-manifest-finalize.md D-A, D-B, D-C ``steps``/``tools``/``seeds`` rows; edge
cases "A failed run", "trace.tsv absent", "-resume", "Preflight error"; lead decisions 5, 8, 9, 12).

A fake ``nextflow`` on ``PATH`` (as tests/unit/test_cli_run_nextflow.py) exits ``$FAKE_NF_RC``,
copies ``$FAKE_NF_TRACE`` to the ``-with-trace`` path and, with ``$FAKE_NF_STEP_EVENTS``, publishes
that file as ``OUT/logs/steps/parse_bakta.events.jsonl`` (what a GIAP_STEP task publishes;
``append_step_events`` copies it into ``OUT/logs/events.jsonl``). The trace rows point at task
directories holding ``versions.yml``. Nothing real is launched; ``$GIAP_ROOT`` and
``$GIAP_DB_DIR`` are fakes from ``helpers.fake_env``.
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
from helpers.fake_env import fake_giap_root, write_setup_manifest
from helpers.manifest_data import validate

from giap.exitcodes import ExitCode

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
DIGEST = "7" * 64
STEP_EVENT = {"ts": "2026-10-22T00:00:00+00:00", "step": "parse_bakta", "level": "WARNING",
              "event": "parse_bakta.run_flag", "data": {"flag": "anti_sd_default"}}  # fmt: skip


#: The tools/databases readers moved to MVP-03c (c06f140); 03c reverts it and drops this mark.
MVP03C = pytest.mark.xfail(reason="MVP-03c: tools/databases readers", strict=True)


@pytest.fixture
def nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "giap_root")))
    db = tmp_path / "db"
    db.mkdir()
    write_setup_manifest(db)
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
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
    return {
        "genome": genome,
        "proteins": proteins,
        "events": events,
        "args": tmp_path / "nf_args.txt",
    }


def _task(tmp_path: Path, name: str, versions: str) -> Path:
    d = tmp_path / "tasks" / name
    d.mkdir(parents=True)
    (d / "versions.yml").write_text(versions)
    (d / ".command.log").write_text("log\n")
    return d


def _trace(tmp_path: Path, rows: list[dict[str, str]]) -> Path:
    path = tmp_path / "trace_src.tsv"
    lines = ["\t".join(COLS)]
    for i, r in enumerate(rows, 1):
        full = {c: "-" for c in COLS} | {"task_id": str(i), "attempt": "1", "exit": "0"} | r
        lines.append("\t".join(full[c] for c in COLS))
    path.write_text("\n".join(lines) + "\n")
    return path


def _ok_trace(tmp_path: Path) -> Path:
    pb = _task(tmp_path, "pb", '"ANNOTATE:GIAP_STEP":\n    giap: 0.1.0\n')
    gm = _task(tmp_path, "gm", '"ANNOTATE:GENOMAD":\n    genomad: 1.12.0\n')
    return _trace(tmp_path, [
        {"name": "ANNOTATE:GIAP_STEP (parse_bakta)", "status": "COMPLETED", "realtime": "5.5s",
         "%cpu": "65.8%", "peak_rss": "269.9 MB", "workdir": str(pb)},
        {"name": "ANNOTATE:GENOMAD (genome)", "status": "COMPLETED", "realtime": "2m",
         "%cpu": "100.0%", "peak_rss": "1 GB", "workdir": str(gm),
         "container": f"docker://quay.io/biocontainers/genomad@sha256:{DIGEST}"},
    ])  # fmt: skip


def _run(nf: dict[str, Path], out: Path, *extra: str) -> tuple[int, str]:
    from giap.cli import cli

    argv = ["run", "--genome", str(nf["genome"]), "--proteins", str(nf["proteins"]),
            "--outdir", str(out), *extra]  # fmt: skip
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


def _events(path: Path) -> list[dict[str, Any]]:
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


@pytest.fixture
def ok(nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_STEP_EVENTS", str(nf["events"]))
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.SUCCESS, log
    return out


# ----------------------------------------------------------------------------- success


def test_success_manifest_validates(ok: Path) -> None:
    m = _manifest(ok)
    validate(m)
    assert (m["status"], m["exit_code"], m["executor"]) == ("completed", 0, "nextflow")
    assert m["profile"] == "farmshare,apptainer"
    assert m["resume"] is False
    assert m["detection_tier"] == json.loads(
        (ok / "preflight" / "preflight.json").read_text())["evidence_tier"]  # fmt: skip


def test_success_steps_from_the_trace(ok: Path) -> None:
    steps = _manifest(ok)["steps"]
    assert [s["name"] for s in steps] == ["parse_bakta", "genomad"]
    pb = steps[0]
    assert (pb["status"], pb["n_tasks"], pb["n_retries"]) == ("completed", 1, 0)
    assert pb["wall_time_s"] == pytest.approx(5.5)
    assert pb["cpu_time_s"] == pytest.approx(3.619)
    assert pb["peak_rss_bytes"] == 283010662
    assert steps[1]["wall_time_s"] == pytest.approx(120.0)


@MVP03C
def test_success_tools_from_versions_yml(ok: Path) -> None:
    tools = _manifest(ok)["tools"]
    assert [(t["name"], t["version"]) for t in tools] == [("genomad", "1.12.0")]  # giap left out
    assert tools[0]["container_digest"] == f"sha256:{DIGEST}"
    assert tools[0]["processes"]


def test_success_seeds_list_executed_py_steps(ok: Path) -> None:
    import zlib

    seeds = _manifest(ok)["seeds"]
    assert seeds["parse_bakta"] == zlib.crc32(b"parse_bakta")
    assert "genomad" not in seeds  # an external step has no stream


def test_success_step_events_reach_flags(ok: Path) -> None:
    m = _manifest(ok)
    assert "anti_sd_default" in m["flags"]
    assert "parse_bakta:parse_bakta.run_flag" in m["warnings"]


@MVP03C
def test_success_databases_and_provenance(ok: Path) -> None:
    import hashlib

    m = _manifest(ok)
    names = [d["name"] for d in m["databases"]]
    assert names and names == sorted(names)
    prov = ok / "work" / "provenance.json"
    assert prov.is_file()
    assert m["provenance_digest"] == hashlib.sha256(prov.read_bytes()).hexdigest()


def test_resume_recorded(nf: dict[str, Path], tmp_path: Path,
                         monkeypatch: pytest.MonkeyPatch) -> None:  # fmt: skip
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    out = tmp_path / "OUT"
    code, log = _run(nf, out, "--resume")
    assert code == ExitCode.SUCCESS, log
    assert _manifest(out)["resume"] is True


# ----------------------------------------------------------------------------- rotation


def test_rerun_holds_each_step_event_once(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_NF_TRACE", str(_ok_trace(tmp_path)))
    monkeypatch.setenv("FAKE_NF_STEP_EVENTS", str(nf["events"]))
    out = tmp_path / "OUT"
    for _ in range(2):
        code, log = _run(nf, out, "--resume")
        assert code == ExitCode.SUCCESS, log
    cur = _events(out / "logs" / "events.jsonl")
    prev = _events(out / "logs" / "events.previous.jsonl")
    assert sum(e.get("event") == "parse_bakta.run_flag" for e in cur) == 1
    assert sum(e.get("event") == "parse_bakta.run_flag" for e in prev) == 1
    assert sum(e.get("step") == "preflight" for e in cur) <= 1


# ----------------------------------------------------------------------------- failure


def test_failed_step_writes_a_failed_manifest(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    pb = _task(tmp_path, "pb", '"ANNOTATE:GIAP_STEP":\n    giap: 0.1.0\n')
    pt = _task(tmp_path, "pt", '"ANNOTATE:PROMOTECH":\n    promotech: 1.0\n')
    trace = _trace(tmp_path, [
        {"name": "ANNOTATE:GIAP_STEP (parse_bakta)", "status": "COMPLETED", "realtime": "1s",
         "workdir": str(pb)},
        {"name": "ANNOTATE:PROMOTECH (chunk_1)", "status": "FAILED", "exit": "137",
         "attempt": "1", "realtime": "3s", "workdir": str(pt)},
        {"name": "ANNOTATE:PROMOTECH (chunk_1)", "status": "FAILED", "exit": "137",
         "attempt": "2", "realtime": "3s", "workdir": str(pt)},
    ])  # fmt: skip
    monkeypatch.setenv("FAKE_NF_TRACE", str(trace))
    monkeypatch.setenv("FAKE_NF_RC", "1")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    assert "step promotech failed" in log, log
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["failed_step"]) == ("failed", 4, "promotech")
    assert m["failed_log"] == f"{pt}/.command.log"
    by = {s["name"]: s for s in m["steps"]}
    assert by["promotech"]["status"] == "failed"
    assert by["promotech"]["n_retries"] == 1
    assert by["parse_bakta"]["status"] == "completed"
    assert "finalize" not in by


def test_failure_without_trace_still_writes_a_manifest(
    nf: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_NF_RC", "1")
    out = tmp_path / "OUT"
    code, log = _run(nf, out)
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    validate(m)
    assert m["status"] == "failed"
    assert m["steps"] == []
    assert m["failed_step"] is None


def test_preflight_error_writes_no_manifest(nf: dict[str, Path], tmp_path: Path) -> None:
    from giap.cli import cli

    out = tmp_path / "OUT"
    res = CliRunner().invoke(cli, ["run", "--genome", str(nf["genome"]), "--outdir", str(out)])
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert not nf["args"].exists(), "nextflow must not start after a preflight error"
    assert not (out / "run_manifest.json").exists()
    assert not (out / "work" / "provenance.json").exists()
