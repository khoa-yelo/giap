"""MVP-03a ``giap run --executor local`` run order, events rotation, exit codes and the manifest
(docs/specs/MVP-03-manifest-finalize.md D-A, D-B, D-C, D-D; edge cases "A failed run", "Preflight
error", "The tail fails after a successful graph", "Re-run of the local executor"; lead decisions
2, 5, 6 and 12; "Test ideas": ``giap _step manifest`` on its own).

Interface assumed:

* ``giap.steps.manifest`` registers the py step ``manifest`` (module ``core``): required input
  ``preflight``; optional ``provenance``, ``events``, ``trace``, ``annotation_meta``,
  ``sensitivity``, ``island_diagnostics``; output ``run_manifest``.
* ``giap.runner.local.TAIL_STEPS == ("manifest",)``; ``plan_local`` never plans a tail step.
* ``giap run`` writes ``OUT/run_manifest.json`` and, before the graph, ``OUT/work/provenance.json``
  (= ``canonical_json(provenance_core(manifest))``, so ``provenance_digest`` = its SHA-256).
* A tail failure after a successful graph prints ``step manifest failed; log:
  OUT/logs/manifest.log`` and exits 4; ``logs/manifest.log`` holds the traceback.

The graph is a toy registry (``tests/unit/conftest.py::registry`` plus a reload of
``giap.steps.manifest``): ``x`` (core, emits the ``c_assumed`` run flag as a WARNING event) →
``y``; ``mob`` (module mobile); ``bactermfinder`` (opt-in); ``term_tool`` (ext, replaced by
``annotation.precomputed.terminators``); ``never`` (needs a root nobody stages). Preflight and the
CLI are real; the database directory is ``helpers.mapping_steps.write_test_db`` (never the real
``$GIAP_DB_DIR``).
"""

from __future__ import annotations

import hashlib
import importlib
import json
import re
import zlib
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.manifest_data import validate
from helpers.mapping_steps import write_test_db

from giap.exitcodes import ExitCode

FLAG_EVENT = {"event": "x.run_flag", "data": {"flag": "c_assumed"}}


# ----------------------------------------------------------------------------- toy graph


def _with_manifest_step(reg: ModuleType) -> None:
    importlib.reload(importlib.import_module("giap.steps.manifest"))
    assert "manifest" in {s.name for s in reg.all_steps()}


def _toy(reg: ModuleType, *, boom: bool = False, bad_event_line: bool = False) -> None:
    def x(ctx: Any) -> None:
        ctx.log.warning("x: run flag c_assumed", extra=FLAG_EVENT)
        out = Path(ctx.outputs["b.txt"])
        out.write_text("b\n")
        if bad_event_line:  # OUT/work/.tmp.x.b.txt -> OUT/logs/events.jsonl
            logs = out.parents[1] / "logs"
            logs.mkdir(parents=True, exist_ok=True)
            with (logs / "events.jsonl").open("a", encoding="utf-8") as fh:
                fh.write("this is not json\n")

    def y(ctx: Any) -> None:
        ctx.outputs["c.txt"].write_text(ctx.inputs["b.txt"].read_text() + "y\n")

    def write(name: str) -> Any:
        def fn(ctx: Any) -> None:
            for p in ctx.outputs.values():
                Path(p).write_text(name)

        return fn

    def explode(ctx: Any) -> None:
        raise RuntimeError("boom on purpose")

    reg.step("x", inputs=[], outputs=["b.txt"])(x)
    reg.step("y", inputs=["b.txt"], outputs=["c.txt"])(y)
    reg.step("mob", inputs=[], outputs=["m.txt"], module="mobile")(write("mob"))
    reg.step("bactermfinder", inputs=[], outputs=["bt.txt"], module="opt_in")(write("bt"))
    reg.external("term_tool", inputs=[], outputs=["term_tool.tsv"], precomputed="terminators")
    reg.step("never", inputs=["absent_root"], outputs=["n.txt"])(write("never"))
    if boom:
        reg.step("boom", inputs=["c.txt"], outputs=["boom.txt"])(explode)
    _with_manifest_step(reg)


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


def _argv(env: dict[str, Path], out: Path, *, peptides: bool = True) -> list[str]:
    argv = ["run", "--genome", str(env["genome"]), "--outdir", str(out), "--executor", "local",
            "--modules", "core", "--seed", "5",
            "--set", f"annotation.precomputed.terminators={env['layer']}"]  # fmt: skip
    if peptides:
        argv += ["--peptides", str(env["peptides"])]
    return argv


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


def _events(path: Path) -> list[dict[str, Any]]:
    return [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


@pytest.fixture
def ok_run(registry: ModuleType, env: dict[str, Path], tmp_path: Path) -> Iterator[Path]:
    _toy(registry)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    yield out


# ----------------------------------------------------------------------------- registry


def test_manifest_step_is_registered(registry: ModuleType) -> None:
    _with_manifest_step(registry)
    s = registry.get("manifest")
    assert (s.kind, s.module) == ("py", "core")
    assert s.inputs == ("preflight",)
    assert {"provenance", "events", "trace", "annotation_meta", "sensitivity",
            "island_diagnostics"} <= set(s.optional_inputs)  # fmt: skip
    assert s.outputs == ("run_manifest",)


def test_tail_steps_are_left_out_of_the_plan(registry: ModuleType) -> None:
    from giap.config import GiapConfig
    from giap.runner import local

    assert local.TAIL_STEPS == ("manifest",)
    _toy(registry)
    cfg = GiapConfig.model_validate(
        {"seed": 5, "modules": ["core"], "annotation": {"precomputed": {"terminators": "/t.tsv"}}}
    )
    plan = local.plan_local(cfg, [s.name for s in registry.all_steps()], {"preflight"})
    assert "manifest" not in plan
    assert plan == ["x", "y"]


# ----------------------------------------------------------------------------- success


def test_success_writes_a_valid_manifest(ok_run: Path) -> None:
    m = _manifest(ok_run)
    validate(m)
    assert (m["status"], m["exit_code"], m["executor"]) == ("completed", 0, "local")
    assert m["failed_step"] is None
    assert m["tools"] == []
    assert m["resume"] is False
    assert m["modules_off"] == ["bgc", "function", "mobile"]


def test_success_skipped_steps_and_reasons(ok_run: Path) -> None:
    skipped = {s["name"]: s["reason"] for s in _manifest(ok_run)["skipped"]}
    assert skipped["mob"] == "module_off:mobile"
    assert skipped["bactermfinder"] == "opt_in_off"
    assert skipped["term_tool"] == "precomputed:terminators"
    assert skipped["never"] == "not_run"
    assert "x" not in skipped and "y" not in skipped


def test_success_steps_in_plan_order_with_local_timing(ok_run: Path) -> None:
    steps = _manifest(ok_run)["steps"]
    assert [s["name"] for s in steps] == ["x", "y"]
    for s in steps:
        assert s["status"] == "completed"
        assert s["wall_time_s"] >= 0 and s["cpu_time_s"] >= 0
        assert "peak_rss_bytes" not in s  # local records no RSS (question 2)


def test_success_seeds(ok_run: Path) -> None:
    seeds = _manifest(ok_run)["seeds"]
    assert seeds["master"] == 5
    assert seeds["x"] == zlib.crc32(b"x") and seeds["y"] == zlib.crc32(b"y")
    assert set(seeds) - {"manifest"} == {"master", "x", "y"}


def test_success_inputs_reuse_preflight_and_hash_the_layer(
    ok_run: Path, env: dict[str, Path]
) -> None:
    m = _manifest(ok_run)
    pf = json.loads((ok_run / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    by = {r["role"]: r for r in m["inputs"]}
    assert len(by) == len(m["inputs"]), "a role appears twice"
    want = {r["role"]: r["sha256"] for r in pf["inputs"]}
    assert by["genome"]["sha256"] == want["genome"]
    assert by["peptides"]["sha256"] == want["peptides"]
    layer = by["precomputed:terminators"]
    assert layer["sha256"] == hashlib.sha256(env["layer"].read_bytes()).hexdigest()
    assert layer["size"] == env["layer"].stat().st_size
    assert {"proteins", "abundance", "annotation", "reference_peptides",
            "replicon_depth"} <= set(m["inputs_absent"])  # fmt: skip
    assert "peptides" not in m["inputs_absent"]
    assert m["detection_tier"] == pf["evidence_tier"]


def test_success_flags_and_warnings_from_events(ok_run: Path) -> None:
    m = _manifest(ok_run)
    assert "c_assumed" in m["flags"]
    assert "x:x.run_flag" in m["warnings"]


def test_success_databases_from_the_test_db(ok_run: Path) -> None:
    dbs = _manifest(ok_run)["databases"]
    assert dbs, "write_test_db lists every requirement"
    assert [d["name"] for d in dbs] == sorted(d["name"] for d in dbs)
    assert all(d["checksum"] is None for d in dbs)  # provisional form


def test_provenance_core_written_and_digest_matches(ok_run: Path) -> None:
    from giap.manifest import canonical_json, provenance_core

    m = _manifest(ok_run)
    path = ok_run / "work" / "provenance.json"
    assert path.is_file()
    assert path.read_text(encoding="utf-8") == canonical_json(provenance_core(m))
    assert m["provenance_digest"] == hashlib.sha256(path.read_bytes()).hexdigest()


def test_manifest_file_is_canonical_json(ok_run: Path) -> None:
    from giap.manifest import canonical_json

    text = (ok_run / "run_manifest.json").read_text(encoding="utf-8")
    assert text == canonical_json(json.loads(text))


def test_two_fresh_runs_have_equal_stable_views(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    from giap.manifest import stable_view

    _toy(registry)
    a, b = tmp_path / "A" / "OUT", tmp_path / "B" / "OUT"
    for out in (a, b):
        code, log = _invoke(_argv(env, out))
        assert code == ExitCode.SUCCESS, log
    ma, mb = _manifest(a), _manifest(b)
    assert ma["outdir"] != mb["outdir"]
    assert stable_view(ma) == stable_view(mb)


# ----------------------------------------------------------------------------- rotation


def test_events_rotated_at_start(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry)
    out = tmp_path / "OUT"
    (out / "logs").mkdir(parents=True)
    (out / "logs" / "events.jsonl").write_text('{"sentinel": true}\n')
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    prev = (out / "logs" / "events.previous.jsonl").read_text(encoding="utf-8")
    assert '"sentinel"' in prev
    cur = _events(out / "logs" / "events.jsonl")
    assert not any("sentinel" in e for e in cur)
    assert sum(e.get("event") == "x.run_flag" for e in cur) == 1


def test_rerun_into_same_out_keeps_one_generation(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry)
    out = tmp_path / "OUT"
    for _ in range(2):
        code, log = _invoke(_argv(env, out))
        assert code == ExitCode.SUCCESS, log
    cur = _events(out / "logs" / "events.jsonl")
    prev = _events(out / "logs" / "events.previous.jsonl")
    assert sum(e.get("event") == "x.run_flag" for e in cur) == 1
    assert sum(e.get("event") == "x.run_flag" for e in prev) == 1
    assert sum(e.get("step") == "preflight" for e in cur) <= 1
    m = _manifest(out)
    assert [s["name"] for s in m["steps"]] == ["x", "y"]  # this run only


# ----------------------------------------------------------------------------- failures


def test_failed_step_still_writes_the_manifest(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry, boom=True)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.STEP_FAILED, log
    m = _manifest(out)
    validate(m)
    assert (m["status"], m["exit_code"], m["failed_step"]) == ("failed", 4, "boom")
    assert str(m["failed_log"]).endswith("boom.log")
    by = {s["name"]: s["status"] for s in m["steps"]}
    assert by["x"] == by["y"] == "completed"
    assert by["boom"] == "failed"
    assert "c_assumed" in m["flags"]  # flags known so far


def test_preflight_error_writes_no_manifest(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out, peptides=False))  # no evidence: preflight exit 1
    assert code == ExitCode.INVALID_INPUT, log
    assert (out / "preflight" / "preflight.json").is_file()
    assert not (out / "run_manifest.json").exists()
    assert not (out / "work" / "provenance.json").exists()


def test_tail_failure_after_success_exits_4(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry, bad_event_line=True)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.STEP_FAILED, log
    assert re.search(r"step manifest failed; log: \S*logs/manifest\.log", log), log
    text = (out / "logs" / "manifest.log").read_text(encoding="utf-8")
    assert "Traceback" in text or "Error" in text, text


def test_tail_failure_never_hides_the_graph_failure(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path
) -> None:
    _toy(registry, boom=True, bad_event_line=True)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.STEP_FAILED, log
    assert "boom" in log, log
    assert (out / "logs" / "boom.log").is_file()


# ----------------------------------------------------------------------------- _step manifest


def test_step_manifest_runs_on_its_own(
    ok_run: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The registry path (``giap _step manifest``) works with only the required input."""
    task = tmp_path / "task"
    task.mkdir()
    monkeypatch.chdir(task)
    target = task / "run_manifest.json"
    code, log = _invoke(["_step", "manifest", "--in",
                         f"preflight={ok_run / 'preflight' / 'preflight.json'}",
                         "--out", f"run_manifest={target}"])  # fmt: skip
    assert code == ExitCode.SUCCESS, log
    m = json.loads(target.read_text(encoding="utf-8"))
    validate(m)
