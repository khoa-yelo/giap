"""M9-16b: ``giap run --executor local`` re-renders ``OUT/report/report.html`` after writing
``run_manifest.json`` (docs/specs/M9-16b-report-after-manifest.md behaviours 1, 2, 5, 6, 7;
examples 1, 3, 5, 6, 7, 8, 9; lead decisions 1, 3, 5).

Interface pinned here (M9-16b test author, beyond the spec text):

* the CLI tail looks ``rerender_report`` up on ``giap.steps.plots`` at call time (a lazy
  ``from giap.steps.plots import rerender_report`` inside ``_run_tail``), so a monkeypatched
  attribute is the one called; it is called once, as ``rerender_report(outdir, cfg)``;
* ``_run_tail`` looks ``run_tail`` up on ``giap.steps.manifest`` at call time (as today);
* the provenance table of the page is ``<table data-table-provenance="run">`` with one
  ``<tr><th>label</th><td>value</td></tr>`` per line (M9-16);
* the KEGG date reaches the page end to end: ``build_manifest`` fills ``databases[]`` from
  ``$GIAP_DB_DIR/db_manifest.json`` (MVP-03c, #236, on main; this branch needs origin/main merged
  in for example 5 to pass), and the items' ``retrieved`` date is what the page prints.

The graph is a toy registry (``tests/unit/conftest.py::registry``): ``x`` -> ``y``, plus the real
``plots`` and ``report`` step functions of ``giap.steps.plots`` registered into it, and a reload
of ``giap.steps.manifest``. Plots are fakes (tests/helpers/report_fixtures.py). Preflight and the
CLI are real; ``$GIAP_DB_DIR`` is a temporary copy of tests/data/db_manifest (never the real one).
"""

from __future__ import annotations

import hashlib
import importlib
import json
import logging
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.fake_env import setup_manifest
from helpers.report_content import write_category_summary
from helpers.report_fixtures import Node, parse, register_fake, swap_registry

from giap.exitcodes import ExitCode

NOT_WRITTEN = "run_manifest.json not written by this GIAP version"
KEGG_RETRIEVED = "2026-10-01"
DIGEST_LABEL = "manifest digest (SHA-256 of run_manifest.json)"
TMP_NAME = ".report.html.tmp"


# ----------------------------------------------------------------------------- toy graph


def _toy(reg: ModuleType, *, boom: bool = False) -> None:
    from giap.steps import plots as P

    def x(ctx: Any) -> None:
        Path(ctx.outputs["b.txt"]).write_text("b\n")

    def y(ctx: Any) -> None:
        ctx.outputs["c.txt"].write_text(ctx.inputs["b.txt"].read_text() + "y\n")

    def explode(ctx: Any) -> None:
        raise RuntimeError("operons raised on purpose")

    reg.step("x", inputs=[], outputs=["b.txt"])(x)
    reg.step("y", inputs=["b.txt"], outputs=["c.txt"])(y)
    if boom:
        reg.step("operons", inputs=["c.txt"], outputs=["operons.txt"])(explode)
    reg.step("plots", inputs=["preflight"], outputs=["plots"])(P.plots_step)
    reg.step("report", inputs=["plots", "preflight"], outputs=["report"])(P.report_step)
    importlib.reload(importlib.import_module("giap.steps.manifest"))
    assert {"plots", "report", "manifest"} <= {s.name for s in reg.all_steps()}


def _db(path: Path) -> Path:
    """A ``giap setup``-form db_manifest.json with a KEGG retrieval date (no contaminants)."""
    data: dict[str, Any] = setup_manifest()
    items = []
    for item in data["items"]:
        if item["name"] == "contaminants":
            continue
        item = dict(item)
        if str(item["name"]).startswith("kegg"):
            item["retrieved"] = KEGG_RETRIEVED
        items.append(item)
    path.mkdir(parents=True, exist_ok=True)
    (path / "db_manifest.json").write_text(json.dumps({**data, "items": items}, indent=2))
    return path


@pytest.fixture
def env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    genome = tmp_path / "in" / "g.fa"
    genome.parent.mkdir()
    genome.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    pep = tmp_path / "in" / "p.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    monkeypatch.setenv("GIAP_DB_DIR", str(_db(tmp_path / "test_db")))
    swap_registry(monkeypatch)
    register_fake("V07")
    return {"genome": genome, "peptides": pep}


def _argv(env: dict[str, Path], out: Path, modules: str = "core") -> list[str]:
    return ["run", "--genome", str(env["genome"]), "--peptides", str(env["peptides"]),
            "--outdir", str(out), "--executor", "local", "--modules", modules, "--seed", "5",
            "--set", "report.formats=[svg]", "--set", "report.dpi=50"]  # fmt: skip


def _invoke(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output + (res.stderr if res.stderr_bytes is not None else "")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _files(out: Path) -> dict[str, str]:
    """SHA-256 of every file under ``OUT`` outside ``logs/`` and ``report/report.html``."""
    keep: dict[str, str] = {}
    for p in sorted(out.rglob("*")):
        rel = p.relative_to(out).as_posix()
        if p.is_file() and not rel.startswith("logs/") and rel != "report/report.html":
            keep[rel] = _sha(p)
    return keep


def _manifest(out: Path) -> dict[str, Any]:
    data = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _prov(page_path: Path) -> dict[str, str]:
    page = parse(page_path.read_text(encoding="utf-8"))
    tables = page.root.by_attr("data-table-provenance", "run")
    assert tables, "no provenance table: the page has no manifest"
    rows: dict[str, str] = {}
    for tr in tables[0].by_tag("tr"):
        cells = [c for c in tr.children if isinstance(c, Node) and c.tag in ("th", "td")]
        if len(cells) == 2:
            rows[cells[0].text().strip()] = cells[1].text().strip()
    return rows


def _section(page_path: Path, sid: str) -> str:
    page = parse(page_path.read_text(encoding="utf-8"))
    return " ".join(page.section(sid).text().split())


@dataclass
class Call:
    args: tuple[Any, ...]
    manifest: bytes | None
    files: dict[str, str]
    events: bytes | None
    records: list[logging.LogRecord] = field(default_factory=list)


class _Capture(logging.Handler):
    def __init__(self) -> None:
        super().__init__(logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


@pytest.fixture
def spy(monkeypatch: pytest.MonkeyPatch) -> list[Call]:
    """Wrap ``giap.steps.plots.rerender_report``: record the call and the state of ``OUT`` then."""
    from giap.steps import plots as P

    calls: list[Call] = []
    real = getattr(P, "rerender_report", None)

    def wrapper(outdir: Path, cfg: Any, *a: Any, **k: Any) -> Any:
        out = Path(outdir)
        man, ev = out / "run_manifest.json", out / "logs" / "events.jsonl"
        call = Call((outdir, cfg, *a), man.read_bytes() if man.is_file() else None, _files(out),
                    ev.read_bytes() if ev.is_file() else None)  # fmt: skip
        calls.append(call)
        assert real is not None
        cap = _Capture()
        giap_log = logging.getLogger("giap")
        old = giap_log.level
        giap_log.addHandler(cap)
        giap_log.setLevel(logging.DEBUG)
        try:
            return real(outdir, cfg, *a, **k)
        finally:
            giap_log.removeHandler(cap)
            giap_log.setLevel(old)
            call.records = cap.records

    monkeypatch.setattr(P, "rerender_report", wrapper, raising=False)
    return calls


@pytest.fixture
def ok_run(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path, spy: list[Call]
) -> Iterator[Path]:
    _toy(registry)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.SUCCESS, log
    yield out


# ----------------------------------------------------------------------------- example 1


def test_example_1_the_published_page_has_the_full_provenance(ok_run: Path) -> None:
    import giap

    m = _manifest(ok_run)
    assert m["status"] == "completed"
    page = ok_run / "report" / "report.html"
    assert NOT_WRITTEN not in _section(page, "provenance"), "the page predates the manifest"
    rows = _prov(page)
    assert rows["GIAP version"] == m["giap_version"] == giap.__version__
    if m.get("git_commit"):
        assert rows["commit"].startswith(m["git_commit"])
    assert "--executor local" in rows["command line"]
    assert rows["seed"] == "5"
    assert rows[DIGEST_LABEL] == _sha(ok_run / "run_manifest.json")
    assert rows["provenance digest"] == m["provenance_digest"]


def test_example_1_the_work_copy_keeps_the_in_graph_page(ok_run: Path) -> None:
    work = ok_run / "work" / "report" / "report.html"
    assert work.is_file()
    assert NOT_WRITTEN in _section(work, "provenance")
    assert work.read_bytes() != (ok_run / "report" / "report.html").read_bytes()


# ----------------------------------------------------------------------------- run order


def test_rerender_runs_once_after_the_manifest_with_the_run_cfg(
    ok_run: Path, spy: list[Call]
) -> None:
    assert len(spy) == 1, "giap run did not call rerender_report in the tail"
    (call,) = spy
    assert Path(call.args[0]) == ok_run
    assert call.args[1].seed == 5, "the run's cfg"
    assert call.manifest is not None, "rerender_report ran before run_manifest.json was written"
    assert json.loads(call.manifest)["status"] == "completed"


def test_nothing_outside_logs_is_written_after_the_rerender(ok_run: Path, spy: list[Call]) -> None:
    assert len(spy) == 1, "giap run did not call rerender_report in the tail"
    (call,) = spy
    assert (ok_run / "run_manifest.json").read_bytes() == call.manifest, "manifest rewritten"
    assert _files(ok_run) == call.files


def test_the_rerender_is_not_a_step_and_emits_no_event(ok_run: Path, spy: list[Call]) -> None:
    assert len(spy) == 1, "giap run did not call rerender_report in the tail"
    ev = ok_run / "logs" / "events.jsonl"
    assert (ev.read_bytes() if ev.is_file() else None) == spy[0].events
    m = _manifest(ok_run)
    steps = {s["name"]: s["status"] for s in m["steps"]}
    assert steps.get("report") == "completed", "steps[] keeps the in-graph report step"
    assert set(steps) == {"x", "y", "plots", "report"}
    assert set(m["seeds"]) - {"manifest"} == {"master", "x", "y", "plots", "report"}


def test_the_rerender_logs_an_info_line_to_the_giap_logger(ok_run: Path, spy: list[Call]) -> None:
    assert len(spy) == 1, "giap run did not call rerender_report in the tail"
    info = [r for r in spy[0].records if r.levelno == logging.INFO and r.name.startswith("giap")]
    assert info, "rerender_report logged no INFO line to the giap logger"


# ----------------------------------------------------------------------------- example 3


def test_example_3_giap_report_right_after_gives_the_same_bytes(ok_run: Path) -> None:
    tail = (ok_run / "report" / "report.html").read_bytes()
    code, log = _invoke(["report", str(ok_run), "--set", "report.formats=[svg]",
                         "--set", "report.dpi=50"])  # fmt: skip
    assert code == ExitCode.SUCCESS, log
    assert (ok_run / "report" / "report.html").read_bytes() == tail


# ----------------------------------------------------------------------------- example 9


def _strings(obj: Any) -> Iterator[str]:
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield str(k)
            yield from _strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _strings(v)
    elif isinstance(obj, str):
        yield obj


def test_example_9_the_manifest_records_nothing_under_report_or_plots(ok_run: Path) -> None:
    for s in _strings(_manifest(ok_run)):
        assert "report.html" not in s, s
        assert not re.search(r"(^|/)plots/", s), s


# ----------------------------------------------------------------------------- example 5


def test_example_5_kegg_retrieval_date_end_to_end(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path, spy: list[Call]
) -> None:
    _toy(registry)
    out = tmp_path / "OUT"
    write_category_summary(out, {"cog": 0.5, "kegg": 0.4})  # OUT is never cleaned
    code, log = _invoke(_argv(env, out, modules="core,function"))
    assert code == ExitCode.SUCCESS, log
    kegg = [d for d in _manifest(out)["databases"] if str(d["name"]).startswith("kegg")]
    assert kegg and {d["retrieved"] for d in kegg} == {KEGG_RETRIEVED}, (
        "build_manifest does not fill databases[] from db_manifest.json (MVP-03c, #236)"
    )
    page = out / "report" / "report.html"
    assert _prov(page)["KEGG retrieval date"] == KEGG_RETRIEVED
    assert f"KEGG retrieval date: {KEGG_RETRIEVED}." in _section(page, "function")


# ----------------------------------------------------------------------------- example 6


def test_example_6_failed_graph_no_rerender(
    registry: ModuleType, env: dict[str, Path], tmp_path: Path, spy: list[Call]
) -> None:
    _toy(registry, boom=True)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.STEP_FAILED, log
    assert _manifest(out)["status"] == "failed"
    assert spy == []
    assert not (out / "logs" / "report_tail.log").exists()


# ----------------------------------------------------------------------------- example 7


def test_example_7_failed_manifest_write_no_rerender(
    registry: ModuleType,
    env: dict[str, Path],
    tmp_path: Path,
    spy: list[Call],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _toy(registry)
    import giap.steps.manifest as M

    def boom(*_a: Any, **_k: Any) -> None:
        raise OSError("manifest write failed on purpose")

    monkeypatch.setattr(M, "run_tail", boom)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert code == ExitCode.STEP_FAILED, log
    assert "step manifest failed" in log
    assert spy == []
    page, work = out / "report" / "report.html", out / "work" / "report" / "report.html"
    assert page.read_bytes() == work.read_bytes()
    assert not (out / "logs" / "report_tail.log").exists()


# ----------------------------------------------------------------------------- example 8


def test_example_8_failed_rerender_exits_4_and_keeps_the_in_graph_page(
    registry: ModuleType,
    env: dict[str, Path],
    tmp_path: Path,
    spy: list[Call],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import giap.report.build as B
    from giap.steps import plots as P

    _toy(registry)
    wrapper = P.rerender_report

    def boom(*_a: Any, **_k: Any) -> str:
        raise RuntimeError("build_report failed on purpose")

    def failing(outdir: Path, cfg: Any) -> Any:  # build_report raises only in the tail call
        monkeypatch.setattr(B, "build_report", boom)
        if hasattr(P, "build_report"):
            monkeypatch.setattr(P, "build_report", boom)
        return wrapper(outdir, cfg)

    monkeypatch.setattr(P, "rerender_report", failing)
    out = tmp_path / "OUT"
    code, log = _invoke(_argv(env, out))
    assert len(spy) == 1, "rerender_report was not called"
    assert code == ExitCode.STEP_FAILED, log
    tail_log = out / "logs" / "report_tail.log"
    assert f"giap run: report re-render failed; log: {tail_log}" in log
    assert tail_log.is_file() and "build_report failed on purpose" in tail_log.read_text()
    page, work = out / "report" / "report.html", out / "work" / "report" / "report.html"
    assert page.read_bytes() == work.read_bytes()
    assert (out / "run_manifest.json").read_bytes() == spy[0].manifest
    assert _manifest(out)["status"] == "completed"
    assert not (out / "report" / TMP_NAME).exists()
