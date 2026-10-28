"""PR #205 review round B6 (reviewer N1): run events reach the ``plots`` and ``report`` steps in
``gippa run --executor local`` as they reach ``gippa report`` (``OUT/logs``), so a plot declaring
``layers=[...]`` skips with ``layer_missing:<layer>`` in both (M9-01 behaviour 5, M9-02
behaviour 1, 2, 8; GIPPA_SPEC.md §6.4 executor equality, §6.7 events).

A swapped step registry holds the real ``plots`` and ``report`` steps and one toy step that
logs the run flag ``layer_missing:fakelayer`` (as ``context_layers`` / ``novel`` do) and writes
``pairs`` (an optional input of ``plots``, so it runs first). A swapped plot registry holds a fake
V11 declaring ``layers=["fakelayer"]`` and a plain fake V07. Inputs as
tests/unit/test_steps_review.py (a toy genome, one peptide, a test ``$GIPPA_DB_DIR``).
"""

from __future__ import annotations

import json
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.report_fixtures import figure, frame, parse, register_fake, swap_registry
from helpers.viz_out_tree import write_schema_table

from gippa.exitcodes import ExitCode

FLAG = "layer_missing:fakelayer"


@pytest.fixture(autouse=True)
def _test_db(tmp_path_factory: pytest.TempPathFactory, monkeypatch: pytest.MonkeyPatch) -> None:
    from helpers.mapping_steps import write_test_db

    db = write_test_db(tmp_path_factory.mktemp("test_db"))
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))


@pytest.fixture
def steps(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Iterator[ModuleType]:
    """The real ``plots`` / ``report`` specs plus a toy flag-emitting step, nothing else."""
    import gippa.steps
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    keep = {n: reg.get(n) for n in ("plots", "report")}
    monkeypatch.setattr(reg, "_REGISTRY", dict(keep))
    empty = tmp_path / "steps_pkg"
    empty.mkdir()
    monkeypatch.setattr(gippa.steps, "__path__", [str(empty)])

    def toy_flags(ctx: Any) -> None:
        data = {"layer": "fakelayer", "flag": FLAG}
        ctx.log.warning("%s: no fakelayer layer", FLAG,
                        extra={"event": "toy_flags.layer_missing", "data": data})  # fmt: skip
        base = {"class": "adjacent", "gene_a": "g1", "gene_b": "g2", "strand": "+", "test": "none"}
        write_schema_table(ctx.outputs["pairs"], "pairs", [{**base, "pair_id": "P1"}])

    reg.step("toy_flags", inputs=[], outputs=["pairs"])(toy_flags)
    yield reg


@pytest.fixture
def plots(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    swap_registry(monkeypatch)
    from gippa.viz import base

    calls: list[str] = []

    @base.plot("V11", "operon_pairs", layers=["fakelayer"])
    def _v11(view: Any, cfg: Any) -> Any:
        calls.append("V11")
        return [base.Instance(suffix=None, figure=figure(), data=frame(), title="fake")]

    register_fake("V07")
    return calls


def _genome(tmp_path: Path) -> Path:
    path = tmp_path / "g.fa"
    path.write_text(">chr1 [topology=circular]\n" + "ACGT" * 50 + "\n")
    return path


def _exit_code(argv: list[str]) -> int:
    from gippa.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    return 0 if code is None else int(code)


def _run(tmp_path: Path, out: Path) -> None:
    genome = _genome(tmp_path)
    pep = tmp_path / "evidence_peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n")
    argv = ["run", "--genome", str(genome), "--peptides", str(pep), "--outdir", str(out),
            "--executor", "local", "--set", "report.formats=[svg]",
            "--set", "report.dpi=50"]  # fmt: skip
    assert _exit_code(argv) == ExitCode.SUCCESS


def _slot(out: Path, pid: str) -> str:
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    return " ".join(page.slot(pid).text().split())


def _skips(out: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    path = out / "logs" / "events.jsonl"
    for line in path.read_text(encoding="utf-8").splitlines():
        rec = json.loads(line)
        if rec.get("event") == "plots.skipped":
            found[rec["data"]["plot"]] = rec["data"]["reason"]
    return found


def test_b6_layer_missing_skip_fires_in_gippa_run_local(
    steps: ModuleType, plots: list[str], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    _run(tmp_path, out)
    assert FLAG in (out / "logs" / "events.jsonl").read_text(encoding="utf-8"), "toy step not run"
    assert list((out / "plots").glob("V07_*.svg")), "the plain fake plot did not render"
    assert plots == [], "V11 was drawn although its layer is missing in this run"
    left = sorted(p.name for p in (out / "plots").iterdir())
    assert not [n for n in left if n.startswith("V11_")], left
    assert _skips(out).get("V11") == FLAG, _skips(out)
    text = _slot(out, "V11")
    assert "not run" in text and "fakelayer" in text, text


def test_b6_gippa_run_and_gippa_report_agree_on_the_layer_missing_skip(
    steps: ModuleType, plots: list[str], tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    _run(tmp_path, out)
    in_run = _slot(out, "V11")
    code = _exit_code(["report", str(out), "--set", "report.formats=[svg]",
                       "--set", "report.dpi=50"])  # fmt: skip
    assert code == ExitCode.SUCCESS
    in_report = _slot(out, "V11")
    assert "fakelayer" in in_report, in_report
    assert in_run == in_report, (in_run, in_report)
