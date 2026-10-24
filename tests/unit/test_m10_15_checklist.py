"""Checklist hooks of ``metrics.py compute`` (M10-15 acceptance criterion 1; examples 1, 5, 6, 9;
edge case "A3 skipped by the M10-08 verdict"; deliverable D1).

Spec: docs/specs/M10-15-pilot-gate.md, "Checklist rows: status and evidence (D1)". ``compute``
gains ``--baseline``, ``--benchmark-adr``, ``--a3-feasibility`` and ``--runbook``; bullet 1 names
the A3 branch, bullet 4 reads ``baseline.json`` (``standard_run.wheat``) and the ADR's ``Status:``
line, every row carries an ``evidence`` string, the ``metrics.json`` key ``evidence`` copies
the runbook's Evidence rows and ``inputs`` hashes the four new files.

Inputs: M10-05's synthetic pilot (``helpers.pilot_metrics_m10``, giapsim SIM100K) plus the
hand-made runbook, baseline, ADR and verdict of ``helpers.pilot_gate_m10_15`` (completeness rule
for a benchmark cell in that module's docstring). Without the new options the M10-05 checklist
is unchanged (its own tests cover that).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import (
    CELLS,
    EVIDENCE_JOBS,
    runbook_text,
    sha256,
    write_adr,
    write_baseline,
    write_feasibility,
    write_text,
)
from helpers.pilot_metrics_m10 import Pilot, World, build_pilot, run_cli, sim_world, write_a3


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_15_sim"))


@pytest.fixture
def pilot(tmp_path: Path, world: World) -> Pilot:
    return build_pilot(tmp_path / "pilot", world)


def _inputs(
    p: Pilot,
    *,
    verdict: str = "pass",
    reason: str | None = None,
    a3_branch: str = "cascadia",
    cells: dict[str, str | None] | None = None,
    adr_status: str = "accepted",
    runbook_status: str = "final",
) -> dict[str, Path]:
    d = p.root / "committed"
    steps = [s for s in EVIDENCE_JOBS if s != "S11" or a3_branch in ("cascadia", "fallback")]
    return {
        "baseline": write_baseline(d / "benchmarks" / "baseline.json", cells),
        "adr": write_adr(d / "docs" / "adr" / "ADR-0099-benchmark-definition.md", adr_status),
        "feasibility": write_feasibility(d / "a3_feasibility.json", verdict, reason),
        "runbook": write_text(
            d / "reproduce.md",
            runbook_text(status=runbook_status, a3_branch=a3_branch, steps=steps),
        ),
    }


def _compute(p: Pilot, out: Path, files: dict[str, Path], *, a3_exit: int | None = 0) -> Any:
    extra = [
        "--baseline",
        str(files["baseline"]),
        "--benchmark-adr",
        str(files["adr"]),
        "--a3-feasibility",
        str(files["feasibility"]),
        "--runbook",
        str(files["runbook"]),
    ]
    if a3_exit is not None:
        extra += ["--a3", str(write_a3(p.root, sha256(p.genome), exit_code=a3_exit))]
    res = run_cli(*p.compute_args(out, *extra), cwd=p.root)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _status(doc: dict[str, Any], i: int) -> str:
    return str(doc["checklist"][i]["status"])


def _evidence(doc: dict[str, Any], i: int) -> str:
    return str(doc["checklist"][i]["evidence"])


# ------------------------------------------------------------------ bullet 1


def test_example_1_bullet_1_cascadia(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    assert _status(doc, 0) == "pass (A3 = Cascadia)"


def test_example_5_bullet_1_fallback(pilot: Pilot, tmp_path: Path) -> None:
    files = _inputs(pilot, verdict="fail", reason="low_genome_match", a3_branch="fallback")
    doc = _compute(pilot, tmp_path / "m", files)
    assert _status(doc, 0) == "pass (A3 = Casanovo fallback on UTI89, different strain)"


def test_example_6_bullet_1_no_a3_run(pilot: Pilot, tmp_path: Path) -> None:
    files = _inputs(pilot, verdict="fail", reason="route_unavailable", a3_branch="fallback")
    doc = _compute(pilot, tmp_path / "m", files, a3_exit=None)
    st = _status(doc, 0)
    assert st.startswith("fail (A3 not run: route_unavailable"), st


def test_bullet_1_runbook_pending_stays_partial(pilot: Pilot, tmp_path: Path) -> None:
    """Edge case: ``a3_branch=pending`` and no A3 run: bullet 1 stays ``partial`` (which
    ``render --final`` rejects), never a ``fail (A3 not run …)``."""
    files = _inputs(pilot, a3_branch="pending", runbook_status="draft")
    doc = _compute(pilot, tmp_path / "m", files, a3_exit=None)
    assert _status(doc, 0).startswith("partial"), _status(doc, 0)


def test_bullet_1_a3_exit_nonzero_names_the_arm(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot), a3_exit=1)
    assert _status(doc, 0) == "fail (A3 exit 1)"


def test_bullet_1_evidence_names_the_runbook_jobs(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    ev = _evidence(doc, 0)
    for step in ("S8", "S9", "S11"):
        assert EVIDENCE_JOBS[step][1] in ev, (step, ev)


# ------------------------------------------------------------------ bullet 4


@pytest.mark.parametrize("status", ["accepted", "accepted-by-khoa"])
def test_bullet_4_four_complete_cells_and_accepted_adr_pass(
    pilot: Pilot, tmp_path: Path, status: str
) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, adr_status=status))
    assert _status(doc, 3) == "pass"


def test_example_9_bullet_4_proposed_adr_is_pending(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, adr_status="proposed"))
    assert _status(doc, 3) == "pending (ADR-M10-D proposed)"


def test_example_9_incomplete_cell_wins_over_the_adr_status(pilot: Pilot, tmp_path: Path) -> None:
    cells: dict[str, str | None] = {c: "complete" for c in CELLS} | {"sim5m_core": "incomplete"}
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, cells=cells, adr_status="proposed"))
    assert _status(doc, 3) == "fail (incomplete: sim5m_core)"


def test_bullet_4_missing_cell_fails(pilot: Pilot, tmp_path: Path) -> None:
    cells: dict[str, str | None] = {c: "complete" for c in CELLS} | {"k12_core": None}
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, cells=cells))
    assert _status(doc, 3) == "fail (incomplete: k12_core)"


def test_bullet_4_no_standard_run_names_every_cell(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, cells={}))
    st = _status(doc, 3)
    assert st.startswith("fail (incomplete:"), st
    assert all(c in st for c in CELLS), st


def test_bullet_4_evidence_names_baseline_and_adr(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    ev = _evidence(doc, 3)
    assert "baseline.json" in ev and "standard_run.wheat" in ev, ev
    assert "ADR-0099" in ev and "accepted" in ev, ev


# ------------------------------------------------------------------ every row: evidence


def test_every_row_carries_an_evidence_string_naming_a_file_and_key(
    pilot: Pilot, tmp_path: Path
) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    assert len(doc["checklist"]) == 5
    for i, row in enumerate(doc["checklist"]):
        ev = str(row.get("evidence") or "")
        assert ev.strip(), f"row {i} has no evidence"
        assert re.search(r"[\w-]+\.(json|md)\b", ev), f"row {i} names no file: {ev!r}"
    assert "te1.share_exact_inframe" in _evidence(doc, 1) and "te1.passed" in _evidence(doc, 1)
    assert EVIDENCE_JOBS["S14"][1] in _evidence(doc, 1)
    assert "m4" in _evidence(doc, 2) and "m5" in _evidence(doc, 2)
    assert EVIDENCE_JOBS["S12"][1] in _evidence(doc, 2)
    assert "te_real" in _evidence(doc, 4)


def test_metrics_json_evidence_copies_the_runbook_rows(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    assert "evidence" in doc
    dumped = json.dumps(doc["evidence"])
    for step in ("S8", "S9", "S10", "S11", "S12", "S13", "S14", "S-bench"):
        assert EVIDENCE_JOBS[step][1] in dumped, step


def test_the_new_inputs_are_hashed(pilot: Pilot, tmp_path: Path) -> None:
    files = _inputs(pilot)
    doc = _compute(pilot, tmp_path / "m", files)
    hashes = {v["sha256"] for v in doc["inputs"].values()}
    for name, path in files.items():
        assert sha256(path) in hashes, name


def test_compute_with_the_new_options_is_byte_identical(pilot: Pilot, tmp_path: Path) -> None:
    files = _inputs(pilot)
    _compute(pilot, tmp_path / "m1", files)
    _compute(pilot, tmp_path / "m2", files)
    a = (tmp_path / "m1" / "metrics.json").read_bytes()
    assert a == (tmp_path / "m2" / "metrics.json").read_bytes()
