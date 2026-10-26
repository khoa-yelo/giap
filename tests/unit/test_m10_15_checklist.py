"""Checklist hooks of ``metrics.py compute`` and the compute -> ``render --final`` path on the real
template (M10-15 acceptance criterion 1; examples 1, 5, 6, 9; edge case "A3 skipped by the
M10-08 verdict"; deliverable D1; review rev-1015 B1, S7, S9(a), N1, N2, N4).

Spec: docs/specs/M10-15-pilot-gate.md, "Checklist rows: status and evidence (D1)", amended by
"Lead decisions after the 15a code review": ``compute`` gains ``--baseline``,
``--benchmark-adr`` and ``--a3-feasibility`` but not ``--runbook`` (B1: ``metrics.json`` never
records the runbook); bullet 1's label and the S8-S14 job IDs are filled by ``render --final
--runbook``. Bullet 4 reads ``baseline.json`` (``standard_run.wheat``) and the ADR's ``Status:``
line, every row carries an ``evidence`` string, and ``inputs`` hashes the new files.

Inputs: M10-05's synthetic pilot (``helpers.pilot_metrics_m10``, giapsim SIM100K) plus the
hand-made runbook, baseline, ADR, verdict and TE files of ``helpers.pilot_gate_m10_15``. The keys
M10-07 (``m8``) and M10-10 (``m6``, ``h4``) will add are injected after ``compute`` where a test
needs a final report (``with_later_hooks``). Without the new options the M10-05 checklist is
unchanged (its own tests cover that).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_metrics_m10_10b import build_a3
from helpers.pilot_gate_m10_15 import (
    CELLS,
    EVIDENCE_JOBS,
    NOTES_OK,
    render_final,
    runbook_text,
    sha256,
    status_row,
    with_later_hooks,
    write_adr,
    write_baseline,
    write_feasibility,
    write_json,
    write_te_inputs,
    write_text,
)
from helpers.pilot_metrics_m10 import (
    TEMPLATE,
    Pilot,
    World,
    build_pilot,
    run_cli,
    sim_world,
)


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_15_sim"))


#: the simulated world of each pilot root, for the A3 inputs (``build_a3``, M10-10b)
_WORLDS: dict[Path, World] = {}


@pytest.fixture
def pilot(tmp_path: Path, world: World) -> Pilot:
    p = build_pilot(tmp_path / "pilot", world)
    _WORLDS[p.root] = world
    return p


def _inputs(
    p: Pilot,
    *,
    verdict: str | None = "pass",
    reason: str | None = None,
    cells: dict[str, str | None] | None = None,
    adr_status: str = "accepted",
) -> dict[str, Path]:
    """The committed inputs of ``compute`` (``verdict=None``: no ``--a3-feasibility``)."""
    d = p.root / "committed"
    files = {
        "baseline": write_baseline(d / "benchmarks" / "baseline.json", cells),
        "benchmark_adr": write_adr(
            d / "docs" / "adr" / "ADR-0099-benchmark-definition.md", adr_status
        ),
    }
    if verdict is not None:
        files["a3_feasibility"] = write_feasibility(d / "a3_feasibility.json", verdict, reason)
    return files


def _run(
    p: Pilot, out: Path, files: dict[str, Path], *, a3_exit: int | None = 0, te: bool = False
) -> Any:
    extra: list[str] = []
    for key, path in files.items():
        extra += [f"--{key.replace('_', '-')}", str(path)]
    if te:
        te_real, te7 = write_te_inputs(p.root / "te")
        extra += ["--te-real", str(te_real), "--te7", str(te7)]
    if a3_exit is not None:  # M10-10b: --a3 comes with its companions (build_a3)
        extra += build_a3(p, _WORLDS[p.root], exit_code=a3_exit).args()
    return run_cli(*p.compute_args(out, *extra), cwd=p.root)


def _compute(
    p: Pilot, out: Path, files: dict[str, Path], *, a3_exit: int | None = 0, te: bool = False
) -> Any:
    res = _run(p, out, files, a3_exit=a3_exit, te=te)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _status(doc: dict[str, Any], i: int) -> str:
    return str(doc["checklist"][i]["status"])


def _evidence(doc: dict[str, Any], i: int) -> str:
    return str(doc["checklist"][i]["evidence"])


# ------------------------------------------------------------------ B1: compute has no runbook


def test_b1_compute_refuses_runbook(pilot: Pilot, tmp_path: Path) -> None:
    """B1: ``--runbook`` is an input of ``render --final`` and of the packet, not of
    ``compute`` (argparse exit 2)."""
    rb = write_text(pilot.root / "committed" / "reproduce.md", runbook_text())
    files = _inputs(pilot) | {"runbook": rb}
    res = _run(pilot, tmp_path / "m", files)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stderr}"
    assert not (tmp_path / "m" / "metrics.json").exists()


def test_b1_metrics_json_never_records_the_runbook(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    assert "runbook" not in doc["inputs"]
    dumped = json.dumps(doc)
    for step, (_, job) in EVIDENCE_JOBS.items():
        assert job not in dumped, step
    assert not (doc["evidence"].get("steps") or {}), doc["evidence"].get("steps")


@pytest.mark.parametrize(("a3_exit", "expected"), [(0, 0), (None, None)])
def test_evidence_records_the_arm_exit_codes(
    pilot: Pilot, tmp_path: Path, a3_exit: int | None, expected: int | None
) -> None:
    """``render --final`` labels bullet 1 from these and the runbook (helper docstring)."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot), a3_exit=a3_exit)
    assert doc["evidence"]["exit_codes"] == {"A1": 0, "A2": 0, "A3": expected}


@pytest.mark.parametrize(
    ("a3_exit", "expected"),
    [(0, "pass"), (None, "partial (A1, A2; A3 pending M10-10)")],
)
def test_compute_bullet_1_keeps_m10_05_strings(
    pilot: Pilot, tmp_path: Path, a3_exit: int | None, expected: str
) -> None:
    """Lead decision (15a tests): without a runbook bullet 1 keeps M10-05's strings."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot), a3_exit=a3_exit)
    assert _status(doc, 0) == expected


def test_a3_run_that_failed_is_a_precondition(pilot: Pilot, tmp_path: Path) -> None:
    """M10-10b: an A3 run with a non-zero exit is refused before compute (exit 2, nothing
    written); it was a computed ``fail`` in 15a."""
    out = tmp_path / "m"
    res = _run(pilot, out, _inputs(pilot), a3_exit=1)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stderr}"
    assert "precondition:" in res.stderr, res.stderr
    assert not (out / "metrics.json").exists()


# ------------------------------------------------------------------ S7: absent feasibility


def test_s7_absent_feasibility_is_a_pending_reason(pilot: Pilot, tmp_path: Path) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, verdict=None))
    feas = doc["evidence"]["a3_feasibility"]
    assert feas.get("verdict") is None
    assert feas.get("reason") == "pending (M10-08 verdict)", feas


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


@pytest.mark.parametrize("status", ["rejected", "superseded"])
def test_n1_rejected_or_superseded_adr_fails(pilot: Pilot, tmp_path: Path, status: str) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, adr_status=status))
    assert _status(doc, 3) == f"fail (ADR-M10-D {status})"


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
    assert "m4" in _evidence(doc, 2) and "m5" in _evidence(doc, 2)
    assert "te_real" in _evidence(doc, 4)


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


# ------------------------------------------------------------------ S9(a): compute, render --final


def _render(
    tmp: Path,
    doc: dict[str, Any],
    *,
    a3_branch: str,
    steps: list[str] | None = None,
    final: bool = True,
) -> tuple[Any, Path]:
    m = write_json(tmp / "r" / "metrics.json", doc)
    rb = write_text(tmp / "r" / "reproduce.md", runbook_text(a3_branch=a3_branch, steps=steps))
    notes = write_text(tmp / "r" / "report_notes.md", NOTES_OK)
    out = tmp / "r" / "out" / "REPORT.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    res = render_final(m, out, notes=notes, template=TEMPLATE, runbook=rb, final=final)
    return res, out


def _final_ok(res: Any, out: Path) -> str:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    text = out.read_text(encoding="utf-8")
    assert "**missing:" not in text and "{{" not in text
    assert not re.search(r"\bNA \(reason: NA\)", text), "a bare null reached the final report"
    assert re.search(r"runbook status \**final", text), "section 9 names the final runbook"
    return text


def _status_cells(text: str) -> list[str]:
    return status_row(text, "GIAP runs end to end")


def test_s9a_example_1_compute_then_render_final_on_the_real_template(
    pilot: Pilot, tmp_path: Path
) -> None:
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot), te=True)
    res, out = _render(tmp_path, with_later_hooks(doc, None), a3_branch="cascadia")
    text = _final_ok(res, out)
    cells = _status_cells(text)
    assert cells[1] == "pass (A3 = Cascadia)", cells
    for step in ("S8", "S9", "S11"):
        assert EVIDENCE_JOBS[step][1] in " ".join(cells[2:]), (step, cells)


def test_s9a_example_6_no_a3_run(pilot: Pilot, tmp_path: Path) -> None:
    """Example 6 (N2: no "; ADR-M10-C declined"), S1: M10-10's ``a3_not_run`` is a finding."""
    files = _inputs(pilot, verdict="fail", reason="route_unavailable")
    doc = _compute(pilot, tmp_path / "m", files, a3_exit=None, te=True)
    steps = [s for s in EVIDENCE_JOBS if s != "S11"]
    res, out = _render(
        tmp_path, with_later_hooks(doc, "a3_not_run"), a3_branch="fallback", steps=steps
    )
    text = _final_ok(res, out)
    assert _status_cells(text)[1] == "fail (A3 not run: route_unavailable)"
    assert "not evaluable (a3_not_run)" in text
    section = re.split(r"^## \d+\. Not evaluated.*$", text, flags=re.M)[1].split("\n## ", 1)[0]
    assert "m6" in section and "h4" in section, section
    verdict = next(ln for ln in text.splitlines() if "M10-08 verdict" in ln)
    assert "fail" in verdict and "route_unavailable" in verdict, verdict


def test_s9a_example_5_fallback(pilot: Pilot, tmp_path: Path) -> None:
    files = _inputs(pilot, verdict="fail", reason="low_genome_match")
    doc = _compute(pilot, tmp_path / "m", files, te=True)
    res, out = _render(tmp_path, with_later_hooks(doc, "different_strain"), a3_branch="fallback")
    text = _final_ok(res, out)
    cells = _status_cells(text)
    assert cells[1] == "pass (A3 = Casanovo fallback on UTI89, different strain)", cells
    assert EVIDENCE_JOBS["S11"][1] in " ".join(cells[2:]), cells
    assert "not evaluable (different_strain)" in text


def test_s7_compute_without_feasibility_is_refused_by_render_final(
    pilot: Pilot, tmp_path: Path
) -> None:
    """S7 (probe P5): no "M10-08 verdict NA (reason: NA)" in a final report."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot, verdict=None), te=True)
    res, out = _render(tmp_path, with_later_hooks(doc, None), a3_branch="cascadia")
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stderr}"
    assert "M10-08" in res.stderr + res.stdout, res.stderr
    assert not out.exists()


def test_compute_alone_is_refused_by_render_final_until_m10_07_and_m10_10(
    pilot: Pilot, tmp_path: Path
) -> None:
    """``m8`` is a pending placeholder until M10-07: a defect, never a finding. Since M10-10b
    ``m6``/``h4`` are computed, so only M10-07 is named."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot), te=True)
    res, out = _render(tmp_path, doc, a3_branch="cascadia")
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stderr}"
    assert "M10-07" in res.stderr + res.stdout, res.stderr
    assert not out.exists()


def test_n4_draft_render_shows_the_pending_placeholders(pilot: Pilot, tmp_path: Path) -> None:
    """N4: the draft of compute's output (no ``--runbook``) renders ``m8`` with its pending reason,
    never a bare ``NA``; ``m6`` and ``h4`` are computed since M10-10b."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    m = write_json(tmp_path / "d" / "metrics.json", doc)
    tpl = write_text(
        tmp_path / "d" / "t.md",
        "M6: {{m6.overlap_a2_a3.novel_orfs.called.jaccard}}\nH4: {{h4.verdict}}\n"
        "M8: {{m8.peptide.fdp_combined}}\n",
    )
    out = tmp_path / "d" / "out.md"
    res = render_final(m, out, template=tpl, final=False)
    assert res.returncode == 0, res.stderr
    lines = dict(ln.split(": ", 1) for ln in out.read_text(encoding="utf-8").splitlines() if ln)
    assert "M10-07" in lines["M8"] and not lines["M8"].startswith("NA"), lines["M8"]
    for key in ("M6", "H4"):
        assert "pending" not in lines[key] and not lines[key].startswith("NA"), (key, lines[key])


def test_draft_render_of_the_real_template_needs_no_runbook(pilot: Pilot, tmp_path: Path) -> None:
    """Draft mode keeps M10-05's behaviour: every template path (``evidence.runbook.*`` too)
    resolves without ``--runbook``."""
    doc = _compute(pilot, tmp_path / "m", _inputs(pilot))
    m = write_json(tmp_path / "d" / "metrics.json", doc)
    out = tmp_path / "d" / "REPORT.md"
    res = render_final(m, out, template=TEMPLATE, final=False)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stderr}"
