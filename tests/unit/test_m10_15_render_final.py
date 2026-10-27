"""``metrics.py render --final`` (M10-15 acceptance criterion 3; examples 2-8; edge cases "A failed
TE", "A3 skipped" (runbook pending -> ``partial``), "A missing metric", "Empty input"; the
byte-identical re-render behind criterion 2; review rev-1015 B1, B2, S1, S2, S6, S9(b), N2).

Spec: docs/specs/M10-15-pilot-gate.md, "The final report", "Checklist rows", "Examples", "Edge
cases" and "Lead decisions after the 15a code review" (binding: ``--runbook`` is an input of
``render --final``, not of ``compute``; bullet 1 and the S8-S14 job IDs are computed at render
time from the final runbook). Inputs are hand-made (``tests/helpers/pilot_gate_m10_15.py``: a
small final-mode template passed with ``--template``, a minimal ``metrics.json`` and a runbook in
M10-14's format); the interface pinned beyond the spec is listed in that module's docstring.
Nothing comes from ``$GIPPA_ROOT``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import (
    DRAFT_BULLET1,
    EVIDENCE_JOBS,
    FIXTURE_TEMPLATE,
    NOTES_OK,
    OTHER_COMMIT,
    diagnosis,
    final_metrics,
    render_final,
    runbook_text,
    status_row,
    write_json,
    write_text,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

#: the fixture template without its TE table (an absent TE status is then not an absent path)
NO_TE_TABLE = (
    FIXTURE_TEMPLATE.split("## 5. TE2-TE5 and TE7")[0] + "## 10. Notes\n\n<!-- notes -->\n"
)


def _setup(
    tmp: Path,
    doc: dict[str, Any],
    notes: str | None = NOTES_OK,
    template: str = FIXTURE_TEMPLATE,
    runbook: str | None = None,
) -> tuple[Path, Path | None, Path, Path | None]:
    m = write_json(tmp / "in" / "metrics.json", doc)
    n = None if notes is None else write_text(tmp / "in" / "report_notes.md", notes)
    t = write_text(tmp / "in" / "report_template.md", template)
    r = None if runbook is None else write_text(tmp / "in" / "reproduce.md", runbook)
    return m, n, t, r


def _final(
    tmp: Path,
    doc: dict[str, Any],
    notes: str | None = NOTES_OK,
    template: str = FIXTURE_TEMPLATE,
    final: bool = True,
    runbook: str | None = "default",
) -> tuple[Any, Path]:
    """Renders ``doc``; ``runbook`` is the runbook text (``"default"``: the final Cascadia
    runbook with every Evidence row; ``None``: no ``--runbook``)."""
    if runbook == "default":
        runbook = runbook_text()
    m, n, t, r = _setup(tmp, doc, notes, template, runbook)
    out = tmp / "out" / "REPORT.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    return render_final(m, out, notes=n, template=t, runbook=r, final=final), out


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _refused(res: Any, out: Path, code: int = 1, *names: str) -> None:
    assert res.returncode == code, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "unrecognized arguments" not in res.stderr, res.stderr  # a usage error is not a check
    assert not out.exists(), "nothing is written on a refused render"
    for name in names:
        assert name in res.stderr + res.stdout, f"{name!r} not named in the message"


def _bullet1(out: Path) -> list[str]:
    return status_row(out.read_text(encoding="utf-8"), "runs")


# ------------------------------------------------------------------ example 1: all pass


def test_example_1_all_pass_renders_exit_0(tmp_path: Path) -> None:
    res, out = _final(tmp_path, final_metrics())
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert _bullet1(out)[1] == "pass (A3 = Cascadia)"
    assert "Submitted for Khoa's review" in text
    assert "**missing:" not in text and "{{" not in text


def test_final_render_is_byte_identical_when_re_rendered(tmp_path: Path) -> None:
    res1, out1 = _final(tmp_path / "a", final_metrics(), notes=NOTES_OK + diagnosis("TE3"))
    res2, out2 = _final(tmp_path / "b", final_metrics(), notes=NOTES_OK + diagnosis("TE3"))
    _ok(res1)
    _ok(res2)
    assert out1.read_bytes() == out2.read_bytes()


# ------------------------------------------------------------------ B1: runbook input


def test_final_without_runbook_is_a_precondition_exit_2(tmp_path: Path) -> None:
    """B1: bullet 1 and the S8-S14 evidence come from the final runbook, so ``--final`` needs
    ``--runbook``."""
    res, out = _final(tmp_path, final_metrics(), runbook=None)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stderr}"
    assert not out.exists()


def test_final_with_a_draft_runbook_exit_1(tmp_path: Path) -> None:
    """B1(b): "runbook status draft" never reaches a final report."""
    res, out = _final(tmp_path, final_metrics(), runbook=runbook_text(status="draft"))
    _refused(res, out, 1, "runbook")


def test_runbook_status_is_rendered_from_the_runbook(tmp_path: Path) -> None:
    res, out = _final(tmp_path, final_metrics())
    _ok(res)
    assert "runbook status final" in out.read_text(encoding="utf-8")


def test_evidence_strings_carry_the_runbook_job_ids(tmp_path: Path) -> None:
    """B1(a): bullet 1 names S8, S9, S11; bullet 2 the S14 job; bullet 3 the S12 job (spec
    "Checklist rows"), read from the runbook at render time (``compute`` has no job ID)."""
    res, out = _final(tmp_path, final_metrics())
    _ok(res)
    report = out.read_text(encoding="utf-8")
    rows = {k: " ".join(status_row(report, k)[2:]) for k in ("runs", "TE1", "compare")}
    for step in ("S8", "S9", "S11"):
        assert EVIDENCE_JOBS[step][1] in rows["runs"], (step, rows["runs"])
    assert EVIDENCE_JOBS["S14"][1] in rows["TE1"], rows["TE1"]
    assert EVIDENCE_JOBS["S12"][1] in rows["compare"], rows["compare"]


@pytest.mark.parametrize("first", [False, True])
def test_b2_superseded_evidence_rows_are_never_chosen(tmp_path: Path, first: bool) -> None:
    """B2: a superseded S8 attempt (another commit, exit 1) is ignored wherever it stands in
    the table; the accepted job is the evidence."""
    rb = runbook_text(superseded=[("S8", "9000001", OTHER_COMMIT)], superseded_first=first)
    res, out = _final(tmp_path, final_metrics(), runbook=rb)
    _ok(res)
    ev = " ".join(_bullet1(out)[2:])
    assert EVIDENCE_JOBS["S8"][1] in ev and "9000001" not in ev, ev


# ------------------------------------------------------------------ bullet 1 at render time (S2)


@pytest.mark.parametrize(
    ("branch", "a3_exit", "verdict", "reason", "expected"),
    [
        ("cascadia", 0, "pass", None, "pass (A3 = Cascadia)"),
        (
            "fallback",
            0,
            "fail",
            "low_genome_match",
            "pass (A3 = Casanovo fallback on UTI89, different strain)",
        ),
        # example 6, with N2: no "; ADR-M10-C declined" (the code cannot know it)
        ("fallback", None, "fail", "route_unavailable", "fail (A3 not run: route_unavailable)"),
        ("cascadia", 1, "pass", None, "fail (A3 exit 1)"),
        # S2: a Cascadia branch against a failed feasibility verdict is not labelled Cascadia
        ("cascadia", 0, "fail", "low_genome_match", "fail (A3 feasibility fail)"),
    ],
)
def test_bullet_1_label_comes_from_the_runbook_branch(
    tmp_path: Path,
    branch: str,
    a3_exit: int | None,
    verdict: str,
    reason: str | None,
    expected: str,
) -> None:
    doc = final_metrics(a3_exit=a3_exit, verdict=verdict, reason=reason)
    steps = [s for s in EVIDENCE_JOBS if s != "S11" or a3_exit is not None]
    res, out = _final(tmp_path, doc, runbook=runbook_text(a3_branch=branch, steps=steps))
    _ok(res)
    assert _bullet1(out)[1] == expected


def test_bullet_1_arm_failure_wins(tmp_path: Path) -> None:
    doc = final_metrics()
    doc["evidence"]["exit_codes"]["A1"] = 3
    res, out = _final(tmp_path, doc)
    _ok(res)
    assert _bullet1(out)[1] == "fail (A1 exit 3)"


@pytest.mark.parametrize(("a3_exit", "verdict"), [(None, "pass"), (0, "pass"), (0, None)])
def test_s2_branch_pending_is_partial_and_refused(
    tmp_path: Path, a3_exit: int | None, verdict: str | None
) -> None:
    """S2 (probe P2): with ``a3_branch=pending`` the branch is never guessed from the A3 exit
    code or the feasibility verdict; bullet 1 is ``partial``, which ``--final`` refuses (edge
    case "A3 skipped": the runbook still at ``pending``)."""
    doc = final_metrics(a3_exit=a3_exit, verdict=verdict)
    res, out = _final(tmp_path, doc, runbook=runbook_text(a3_branch="pending"))
    _refused(res, out, 1, "checklist.0.status", "partial")


def test_bullet_1_status_in_metrics_json_is_relabelled_not_trusted(tmp_path: Path) -> None:
    """The checklist row 0 that ``compute`` wrote (M10-05's string) is replaced at render time:
    a hand-edited ``pass (A3 = Cascadia)`` does not survive a pending runbook."""
    doc = final_metrics(statuses=["pass (A3 = Cascadia)", "pass", "pass", "pass", "pass"])
    res, out = _final(tmp_path, doc, runbook=runbook_text(a3_branch="pending"))
    _refused(res, out, 1)


# ------------------------------------------------------------------ pending statuses


@pytest.mark.parametrize(
    ("bullet", "status"),
    [
        (3, "pending (ADR-M10-D proposed)"),  # example 9, second half
        (3, "pending (M10-12)"),
        (4, "pending (M10-06)"),
    ],
)
def test_pending_or_partial_checklist_status_exit_1(
    tmp_path: Path, bullet: int, status: str
) -> None:
    statuses = ["pass", "pass", "pass", "pass", "pass"]
    statuses[bullet] = status
    res, out = _final(tmp_path, final_metrics(statuses=statuses))
    _refused(res, out, 1, f"checklist.{bullet}.status")


@pytest.mark.parametrize("te", ["te2", "te5"])
def test_pending_te_status_exit_1(tmp_path: Path, te: str) -> None:
    res, out = _final(tmp_path, final_metrics(te={te: "pending (M10-06)"}))
    _refused(res, out, 1)


def test_pending_te7_status_exit_1(tmp_path: Path) -> None:
    res, out = _final(tmp_path, final_metrics(te7="pending (M10-06)"))
    _refused(res, out, 1)


def _te_placeholder_shape(doc: dict[str, Any]) -> dict[str, Any]:
    """``te_real``/``te7`` as ``te_placeholders()`` leaves them (no ``schema`` key): the
    diagnosis rule ignores them, so only ``final_problems`` can refuse (S9(b))."""
    p = {"status": "pending (M10-06)"}
    doc["te_real"] = {k: dict(p) for k in ("te2", "te3", "te4", "te5")}
    doc["te7"] = {"gippa": dict(p), "diann": dict(p), "end_to_end": None}
    return doc


def test_s9b_te_placeholders_are_refused_by_final_problems(tmp_path: Path) -> None:
    res, out = _final(tmp_path, _te_placeholder_shape(final_metrics()))
    _refused(res, out, 1, "te_real.te2.status", "te7.gippa.status")


def test_s9b_te7_placeholder_alone_is_refused(tmp_path: Path) -> None:
    """Mutant "TE7 not checked in final_problems": TE2-TE5 pass, only TE7 is pending."""
    doc = final_metrics()
    doc["te7"] = {"gippa": {"status": "pending (M10-06)"}, "diann": {"status": "pending (M10-06)"}}
    res, out = _final(tmp_path, doc)
    _refused(res, out, 1, "te7.gippa.status")


@pytest.mark.parametrize(
    ("key", "path"), [("te5", "te_real.te5.status"), ("te7", "te7.gippa.status")]
)
def test_s9b_absent_te_status_is_refused(tmp_path: Path, key: str, path: str) -> None:
    """Mutant "absent TE status accepted": the status key is gone, the template does not name
    it and the notes hold the diagnosis, so only ``final_problems`` sees it."""
    doc = final_metrics()
    if key == "te7":
        doc["te7"] = {"gippa": {}, "diann": {"status": "pass"}}
    else:
        doc["te_real"][key] = {}
    notes = NOTES_OK + "\n" + diagnosis(key.upper())
    res, out = _final(tmp_path, doc, notes=notes, template=NO_TE_TABLE)
    _refused(res, out, 1, path)


def test_draft_mode_keeps_m10_05_behaviour_on_pending(tmp_path: Path) -> None:
    """Criterion 3: without ``--final`` (and without ``--runbook``) a pending status is
    rendered (M10-05's draft)."""
    doc = final_metrics(
        statuses=[
            "partial (A1, A2; A3 pending M10-10)",
            "pass",
            "pass",
            "pending (M10-12)",
            "pending (M10-06)",
        ],
        a3_exit=None,
    )
    doc["evidence"]["runbook"] = {"reason": "pending (render --final --runbook)"}
    res, out = _final(tmp_path, doc, final=False, runbook=None)
    _ok(res)
    assert "pending (M10-12)" in out.read_text(encoding="utf-8")


# ------------------------------------------------------------------ S6: check 4 in final_problems


@pytest.mark.parametrize("value", ["", "   "])
def test_s6_empty_evidence_string_exit_1(tmp_path: Path, value: str) -> None:
    """Check 4 of the packet, folded into ``final_problems`` (render enforces it)."""
    doc = final_metrics()
    doc["checklist"][3]["evidence"] = value
    res, out = _final(tmp_path, doc)
    _refused(res, out, 1, "evidence")


# ------------------------------------------------------------------ failures are reported


def test_example_4_te1_fail_is_rendered_not_refused(tmp_path: Path) -> None:
    """A ``fail`` is a finding: the final report is built (Khoa reviews failures too)."""
    doc = final_metrics(statuses=["pass", "fail", "pass", "pass", "pass"])
    res, out = _final(tmp_path, doc)
    _ok(res)
    assert "| TE1 | fail" in out.read_text(encoding="utf-8")


def test_example_6_bullet_1_fail_is_rendered(tmp_path: Path) -> None:
    doc = final_metrics(a3_exit=None, verdict="fail", reason="route_unavailable")
    rb = runbook_text(a3_branch="fallback", steps=[s for s in EVIDENCE_JOBS if s != "S11"])
    res, out = _final(tmp_path, doc, runbook=rb)
    _ok(res)
    assert _bullet1(out)[1] == "fail (A3 not run: route_unavailable)"


# ------------------------------------------------------------------ examples 2 and 3: diagnoses


def test_example_2_te3_fail_with_diagnosis_exit_0(tmp_path: Path) -> None:
    doc = final_metrics(
        statuses=["pass", "pass", "pass", "pass", "reported (1 not passed, explained)"],
        te={"te3": "fail"},
    )
    res, out = _final(tmp_path, doc, notes=NOTES_OK + "\n" + diagnosis("TE3", "biology"))
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "## Diagnosis: TE3" in text and "biology" in text
    assert "reported (1 not passed, explained)" in text


def test_example_3_te3_fail_without_diagnosis_exit_1_naming_te3(tmp_path: Path) -> None:
    doc = final_metrics(
        statuses=["pass", "pass", "pass", "pass", "reported (1 not passed, explained)"],
        te={"te3": "fail"},
    )
    res, out = _final(tmp_path, doc, notes=NOTES_OK)
    _refused(res, out, 1, "TE3")


def test_diagnosis_without_a_cause_class_exit_1(tmp_path: Path) -> None:
    doc = final_metrics(te={"te4": "fail"})
    notes = NOTES_OK + "\n## Diagnosis: TE4\n\nWe looked at it and it is fine.\n"
    res, out = _final(tmp_path, doc, notes=notes)
    _refused(res, out, 1, "TE4")


def test_te7_fail_needs_its_diagnosis(tmp_path: Path) -> None:
    """Edge case "A failed TE": TE7 is not a bullet, but a fail is shown with its diagnosis."""
    res, out = _final(tmp_path / "a", final_metrics(te7="fail"), notes=NOTES_OK)
    _refused(res, out, 1, "TE7")
    res, out = _final(
        tmp_path / "b",
        final_metrics(te7="fail"),
        notes=NOTES_OK + "\n" + diagnosis("TE7", "pipeline_bug"),
    )
    _ok(res)
    assert "## Diagnosis: TE7" in out.read_text(encoding="utf-8")


def test_absent_notes_allowed_only_when_every_te_passes(tmp_path: Path) -> None:
    """Edge case "Empty input": an absent or empty notes file only when every TE passes."""
    res, _ = _final(tmp_path / "a", final_metrics(), notes=None)
    _ok(res)
    res, out = _final(tmp_path / "b", final_metrics(te={"te2": "fail"}), notes=None)
    _refused(res, out, 1, "TE2")
    res, out = _final(tmp_path / "c", final_metrics(te={"te2": "fail"}), notes="")
    _refused(res, out, 1, "TE2")


# ------------------------------------------------------------------ examples 7 and 8: missing keys


def test_example_7_absent_key_is_a_defect_exit_1(tmp_path: Path) -> None:
    template = FIXTURE_TEMPLATE + "\nOverlap: {{m6.A2_A3.jaccard}}\n"
    res, out = _final(tmp_path, final_metrics(), template=template)
    _refused(res, out, 1, "m6.A2_A3.jaccard")


def test_missing_marker_in_the_notes_exit_1(tmp_path: Path) -> None:
    """Check 3 of the packet, enforced at render: no ``**missing:`` string in the report."""
    res, out = _final(tmp_path, final_metrics(), notes=NOTES_OK + "\n**missing: m6.x**\n")
    _refused(res, out, 1)


def test_reason_missing_file_is_a_defect_exit_1(tmp_path: Path) -> None:
    """A metric whose table was absent (M10-05's ``reason = "missing:<file>"``) renders as
    ``**missing:`` in draft mode (M10-05 B2); the final report refuses it, naming the file."""
    doc = final_metrics()
    doc["m4"] = {"reason": "missing:compare_abundance.tsv"}
    res, out = _final(tmp_path / "final", doc)
    _refused(res, out, 1, "compare_abundance.tsv")
    res, out = _final(tmp_path / "draft", doc, final=False, runbook=None)
    _ok(res)
    assert "**missing:" in out.read_text(encoding="utf-8")


def test_example_8_null_value_with_reason_is_a_finding_exit_0(tmp_path: Path) -> None:
    doc = final_metrics()
    doc["m8"] = {
        "peptide": {
            "fdp_combined": None,
            "fdp_lower": None,
            "assessment": "not_evaluable",
            "reason": "empty_discoveries",
        }
    }
    template = FIXTURE_TEMPLATE + (
        "\n## 8. Metric 8\n\nFDP: {{m8.peptide.fdp_combined}}; {{m8.peptide.assessment}}\n"
        "Universe: {{m8.peptide.ci.low}}\n"
    )
    res, out = _final(tmp_path, doc, template=template)
    _ok(res)
    text = out.read_text(encoding="utf-8")
    fdp = next(ln for ln in text.splitlines() if ln.startswith("FDP:"))
    assert "empty_discoveries" in fdp, fdp  # "not evaluable (empty_discoveries)", not "NA"
    assert "**missing:" not in text


def test_example_5_m6_null_with_reason_exit_0(tmp_path: Path) -> None:
    doc = final_metrics(verdict="fail", reason="low_genome_match")
    doc["m6"] = {"value": None, "reason": "different_strain"}
    doc["h4"] = {"verdict": "unestimated", "reason": "different_strain"}
    template = (
        FIXTURE_TEMPLATE
        + "\nM6: {{m6.overlap_a2_a3.novel_orfs.called.jaccard}}; H4: {{h4.verdict}}\n"
    )
    res, out = _final(tmp_path, doc, template=template, runbook=runbook_text(a3_branch="fallback"))
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "different_strain" in text and "unestimated" in text
    assert "**missing:" not in text
    assert _bullet1(out)[1] == "pass (A3 = Casanovo fallback on UTI89, different strain)"


def test_s1_a3_not_run_is_a_finding_listed_under_not_evaluated(tmp_path: Path) -> None:
    """S1 (probe P1b): M10-10's ``reason: "a3_not_run"`` (no ``--a3``) is a finding: the
    placeholders render ``not evaluable (a3_not_run)``, "Not evaluated" lists ``m6`` and
    ``h4``, and the final report is built (example 6)."""
    doc = final_metrics(a3_exit=None, verdict="fail", reason="route_unavailable")
    doc["m6"] = {"reason": "a3_not_run"}
    doc["h4"] = {"reason": "a3_not_run"}
    template = FIXTURE_TEMPLATE + (
        "\nM6: {{m6.overlap_a2_a3.novel_orfs.called.jaccard}}; H4: {{h4.verdict}}\n"
    )
    rb = runbook_text(a3_branch="fallback", steps=[s for s in EVIDENCE_JOBS if s != "S11"])
    res, out = _final(tmp_path, doc, template=template, runbook=rb)
    _ok(res)
    text = out.read_text(encoding="utf-8")
    m6 = next(ln for ln in text.splitlines() if ln.startswith("M6:"))
    assert m6.count("not evaluable (a3_not_run)") == 2, m6
    section = text.split("## 7. Not evaluated", 1)[1].split("\n## ", 1)[0]
    assert "m6" in section and "h4" in section and "a3_not_run" in section, section


def test_metrics_without_checklist_is_a_precondition_exit_2(tmp_path: Path) -> None:
    doc = final_metrics()
    del doc["checklist"]
    template = FIXTURE_TEMPLATE.split("## 1. Status")[0] + "Commit {{gippa_commit}}.\n"
    res, out = _final(tmp_path, doc, template=template)
    _refused(res, out, 2)


def test_out_inside_the_checkout_exit_2(tmp_path: Path) -> None:
    """M10-05's rule holds in final mode: the report is rendered outside the checkout."""
    from helpers.pilot_metrics_m10 import REPO

    m, n, t, r = _setup(tmp_path, final_metrics(), runbook=runbook_text())
    out = REPO / "validation" / "pilot" / "metrics" / "REPORT_m10_15_test.md"
    res = render_final(m, out, notes=n, template=t, runbook=r)
    try:
        assert res.returncode == 2, res.stderr
        assert "unrecognized arguments" not in res.stderr, res.stderr
        assert not out.exists()
    finally:
        out.unlink(missing_ok=True)


# ------------------------------------------------------------------ N4: draft pending markers


def test_n4_draft_mode_renders_the_pending_placeholders_as_pending(tmp_path: Path) -> None:
    """N4: before M10-07/M10-10 exist, ``m6``/``h4``/``m8`` render their pending reason in the
    draft, not a bare ``NA``."""
    level = {"fdp_combined": None, "fdp_lower": None, "assessment": None}
    doc = final_metrics()
    doc["m6"] = {
        "overlap_a2_a3": {"novel_orfs": {"called": {"jaccard": None}}},
        "reason": "pending (M10-10)",
    }
    doc["h4"] = {"verdict": None, "reason": "pending (M10-10)"}
    doc["m8"] = {"peptide": dict(level), "reason": "pending (M10-07)"}
    doc["evidence"]["runbook"] = {"reason": "pending (render --final --runbook)"}
    template = FIXTURE_TEMPLATE + (
        "\nM6: {{m6.overlap_a2_a3.novel_orfs.called.jaccard}}\nH4: {{h4.verdict}}\n"
        "M8: {{m8.peptide.fdp_combined}}\n"
    )
    res, out = _final(tmp_path, doc, template=template, final=False, runbook=None)
    _ok(res)
    lines = {
        ln[:3]: ln
        for ln in out.read_text(encoding="utf-8").splitlines()
        if ln[:3] in ("M6:", "H4:", "M8:")
    }
    assert "pending (M10-10)" in lines["M6:"] and "NA" not in lines["M6:"], lines["M6:"]
    assert "pending (M10-10)" in lines["H4:"] and "NA" not in lines["H4:"], lines["H4:"]
    assert "pending (M10-07)" in lines["M8:"] and "NA" not in lines["M8:"], lines["M8:"]


# ------------------------------------------------------------------ property: determinism


STATUS = st.sampled_from(
    [
        "pass",
        "fail",
        "fail (no data)",
        "fail (incomplete: sim5m_core)",
        "reported (2 not passed, explained)",
        "pending (ADR-M10-D proposed)",
        "pending (M10-06)",
        "partial (A1, A2; A3 pending M10-10)",
    ]
)
TE_STATUS = st.sampled_from(["pass", "fail"])
#: (runbook a3_branch, A3 exit code, refused): bullet 1 is ``partial`` only on a pending branch
BULLET1 = st.sampled_from(
    [
        ("cascadia", 0, False),
        ("fallback", None, False),
        ("cascadia", 2, False),
        ("pending", None, True),
    ]
)


@settings(
    max_examples=10, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    statuses=st.lists(STATUS, min_size=4, max_size=4),
    bullet1=BULLET1,
    te=st.fixed_dictionaries({k: TE_STATUS for k in ("te2", "te3", "te4", "te5")}),
    rho=st.floats(-1, 1, allow_nan=False),
)
def test_property_final_render_deterministic_and_refuses_exactly_pending(
    tmp_path_factory: pytest.TempPathFactory,
    statuses: list[str],
    bullet1: tuple[str, int | None, bool],
    te: dict[str, str],
    rho: float,
) -> None:
    """Render is a pure function of its inputs (byte-identical twice), and in final mode it
    refuses exactly when a checklist status (bullet 1 relabelled from the runbook) starts with
    ``pending`` or ``partial`` (every failed TE has its diagnosis here)."""
    branch, a3_exit, pending_branch = bullet1
    doc = final_metrics(
        statuses=[DRAFT_BULLET1.get(a3_exit, "fail"), *statuses],
        te=te,
        rho=rho,
        a3_exit=a3_exit,
        verdict="pass" if a3_exit is not None else "fail",
        reason=None if a3_exit is not None else "route_unavailable",
    )
    steps = [s for s in EVIDENCE_JOBS if s != "S11" or a3_exit is not None]
    rb = runbook_text(a3_branch=branch, steps=steps)
    notes = NOTES_OK + "".join("\n" + diagnosis(f"TE{k[2:]}") for k, v in te.items() if v != "pass")
    tmp = tmp_path_factory.mktemp("prop")
    r1, o1 = _final(tmp / "a", doc, notes=notes, runbook=rb)
    r2, o2 = _final(tmp / "b", doc, notes=notes, runbook=rb)
    refuse = pending_branch or any(s.startswith(("pending", "partial")) for s in statuses)
    assert r1.returncode == r2.returncode == (1 if refuse else 0), r1.stderr
    if refuse:
        assert not o1.exists() and not o2.exists()
    else:
        text = o1.read_text(encoding="utf-8")
        assert o1.read_bytes() == o2.read_bytes()
        assert "**missing:" not in text and "{{" not in text
        assert json.loads((tmp / "a" / "in" / "metrics.json").read_text())["checklist"]
