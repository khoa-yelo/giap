"""``metrics.py render --final`` (M10-15 acceptance criterion 3; examples 2, 3, 7, 8; edge cases
"A failed TE", "A3 skipped" (runbook pending -> ``partial``), "A missing metric", "Empty input";
the byte-identical re-render behind criterion 2).

Spec: docs/specs/M10-15-pilot-gate.md, "The final report", "Checklist rows", "Examples", "Edge
cases". Inputs are hand-made (``tests/helpers/pilot_gate_m10_15.py``: a small final-mode
template passed with ``--template`` and a minimal ``metrics.json``); the interface pinned beyond
the spec is listed in that module's docstring. Nothing comes from ``$GIAP_ROOT``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import (
    FIXTURE_TEMPLATE,
    NOTES_OK,
    diagnosis,
    final_metrics,
    render_final,
    write_json,
    write_text,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def _setup(
    tmp: Path, doc: dict[str, Any], notes: str | None = NOTES_OK, template: str = FIXTURE_TEMPLATE
) -> tuple[Path, Path | None, Path]:
    m = write_json(tmp / "in" / "metrics.json", doc)
    n = None if notes is None else write_text(tmp / "in" / "report_notes.md", notes)
    t = write_text(tmp / "in" / "report_template.md", template)
    return m, n, t


def _final(
    tmp: Path,
    doc: dict[str, Any],
    notes: str | None = NOTES_OK,
    template: str = FIXTURE_TEMPLATE,
    final: bool = True,
) -> tuple[Any, Path]:
    m, n, t = _setup(tmp, doc, notes, template)
    out = tmp / "out" / "REPORT.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    return render_final(m, out, notes=n, template=t, final=final), out


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _refused(res: Any, out: Path, code: int = 1, *names: str) -> None:
    assert res.returncode == code, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "unrecognized arguments" not in res.stderr, res.stderr  # a usage error is not a check
    assert not out.exists(), "nothing is written on a refused render"
    for name in names:
        assert name in res.stderr + res.stdout, f"{name!r} not named in the message"


# ------------------------------------------------------------------ example 1: all pass


def test_example_1_all_pass_renders_exit_0(tmp_path: Path) -> None:
    res, out = _final(tmp_path, final_metrics())
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "pass (A3 = Cascadia)" in text
    assert "Submitted for Khoa's review" in text
    assert "**missing:" not in text and "{{" not in text


def test_final_render_is_byte_identical_when_re_rendered(tmp_path: Path) -> None:
    res1, out1 = _final(tmp_path / "a", final_metrics(), notes=NOTES_OK + diagnosis("TE3"))
    res2, out2 = _final(tmp_path / "b", final_metrics(), notes=NOTES_OK + diagnosis("TE3"))
    _ok(res1)
    _ok(res2)
    assert out1.read_bytes() == out2.read_bytes()


# ------------------------------------------------------------------ pending statuses


@pytest.mark.parametrize(
    ("bullet", "status"),
    [
        (3, "pending (ADR-M10-D proposed)"),  # example 9, second half
        (3, "pending (M10-12)"),
        (4, "pending (M10-06)"),
        (0, "partial (A1, A2; A3 pending M10-10)"),  # edge case: runbook a3_branch=pending
    ],
)
def test_pending_or_partial_checklist_status_exit_1(
    tmp_path: Path, bullet: int, status: str
) -> None:
    statuses = ["pass (A3 = Cascadia)", "pass", "pass", "pass", "pass"]
    statuses[bullet] = status
    res, out = _final(tmp_path, final_metrics(statuses=statuses))
    _refused(res, out, 1)


@pytest.mark.parametrize("te", ["te2", "te5"])
def test_pending_te_status_exit_1(tmp_path: Path, te: str) -> None:
    res, out = _final(tmp_path, final_metrics(te={te: "pending (M10-06)"}))
    _refused(res, out, 1)


def test_pending_te7_status_exit_1(tmp_path: Path) -> None:
    res, out = _final(tmp_path, final_metrics(te7="pending (M10-06)"))
    _refused(res, out, 1)


def test_draft_mode_keeps_m10_05_behaviour_on_pending(tmp_path: Path) -> None:
    """Criterion 3: without ``--final`` a pending status is rendered (M10-05's draft)."""
    doc = final_metrics(
        statuses=[
            "partial (A1, A2; A3 pending M10-10)",
            "pass",
            "pass",
            "pending (M10-12)",
            "pending (M10-06)",
        ]
    )
    res, out = _final(tmp_path, doc, final=False)
    _ok(res)
    assert "pending (M10-12)" in out.read_text(encoding="utf-8")


# ------------------------------------------------------------------ failures are reported


def test_example_4_te1_fail_is_rendered_not_refused(tmp_path: Path) -> None:
    """A ``fail`` is a finding: the final report is built (Khoa reviews failures too)."""
    doc = final_metrics(statuses=["pass (A3 = Cascadia)", "fail", "pass", "pass", "pass"])
    res, out = _final(tmp_path, doc)
    _ok(res)
    assert "| TE1 | fail" in out.read_text(encoding="utf-8")


def test_example_6_bullet_1_fail_is_rendered(tmp_path: Path) -> None:
    doc = final_metrics(
        statuses=[
            "fail (A3 not run: route_unavailable; ADR-M10-C declined)",
            "pass",
            "pass",
            "pass",
            "pass",
        ]
    )
    res, out = _final(tmp_path, doc)
    _ok(res)
    assert "fail (A3 not run: route_unavailable" in out.read_text(encoding="utf-8")


# ------------------------------------------------------------------ examples 2 and 3: diagnoses


def test_example_2_te3_fail_with_diagnosis_exit_0(tmp_path: Path) -> None:
    doc = final_metrics(
        statuses=[
            "pass (A3 = Cascadia)",
            "pass",
            "pass",
            "pass",
            "reported (1 not passed, explained)",
        ],
        te={"te3": "fail"},
    )
    res, out = _final(tmp_path, doc, notes=NOTES_OK + "\n" + diagnosis("TE3", "biology"))
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "## Diagnosis: TE3" in text and "biology" in text
    assert "reported (1 not passed, explained)" in text


def test_example_3_te3_fail_without_diagnosis_exit_1_naming_te3(tmp_path: Path) -> None:
    doc = final_metrics(
        statuses=[
            "pass (A3 = Cascadia)",
            "pass",
            "pass",
            "pass",
            "reported (1 not passed, explained)",
        ],
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
    ``**missing:`` in draft mode (M10-05 B2); the final report refuses it."""
    doc = final_metrics()
    doc["m4"] = {"reason": "missing:compare_abundance.tsv"}
    res, out = _final(tmp_path / "final", doc)
    _refused(res, out, 1)
    res, out = _final(tmp_path / "draft", doc, final=False)
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
    doc = final_metrics(
        statuses=[
            "pass (A3 = Casanovo fallback on UTI89, different strain)",
            "pass",
            "pass",
            "pass",
            "pass",
        ]
    )
    doc["m6"] = {"value": None, "reason": "different_genome"}
    doc["h4"] = {"verdict": "unestimated", "reason": "different_genome"}
    template = FIXTURE_TEMPLATE + "\nM6: {{m6.A2_A3.jaccard}}; H4: {{h4.verdict}}\n"
    res, out = _final(tmp_path, doc, template=template)
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "different_genome" in text and "unestimated" in text
    assert "**missing:" not in text


def test_metrics_without_checklist_is_a_precondition_exit_2(tmp_path: Path) -> None:
    doc = final_metrics()
    del doc["checklist"]
    template = FIXTURE_TEMPLATE.split("## 1. Status")[0] + "Commit {{giap_commit}}.\n"
    res, out = _final(tmp_path, doc, template=template)
    _refused(res, out, 2)


def test_out_inside_the_checkout_exit_2(tmp_path: Path) -> None:
    """M10-05's rule holds in final mode: the report is rendered outside the checkout."""
    from helpers.pilot_metrics_m10 import REPO

    m, n, t = _setup(tmp_path, final_metrics())
    out = REPO / "validation" / "pilot" / "metrics" / "REPORT_m10_15_test.md"
    res = render_final(m, out, notes=n, template=t)
    try:
        assert res.returncode == 2, res.stderr
        assert "unrecognized arguments" not in res.stderr, res.stderr
        assert not out.exists()
    finally:
        out.unlink(missing_ok=True)


# ------------------------------------------------------------------ property: determinism


STATUS = st.sampled_from(
    [
        "pass",
        "fail",
        "fail (no data)",
        "pass (A3 = Cascadia)",
        "pass (A3 = Casanovo fallback on UTI89, different strain)",
        "fail (A3 not run: route_unavailable)",
        "fail (incomplete: sim5m_core)",
        "reported (2 not passed, explained)",
        "pending (ADR-M10-D proposed)",
        "pending (M10-06)",
        "partial (A1, A2; A3 pending M10-10)",
    ]
)
TE_STATUS = st.sampled_from(["pass", "fail"])


@settings(
    max_examples=10, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    statuses=st.lists(STATUS, min_size=5, max_size=5),
    te=st.fixed_dictionaries({k: TE_STATUS for k in ("te2", "te3", "te4", "te5")}),
    rho=st.floats(-1, 1, allow_nan=False),
)
def test_property_final_render_deterministic_and_refuses_exactly_pending(
    tmp_path_factory: pytest.TempPathFactory,
    statuses: list[str],
    te: dict[str, str],
    rho: float,
) -> None:
    """Render is a pure function of its inputs (byte-identical twice), and in final mode it
    refuses exactly when a checklist status starts with ``pending`` or ``partial`` (every failed
    TE has its diagnosis here)."""
    doc = final_metrics(statuses=statuses, te=te, rho=rho)
    notes = NOTES_OK + "".join("\n" + diagnosis(f"TE{k[2:]}") for k, v in te.items() if v != "pass")
    tmp = tmp_path_factory.mktemp("prop")
    r1, o1 = _final(tmp / "a", doc, notes=notes)
    r2, o2 = _final(tmp / "b", doc, notes=notes)
    refuse = any(s.startswith(("pending", "partial")) for s in statuses)
    assert r1.returncode == r2.returncode == (1 if refuse else 0), r1.stderr
    if refuse:
        assert not o1.exists() and not o2.exists()
    else:
        text = o1.read_text(encoding="utf-8")
        assert o1.read_bytes() == o2.read_bytes()
        assert "**missing:" not in text and "{{" not in text
        assert json.loads((tmp / "a" / "in" / "metrics.json").read_text())["checklist"]
