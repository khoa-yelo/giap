"""``validation/pilot/review_packet.py`` (M10-15 acceptance criterion 5; checks 1-6 of "The review
packet"; examples 1, 4, 6, 10 (round r2 with a diff); edge case "``metrics.json`` from another
commit").

Spec: docs/specs/M10-15-pilot-gate.md, "The review packet (``local/review/m10-pilot/r<k>/``,
D3)". Everything runs in a ``tmp_path`` sandbox built by ``helpers.pilot_gate_m10_15``: a
checkout-like ``repo/`` (``--repo``), a ``root/`` standing in for ``$GIAP_ROOT`` (also exported to
the subprocess) and ``out/m10-pilot`` (``--out``). The command-line options beyond the spec's four
are pinned in that module's docstring. The committed ``REPORT.md`` of the sandbox is rendered by
``metrics.py render --final``.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import (
    F7,
    OTHER_COMMIT,
    Sandbox,
    build_sandbox,
    packet_files,
    runbook_text,
    tree,
)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _files(round_dir: Path) -> set[str]:
    return {str(p.relative_to(round_dir)) for p in round_dir.rglob("*") if p.is_file()}


def _refused(sb: Sandbox, res: Any, code: int, round_k: int = 1) -> None:
    assert res.returncode == code, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "unrecognized arguments" not in res.stderr, res.stderr  # a usage error is not a check
    assert not (sb.out / f"r{round_k}").exists(), "nothing is written on a failed check"


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    return build_sandbox(tmp_path)


# ------------------------------------------------------------------ example 1: round r1


def test_example_1_writes_exactly_the_packet_files(sb: Sandbox) -> None:
    _ok(sb.run(1))
    r1 = sb.out / "r1"
    assert _files(r1) == packet_files(sb)


def test_manifest_lists_every_file_sorted_with_matching_hashes(sb: Sandbox) -> None:
    _ok(sb.run(1))
    r1 = sb.out / "r1"
    lines = (r1 / "MANIFEST.sha256").read_text(encoding="utf-8").splitlines()
    entries = [ln.split("  ", 1) for ln in lines if ln.strip()]
    paths = [p for _, p in entries]
    assert paths == sorted(paths), "MANIFEST.sha256 is sorted by path"
    assert set(paths) == packet_files(sb) - {"MANIFEST.sha256"}
    for digest, path in entries:
        assert digest == hashlib.sha256((r1 / path).read_bytes()).hexdigest(), path


def test_copies_are_byte_identical_to_their_sources(sb: Sandbox) -> None:
    _ok(sb.run(1))
    r1 = sb.out / "r1"
    pairs = {
        "REPORT.md": sb.report,
        "metrics.json": sb.metrics,
        "report_notes.md": sb.notes,
        "reproduce.md": sb.runbook,
        f"verify-{F7}.json": sb.verify,
        "benchmark/summary.md": sb.bench_summary,
        "a3_feasibility.json": sb.repo / "validation/pilot/metrics/a3_feasibility.json",
        "runs/A1/report.html": sb.a1 / "report" / "report.html",
        "runs/A2/report.html": sb.a2 / "report" / "report.html",
        "runs/compare.html": sb.cmp / "report" / "compare.html",
    }
    assert sb.a3 is not None
    pairs["runs/A3/report.html"] = sb.a3 / "report" / "report.html"
    for rel, src in pairs.items():
        assert (r1 / rel).read_bytes() == src.read_bytes(), rel
    for adr in sb.adrs:
        assert (r1 / "adrs" / adr.name).read_bytes() == adr.read_bytes()


def test_index_holds_the_status_table_and_the_adr_status_lines(sb: Sandbox) -> None:
    _ok(sb.run(1))
    index = (sb.out / "r1" / "index.md").read_text(encoding="utf-8")
    report = sb.report.read_text(encoding="utf-8")
    section1 = report.split("## 1. Status", 1)[1].split("\n## ", 1)[0]
    table = [ln for ln in section1.splitlines() if ln.startswith("|")]
    assert table
    for ln in table:
        assert ln in index, f"status table line not copied: {ln}"
    for adr, status in zip(sb.adrs, ("proposed", "accepted"), strict=True):
        line = next(ln for ln in index.splitlines() if adr.stem in ln or adr.name in ln)
        assert status in line, line


def test_never_writes_outside_out(sb: Sandbox) -> None:
    before = tree(sb.base, exclude=sb.out)
    _ok(sb.run(1))
    assert tree(sb.base, exclude=sb.out) == before


def test_example_6_no_a3_run_packet_without_a3_html(tmp_path: Path) -> None:
    """Check 5: A3's ``report.html`` only when bullet 1 reports an A3 run; a failing bullet
    still gets a packet (Khoa reviews failures too)."""
    sb = build_sandbox(tmp_path, a3=False)
    _ok(sb.run(1))
    assert _files(sb.out / "r1") == packet_files(sb, a3=False)


def test_example_4_failing_bullet_still_builds_the_packet(tmp_path: Path) -> None:
    sb = build_sandbox(tmp_path, statuses=["pass (A3 = Cascadia)", "fail", "pass", "pass", "pass"])
    _ok(sb.run(1))
    assert (sb.out / "r1" / "MANIFEST.sha256").is_file()


# ------------------------------------------------------------------ check 1: re-render


def test_check_1_report_differs_from_the_re_render_exit_1(sb: Sandbox) -> None:
    sb.report.write_text(
        sb.report.read_text(encoding="utf-8") + "\nEdited by hand.\n", encoding="utf-8"
    )
    _refused(sb, sb.run(1), 1)


def test_check_1_metrics_changed_after_rendering_exit_1(sb: Sandbox) -> None:
    doc = json.loads(sb.metrics.read_text(encoding="utf-8"))
    doc["m4"]["spearman_rho"] = 0.5
    sb.metrics.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    _refused(sb, sb.run(1), 1)


# ------------------------------------------------------------------ check 2: commits


def test_check_2_giap_commit_differs_from_freeze_commit_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(freeze_commit=OTHER_COMMIT), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


def test_check_2_evidence_row_from_another_commit_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(commits={"S9": OTHER_COMMIT}), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


def test_check_2_runbook_not_final_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(status="draft"), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


# ------------------------------------------------------------------ checks 3 and 4


def test_check_3_missing_marker_in_report_exit_1(sb: Sandbox) -> None:
    text = sb.report.read_text(encoding="utf-8")
    sb.report.write_text(
        text.replace("Submitted", "**missing: m6.A2_A3.jaccard** Submitted"), encoding="utf-8"
    )
    _refused(sb, sb.run(1), 1)


def test_check_3_pending_status_in_report_exit_1(sb: Sandbox) -> None:
    text = sb.report.read_text(encoding="utf-8")
    sb.report.write_text(
        text.replace("| runtimes | pass", "| runtimes | pending (M10-12)"), encoding="utf-8"
    )
    _refused(sb, sb.run(1), 1)


def test_check_4_empty_evidence_exit_1(sb: Sandbox) -> None:
    doc = json.loads(sb.metrics.read_text(encoding="utf-8"))
    doc["checklist"][3]["evidence"] = ""
    sb.metrics.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    _refused(sb, sb.run(1), 1)


# ------------------------------------------------------------------ check 5: inputs exist


@pytest.mark.parametrize("which", ["compare", "a3_html", "verify", "bench"])
def test_check_5_missing_input_exit_2(sb: Sandbox, which: str) -> None:
    assert sb.a3 is not None
    path = {
        "compare": sb.cmp / "report" / "compare.html",
        "a3_html": sb.a3 / "report" / "report.html",
        "verify": sb.verify,
        "bench": sb.bench_summary,
    }[which]
    path.unlink()
    _refused(sb, sb.run(1), 2)


def test_missing_metrics_dir_exit_2(sb: Sandbox) -> None:
    import shutil

    shutil.rmtree(sb.metrics_dir)
    _refused(sb, sb.run(1), 2)


# ------------------------------------------------------------------ check 6: rounds


def test_check_6_existing_round_is_never_overwritten(sb: Sandbox) -> None:
    _ok(sb.run(1))
    before = tree(sb.out / "r1")
    res = sb.run(1)
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stderr}"
    assert tree(sb.out / "r1") == before


def test_example_10_round_2_with_the_diff_against_round_1(sb: Sandbox) -> None:
    _ok(sb.run(1))
    r1 = tree(sb.out / "r1")
    sb.notes.write_text("The H1 paragraph was reworded at Khoa's request.\n", encoding="utf-8")
    sb.rerender()
    _ok(sb.run(2))
    assert tree(sb.out / "r1") == r1, "round r1 is left as it was"
    index = (sb.out / "r2" / "index.md").read_text(encoding="utf-8")
    assert "The H1 paragraph was reworded at Khoa's request." in index
    assert "The pilot ran at the freeze commit." in index  # the removed line, in the diff
    assert _files(sb.out / "r2") == packet_files(sb)
