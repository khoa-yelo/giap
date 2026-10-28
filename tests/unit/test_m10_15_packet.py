"""``validation/pilot/review_packet.py`` (M10-15 acceptance criterion 5; checks 1-6 of "The review
packet"; examples 1, 4, 6, 10 (round r2 with a diff); edge case "``metrics.json`` from another
commit"; review rev-1015 B1, B2, S3-S6, S9(c), S9(d), N6, N7 and the mutation survivors: check 6's
own message, the freeze-directory comparison, the ``fail (A3 exit n)`` branch).

Spec: docs/specs/M10-15-pilot-gate.md, "The review packet (``local/review/m10-pilot/r<k>/``,
D3)". Everything runs in a ``tmp_path`` sandbox built by ``helpers.pilot_gate_m10_15``: a
checkout-like ``repo/`` (``--repo``), a ``root/`` standing in for ``$GIPPA_ROOT`` (also exported to
the subprocess) and ``out/m10-pilot`` (``--out``). The command-line options beyond the spec's four
are pinned in that module's docstring. The committed ``REPORT.md`` of the sandbox is rendered by
``metrics.py render --final``.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_gate_m10_15 import (
    EVIDENCE_JOBS,
    F7,
    OTHER_COMMIT,
    REPO,
    Sandbox,
    build_sandbox,
    packet_files,
    run_packet,
    runbook_text,
    status_row,
    tree,
)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _files(round_dir: Path) -> set[str]:
    return {str(p.relative_to(round_dir)) for p in round_dir.rglob("*") if p.is_file()}


def _refused(sb: Sandbox, res: Any, code: int, round_k: int = 1, *names: str) -> None:
    assert res.returncode == code, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert "unrecognized arguments" not in res.stderr, res.stderr  # a usage error is not a check
    assert "Traceback" not in res.stderr, res.stderr  # S4: only the documented exit codes
    assert not (sb.out / f"r{round_k}").exists(), "nothing is written on a failed check"
    assert not (sb.out / f".r{round_k}.tmp").exists(), "S4: the temporary round is removed"
    for name in names:
        assert name in res.stderr + res.stdout, f"{name!r} not named in the message"


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
    index = (sb.out / "r1" / "index.md").read_text(encoding="utf-8")
    assert "fail (A3 not run: route_unavailable)" in index


def test_example_4_failing_bullet_still_builds_the_packet(tmp_path: Path) -> None:
    sb = build_sandbox(tmp_path, statuses=["pass", "fail", "pass", "pass", "pass"])
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


def test_check_2_gippa_commit_differs_from_freeze_commit_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(freeze_commit=OTHER_COMMIT), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


def test_check_2_evidence_row_from_another_commit_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(commits={"S9": OTHER_COMMIT}), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


def test_check_2_runbook_not_final_exit_1(sb: Sandbox) -> None:
    sb.runbook.write_text(runbook_text(status="draft"), encoding="utf-8")
    _refused(sb, sb.run(1), 1)


@pytest.mark.parametrize("first", [False, True])
def test_b2_superseded_evidence_row_from_another_commit_is_ignored(
    tmp_path: Path, first: bool
) -> None:
    """B2 (probe P8): M10-14 keeps superseded attempts in the Evidence table; check 2 skips
    them, and the report names the accepted job."""
    rb = runbook_text(superseded=[("S8", "9000001", OTHER_COMMIT)], superseded_first=first)
    sb = build_sandbox(tmp_path, runbook=rb)
    _ok(sb.run(1))
    index = (sb.out / "r1" / "index.md").read_text(encoding="utf-8")
    assert EVIDENCE_JOBS["S8"][1] in index and "9000001" not in index


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
    """S6: check 4 lives in ``final_problems``; the committed and freeze ``metrics.json`` agree
    here, so only the re-render can refuse, and its message says why (S5)."""
    doc = json.loads(sb.metrics.read_text(encoding="utf-8"))
    doc["checklist"][3]["evidence"] = ""
    sb.write_metrics(doc)
    _refused(sb, sb.run(1), 1, 1, "evidence")


def test_s5_check_1_message_carries_the_refusal_reason(sb: Sandbox) -> None:
    doc = json.loads(sb.metrics.read_text(encoding="utf-8"))
    doc["checklist"][3]["status"] = "pending (M10-12)"
    sb.write_metrics(doc)
    _refused(sb, sb.run(1), 1, 1, "1:", "pending (M10-12)")


def test_s6_freeze_metrics_differs_from_the_committed_one_exit_1(sb: Sandbox) -> None:
    """Mutation survivor: the committed ``metrics.json`` and ``REPORT.md`` agree, but the
    freeze directory's copy holds another value."""
    doc = json.loads((sb.metrics_dir / "metrics.json").read_text(encoding="utf-8"))
    doc["m4"]["spearman_rho"] = 0.5
    (sb.metrics_dir / "metrics.json").write_text(
        json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    _refused(sb, sb.run(1), 1, 1, "1:")


def test_n7_freeze_metrics_compared_as_bytes(sb: Sandbox) -> None:
    """N7: the same text with CRLF line ends is not the same file."""
    p = sb.metrics_dir / "metrics.json"
    p.write_bytes(p.read_bytes().replace(b"\n", b"\r\n"))
    _refused(sb, sb.run(1), 1, 1, "1:")


def test_n7_report_compared_as_bytes(sb: Sandbox) -> None:
    """N7: check 1 compares bytes; ``read_text`` would hide CRLF line ends."""
    sb.report.write_bytes(sb.report.read_bytes().replace(b"\n", b"\r\n"))
    _refused(sb, sb.run(1), 1, 1, "1:")


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
    """S4 / mutation survivor: the round-exists check fires with its own message (not an
    ``OSError`` from the rename), and no ``.r1.tmp`` is left."""
    _ok(sb.run(1))
    before = tree(sb.out / "r1")
    res = sb.run(1)
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stderr}"
    assert "6:" in res.stderr and "Traceback" not in res.stderr, res.stderr
    assert tree(sb.out / "r1") == before
    assert not (sb.out / ".r1.tmp").exists()


def test_check_6_an_empty_existing_round_is_refused_too(sb: Sandbox) -> None:
    """``rename`` onto an empty directory succeeds on POSIX: only the explicit check stops it."""
    (sb.out / "r1").mkdir(parents=True)
    res = sb.run(1)
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stderr}"
    assert "6:" in res.stderr, res.stderr
    assert list((sb.out / "r1").iterdir()) == []


def test_s4_copy_error_exits_2_and_cleans_up(sb: Sandbox) -> None:
    """S4: a source that cannot be read: exit 2, no traceback, no round, no ``.r1.tmp``."""
    if os.geteuid() == 0:
        pytest.skip("root reads unreadable files")
    mode = sb.bench_summary.stat().st_mode
    sb.bench_summary.chmod(0)
    try:
        _refused(sb, sb.run(1), 2)
    finally:
        sb.bench_summary.chmod(stat.S_IMODE(mode))


# ------------------------------------------------------------------ S3: where --out may point


def test_s3_out_inside_the_checkout_exit_2(tmp_path: Path) -> None:
    """Probe P9: the per-gene HTML must never be one ``git add`` away from a commit."""
    sb = build_sandbox(tmp_path, git=True)
    sb.out = sb.repo / "validation" / "pilot" / "packet"
    before = tree(sb.repo)
    _refused(sb, sb.run(1), 2)
    assert tree(sb.repo) == before


def test_s3_out_inside_a_non_git_repo_dir_exit_2(sb: Sandbox) -> None:
    sb.out = sb.repo / "review"
    _refused(sb, sb.run(1), 2)


def test_s3_out_under_a_git_ignored_path_is_allowed(tmp_path: Path) -> None:
    sb = build_sandbox(tmp_path, git=True)
    sb.out = sb.repo / "local" / "review" / "m10-pilot"
    _ok(sb.run(1))
    assert _files(sb.out / "r1") == packet_files(sb)


# ------------------------------------------------------------------ B1 / S9(c): bullet 1 at render


def test_s9c_a3_exit_nonzero_requires_the_a3_html(tmp_path: Path) -> None:
    """``fail (A3 exit n)`` reports an A3 run, so its ``report.html`` is a packet input."""
    sb = build_sandbox(tmp_path / "a", a3_exit=1)
    assert status_row(sb.report.read_text(encoding="utf-8"), "runs")[1] == "fail (A3 exit 1)"
    _ok(sb.run(1))
    assert "runs/A3/report.html" in _files(sb.out / "r1")
    sb = build_sandbox(tmp_path / "b", a3_exit=1)
    assert sb.a3 is not None
    (sb.a3 / "report" / "report.html").unlink()
    _refused(sb, sb.run(1), 2)


def test_b1_check_1_re_renders_with_the_given_runbook(sb: Sandbox) -> None:
    """The committed report was rendered from the Cascadia runbook; a runbook whose branch is
    now ``fallback`` re-labels bullet 1, so the re-render differs."""
    sb.runbook.write_text(runbook_text(a3_branch="fallback"), encoding="utf-8")
    _refused(sb, sb.run(1), 1, 1, "1:")


# ------------------------------------------------------------------ N6: no bytecode in the repo


def test_n6_the_packet_writes_no_pycache(sb: Sandbox, tmp_path: Path) -> None:
    """N6: run a copy of ``validation/`` (without ``__pycache__``) with bytecode writing
    enabled in the environment; the script must not write any."""
    copy = tmp_path / "code"
    shutil.copytree(
        REPO / "validation",
        copy / "validation",
        ignore=shutil.ignore_patterns("__pycache__", "db_a2", "inputs"),
    )
    pythonpath = os.pathsep.join([str(copy), str(REPO / "src"), str(REPO / "tests" / "sim")])
    res = run_packet(
        *sb.args(1),
        root=sb.root,
        cwd=sb.base,
        script=copy / "validation" / "pilot" / "review_packet.py",
        pythonpath=pythonpath,
    )
    _ok(res)
    caches = sorted(str(p.relative_to(copy)) for p in copy.rglob("__pycache__"))
    assert caches == [], caches


def test_index_status_table_comes_from_the_re_rendered_report(sb: Sandbox) -> None:
    """``index.md`` copies section 1 of ``REPORT.md``: bullet 1 there is the render-time label,
    not the M10-05 string in ``metrics.json``."""
    _ok(sb.run(1))
    index = (sb.out / "r1" / "index.md").read_text(encoding="utf-8")
    assert "pass (A3 = Cascadia)" in index
