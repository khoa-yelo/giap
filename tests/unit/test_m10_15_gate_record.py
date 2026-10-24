"""Shape of the M10 gate record in ``PROGRESS.md`` (M10-15 acceptance criterion 7; examples 1,
10, 11; "The gate record (D7) and Khoa's review").

When a ``### Gate M10`` block exists, its status line is one of ``Gate M10: awaiting Khoa's
review`` / ``Gate M10: not met: <bullets>`` / ``Gate M10: met``; a ``met`` line needs an
Answered-by-Khoa entry ``- <date>: **M10 pilot report** (REPORT.md sha256 <64 hex>, packet
r<k>): "<quote>" → verdict accept`` whose sha256 equals the sha256 of the committed
``validation/pilot/REPORT.md``. This test checks the shape only: it cannot tell who wrote the
quote (example 11), which is why the recording rule binds the lead.

The checker is defined here and tested on hand-made PROGRESS texts; in 15a ``PROGRESS.md`` has no
gate block yet, so the check on the real file passes vacuously until 15b.
"""

from __future__ import annotations

import hashlib
import re

import pytest
from helpers.pilot_metrics_m10 import REPO

PROGRESS = REPO / "PROGRESS.md"
REPORT = REPO / "validation" / "pilot" / "REPORT.md"

GATE_HEAD = re.compile(r"^### Gate M10\b.*$", re.M)
STATUS = re.compile(r"Gate M10: (awaiting Khoa's review|not met: \S.*|met)\s*$")
STATUS_ANY = re.compile(r"Gate M10:")
ENTRY = re.compile(
    r"^- (\d{4}-\d{2}-\d{2})[^:]*: \*\*M10 pilot report\*\* \(REPORT\.md sha256 ([0-9a-f]{64}), "
    r"packet r(\d+)\): \"([^\"]+)\" → verdict (accept|changes|reject)\s*$",
    re.M,
)


def _section(text: str, heading: str) -> str:
    m = re.search(rf"^## {re.escape(heading)}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def gate_problems(progress: str, report_sha: str | None) -> list[str]:
    """Shape problems of the gate record (empty: fine or no gate block)."""
    m = GATE_HEAD.search(progress)
    if m is None:
        return []
    nxt = re.search(r"^#{2,3} ", progress[m.end() :], re.M)
    block = progress[m.end() : m.end() + nxt.start() if nxt else len(progress)]
    lines = [ln.strip().strip("*").strip() for ln in block.splitlines() if STATUS_ANY.search(ln)]
    if len(lines) != 1:
        return [f"expected one 'Gate M10:' status line, found {len(lines)}"]
    sm = STATUS.search(lines[0])
    if sm is None:
        return [f"status line not one of the three forms: {lines[0]!r}"]
    if sm.group(1) != "met":
        return []
    accepted = [
        e for e in ENTRY.finditer(_section(progress, "Answered by Khoa")) if e.group(5) == "accept"
    ]
    if not accepted:
        return ["'met' without an Answered-by-Khoa entry with verdict accept"]
    if report_sha is None or all(e.group(2) != report_sha for e in accepted):
        return ["'met' but no accept entry names the committed REPORT.md sha256"]
    return []


SHA = "ab" * 32


def _progress(status: str, entry: str | None = None) -> str:
    answered = f"\n{entry}\n" if entry else "\n- 2026-10-14: unrelated answer.\n"
    return (
        "# PROGRESS\n\n## Done (with evidence)\n\n### Gate M10 (M10-15, 2026-10-30)\n\n"
        f"{status}\n\n| bullet | status | evidence |\n|---|---|---|\n\n### Something else\n\n"
        f"## Answered by Khoa\n{answered}\n## Questions for Khoa\n\n- none\n"
    )


def _entry(verdict: str = "accept", sha: str = SHA, quote: str = "Report reads fine.") -> str:
    return (
        f'- 2026-10-31: **M10 pilot report** (REPORT.md sha256 {sha}, packet r1): "{quote}" '
        f"→ verdict {verdict}"
    )


@pytest.mark.parametrize(
    "status",
    ["Gate M10: awaiting Khoa's review", "Gate M10: not met: TE1", "Gate M10: not met: bullet 1"],
)
def test_checker_accepts_the_open_forms(status: str) -> None:
    assert gate_problems(_progress(status), SHA) == []


def test_checker_accepts_met_with_a_matching_accept_entry() -> None:
    assert gate_problems(_progress("Gate M10: met", _entry()), SHA) == []


@pytest.mark.parametrize(
    ("status", "entry", "sha"),
    [
        ("Gate M10: done", None, SHA),  # not one of the three forms
        ("Gate M10: met", None, SHA),  # example 11: no transcribed answer
        ("Gate M10: met", _entry("changes"), SHA),  # example 10: changes is not accept
        ("Gate M10: met", _entry(sha="cd" * 32), SHA),  # a verdict on another report
        ("Gate M10: met", _entry().replace('"Report reads fine."', "Report reads fine."), SHA),
        ("Gate M10: met", _entry(), None),  # no committed report to match
    ],
)
def test_checker_rejects_bad_records(status: str, entry: str | None, sha: str | None) -> None:
    assert gate_problems(_progress(status, entry), sha), (status, entry)


def test_checker_ignores_a_progress_without_gate_block() -> None:
    assert gate_problems("# PROGRESS\n\n## Done\n\n- nothing\n", None) == []


def test_two_status_lines_are_rejected() -> None:
    text = _progress("Gate M10: awaiting Khoa's review\nGate M10: met", _entry())
    assert gate_problems(text, SHA)


def test_progress_gate_record_shape() -> None:
    text = PROGRESS.read_text(encoding="utf-8")
    sha = hashlib.sha256(REPORT.read_bytes()).hexdigest() if REPORT.is_file() else None
    assert gate_problems(text, sha) == []
