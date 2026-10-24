"""Shape and consistency of the M10 gate record in ``PROGRESS.md`` (M10-15 acceptance criterion 7;
examples 1, 6, 10, 11; "The gate record (D7) and Khoa's review"; review rev-1015 S8).

Every ``### Gate M10`` block has one status line: ``Gate M10: awaiting Khoa's review`` /
``Gate M10: not met: <bullets>`` / ``Gate M10: met``.

* ``awaiting`` and ``met`` need the committed ``REPORT.md`` and bullets 1-4 ``pass`` and bullet 5
  ``pass`` or ``reported (...)`` in its section-1 table, unless the block cites a waiver
  ``Waiver: ADR-NNNN`` whose ``Status:`` is ``accepted-by-khoa`` (L3, R3, R13).
* ``met`` needs an Answered-by-Khoa entry ``- <date>: **M10 pilot report** (REPORT.md sha256
  <64 hex>, packet r<k>): "<quote>" → verdict accept`` whose sha256 equals the committed
  ``REPORT.md``'s. The quote may hold escaped quotes (``\\"``) and line breaks. The verdict on a
  sha256 is the one of its latest-dated entry (PROGRESS lists newest first); on a tie a
  ``changes`` or ``reject`` wins over ``accept``.

The test checks the record only: it cannot tell who wrote the quote (example 11), which is why
the recording rule binds the lead. The checker is defined here and tested on hand-made texts; in
15a ``PROGRESS.md`` has no gate block, so the check on the real file passes vacuously until 15b.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Callable

import pytest
from helpers.pilot_metrics_m10 import REPO

PROGRESS = REPO / "PROGRESS.md"
REPORT = REPO / "validation" / "pilot" / "REPORT.md"
ADR_DIR = REPO / "docs" / "adr"

GATE_HEAD = re.compile(r"^### Gate M10\b.*$", re.M)
STATUS = re.compile(r"Gate M10: (awaiting Khoa's review|not met: \S.*|met)\s*$")
STATUS_ANY = re.compile(r"Gate M10:")
ENTRY = re.compile(
    r"^- (\d{4}-\d{2}-\d{2})[^:\n]*: \*\*M10 pilot report\*\* \(REPORT\.md sha256 ([0-9a-f]{64}), "
    r"packet r(\d+)\): \"((?:[^\"\\]|\\.)+)\" → verdict (accept|changes|reject)[ \t]*$",
    re.M,
)
WAIVER = re.compile(r"Waiver: (ADR-\d{4})\b")
ADR_STATUS = re.compile(r"^-?\s*Status:\s*([\w-]+)", re.M)


def _section(text: str, heading: str) -> str:
    m = re.search(rf"^## {re.escape(heading)}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    return m.group(1) if m else ""


def _blocks(progress: str) -> list[str]:
    out = []
    for m in GATE_HEAD.finditer(progress):
        nxt = re.search(r"^#{2,3} ", progress[m.end() :], re.M)
        out.append(progress[m.end() : m.end() + nxt.start() if nxt else len(progress)])
    return out


def failing_bullets(report: str) -> list[int]:
    """1-based bullets of ``REPORT.md`` section 1 that are not satisfied (§15.6)."""
    section = report.split("## 1. Status", 1)[-1].split("\n## ", 1)[0]
    rows = [ln for ln in section.splitlines() if ln.startswith("|")][2:7]
    if len(rows) != 5:
        return [1, 2, 3, 4, 5]
    status = [[c.strip() for c in r.strip().strip("|").split("|")][1] for r in rows]
    bad = [i + 1 for i in range(4) if not status[i].startswith("pass")]
    return bad + ([] if status[4].startswith(("pass", "reported (")) else [5])


def verdict_of(progress: str, sha: str) -> str | None:
    """The verdict on ``sha``: its latest-dated entry; on a tie a non-accept wins."""
    entries = [e for e in ENTRY.finditer(_section(progress, "Answered by Khoa")) if e[2] == sha]
    if not entries:
        return None
    last = max(e[1] for e in entries)
    verdicts = {e[5] for e in entries if e[1] == last}
    return "accept" if verdicts == {"accept"} else sorted(verdicts - {"accept"})[0]


def gate_problems(
    progress: str, report: str | None, adr_status: Callable[[str], str | None] = lambda _: None
) -> list[str]:
    """Problems of every gate block (empty: fine, or no gate block)."""
    out: list[str] = []
    sha = None if report is None else hashlib.sha256(report.encode("utf-8")).hexdigest()
    for block in _blocks(progress):
        lines = [
            ln.strip().strip("*").strip() for ln in block.splitlines() if STATUS_ANY.search(ln)
        ]
        if len(lines) != 1:
            out.append(f"expected one 'Gate M10:' status line, found {len(lines)}")
            continue
        sm = STATUS.search(lines[0])
        if sm is None:
            out.append(f"status line not one of the three forms: {lines[0]!r}")
            continue
        if sm.group(1).startswith("not met"):
            continue
        if report is None:
            out.append(f"{sm.group(1)!r} without a committed REPORT.md")
            continue
        waived = any(adr_status(a) == "accepted-by-khoa" for a in WAIVER.findall(block))
        if (bad := failing_bullets(report)) and not waived:
            out.append(f"{sm.group(1)!r} while bullets {bad} fail and no Khoa-accepted waiver")
        if sm.group(1) == "met" and verdict_of(progress, str(sha)) != "accept":
            out.append("'met' but the latest Answered-by-Khoa verdict on REPORT.md is not accept")
    return out


def repo_adr_status(adr: str) -> str | None:
    files = sorted(ADR_DIR.glob(f"{adr}-*.md"))
    m = ADR_STATUS.search(files[0].read_text(encoding="utf-8")) if files else None
    return m.group(1).lower() if m else None


# ------------------------------------------------------------------ hand-made records


def _report(
    statuses: tuple[str, ...] = ("pass (A3 = Cascadia)", "pass", "pass", "pass", "pass"),
) -> str:
    rows = "".join(f"| bullet {i + 1} | {s} | evidence |\n" for i, s in enumerate(statuses))
    head = "| bullet | status | evidence |\n|---|---|---|\n"
    return f"# Report\n\n## 1. Status\n\n{head}{rows}\n## 2. Data\n"


REPORT_OK = _report()
SHA = hashlib.sha256(REPORT_OK.encode("utf-8")).hexdigest()
REPORT_NO_A3 = _report(("fail (A3 not run: route_unavailable)", "pass", "pass", "pass", "pass"))
SHA_NO_A3 = hashlib.sha256(REPORT_NO_A3.encode("utf-8")).hexdigest()


def _progress(status: str, entry: str | None = None, extra: str = "") -> str:
    answered = f"\n{entry}\n" if entry else "\n- 2026-10-14: unrelated answer.\n"
    return (
        "# PROGRESS\n\n## Done (with evidence)\n\n### Gate M10 (M10-15, 2026-10-30)\n\n"
        f"{status}\n{extra}\n| bullet | status | evidence |\n|---|---|---|\n\n"
        "### Something else\n\n"
        f"## Answered by Khoa\n{answered}\n## Questions for Khoa\n\n- none\n"
    )


def _entry(
    verdict: str = "accept",
    sha: str = SHA,
    quote: str = "Report reads fine.",
    date: str = "2026-10-31",
) -> str:
    return (
        f'- {date}: **M10 pilot report** (REPORT.md sha256 {sha}, packet r1): "{quote}" '
        f"→ verdict {verdict}"
    )


@pytest.mark.parametrize(
    "status",
    ["Gate M10: awaiting Khoa's review", "Gate M10: not met: TE1", "Gate M10: not met: bullet 1"],
)
def test_checker_accepts_the_open_forms(status: str) -> None:
    assert gate_problems(_progress(status), REPORT_OK) == []


def test_checker_accepts_met_with_a_matching_accept_entry() -> None:
    assert gate_problems(_progress("Gate M10: met", _entry()), REPORT_OK) == []


@pytest.mark.parametrize(
    ("status", "entry", "report"),
    [
        ("Gate M10: done", None, REPORT_OK),  # not one of the three forms
        ("Gate M10: met", None, REPORT_OK),  # example 11: no transcribed answer
        ("Gate M10: met", _entry("changes"), REPORT_OK),  # example 10: changes is not accept
        ("Gate M10: met", _entry(sha="cd" * 32), REPORT_OK),  # a verdict on another report
        (
            "Gate M10: met",
            _entry().replace('"Report reads fine."', "Report reads fine."),
            REPORT_OK,
        ),
        ("Gate M10: met", _entry(), None),  # no committed report to match
        ("Gate M10: awaiting Khoa's review", None, None),
    ],
)
def test_checker_rejects_bad_records(status: str, entry: str | None, report: str | None) -> None:
    assert gate_problems(_progress(status, entry), report), (status, entry)


def test_checker_ignores_a_progress_without_gate_block() -> None:
    assert gate_problems("# PROGRESS\n\n## Done\n\n- nothing\n", None) == []


def test_two_status_lines_are_rejected() -> None:
    text = _progress("Gate M10: awaiting Khoa's review\nGate M10: met", _entry())
    assert gate_problems(text, REPORT_OK)


# ------------------------------------------------------------------ S8: consistency with the report


@pytest.mark.parametrize("status", ["Gate M10: awaiting Khoa's review", "Gate M10: met"])
def test_s8_awaiting_or_met_while_a_bullet_fails_is_rejected(status: str) -> None:
    entry = _entry(sha=SHA_NO_A3)
    assert gate_problems(_progress(status, entry), REPORT_NO_A3)


@pytest.mark.parametrize(
    ("statuses", "ok"),
    [
        (
            ("pass (A3 = Cascadia)", "pass", "pass", "pass", "reported (1 not passed, explained)"),
            True,
        ),
        (("pass (A3 = Cascadia)", "fail", "pass", "pass", "pass"), False),  # example 4
        (("pass (A3 = Cascadia)", "pass", "pass", "fail (incomplete: sim5m_core)", "pass"), False),
        (("pass (A3 = Cascadia)", "pass", "pass", "pass", "fail"), False),
    ],
)
def test_s8_bullet_rules(statuses: tuple[str, ...], ok: bool) -> None:
    problems = gate_problems(_progress("Gate M10: awaiting Khoa's review"), _report(statuses))
    assert (problems == []) is ok, problems


def test_s8_not_met_while_a_bullet_fails_is_fine() -> None:
    assert gate_problems(_progress("Gate M10: not met: bullet 1"), REPORT_NO_A3) == []


@pytest.mark.parametrize(
    ("adr", "ok"), [("accepted-by-khoa", True), ("proposed", False), ("accepted", False)]
)
def test_s8_a_khoa_accepted_waiver_allows_awaiting_and_met(adr: str, ok: bool) -> None:
    """L3: example 6 can read ``met`` only when Khoa accepted both the waiver and the report."""
    text = _progress("Gate M10: met", _entry(sha=SHA_NO_A3), extra="\nWaiver: ADR-0123\n")
    problems = gate_problems(text, REPORT_NO_A3, {"ADR-0123": adr}.get)
    assert (problems == []) is ok, problems


# ------------------------------------------------------------------ S8: quotes


@pytest.mark.parametrize(
    "quote",
    [
        'Reads fine; reword the \\"H1\\" paragraph later.',
        "Reads fine.\n  Second paragraph of the answer, on its own line.",
    ],
)
def test_s8_quote_with_escaped_quote_or_line_break_is_recorded(quote: str) -> None:
    assert gate_problems(_progress("Gate M10: met", _entry(quote=quote)), REPORT_OK) == []


def test_s8_an_unescaped_quote_inside_the_quote_is_rejected() -> None:
    entry = _entry(quote='Reads fine; the "H1" paragraph.')
    assert gate_problems(_progress("Gate M10: met", entry), REPORT_OK)


# ------------------------------------------------------------------ S8: a later verdict wins


@pytest.mark.parametrize("later", ["changes", "reject"])
@pytest.mark.parametrize("newest_first", [True, False])
def test_s8_a_later_changes_or_reject_overrides_an_earlier_accept(
    later: str, newest_first: bool
) -> None:
    old, new = _entry(date="2026-10-31"), _entry(later, date="2026-11-02")
    entries = f"{new}\n{old}" if newest_first else f"{old}\n{new}"
    assert gate_problems(_progress("Gate M10: met", entries), REPORT_OK)


def test_s8_a_later_accept_overrides_an_earlier_changes() -> None:
    entries = f"{_entry(date='2026-11-02')}\n{_entry('changes', date='2026-10-31')}"
    assert gate_problems(_progress("Gate M10: met", entries), REPORT_OK) == []


def test_s8_same_day_changes_beats_accept() -> None:
    entries = f"{_entry()}\n{_entry('changes')}"
    assert gate_problems(_progress("Gate M10: met", entries), REPORT_OK)


# ------------------------------------------------------------------ S8: every block


def test_s8_every_gate_block_is_checked() -> None:
    text = _progress("Gate M10: awaiting Khoa's review", _entry())
    second = "### Gate M10 (M10-15, round r2)\n\nGate M10: done\n\n"
    text = text.replace("### Something else", second + "### Something else")
    assert gate_problems(text, REPORT_OK)


# ------------------------------------------------------------------ the real PROGRESS.md


def test_progress_gate_record_shape() -> None:
    text = PROGRESS.read_text(encoding="utf-8")
    report = REPORT.read_text(encoding="utf-8") if REPORT.is_file() else None
    assert gate_problems(text, report, repo_adr_status) == []
