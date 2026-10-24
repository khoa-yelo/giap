"""Build the M10 pilot review packet ``OUT/r<k>/`` (docs/specs/M10-15-pilot-gate.md, D3).

    python validation/pilot/review_packet.py --metrics-dir $GIAP_ROOT/pilot/metrics/freeze-<F7> \\
        --runbook validation/pilot/reproduce.md --round K --out local/review/m10-pilot \\
        [--repo .] --a1 R/out_a1 --a2 R/out_a2 [--a3 R/out_a3] --cmp R/cmp \\
        --verify verify-<F7>.json --bench-summary $GIAP_ROOT/bench/<id>/summary.md \\
        --adr docs/adr/ADR-0051-....md [--adr ...]

It copies files and hashes them; nothing is written before checks 1-6 pass. Exit codes: 0 the
round was written; 1 a failed check (re-render differs, commits differ, a ``**missing:`` or
``pending`` status in the report, an empty evidence string, the round exists); 2 a missing
input. The packet holds per-gene HTML, so it lives under the gitignored ``local/``.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import re
import shutil
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from validation.pilot.metrics import (
    ADR_STATUS,
    Inconsistent,
    Precondition,
    _json,
    _load,
    render_report,
)
from validation.pilot.reproduce_check import parse as parse_runbook

#: committed files, relative to ``--repo``
COMMITTED = {
    "REPORT.md": "validation/pilot/REPORT.md",
    "metrics.json": "validation/pilot/metrics/metrics.json",
    "report_notes.md": "validation/pilot/report_notes.md",
    "a3_feasibility.json": "validation/pilot/metrics/a3_feasibility.json",
}
TEMPLATE = "validation/pilot/report_template.md"
PROGRESS = "PROGRESS.md"
FREEZE_STEPS = re.compile(r"S(8|9|1[0-4])")
PENDING_CELL = re.compile(r"\|\s*(pending|partial)\b")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _text(path: Path) -> str:
    return _load(path, lambda p: p.read_text(encoding="utf-8"))


def a3_reported(checklist: Sequence[Mapping[str, Any]]) -> bool:
    """Bullet 1 reports an A3 run: ``pass (A3 = …)`` or ``fail (A3 exit <n>)``."""
    status = str(checklist[0].get("status", "")) if checklist else ""
    return status.startswith("pass") or status.startswith("fail (A3 exit")


def sources(a: argparse.Namespace, metrics: Mapping[str, Any]) -> dict[str, Path]:
    """Packet path -> source file (check 5's table)."""
    files = {name: a.repo / rel for name, rel in COMMITTED.items()}
    files |= {
        "runs/A1/report.html": a.a1 / "report" / "report.html",
        "runs/A2/report.html": a.a2 / "report" / "report.html",
        "runs/compare.html": a.cmp / "report" / "compare.html",
        "reproduce.md": a.runbook,
        a.verify.name: a.verify,
        "benchmark/summary.md": a.bench_summary,
    }
    if a3_reported(metrics.get("checklist") or []):
        if a.a3 is None:
            raise Precondition("bullet 1 reports an A3 run but --a3 is not given")
        files["runs/A3/report.html"] = a.a3 / "report" / "report.html"
    files |= {f"adrs/{p.name}": p for p in a.adr}
    return files


def check_problems(
    a: argparse.Namespace, metrics: Mapping[str, Any], report: str, rerendered: str
) -> list[str]:
    """Checks 1-4 and 6 (every one a failure, exit 1)."""
    out = []
    if rerendered != report:
        out.append("1: render --final does not reproduce the committed REPORT.md")
    if _text(a.metrics_dir / "metrics.json") != _text(a.repo / COMMITTED["metrics.json"]):
        out.append("1: the committed metrics.json differs from the freeze metrics directory's")
    rb = parse_runbook(_text(a.runbook))
    if rb.status != "final":
        out.append(f"2: the runbook is status={rb.status or 'absent'}, not final")
    if metrics.get("giap_commit") != rb.freeze_commit:
        out.append(
            f"2: giap_commit {metrics.get('giap_commit')} != freeze_commit {rb.freeze_commit}"
        )
    for row in rb.evidence:
        if FREEZE_STEPS.fullmatch(row.get("step", "")) and row.get("commit") != rb.freeze_commit:
            out.append(f"2: Evidence {row.get('step')} ran at {row.get('commit')}")
    if "**missing:" in report:
        out.append("3: REPORT.md holds a '**missing:' marker")
    if PENDING_CELL.search(report):
        out.append("3: REPORT.md holds a pending status")
    for i, row in enumerate(metrics.get("checklist") or []):
        if not str(row.get("evidence") or "").strip():
            out.append(f"4: checklist row {i} has no evidence")
    if (a.out / f"r{a.round}").exists():
        out.append(f"6: round r{a.round} exists; rounds are never overwritten")
    return out


def status_table(report: str) -> list[str]:
    section = report.split("## 1. Status", 1)[-1].split("\n## ", 1)[0]
    return [ln for ln in section.splitlines() if ln.startswith("|")]


def index_md(a: argparse.Namespace, report: str, progress: str) -> str:
    """The round's index: status table, the diff against the previous round, the open M10
    questions for Khoa, the ADR status lines and how to answer."""
    prev = sorted(
        (int(p.name[1:]) for p in a.out.glob("r*") if re.fullmatch(r"r\d+", p.name)),
    )
    prev = [k for k in prev if k < a.round]
    if prev and (old := a.out / f"r{prev[-1]}" / "REPORT.md").is_file():
        diff = list(
            difflib.unified_diff(
                old.read_text(encoding="utf-8").splitlines(),
                report.splitlines(),
                f"r{prev[-1]}/REPORT.md",
                f"r{a.round}/REPORT.md",
                lineterm="",
            )
        )
        changed = ["```diff", *diff, "```"] if diff else ["REPORT.md is unchanged."]
    else:
        changed = ["First round."]
    questions = progress.split("## Questions for Khoa", 1)[-1].split("\n## ", 1)[0]
    open_m10 = [ln for ln in questions.splitlines() if ln.startswith(("- ", "|")) and "M10" in ln]
    adrs = []
    for p in a.adr:
        m = ADR_STATUS.search(p.read_text(encoding="utf-8"))
        adrs.append(f"- `{p.name}`: Status: {m.group(1) if m else 'absent'}")
    return "\n".join(
        [
            f"# M10 pilot report review, round r{a.round}",
            "",
            f"REPORT.md sha256 `{hashlib.sha256(report.encode('utf-8')).hexdigest()}`.",
            "",
            "## Status (REPORT.md section 1)",
            "",
            *status_table(report),
            "",
            "## What changed since the previous round",
            "",
            *changed,
            "",
            "## Questions for Khoa",
            "",
            "- K1: verdict on `REPORT.md` (sha256 above): accept, accept with changes, or reject"
            " (§14.8 item 4).",
            *(open_m10 or ["- No other open M10 question in PROGRESS.md."]),
            "",
            "## ADRs in this packet",
            "",
            *adrs,
            "",
            "## How to answer",
            "",
            "Answer in a session (or edit PROGRESS.md yourself) with accept, accept with changes,"
            " or reject, and any reasons. Only your own message counts as the review; it is"
            " transcribed verbatim with the sha256 above.",
            "",
        ]
    )


def build(a: argparse.Namespace) -> Path:
    """Checks, then writes ``OUT/r<k>/`` (via a temporary directory renamed at the end)."""
    metrics = _json(a.repo / COMMITTED["metrics.json"])
    files = sources(a, metrics)
    needed = [a.metrics_dir / "metrics.json", a.repo / TEMPLATE, a.repo / PROGRESS]
    for p in [*needed, *files.values()]:
        if not p.is_file():
            raise Precondition(f"missing input: {p}")
    report = _text(files["REPORT.md"])
    try:
        rerendered = render_report(
            _text(a.repo / TEMPLATE), metrics, _text(files["report_notes.md"]), final=True
        )
    except Inconsistent as err:
        rerendered = f"<refused: {err}>"
    if problems := check_problems(a, metrics, report, rerendered):
        raise Inconsistent("; ".join(problems))
    a.out.mkdir(parents=True, exist_ok=True)
    tmp = a.out / f".r{a.round}.tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    for rel, src in files.items():
        (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, tmp / rel)
    (tmp / "index.md").write_text(index_md(a, report, _text(a.repo / PROGRESS)), encoding="utf-8")
    rels = sorted(str(p.relative_to(tmp)) for p in tmp.rglob("*") if p.is_file())
    manifest = "".join(f"{sha256(tmp / r)}  {r}\n" for r in rels)
    (tmp / "MANIFEST.sha256").write_text(manifest, encoding="utf-8")
    final = a.out / f"r{a.round}"
    tmp.rename(final)
    return Path(final)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    for name in ("metrics-dir", "runbook", "out", "a1", "a2", "cmp", "verify", "bench-summary"):
        ap.add_argument(f"--{name}", type=Path, required=True)
    ap.add_argument("--round", type=int, required=True)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    ap.add_argument("--a3", type=Path, default=None)
    ap.add_argument("--adr", type=Path, action="append", default=[], required=True)
    a = ap.parse_args(argv)
    if a.round < 1:
        ap.error("--round must be >= 1")
    try:
        out = build(a)
    except Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except Inconsistent as err:
        print(f"check failed: {err}", file=sys.stderr)
        return 1
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
