"""Build the M10 pilot review packet ``OUT/r<k>/`` (docs/specs/M10-15-pilot-gate.md, D3).

    python validation/pilot/review_packet.py --metrics-dir $GIPPA_ROOT/pilot/metrics/freeze-<F7> \\
        --runbook validation/pilot/reproduce.md --round K --out local/review/m10-pilot \\
        [--repo .] --a1 R/out_a1 --a2 R/out_a2 [--a3 R/out_a3] --cmp R/cmp \\
        --verify verify-<F7>.json --bench-summary $GIPPA_ROOT/bench/<id>/summary.md \\
        --adr docs/adr/ADR-0051-....md [--adr ...]

It copies files and hashes them; nothing is written before checks 1-6 pass. Exit codes: 0 the
round was written; 1 a failed check (the ``render --final --runbook`` re-render refuses or
differs, the freeze ``metrics.json`` differs, commits differ, the round exists); 2 a missing or
unreadable input, or ``--out`` inside the checkout and not git-ignored. The packet holds per-gene
HTML, so it lives under the gitignored ``local/``.
"""

from __future__ import annotations

import sys

sys.dont_write_bytecode = True  # N6: no __pycache__ in the checkout

import argparse  # noqa: E402
import difflib  # noqa: E402
import hashlib  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
from collections.abc import Sequence  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

from validation.pilot.metrics import (  # noqa: E402
    Inconsistent,
    Precondition,
    _json,
    _load,
    accepted_evidence,
    adr_status,
    render_report,
    with_runbook,
)
from validation.pilot.reproduce_check import parse as parse_runbook  # noqa: E402

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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _text(path: Path) -> str:
    return _load(path, lambda p: p.read_text(encoding="utf-8"))


def _bytes(path: Path) -> bytes:
    return _load(path, lambda p: p.read_bytes())


def check_out(out: Path, repo: Path) -> None:
    """S3: ``--out`` lies outside the checkout or under a path the checkout git-ignores."""
    o, r = out.resolve(), repo.resolve()
    if o != r and r not in o.parents:
        return
    res = subprocess.run(
        ["git", "-C", str(r), "check-ignore", "-q", str(o)], capture_output=True, check=False
    )
    if res.returncode != 0:
        raise Precondition(f"--out {out} is inside {repo} and not git-ignored")


def sources(a: argparse.Namespace, bullet1: str) -> dict[str, Path]:
    """Packet path -> source file (check 5's table); A3's HTML when bullet 1 reports a run."""
    files = {name: a.repo / rel for name, rel in COMMITTED.items()}
    files |= {
        "runs/A1/report.html": a.a1 / "report" / "report.html",
        "runs/A2/report.html": a.a2 / "report" / "report.html",
        "runs/compare.html": a.cmp / "report" / "compare.html",
        "reproduce.md": a.runbook,
        a.verify.name: a.verify,
        "benchmark/summary.md": a.bench_summary,
    }
    if bullet1.startswith(("pass", "fail (A3 exit")):
        if a.a3 is None:
            raise Precondition("bullet 1 reports an A3 run but --a3 is not given")
        files["runs/A3/report.html"] = a.a3 / "report" / "report.html"
    return files | {f"adrs/{p.name}": p for p in a.adr}


def check_problems(a: argparse.Namespace, metrics: dict[str, Any], rb: Any) -> list[str]:
    """Checks 1, 2 and 6 (each a failure, exit 1); 3 and 4 are part of ``render --final``."""
    out = []
    try:
        text = render_report(
            _text(a.repo / TEMPLATE), metrics, _text(a.repo / COMMITTED["report_notes.md"]),
            final=True, runbook=rb,
        )  # fmt: skip
        if text.encode("utf-8") != _bytes(a.repo / COMMITTED["REPORT.md"]):
            out.append("1: render --final does not reproduce the committed REPORT.md")
    except Inconsistent as err:
        out.append(f"1: render --final refuses the committed inputs: {err}")
    if _bytes(a.metrics_dir / "metrics.json") != _bytes(a.repo / COMMITTED["metrics.json"]):
        out.append("1: the freeze directory's metrics.json differs from the committed one")
    if rb.status != "final":
        out.append(f"2: the runbook is status={rb.status or 'absent'}, not final")
    if metrics.get("gippa_commit") != rb.freeze_commit:
        out.append(f"2: gippa_commit {metrics.get('gippa_commit')} != freeze {rb.freeze_commit}")
    for step, row in accepted_evidence(rb).items():
        if FREEZE_STEPS.fullmatch(step) and row.get("commit") != rb.freeze_commit:
            out.append(f"2: Evidence {step} ran at {row.get('commit')}")
    if (a.out / f"r{a.round}").exists():
        out.append(f"6: round r{a.round} exists; rounds are never overwritten")
    return out


def index_md(a: argparse.Namespace, report: str, progress: str) -> str:
    """The round's index: status table, the diff against the previous round, the open M10
    questions for Khoa, the ADR status lines and how to answer."""
    rounds = (int(p.name[1:]) for p in a.out.glob("r*") if re.fullmatch(r"r\d+", p.name))
    prev = max((k for k in rounds if k < a.round), default=None)
    changed = ["First round."]
    if prev is not None and (old := a.out / f"r{prev}" / "REPORT.md").is_file():
        diff = difflib.unified_diff(
            old.read_text(encoding="utf-8").splitlines(), report.splitlines(),
            f"r{prev}/REPORT.md", f"r{a.round}/REPORT.md", lineterm="",
        )  # fmt: skip
        changed = ["```diff", *diff, "```"]
    section = report.split("## 1. Status", 1)[-1].split("\n## ", 1)[0]
    questions = progress.split("## Questions for Khoa", 1)[-1].split("\n## ", 1)[0]
    open_m10 = [ln for ln in questions.splitlines() if ln.startswith(("- ", "|")) and "M10" in ln]
    adrs = [f"- `{p.name}`: Status: {adr_status(_text(p)) or 'absent'}" for p in a.adr]
    return "\n".join(
        [
            f"# M10 pilot report review, round r{a.round}",
            "",
            f"REPORT.md sha256 `{hashlib.sha256(report.encode('utf-8')).hexdigest()}`.",
            "",
            "## Status (REPORT.md section 1)",
            "",
            *(ln for ln in section.splitlines() if ln.startswith("|")),
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
    """Checks, then writes ``OUT/r<k>/`` through ``OUT/.r<k>.tmp`` (removed on any error)."""
    check_out(a.out, a.repo)
    metrics = _json(a.repo / COMMITTED["metrics.json"])
    rb = parse_runbook(_text(a.runbook))
    files = sources(a, str(with_runbook(metrics, rb)["checklist"][0]["status"]))
    for p in [
        a.metrics_dir / "metrics.json",
        a.repo / TEMPLATE,
        a.repo / PROGRESS,
        *files.values(),
    ]:
        if not p.is_file():
            raise Precondition(f"missing input: {p}")
    if problems := check_problems(a, metrics, rb):
        raise Inconsistent("; ".join(problems))
    tmp, final = a.out / f".r{a.round}.tmp", a.out / f"r{a.round}"
    try:
        shutil.rmtree(tmp, ignore_errors=True)
        for rel, src in files.items():
            (tmp / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, tmp / rel)
        report = _text(files["REPORT.md"])
        (tmp / "index.md").write_text(index_md(a, report, _text(a.repo / PROGRESS)), "utf-8")
        rels = sorted(str(p.relative_to(tmp)) for p in tmp.rglob("*") if p.is_file())
        manifest = "".join(f"{sha256(tmp / r)}  {r}\n" for r in rels)
        (tmp / "MANIFEST.sha256").write_text(manifest, encoding="utf-8")
        tmp.rename(final)
    except (OSError, Precondition) as err:
        shutil.rmtree(tmp, ignore_errors=True)
        raise Precondition(f"writing {final}: {err}") from err
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
