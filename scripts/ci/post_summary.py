#!/usr/bin/env python3
"""Post a nightly array's summary to the "Nightly (FarmShare)" issue (GIAP_SPEC.md §3.10, §12).

Usage: ``post_summary.py --job ARRAY_JOB [--tsv local/logs/nightly.tsv] [--expect 4] [--dry-run]``

``nightly.sbatch`` appends one line per sub-task to ``nightly.tsv`` (M2-21 lead decision 4:
``date  job_task  commit  exit  sub-task: summary``) and then calls this script. Only once every
expected sub-task of the array job has a line is the comment posted (``gh issue comment``); the
issue is created once if missing. A lock plus ``<tsv>.posted`` make the post happen once per job
even when two sub-tasks finish together. ``--dry-run`` prints the comment and touches nothing.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

TITLE = "Nightly (FarmShare)"


@dataclass(frozen=True)
class Line:
    """One ``nightly.tsv`` row."""

    date: str
    job: str
    task: str
    commit: str
    exit: str
    summary: str


def parse(text: str, job: str) -> list[Line]:
    """Rows of array job ``job`` (``<job>_<task>``); the last row per task wins."""
    rows: dict[str, Line] = {}
    for raw in text.splitlines():
        cells = raw.split("\t")
        if len(cells) < 4:
            continue
        jid, _, task = cells[1].partition("_")
        if jid != job:
            continue
        summary = cells[4] if len(cells) > 4 else ""
        rows[task] = Line(cells[0], jid, task, cells[2], cells[3], summary)
    return [rows[k] for k in sorted(rows, key=lambda t: (len(t), t))]


def render(lines: list[Line], job: str) -> str:
    """Markdown comment body for one array job."""
    ok = all(ln.exit == "0" for ln in lines)
    commits = sorted({ln.commit for ln in lines})
    head = f"**{'green' if ok else 'RED'}**: nightly array job {job} on `{', '.join(commits)}`"
    out = [head, "", "| task | exit | sub-task: summary |", "|---|---|---|"]
    out += [f"| {ln.task} | {ln.exit} | {ln.summary.replace('|', '/')} |" for ln in lines]
    out += [
        "",
        f"Logs: `local/logs/giap-nightly.{job}_<task>.out`; rows in `local/logs/nightly.tsv`.",
    ]
    return "\n".join(out) + "\n"


def _gh(*args: str, stdin: str | None = None) -> str:
    res = subprocess.run(["gh", *args], input=stdin, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"gh {' '.join(args[:2])} failed: {res.stderr.strip()}")
    return res.stdout


def issue_number(title: str, body: str) -> str:
    """Number of the open issue titled ``title``; created (once) when missing."""
    query = f'"{title}" in:title'
    args = ["issue", "list", "--state", "open", "--search", query, "--json", "number,title"]
    found = json.loads(_gh(*args, "--limit", "50"))
    for item in found:
        if item["title"] == title:
            return str(item["number"])
    url = _gh("issue", "create", "--title", title, "--body", body).strip()
    return url.rstrip("/").rsplit("/", 1)[-1]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--job", required=True, help="Slurm array job ID")
    ap.add_argument("--tsv", type=Path, default=Path("local/logs/nightly.tsv"))
    ap.add_argument("--expect", type=int, default=1, help="sub-tasks to wait for")
    ap.add_argument("--title", default=TITLE)
    ap.add_argument("--dry-run", action="store_true", help="print the comment, post nothing")
    args = ap.parse_args(argv)
    if not args.tsv.is_file():
        print(f"post_summary: {args.tsv} does not exist", file=sys.stderr)
        return 2
    lines = parse(args.tsv.read_text(), args.job)
    if len(lines) < args.expect:
        print(f"post_summary: {len(lines)}/{args.expect} sub-tasks of {args.job} done; not yet")
        return 0
    body = render(lines, args.job)
    if args.dry_run:
        print(body, end="")
        return 0
    posted = args.tsv.with_name(args.tsv.name + ".posted")
    with open(args.tsv.with_name(args.tsv.name + ".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        done = posted.read_text().split() if posted.is_file() else []
        if args.job in done:
            print(f"post_summary: {args.job} already posted")
            return 0
        try:
            num = issue_number(args.title, "Daily results of scripts/slurm/nightly.sbatch.")
            _gh("issue", "comment", num, "--body-file", "-", stdin=body)
        except (OSError, RuntimeError, ValueError) as exc:
            print(f"post_summary: {exc}", file=sys.stderr)
            return 1
        with posted.open("a") as fh:
            fh.write(args.job + "\n")
    print(f"post_summary: posted {args.job} to issue #{num}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
