#!/usr/bin/env python3
"""Milestone progress bars for PROGRESS.md (done / running / pending per issue).

An issue is an ``Mx-NN`` spec in docs/specs/. It is *done* when a commit subject on ``origin/main``
names it (squash merges are titled ``Mx-NN: …`` or ``feat(Mx-NN): …``), *running* when a commit
on another ``origin/*`` branch that is not yet on main names it, and *pending* otherwise. Branches
whose tip is the head of a closed or merged PR are leftovers and do not count as running. A split
issue's parts (``M8-06a``, ``M8-06b``) count for ``M8-06``, done once any part merged.
Milestones without issue specs yet (M5-M11) show as not designed.
M0 has no issue specs and is done (PR #1).

Usage: ``python scripts/dev/progress_bars.py [--write]`` (``git fetch`` first). With ``--write`` the
block between the ``<!-- progress:start -->`` and ``<!-- progress:end -->`` markers in PROGRESS.md
is replaced; without it the block is printed.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ISSUE = re.compile(r"\bM(\d+)-(\d{2})[a-z]?\b")
START, END = "<!-- progress:start -->", "<!-- progress:end -->"
MILESTONES = {  # GIAP_SPEC.md §13 milestone table
    0: "Scaffold",
    1: "Inputs and preflight",
    2: "Annotation wrappers",
    3: "Mapping and decoy FDR",
    4: "Gene table and region tests",
    5: "Autocorrelation and islands",
    6: "Operons",
    7: "Variants and novel ORFs",
    8: "Function and determinants",
    9: "Visualisation and report",
    10: "Pilot validation",
    11: "Release 1.0",
}
DONE, RUNNING, PENDING = "🟩", "🟨", "⬜"


def git(*args: str) -> str:
    proc = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True, check=True)
    return proc.stdout


def closed_heads() -> set[str]:
    """Head commits of closed and merged PRs (gh); empty when gh is unavailable."""
    try:
        out = subprocess.run(
            ["gh", "pr", "list", "--state", "closed", "--limit", "1000", "--json", "headRefOid"],
            cwd=REPO,
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return set()
    return {pr["headRefOid"] for pr in json.loads(out)}


def ids(text: str) -> set[str]:
    return {f"M{m}-{n}" for m, n in ISSUE.findall(text)}


def collect() -> tuple[dict[int, list[str]], set[str], set[str]]:
    issues: dict[int, list[str]] = {}
    for spec in sorted((REPO / "docs/specs").glob("M*-[0-9][0-9]-*.md")):
        m = ISSUE.match(spec.name)
        if m:
            issues.setdefault(int(m.group(1)), []).append(f"M{m.group(1)}-{m.group(2)}")
    # Squash merges are titled "Mx-NN: … (#PR)" or "feat(Mx-NN): … (#PR)"; other subjects
    # (docs/test/fix follow-ups, partial merges such as "Mx-NN part B: …") merely mention issues.
    merged = re.compile(r"^(?:feat\()?M(\d+)-(\d{2})[a-z]?\)?: .*\(#\d+\)$")
    subjects = git("log", "origin/main", "--format=%s").splitlines()
    done = {f"M{m.group(1)}-{m.group(2)}" for s in subjects if (m := merged.match(s))}
    running: set[str] = set()
    stale = closed_heads()
    for line in git(
        "for-each-ref", "--format=%(refname:short) %(objectname)", "refs/remotes/origin"
    ).splitlines():
        ref, sha = line.split()
        if ref in ("origin/main", "origin/HEAD", "origin") or sha in stale:
            continue
        # A branch works on the issue its unmerged commits name most often.
        counts: dict[str, int] = {}
        for subject in git("log", f"origin/main..{ref}", "--format=%s").splitlines():
            for i in ids(subject) - done:
                counts[i] = counts.get(i, 0) + 1
        if counts:
            running.add(max(sorted(counts), key=lambda i: counts[i]))
    return issues, done, running


def render() -> str:
    issues, done, running = collect()
    rows = ["| Milestone | Progress | Done | Running | Pending |", "|---|---|---|---|---|"]
    tot = [0, 0, 0]
    for m, name in MILESTONES.items():
        ids_m = sorted(set(issues.get(m, [])), key=lambda s: int(s.split("-")[1]))
        if m == 0:
            rows.append(f"| M0 {name} | {DONE} | 1 | 0 | 0 |")
            tot[0] += 1
            continue
        if not ids_m:
            rows.append(f"| M{m} {name} | _not designed yet_ | - | - | - |")
            continue
        d = [i for i in ids_m if i in done]
        r = [i for i in ids_m if i in running]
        p = [i for i in ids_m if i not in done and i not in running]
        bar = DONE * len(d) + RUNNING * len(r) + PENDING * len(p)
        pct = round(100 * len(d) / len(ids_m))
        rows.append(f"| M{m} {name} | {bar} {pct} % | {len(d)} | {len(r)} | {len(p)} |")
        tot = [tot[0] + len(d), tot[1] + len(r), tot[2] + len(p)]
    n = sum(tot)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    return "\n".join(
        [
            START,
            f"{DONE} done (merged) · {RUNNING} running (open branch) · {PENDING} pending; "
            "designed issues: "
            f"{tot[0]}/{n} done ({round(100 * tot[0] / n)} %), {tot[1]} running · "
            f"generated {stamp} by `scripts/dev/progress_bars.py`",
            "",
            *rows,
            END,
        ]
    )


def main() -> int:
    block = render()
    if "--write" not in sys.argv:
        print(block)
        return 0
    path = REPO / "PROGRESS.md"
    text = path.read_text(encoding="utf-8")
    if START in text and END in text:
        text = text[: text.index(START)] + block + text[text.index(END) + len(END) :]
    else:  # first time: insert before "## Current milestone"
        head = "## Milestone progress\n" + block + "\n\n## Current milestone"
        text = text.replace("## Current milestone", head, 1)
    path.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
