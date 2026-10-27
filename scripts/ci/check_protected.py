#!/usr/bin/env python3
"""Fail when a protected file changes without Khoa's approval (GIPPA_SPEC.md R3, §14.6, A.4).

Creating a protected file is allowed. Modifying, deleting or renaming one needs a commit trailer
``Gate-Change-Approved: ADR-NNNN`` naming an ADR whose status on the *base* branch is
``accepted-by-khoa`` (an ADR added in the same PR cannot approve it).
For src/gippa/config/defaults.yaml, a diff that only adds keys passes with any ADR reference.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
import sys
from pathlib import Path

ADDITIVE_OK = "src/gippa/config/defaults.yaml"


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def protected_globs(base: str) -> list[str]:
    try:
        text = git("show", f"{base}:ci/protected.txt")
    except subprocess.CalledProcessError:
        text = Path("ci/protected.txt").read_text()
    return [ln.split("#")[0].strip() for ln in text.splitlines() if ln.split("#")[0].strip()]


def adr_status(base: str, adr: str) -> str | None:
    """Status of an ADR as it exists on the base branch, so a PR cannot approve itself."""
    names = git("ls-tree", "--name-only", f"{base}:docs/adr").splitlines()
    for name in names:
        if name.startswith(f"{adr}-") and name.endswith(".md"):
            text = git("show", f"{base}:docs/adr/{name}")
            m = re.search(r"^- Status:\s*(\S+)", text, re.M)
            if m:
                return m.group(1)
    return None


def only_adds_lines(base: str, path: str) -> bool:
    diff = git("diff", "--unified=0", f"{base}...HEAD", "--", path)
    removed = [ln for ln in diff.splitlines() if ln.startswith("-") and not ln.startswith("---")]
    return not removed


def main() -> int:
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    globs = protected_globs(base)
    changed = git("diff", "--name-status", "-M", f"{base}...HEAD").splitlines()
    modified: list[str] = []
    for ln in changed:
        fields = ln.split("\t")
        if fields[0][:1] in {"M", "D", "T"}:
            modified.append(fields[1])
        elif fields[0][:1] in {"R", "C"}:  # both the old and the new path count
            modified.extend(fields[1:] if fields[0][:1] == "R" else fields[2:])
    hits = [p for p in modified if any(fnmatch.fnmatch(p, g) for g in globs)]
    if not hits:
        print("check_protected: no protected file modified")
        return 0
    log = git("log", "--format=%B", f"{base}..HEAD")
    adrs = re.findall(r"^Gate-Change-Approved:\s*(ADR-\d{4})", log, re.M)
    approved = [a for a in adrs if adr_status(base, a) == "accepted-by-khoa"]
    failed = []
    for path in hits:
        if approved:
            continue
        if path == ADDITIVE_OK and adrs and only_adds_lines(base, path):
            continue
        failed.append(path)
    for path in failed:
        print(f"protected file modified without an accepted-by-khoa ADR: {path}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
