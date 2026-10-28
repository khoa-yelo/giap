#!/usr/bin/env python3
"""Require a changelog fragment when a PR changes code (GIPPA_SPEC.md §14.2, M11-12).

``check_changelog.py BASE``: a PR whose ``BASE...HEAD`` diff touches ``src/`` or ``workflow/``
(``workflow/tests/`` excluded) must add a non-empty ``changelog.d/<PR>.<type>.md`` or
``changelog.d/+<slug>.<type>.md`` (M11-02). Exempt: the label ``no-changelog`` (from
``$GITHUB_EVENT_PATH``) or every non-merge commit subject in ``BASE..HEAD`` of type ``test``,
``ci``, ``chore`` or ``refactor``. Exit 0 on pass or exemption, 1 otherwise.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path

CODE = ("src/", "workflow/")
NOT_CODE = ("workflow/tests/",)
FRAGMENT = re.compile(r"^changelog\.d/(\+[a-z0-9-]+|\d+)\.(feature|bugfix|doc|removal|misc)\.md$")
EXEMPT_SUBJECT = re.compile(r"^(test|ci|chore|refactor)(\([^)]*\))?: ")
LABEL = "no-changelog"


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def diff_entries(base: str) -> list[tuple[str, str]]:
    """(status letter, path) per changed file; a rename or copy is listed at its new path."""
    out = git("diff", "--name-status", "-M", f"{base}...HEAD")
    entries = []
    for line in out.splitlines():
        fields = line.split("\t")
        if len(fields) >= 2:
            entries.append((fields[0][:1], fields[-1]))
    return entries


def is_code(path: str) -> bool:
    return path.startswith(CODE) and not path.startswith(NOT_CODE)


def has_exempt_label(event_path: str | None) -> bool:
    if not event_path:
        return False
    try:
        event = json.loads(Path(event_path).read_text())
    except (OSError, ValueError):
        return False
    pr = event.get("pull_request") if isinstance(event, dict) else None
    labels = pr.get("labels", []) if isinstance(pr, dict) else []
    return any(isinstance(lb, dict) and lb.get("name") == LABEL for lb in labels)


def maintenance_only(subjects: list[str]) -> bool:
    return bool(subjects) and all(EXEMPT_SUBJECT.match(s) for s in subjects)


def main() -> int:
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    entries = diff_entries(base)
    code = sorted({p for _, p in entries if is_code(p)})
    if not code:
        print("ok: no src/ or workflow/ changes")
        return 0
    if has_exempt_label(os.environ.get("GITHUB_EVENT_PATH")):
        print(f"exempt: label {LABEL}")
        return 0
    subjects = git("log", "--no-merges", "--format=%s", f"{base}..HEAD").splitlines()
    if maintenance_only(subjects):
        print("exempt: refactor/test/ci/chore only")
        return 0
    added = [p for s, p in entries if s in {"A", "R", "C"} and FRAGMENT.match(p)]
    filled = [p for p in added if git("show", f"HEAD:{p}").strip()]
    if filled:
        print(f"ok: fragment {', '.join(filled)}")
        return 0
    if added:
        print(f"empty fragment: {', '.join(added)}; write a one-line entry")
        return 1
    print(
        f"changelog fragment missing: changed {', '.join(code)} ; "
        f"add changelog.d/<PR>.<type>.md (or +<slug>.<type>.md), or label {LABEL} and push "
        "again (labels are read from the push event; a re-run reuses the old labels)"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
