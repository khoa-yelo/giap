#!/usr/bin/env python3
"""Enforce the separation of test and implementation roles (GIAP_SPEC.md R2, §14.6).

Fails if, within the PR, one Agent-Role changes both src/ (or workflow/) and test files,
an implementer commit modifies an existing test file, or a test-author commit touches
src/ or workflow/. Commits without an Agent-Role trailer fail too.
"""

from __future__ import annotations

import re
import subprocess
import sys
from collections import defaultdict

CODE = ("src/", "workflow/")
TESTS = ("tests/",)
EXEMPT_TESTS = ("tests/data/recorded/",)  # implementers may add recorded fixtures


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def main() -> int:
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    shas = git("rev-list", "--no-merges", f"{base}..HEAD").split()
    touched: dict[str, set[str]] = defaultdict(set)
    errors = []
    for sha in shas:
        body = git("log", "-1", "--format=%B", sha)
        m = re.search(r"^Agent-Role:\s*(\S+)", body, re.M)
        if not m:
            errors.append(f"{sha[:8]}: missing Agent-Role trailer")
            continue
        role = m.group(1)
        status = git("diff-tree", "--no-commit-id", "--name-status", "-r", sha).splitlines()
        for line in status:
            kind, path = line.split("\t")[0], line.split("\t")[-1]
            is_code = path.startswith(CODE)
            is_test = path.startswith(TESTS) and not path.startswith(EXEMPT_TESTS)
            if is_code:
                touched[role].add("code")
            if is_test:
                touched[role].add("tests")
            if role == "implementer" and is_test and kind.startswith("M"):
                errors.append(f"{sha[:8]}: implementer modified existing test {path}")
            if role == "test-author" and is_code:
                errors.append(f"{sha[:8]}: test-author touched {path}")
    for role, kinds in touched.items():
        if kinds == {"code", "tests"} and role not in {"lead"}:
            errors.append(f"role {role} changed both code and tests in this PR")
    for e in errors:
        print(e)
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
