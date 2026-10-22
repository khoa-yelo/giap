#!/usr/bin/env python3
"""Fail if a diff changes tests/data/golden/r/ without updating its PROVENANCE.md (M8-12).

Regeneration (scripts/slurm/golden_r.sbatch or the golden-r workflow) is the only way the R
golden data change, and it rewrites PROVENANCE.md.  Usage: check_golden_provenance.py [BASE]
(default origin/main); compares BASE...HEAD.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Iterable

GOLDEN_R = "tests/data/golden/r/"
PROVENANCE = GOLDEN_R + "PROVENANCE.md"


def violations(changed: Iterable[str]) -> list[str]:
    paths = list(changed)
    touched = [p for p in paths if p.startswith(GOLDEN_R) and p != PROVENANCE]
    if not touched or PROVENANCE in paths:
        return []
    return [f"{p}: changed without updating {PROVENANCE}" for p in touched]


def main() -> int:
    base = sys.argv[1] if len(sys.argv) > 1 else "origin/main"
    out = subprocess.run(
        ["git", "diff", "--name-only", f"{base}...HEAD"], capture_output=True, text=True, check=True
    ).stdout
    bad = violations(line for line in out.splitlines() if line)
    for line in bad:
        print(line)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
