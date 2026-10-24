#!/usr/bin/env python3
"""Fail if a diff changes a golden directory without updating its PROVENANCE.md (M8-12, M8-12b).

The directories are tests/data/golden/r/ and tests/data/golden/wallenius/. Regeneration
(scripts/slurm/golden_r.sbatch, the golden-r workflow, scripts/slurm/golden_wallenius.sbatch) is
the only way their data change, and it rewrites that directory's PROVENANCE.md.
Usage: check_golden_provenance.py [BASE] (default origin/main); compares BASE...HEAD.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Iterable

GOLDEN_R = ("tests/data/golden/r/", "tests/data/golden/wallenius/")


def violations(changed: Iterable[str]) -> list[str]:
    paths = list(changed)
    bad: list[str] = []
    for root in GOLDEN_R:
        prov = root + "PROVENANCE.md"
        if prov in paths:
            continue
        touched = [p for p in paths if p.startswith(root) and p != prov]
        bad += [f"{p}: changed without updating {prov}" for p in touched]
    return bad


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
