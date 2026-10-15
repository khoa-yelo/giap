#!/usr/bin/env python3
"""Fail when two test modules share a basename (pytest rejects the second one).

Test directories have no ``__init__.py``, so ``tests/unit/a/test_x.py`` and
``tests/unit/b/test_x.py`` collide at collection ("import file mismatch").
Usage: ``python scripts/ci/check_test_names.py``.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEARCH = ("tests", "benchmarks")


def duplicates() -> dict[str, list[str]]:
    seen: dict[str, list[str]] = defaultdict(list)
    for top in SEARCH:
        for path in sorted((ROOT / top).rglob("test_*.py")):
            if (path.parent / "__init__.py").exists():
                continue  # inside a package: its module name is qualified
            seen[path.name].append(str(path.relative_to(ROOT)))
    return {name: paths for name, paths in seen.items() if len(paths) > 1}


def main() -> int:
    dups = duplicates()
    for name, paths in sorted(dups.items()):
        print(f"duplicate test module name {name}: {', '.join(paths)}", file=sys.stderr)
    return 1 if dups else 0


if __name__ == "__main__":
    raise SystemExit(main())
