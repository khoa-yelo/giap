"""Print the test files of CI shard I of N (``python3 scripts/ci/shard_tests.py I N``).

The pytest files under tests/ (default ``python_files``: ``test_*.py`` and ``*_test.py``), sorted,
dealt round-robin, so the shards partition the suite and the calibration files spread across them.
Stdlib only.
"""

from __future__ import annotations

import sys
from pathlib import Path


def test_files(root: Path) -> list[str]:
    tests = root / "tests"
    found = {p for pat in ("test_*.py", "*_test.py") for p in tests.rglob(pat)}
    return sorted(p.relative_to(root).as_posix() for p in found if "__pycache__" not in p.parts)


def shard(files: list[str], index: int, total: int) -> list[str]:
    if not 1 <= index <= total:
        raise ValueError(f"shard {index} of {total}")
    return files[index - 1 :: total]


def main(argv: list[str]) -> int:
    index, total = int(argv[0]), int(argv[1])
    root = Path(argv[2]) if len(argv) > 2 else Path.cwd()
    print("\n".join(shard(test_files(root), index, total)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
