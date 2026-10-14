#!/usr/bin/env python3
"""Fail if src/ imports a GPL, AGPL or unlicensed tool (GIAP_SPEC.md R5, §14.6)."""

from __future__ import annotations

import ast
import sys
from pathlib import Path

FORBIDDEN = {
    "bakta",
    "ostir",
    "promoter_calculator",
    "promotercalculator",
    "antismash",
    "promotech",
    "uniop",
    "transterm",
    "genomad",
}


def violations(root: Path) -> list[str]:
    found = []
    for path in sorted(root.rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module]
            for name in names:
                if name.split(".")[0].lower() in FORBIDDEN:
                    found.append(f"{path}:{node.lineno}: import of {name}")
    return found


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "src")
    bad = violations(root)
    for line in bad:
        print(line)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
