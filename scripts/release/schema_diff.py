#!/usr/bin/env python3
"""Schema diff test (M11-14h; docs/specs/M11-14-container-install.md §4; ADR-0060 §3, §6).

``check(root)`` compares the output schema of the working tree of ``root`` with its copy at the
newest ``v*`` tag reachable from HEAD, read with ``git show <tag>:<path>`` (no network, no
checkout):

- the ``fields[].name`` of every ``src/giap/schemas/tables/*.schema.json`` (a table file that
  disappears removes its columns, a new one adds them);
- the top-level ``properties`` keys of ``src/giap/schemas/run_manifest.schema.json`` (table
  ``run_manifest``).

A removed column (a rename is a removal plus an addition) is a breaking change: it needs a larger
major of both ``__version__`` (``src/giap/__init__.py``) and ``output_schema_version``
(``src/giap/schemas/datapackage.json``); an added one needs a larger ``output_schema_version``
(minor, or major). A shallow clone is a problem (``fetch tags: ...``),
checked first; no reachable ``v*`` tag raises ``NoPreviousRelease``. Exit codes of the command
line: 0 compatible or no previous release, 1 a problem.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path

TABLES = "src/giap/schemas/tables"
MANIFEST = "src/giap/schemas/run_manifest.schema.json"
DATAPACKAGE = "src/giap/schemas/datapackage.json"
INIT = "src/giap/__init__.py"
FETCH_TAGS = "fetch tags: actions/checkout with fetch-depth: 0 and fetch-tags: true"

_VERSION = re.compile(r"""^__version__\s*=\s*["']([^"']+)["']""", re.MULTILINE)
_NUMS = re.compile(r"^\s*v?(\d+)(?:\.(\d+))?")


class NoPreviousRelease(LookupError):
    """No ``v*`` tag is reachable from HEAD (before the first release)."""


def _git(root: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", "-C", str(root), *args], check=False, capture_output=True, text=True
    )


def previous_tag(root: Path) -> str:
    """The newest ``v*`` tag reachable from HEAD; ``NoPreviousRelease`` when there is none."""
    res = _git(root, "describe", "--tags", "--abbrev=0", "--match", "v*", "HEAD")
    tag = res.stdout.strip()
    if res.returncode != 0 or not tag:
        raise NoPreviousRelease("no previous release (no v* tag reachable from HEAD)")
    return tag


def _numbers(text: str | None) -> tuple[int, int]:
    m = _NUMS.match(text or "")
    return (int(m.group(1)), int(m.group(2) or 0)) if m else (0, 0)


def _columns(doc: Mapping[str, object]) -> list[str]:
    fields = doc.get("fields")
    if not isinstance(fields, list):
        return []
    return [str(f["name"]) for f in fields if isinstance(f, dict) and "name" in f]


def _manifest_keys(doc: Mapping[str, object]) -> list[str]:
    props = doc.get("properties")
    return list(props) if isinstance(props, dict) else []


class _Tree:
    """The files ``check`` reads, from the working tree (``tag=None``) or a tag."""

    def __init__(self, root: Path, tag: str | None) -> None:
        self.root, self.tag = root, tag

    def text(self, rel: str) -> str | None:
        if self.tag is None:
            path = self.root / rel
            return path.read_text(encoding="utf-8") if path.is_file() else None
        res = _git(self.root, "show", f"{self.tag}:{rel}")
        return res.stdout if res.returncode == 0 else None

    def json(self, rel: str) -> dict[str, object]:
        text = self.text(rel)
        doc = json.loads(text) if text else {}
        return doc if isinstance(doc, dict) else {}

    def table_files(self) -> list[str]:
        if self.tag is None:
            return sorted(p.name for p in (self.root / TABLES).glob("*.schema.json"))
        res = _git(self.root, "ls-tree", "--name-only", f"{self.tag}:{TABLES}")
        names = res.stdout.split() if res.returncode == 0 else []
        return sorted(n for n in names if n.endswith(".schema.json"))

    def schema(self) -> dict[str, list[str]]:
        """Table (file stem, or ``run_manifest``) → its column names."""
        out = {
            name.removesuffix(".schema.json"): _columns(self.json(f"{TABLES}/{name}"))
            for name in self.table_files()
        }
        out["run_manifest"] = _manifest_keys(self.json(MANIFEST))
        return out

    def package_version(self) -> str | None:
        m = _VERSION.search(self.text(INIT) or "")
        return m.group(1) if m else None

    def output_schema_version(self) -> str | None:
        value = self.json(DATAPACKAGE).get("output_schema_version")
        return value if isinstance(value, str) else None


def diff(
    old: Mapping[str, list[str]], new: Mapping[str, list[str]]
) -> tuple[list[tuple[str, str]], list[tuple[str, str]]]:
    """``(removed, added)`` ``(table, column)`` pairs of ``new`` against ``old``; pure."""
    removed = [(t, c) for t, cols in old.items() for c in cols if c not in new.get(t, [])]
    added = [(t, c) for t, cols in new.items() for c in cols if c not in old.get(t, [])]
    return removed, added


def check(root: Path) -> list[str]:
    """The problems of ``root``'s working tree against the newest reachable ``v*`` tag."""
    if _git(root, "rev-parse", "--is-shallow-repository").stdout.strip() == "true":
        return [f"shallow clone, the previous release cannot be read; {FETCH_TAGS}"]
    tag = previous_tag(root)
    before, now = _Tree(root, tag), _Tree(root, None)
    removed, added = diff(before.schema(), now.schema())
    problems = []
    old_pkg, new_pkg = before.package_version(), now.package_version()
    old_out, new_out = before.output_schema_version(), now.output_schema_version()
    majors = (
        _numbers(new_pkg)[0] > _numbers(old_pkg)[0],
        _numbers(new_out)[0] > _numbers(old_out)[0],
    )
    if removed and not all(majors):
        problems += [
            f"{t}: column {c} removed or renamed since {tag} without a major bump of both "
            f"__version__ ({old_pkg} -> {new_pkg}) and output_schema_version "
            f"({old_out} -> {new_out})"
            for t, c in removed
        ]
    if added and _numbers(new_out) <= _numbers(old_out):
        problems += [
            f"{t}: column {c} added since {tag} without a minor bump of output_schema_version "
            f"({old_out} -> {new_out})"
            for t, c in added
        ]
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    try:
        problems = check(args.repo)
    except NoPreviousRelease as exc:
        print(f"schema_diff: skipped: {exc}")
        return 0
    for p in problems:
        print(f"schema_diff: {p}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
