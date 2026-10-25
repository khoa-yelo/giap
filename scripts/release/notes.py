"""Extract the release notes for one version from CHANGELOG.md (M11-04 build item 2).

Usage: python scripts/release/notes.py VERSION [--draft]

Run from the repository root: reads ``./CHANGELOG.md`` and writes ``./notes.md`` with the body of
the ``## GIAP <version> (...)`` section. An ``rc`` version reads its base version's section,
prefixed "Release candidate". With ``--draft`` a missing section falls back to
``towncrier build --draft`` (nothing is written but notes.md), headed ``DRAFT``. Exit 0 on
success, 1 when the section is missing (without ``--draft``) or empty.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

RC = re.compile(r"^(\d+\.\d+\.\d+)rc\d+$")
EMPTY = "No significant changes."


def base_version(version: str) -> str:
    """``1.0.0rc1`` -> ``1.0.0``; any other version is returned unchanged."""
    m = RC.match(version)
    return m.group(1) if m else version


def extract_section(changelog: str, version: str) -> str | None:
    """The body of ``## GIAP <version> (...)`` up to the next ``## `` heading; None if absent."""
    head = re.compile(rf"^## GIAP {re.escape(version)} \(")
    body: list[str] | None = None
    for line in changelog.splitlines():
        if body is None:
            if head.match(line):
                body = []
        elif line.startswith("## "):
            break
        else:
            body.append(line)
    return None if body is None else "\n".join(body).strip()


def towncrier_draft(version: str) -> str:
    """The towncrier draft for ``version`` rendered in the current directory (writes nothing)."""
    out = subprocess.run(
        [sys.executable, "-m", "towncrier", "build", "--draft", "--version", version],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return extract_section(out, version) or ""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("version")
    ap.add_argument("--draft", action="store_true", help="fall back to the towncrier draft")
    args = ap.parse_args(argv)
    version: str = args.version
    base = base_version(version)
    changelog = Path("CHANGELOG.md")
    text = changelog.read_text(encoding="utf-8") if changelog.is_file() else ""
    body = extract_section(text, base)
    header = ""
    if body is None:
        if not args.draft:
            print(f"no '## GIAP {base} (...)' section in CHANGELOG.md", file=sys.stderr)
            return 1
        body = towncrier_draft(base)
        header = f"DRAFT: towncrier draft for GIAP {base}, not yet in CHANGELOG.md\n\n"
    if not body or body == EMPTY:
        print(f"release notes for GIAP {base} are empty", file=sys.stderr)
        return 1
    if base != version:
        header += f"Release candidate {version} of GIAP {base}.\n\n"
    Path("notes.md").write_text(header + body + "\n", encoding="utf-8")
    print(f"notes.md: GIAP {version}{' (draft)' if header.startswith('DRAFT') else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
