#!/usr/bin/env python3
"""Fail if a shell script uses a construct that bash 3.2 (macOS /bin/bash) cannot run.

PR #255 failed on the macOS CI runner over ``declare -A``; this guard catches such constructs
before CI does. Flagged: associative arrays (``declare/local/typeset -A``), ``mapfile`` and
``readarray``, case modification (``${x,,} ${x^^} ${x,} ${x^}``), ``|&``, the case terminators
``;&`` and ``;;&``, ``coproc``, ``${x@Q}``-style transforms, negative array indices
``${a[-1]}`` and ``shopt -s globstar``.

Comment lines, trailing comments and here-document bodies are ignored. A file that runs only on
Linux opts out with a comment line ``# bash32: exempt (<reason>)``; exempt files are printed
and do not fail the check.

Usage: check_bash32.py [PATH ...] (default ``scripts`` and ``workflow/bin`` when it exists).
Directories are scanned recursively for ``*.sh``, ``*.sbatch`` and files whose shebang names
bash; files named explicitly are always scanned. Prints ``path:line: construct`` per finding
and exits 1 when there is any. Standard library only.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

DEFAULT_PATHS = ("scripts", "workflow/bin")
SUFFIXES = {".sh", ".sbatch"}
SHEBANG_RE = re.compile(r"^#!\S*(?:/|\s)(?:env\s+)?bash(?:\s|$)|^#!\S*/bash$")
EXEMPT_RE = re.compile(r"^\s*#\s*bash32:\s*exempt\s*\((?P<reason>[^)]*)\)")
HEREDOC_RE = re.compile(r"(?<!<)<<(?P<dash>-?)\s*(?P<q>['\"]?)(?P<word>[A-Za-z_]\w*)(?P=q)")

_CMD = r"(?:^|(?<=[\s;&|(`]))"  # start of a word in command position (roughly)
_NAME = r"(?:[A-Za-z_]\w*|[0-9]+|[@*])"
_INDEX = r"(?:\[[^\]]*\])?"

CONSTRUCTS: tuple[tuple[str, re.Pattern[str]], ...] = (
    (
        "associative array (declare -A)",
        re.compile(_CMD + r"(?:declare|local|typeset)\s+(?:-[A-Za-z]+\s+)*-[A-Za-z]*A[A-Za-z]*\b"),
    ),
    ("mapfile/readarray", re.compile(_CMD + r"(?:mapfile|readarray)(?=\s|$)")),
    ("case modification (${x,,} ${x^^})", re.compile(r"\$\{!?" + _NAME + _INDEX + r"(?:\^|,)")),
    ("pipe with stderr (|&)", re.compile(r"(?<!\|)\|&")),
    ("case continue (;;&)", re.compile(r";;&")),
    ("case fallthrough (;&)", re.compile(r"(?<!;);&")),
    ("coproc", re.compile(_CMD + r"coproc(?=\s|$)")),
    ("parameter transform (${x@Q})", re.compile(r"\$\{" + _NAME + _INDEX + r"@[QEPAKaULuk]\}")),
    ("negative array index (${a[-1]})", re.compile(r"\$\{[A-Za-z_]\w*\[\s*-\s*[0-9]")),
    ("shopt -s globstar", re.compile(_CMD + r"shopt\s+(?:-\w+\s+)*-s\b[^;&|]*\bglobstar\b")),
)


def strip_comment(line: str) -> str:
    """Return ``line`` without its comment; ``#`` inside quotes or words does not start one."""
    single = double = False
    i = 0
    while i < len(line):
        c = line[i]
        if single:
            single = c != "'"
        elif c == "\\":
            i += 1
        elif double:
            double = c != '"'
        elif c == "'":
            single = True
        elif c == '"':
            double = True
        elif c == "#" and (i == 0 or line[i - 1] in " \t;&|()"):
            return line[:i]
        i += 1
    return line


def scan_text(text: str) -> list[tuple[int, str]]:
    """Return ``(1-based line, construct)`` for every bash-4 construct in one script's text."""
    found: list[tuple[int, str]] = []
    heredocs: list[tuple[str, bool]] = []  # pending terminators: (word, tabs stripped)
    for lineno, line in enumerate(text.splitlines(), 1):
        if heredocs:
            word, dash = heredocs[0]
            if (line.lstrip("\t") if dash else line) == word:
                heredocs.pop(0)
            continue
        code = strip_comment(line)
        for label, pattern in CONSTRUCTS:
            if pattern.search(code):
                found.append((lineno, label))
        heredocs.extend((m["word"], bool(m["dash"])) for m in HEREDOC_RE.finditer(code))
    return found


def exemption(text: str) -> str | None:
    """Return the reason of a ``# bash32: exempt (<reason>)`` comment line, if any."""
    for line in text.splitlines():
        m = EXEMPT_RE.match(line)
        if m and m["reason"].strip():
            return m["reason"].strip()
    return None


def is_bash_script(path: Path) -> bool:
    if path.suffix in SUFFIXES:
        return True
    try:
        with path.open("rb") as fh:
            first = fh.readline(256).decode("utf-8", errors="replace").rstrip("\r\n")
    except OSError:
        return False
    return bool(SHEBANG_RE.match(first))


def candidates(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for root in paths:
        if root.is_dir():
            files.extend(p for p in sorted(root.rglob("*")) if p.is_file() and is_bash_script(p))
        elif root.is_file():
            files.append(root)
    return files


def main(argv: list[str]) -> int:
    if argv:
        paths = [Path(a) for a in argv]
    else:
        paths = [Path(p) for p in DEFAULT_PATHS if Path(p).exists()]
    bad = False
    for path in candidates(paths):
        text = path.read_text(errors="replace")
        reason = exemption(text)
        if reason is not None:
            print(f"{path}: bash32 exempt ({reason})")
            continue
        for lineno, label in scan_text(text):
            print(f"{path}:{lineno}: {label}")
            bad = True
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
