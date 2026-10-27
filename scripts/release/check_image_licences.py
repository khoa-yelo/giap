#!/usr/bin/env python3
"""Fail on a Python distribution whose licence is not allowed (M11-15 item 8a, slice 15g).

Runs inside a built gippa-step image (standard library only)::

    python check_image_licences.py --allow licences.allow SITE_PACKAGES [...]

Reads every ``*.dist-info/METADATA`` under the given directories. A distribution's licence is its
``License-Expression`` (PEP 639), else its ``License`` field when that is a one-line licence
expression, else all its ``License ::`` classifiers (each must be allowed). In an expression
``OR`` needs one allowed operand, ``AND`` every operand; ``X WITH Y`` is matched as a whole. A
distribution with no licence information is unlicensed and fails. Exit 0 when every distribution
passes, 1 otherwise (each offender named), 2 on a usage error.
"""

from __future__ import annotations

import argparse
import re
import sys
from email.parser import HeaderParser
from pathlib import Path

TOKEN = re.compile(r"\s*(\(|\)|[A-Za-z0-9][A-Za-z0-9.+\-:]*)")
CLASSIFIER = "License ::"


def read_allow(path: Path) -> frozenset[str]:
    """Allowlist entries (lower-cased), ``#`` comments and blank lines dropped."""
    lines = (ln.split("#", 1)[0].strip() for ln in path.read_text(encoding="utf-8").splitlines())
    return frozenset(ln.lower() for ln in lines if ln)


def tokenize(expr: str) -> list[str] | None:
    """SPDX-like tokens, or ``None`` when ``expr`` is not a licence expression."""
    out: list[str] = []
    pos, text = 0, expr.strip()
    while pos < len(text):
        m = TOKEN.match(text, pos)
        if m is None:
            return None
        out.append(m.group(1))
        pos = m.end()
        while pos < len(text) and text[pos].isspace():
            pos += 1
    return out or None


def evaluate(tokens: list[str], allow: frozenset[str]) -> bool | None:
    """Whether the expression is allowed; ``None`` when it does not parse."""
    pos = 0

    def peek() -> str | None:
        return tokens[pos].upper() if pos < len(tokens) else None

    def atom() -> bool | None:
        nonlocal pos
        if peek() == "(":
            pos += 1
            value = disjunction()
            if value is None or peek() != ")":
                return None
            pos += 1
            return value
        if peek() in {None, ")", "AND", "OR", "WITH"}:
            return None
        ident = tokens[pos]
        pos += 1
        if peek() == "WITH":
            if pos + 1 >= len(tokens):
                return None
            ident = f"{ident} WITH {tokens[pos + 1]}"
            pos += 2
        return ident.lower() in allow

    def conjunction() -> bool | None:
        nonlocal pos
        value = atom()
        while value is not None and peek() == "AND":
            pos += 1
            right = atom()
            value = None if right is None else value and right
        return value

    def disjunction() -> bool | None:
        nonlocal pos
        value = conjunction()
        while value is not None and peek() == "OR":
            pos += 1
            right = conjunction()
            value = None if right is None else value or right
        return value

    value = disjunction()
    return value if pos == len(tokens) else None


def verdict(fields: dict[str, list[str]], allow: frozenset[str]) -> tuple[bool, str]:
    """``(allowed, licence as read)`` for one distribution's METADATA fields."""
    for key in ("License-Expression", "License"):
        value = next((v.strip() for v in fields.get(key, []) if v.strip()), "")
        if not value or "\n" in value:
            continue
        tokens = tokenize(value)
        ok = evaluate(tokens, allow) if tokens else None
        if ok is not None:
            return ok, value
        if key == "License-Expression":
            return False, f"{value} (unparsable expression)"
    classifiers = [c.strip() for c in fields.get("Classifier", []) if c.startswith(CLASSIFIER)]
    if classifiers:
        return all(c.lower() in allow for c in classifiers), "; ".join(classifiers)
    seen = next((v.strip() for v in fields.get("License", []) if v.strip()), "")
    if seen:
        return False, f"no licence expression or classifier (License: {seen[:60]!r})"
    return False, "no licence information"


def metadata_fields(path: Path) -> dict[str, list[str]]:
    """METADATA header fields by title-cased name (a field may repeat)."""
    msg = HeaderParser().parsestr(path.read_text(encoding="utf-8", errors="replace"))
    out: dict[str, list[str]] = {}
    for key, value in msg.items():
        out.setdefault(key.title(), []).append(str(value))
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--allow", type=Path, required=True, help="the licences.allow file")
    ap.add_argument("sites", type=Path, nargs="+", help="site-packages directories")
    args = ap.parse_args(argv)
    allow = read_allow(args.allow)
    checked, bad = 0, []
    for site in args.sites:
        for meta in sorted(site.glob("*.dist-info/METADATA")):
            fields = metadata_fields(meta)
            name = (fields.get("Name") or [meta.parent.name])[0]
            version = (fields.get("Version") or ["?"])[0]
            ok, licence = verdict(fields, allow)
            checked += 1
            if not ok:
                bad.append(f"{name} {version}: {licence}")
    for line in bad:
        print(f"check_image_licences: not allowed: {line}", file=sys.stderr)
    print(f"check_image_licences: {checked} distributions checked, {len(bad)} not allowed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
