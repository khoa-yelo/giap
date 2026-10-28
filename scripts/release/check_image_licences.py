#!/usr/bin/env python3
"""Fail on a Python distribution whose licence is not allowed (M11-15 item 8a, slice 15g).

Runs inside a built gippa-step image (standard library only)::

    python check_image_licences.py --allow licences.allow [--require NAME ...] SITE_PACKAGES [...]

Reads every ``*.dist-info/METADATA`` under the given directories. ``License-Expression``
(PEP 639) is authoritative. Without it every signal present must pass: a one-line ``License``
that parses as an expression and every ``License ::`` classifier (an unparsable one-line
``License`` is ignored next to classifiers and fails alone). ``License`` text naming GPL, General
Public, Affero or non-commercial fails unless the allowlist overrides the distribution by name
(``name:<dist>``). In an expression ``OR`` needs one allowed operand, ``AND`` every operand;
``X WITH Y`` is matched as a whole. A distribution with no licence information fails. Exit 0 when
every distribution passes and at least one was checked (plus each ``--require NAME``), 1
otherwise (each offender named), 2 on a usage error (a site that is not a directory).
"""

from __future__ import annotations

import argparse
import re
import sys
from email.parser import HeaderParser
from pathlib import Path

TOKEN = re.compile(r"\s*(\(|\)|[A-Za-z0-9][A-Za-z0-9.+\-:]*)")
CLASSIFIER = "License ::"
PROSE = re.compile(r"GPL|General Public|Affero|non-?commercial", re.IGNORECASE)
OVERRIDE = "name:"


def read_allow(path: Path) -> tuple[frozenset[str], frozenset[str]]:
    """``(licence entries, by-name overrides)``, lower-cased; ``#`` comments and blanks dropped.

    A ``name:<distribution>`` line overrides the GPL-prose check for that distribution only
    (PEP 503 normalised name); its licence signals are still checked.
    """
    lines = (ln.split("#", 1)[0].strip() for ln in path.read_text(encoding="utf-8").splitlines())
    entries = [ln.lower() for ln in lines if ln]
    names = {normalise(e[len(OVERRIDE) :]) for e in entries if e.startswith(OVERRIDE)}
    return frozenset(e for e in entries if not e.startswith(OVERRIDE)), frozenset(names)


def normalise(name: str) -> str:
    """PEP 503 normalised distribution name."""
    return re.sub(r"[-_.]+", "-", name.strip()).lower()


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


def verdict(
    fields: dict[str, list[str]], allow: frozenset[str], prose_ok: bool = False
) -> tuple[bool, str]:
    """``(allowed, licence as read)`` for one distribution's METADATA fields.

    ``prose_ok``: the distribution is overridden by name for GPL-like words in its License text.
    """
    expr = next((v.strip() for v in fields.get("License-Expression", []) if v.strip()), "")
    if expr:  # PEP 639: authoritative, nothing else is read
        tokens = tokenize(expr) if "\n" not in expr else None
        ok = evaluate(tokens, allow) if tokens else None
        return (False, f"{expr} (unparsable expression)") if ok is None else (ok, expr)
    licence = next((v.strip() for v in fields.get("License", []) if v.strip()), "")
    classifiers = [c.strip() for c in fields.get("Classifier", []) if c.startswith(CLASSIFIER)]
    if licence and not prose_ok and (m := PROSE.search(licence)):
        return False, f"License text names {m.group(0)!r} (License: {licence[:60]!r})"
    signals = 0
    if licence and "\n" not in licence:  # one line: a signal when it parses as an expression
        tokens = tokenize(licence)
        ok = evaluate(tokens, allow) if tokens else None
        if ok is False:
            return False, licence
        if ok is None and not classifiers:
            return False, f"{licence} (unparsable License and no classifier)"
        signals += ok is True
    if classifiers:
        bad = [c for c in classifiers if c.lower() not in allow]
        if bad:
            return False, "; ".join(bad)
        signals += 1
    if signals:
        return True, "; ".join([licence, *classifiers] if "\n" not in licence else classifiers)
    if licence:
        return False, f"no licence expression or classifier (License: {licence[:60]!r})"
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
    ap.add_argument(
        "--require",
        action="append",
        default=[],
        metavar="NAME",
        help="a distribution that must be among the checked ones (repeatable)",
    )
    ap.add_argument("sites", type=Path, nargs="+", help="site-packages directories")
    args = ap.parse_args(argv)
    missing = [str(s) for s in args.sites if not s.is_dir()]
    if missing:
        ap.error(f"not a directory: {', '.join(missing)}")
    allow, overrides = read_allow(args.allow)
    checked, bad, seen = 0, [], set()
    for site in args.sites:
        for meta in sorted(site.glob("*.dist-info/METADATA")):
            fields = metadata_fields(meta)
            name = (fields.get("Name") or [meta.parent.name])[0]
            version = (fields.get("Version") or ["?"])[0]
            ok, licence = verdict(fields, allow, prose_ok=normalise(name) in overrides)
            checked += 1
            seen.add(normalise(name))
            if not ok:
                bad.append(f"{name} {version}: {licence}")
    bad += [f"{req}: required but not found" for req in args.require if normalise(req) not in seen]
    if not checked:
        bad.append("no distribution found under the given site directories")
    for line in bad:
        print(f"check_image_licences: not allowed: {line}", file=sys.stderr)
    print(f"check_image_licences: {checked} distributions checked, {len(bad)} not allowed")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
