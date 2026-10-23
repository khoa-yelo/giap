"""Wording and error-estimate checks on the report (M9-16 behaviours 7, 8; GIAP_SPEC.md §1 rule
6, §10.3 "Wording"; M9-index decision 25). Pure; standard library only. A unit test runs them on
the templates, snippets and a rendered page; the ``report`` step never fails the run on them.
"""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence
from html.parser import HTMLParser
from typing import Any

__all__ = ["unlabelled_numbers", "wording"]

#: "detected", never "expressed": any case, at a word start, with the word after it
_EXPRESS = re.compile(r"\bexpress\w*(?:\W+\w+)?", re.IGNORECASE)
#: quoted citations, whose text is skipped
_QUOTES = frozenset({"q", "blockquote", "cite"})


class _Text(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.depth = 0
        self.parts: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.depth += tag in _QUOTES

    def handle_endtag(self, tag: str) -> None:
        if tag in _QUOTES and self.depth:
            self.depth -= 1

    def handle_data(self, data: str) -> None:
        if not self.depth:
            self.parts.append(data)


def wording(html: str) -> list[str]:
    """One text fragment per "express…" hit outside quoted citations, in document order."""
    parser = _Text()
    parser.feed(html)
    parser.close()
    return [m.group(0) for part in parser.parts for m in _EXPRESS.finditer(part)]


def unlabelled_numbers(cards: Sequence[Mapping[str, Any]]) -> list[str]:
    """The names of the cards whose value has no error estimate (nor ``unestimated``)."""
    return [str(c.get("name", "")) for c in cards if not str(c.get("estimate") or "").strip()]
