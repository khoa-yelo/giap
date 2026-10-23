"""Link and offline checks on a rendered report page (M9-13; GIAP_SPEC.md §10.3).

``relative_links`` lists the relative URL references of a page, ready to resolve against the
report directory; ``external_loads`` lists the URLs a browser would fetch from outside the report
at view time. Both are pure, use the standard library only, and return each distinct value once in
first-occurrence document order.
"""

from __future__ import annotations

import re
from html.parser import HTMLParser
from urllib.parse import unquote

__all__ = ["external_loads", "relative_links"]

#: elements whose ``href``/``xlink:href`` is fetched when the page opens (``src`` always is)
_HREF_LOADERS = frozenset({"link", "script", "img", "iframe", "image", "use"})
_HREF_ATTRS = ("href", "xlink:href")
_SCHEME = re.compile(r"^([a-zA-Z][a-zA-Z0-9+.\-]*):")
# CSS references in document order: ``@import "x"`` / ``@import 'x'`` and ``url(...)`` (which also
# covers ``@import url(...)``)
_CSS_REF = re.compile(
    r"""@import\s+(?:"([^"]*)"|'([^']*)')"""
    r"""|url\(\s*(?:"([^"]*)"|'([^']*)'|([^)\s]*))\s*\)""",
    re.IGNORECASE,
)


def _clean(value: str) -> str:
    return value.strip().strip("\"'").strip()


def _css_refs(css: str) -> list[str]:
    return [_clean(next(g for g in m.groups() if g is not None)) for m in _CSS_REF.finditer(css)]


class _Refs(HTMLParser):
    """Collects (url, is_load) pairs in document order; comments and script bodies are ignored."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.refs: list[tuple[str, bool]] = []
        self._in_style = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        for name, value in attrs:
            if value is None:
                continue
            if name == "src":
                self.refs.append((_clean(value), True))
            elif name in _HREF_ATTRS:
                self.refs.append((_clean(value), tag in _HREF_LOADERS))
            elif name == "style":
                self.refs.extend((u, True) for u in _css_refs(value))
        if tag == "style":
            self._in_style += 1

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        self.handle_starttag(tag, attrs)
        if tag == "style":
            self._in_style -= 1

    def handle_endtag(self, tag: str) -> None:
        if tag == "style" and self._in_style:
            self._in_style -= 1

    def handle_data(self, data: str) -> None:
        if self._in_style:
            self.refs.extend((u, True) for u in _css_refs(data))


def _references(html: str) -> list[tuple[str, bool]]:
    parser = _Refs()
    parser.feed(html)
    parser.close()
    return parser.refs


def _is_external(url: str) -> bool:
    if url.startswith("//"):
        return True
    m = _SCHEME.match(url)
    return m is not None and m.group(1).lower() != "data"


def _is_relative(url: str) -> bool:
    return bool(url) and not url.startswith(("#", "//")) and _SCHEME.match(url) is None


def _unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(values))


def relative_links(html: str) -> list[str]:
    """Relative URL references of ``html`` (``href``, ``src``, ``xlink:href`` on any element and
    CSS ``url(...)``/``@import``), without ``#fragment`` or ``?query`` and percent-decoded."""
    out: list[str] = []
    for url, _ in _references(html):
        if not _is_relative(url):
            continue
        path = unquote(url.split("#", 1)[0].split("?", 1)[0])
        if path:
            out.append(path)
    return _unique(out)


def external_loads(html: str) -> list[str]:
    """URLs fetched from outside the report when ``html`` is opened: absolute (any scheme but
    ``data:``), protocol-relative or ``file:`` references in ``src``, in ``href`` of loading
    elements, and in CSS ``url(...)``/``@import``; values as written."""
    return _unique([u for u, load in _references(html) if load and _is_external(u)])
