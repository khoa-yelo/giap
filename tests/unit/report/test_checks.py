"""``giap.report.checks``: the link checker and the offline (no external load) detector on crafted
HTML (M9-13 work item 2 and acceptance criterion 2; GIAP_SPEC.md §10.3 tests "A link checker finds
no broken relative links. A test asserts that no external URLs are loaded at view time").

Interface pinned here (M9-13 test author, beyond the spec text):

- ``relative_links(html: str) -> list[str]``: every relative URL reference of the page, from any
  element's ``href``, ``src`` or ``xlink:href`` attribute (``<a>`` included: a link to
  ``../tables/islands.tsv`` must resolve) and from CSS ``url(...)`` in ``<style>`` blocks and
  ``style=""`` attributes. Each value is returned with surrounding whitespace and quotes removed,
  its ``#fragment`` and ``?query`` removed and percent-escapes decoded, so a caller can resolve it
  against the report directory directly. Not relative, so never returned: fragment-only
  references (``#V01``, ``url(#clip1)``), ``data:`` URIs, any URL with a scheme (``https:``,
  ``file:``, ``mailto:``, ``javascript:``) and protocol-relative ``//host/...`` URLs; empty values.
- ``external_loads(html: str) -> list[str]``: every URL the browser would fetch from outside the
  report when the page is opened: an absolute (any scheme but ``data:``), protocol-relative
  (``//cdn...``) or ``file:`` URL in any element's ``src``; in ``href``/``xlink:href`` of
  ``<link>``, ``<script>``, ``<img>``, ``<iframe>`` and the SVG ``<image>``/``<use>`` elements; in
  CSS ``url(...)`` (``<style>`` blocks and ``style`` attributes) or CSS ``@import`` (both the
  string and the ``url(...)`` forms). Values are returned as written, stripped of surrounding
  whitespace and quotes. Not loads: ``<a href>`` (a navigation the reader chooses, M9-16), plain
  text that mentions a URL, HTML comments, and the text of inline ``<script>`` bodies (inlined
  Plotly.js, M9-14, holds URL strings that are never fetched at view time).
- Both functions return each distinct value once, in first-occurrence document order, and parse
  with the standard library only; tag and attribute names and schemes are case-insensitive.
"""

from __future__ import annotations

from types import ModuleType

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st


def _checks() -> ModuleType:
    from giap.report import checks

    return checks


def _page(body: str, head: str = "") -> str:
    return (
        '<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8"><title>t</title>'
        f"{head}</head><body>{body}</body></html>\n"
    )


# ---------------------------------------------------------------- external_loads


def test_clean_page_has_no_external_loads() -> None:
    html = _page('<section id="summary"><p>12 of 40 genes detected</p></section>')
    assert _checks().external_loads(html) == []


@pytest.mark.parametrize(
    "snippet, url",
    [
        ('<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>',
         "https://cdn.plot.ly/plotly-2.35.2.min.js"),
        ('<img src="http://example.org/logo.png" alt="logo">', "http://example.org/logo.png"),
        ('<iframe src="https://example.org/embed"></iframe>', "https://example.org/embed"),
        ('<link rel="stylesheet" href="https://fonts.googleapis.com/css?family=Inter">',
         "https://fonts.googleapis.com/css?family=Inter"),
        ('<link rel="icon" href="http://example.org/favicon.ico">',
         "http://example.org/favicon.ico"),
        ('<video src="https://example.org/a.mp4"></video>', "https://example.org/a.mp4"),
    ],
)  # fmt: skip
def test_absolute_src_and_link_href_are_loads(snippet: str, url: str) -> None:
    assert _checks().external_loads(_page(snippet)) == [url]


@pytest.mark.parametrize(
    "snippet, url",
    [
        ('<script src="//cdn.jsdelivr.net/npm/x.js"></script>', "//cdn.jsdelivr.net/npm/x.js"),
        ('<img src="//cdn.example.org/a.png">', "//cdn.example.org/a.png"),
        ('<link rel="stylesheet" href="//cdn.example.org/s.css">', "//cdn.example.org/s.css"),
    ],
)
def test_protocol_relative_urls_are_loads(snippet: str, url: str) -> None:
    assert _checks().external_loads(_page(snippet)) == [url]


@pytest.mark.parametrize(
    "snippet, url",
    [
        ('<img src="file:///home/u/OUT/plots/V01_circos.png">',
         "file:///home/u/OUT/plots/V01_circos.png"),
        ('<link rel="stylesheet" href="file:///tmp/style.css">', "file:///tmp/style.css"),
        ('<script src="FILE:///tmp/x.js"></script>', "FILE:///tmp/x.js"),
    ],
)  # fmt: skip
def test_file_urls_are_loads(snippet: str, url: str) -> None:
    assert _checks().external_loads(_page(snippet)) == [url]


def test_css_url_in_style_block_is_a_load() -> None:
    head = "<style>body { background: url(https://example.org/bg.png) no-repeat; }</style>"
    assert _checks().external_loads(_page("<p>x</p>", head)) == ["https://example.org/bg.png"]


@pytest.mark.parametrize("quote", ["", "'", '"'])
def test_css_url_quoting_is_stripped(quote: str) -> None:
    head = f"<style>.a {{ background: url( {quote}http://example.org/a.png{quote} ); }}</style>"
    assert _checks().external_loads(_page("<p>x</p>", head)) == ["http://example.org/a.png"]


def test_css_url_in_style_attribute_is_a_load() -> None:
    body = "<div style=\"background-image: url('//cdn.example.org/tile.png')\">x</div>"
    assert _checks().external_loads(_page(body)) == ["//cdn.example.org/tile.png"]


@pytest.mark.parametrize(
    "rule, url",
    [
        ('@import "https://fonts.googleapis.com/css2?family=Inter";',
         "https://fonts.googleapis.com/css2?family=Inter"),
        ("@import url(https://example.org/theme.css);", "https://example.org/theme.css"),
        ("@import url('//cdn.example.org/t.css') screen;", "//cdn.example.org/t.css"),
        ("@IMPORT 'file:///tmp/t.css';", "file:///tmp/t.css"),
    ],
)  # fmt: skip
def test_css_import_is_a_load(rule: str, url: str) -> None:
    head = f"<style>{rule}\nbody {{ margin: 0; }}</style>"
    assert _checks().external_loads(_page("<p>x</p>", head)) == [url]


def test_svg_image_and_use_with_external_href_are_loads() -> None:
    body = (
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">'
        '<image xlink:href="https://example.org/tile.png" width="5" height="5"/>'
        '<use href="https://example.org/sprites.svg#dot"/>'
        "</svg>"
    )
    assert _checks().external_loads(_page(body)) == [
        "https://example.org/tile.png",
        "https://example.org/sprites.svg#dot",
    ]


def test_svg_namespace_declarations_are_not_loads() -> None:
    body = (
        '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
        'version="1.1"><g id="V01-g1"><path d="M0 0L1 1"/></g></svg>'
    )
    assert _checks().external_loads(_page(body)) == []


@pytest.mark.parametrize(
    "snippet",
    [
        '<img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAAB">',
        '<link rel="icon" href="data:image/svg+xml,%3Csvg%3E%3C/svg%3E">',
        "<style>.a { background: url(data:image/gif;base64,R0lGODlhAQABAAAAACw=); }</style>",
        '<img src="../plots/V01_circos.svg">',
        '<link rel="stylesheet" href="assets/report.css">',
        '<script src="./assets/table.js"></script>',
        "<style>.a { background: url(../plots/V02_abundance_landscape.png); }</style>",
        "<style>@import 'assets/print.css';</style>",
        '<svg><use xlink:href="#V18-marker"/><path clip-path="url(#V18-clip)"/></svg>',
        '<iframe src="../plots/V03_locus_ISL0001.svg"></iframe>',
    ],
)
def test_data_uris_relative_paths_and_fragments_are_not_loads(snippet: str) -> None:
    assert _checks().external_loads(_page(snippet)) == []


@pytest.mark.parametrize(
    "snippet",
    [
        '<a href="https://github.com/khoa/giap/blob/main/docs/methods/detection.md">methods</a>',
        '<a href="//example.org/x">x</a>',
        '<a href="file:///tmp/x.tsv">x</a>',
        "<p>Data from https://www.uniprot.org/ retrieved 2026-10-01.</p>",
        '<!-- <img src="https://example.org/tracker.gif"> -->',
        "<script>var u = 'https://cdn.plot.ly/x.js'; /* never fetched */</script>",
    ],
)
def test_anchors_text_comments_and_inline_script_bodies_are_not_loads(snippet: str) -> None:
    assert _checks().external_loads(_page(snippet)) == []


def test_case_insensitive_tags_attributes_and_schemes() -> None:
    html = _page('<IMG SRC="HTTPS://Example.org/A.png"><SCRIPT Src="Http://x.org/b.js"></SCRIPT>')
    assert _checks().external_loads(html) == ["HTTPS://Example.org/A.png", "Http://x.org/b.js"]


def test_loads_in_document_order_each_once() -> None:
    head = (
        '<link rel="stylesheet" href="https://a.example/s.css">'
        "<style>@import url(https://b.example/t.css);</style>"
    )
    body = (
        '<img src="//c.example/1.png"><img src="https://a.example/s.css">'
        '<img src="//c.example/1.png"><img src="../plots/V05_mapping_sankey.svg">'
    )
    assert _checks().external_loads(_page(body, head)) == [
        "https://a.example/s.css",
        "https://b.example/t.css",
        "//c.example/1.png",
    ]


def test_whitespace_around_url_is_stripped() -> None:
    html = _page('<img src="  https://example.org/a.png  ">')
    assert _checks().external_loads(html) == ["https://example.org/a.png"]


# ---------------------------------------------------------------- relative_links


def test_relative_links_from_anchors_images_links_and_scripts() -> None:
    head = '<link rel="stylesheet" href="assets/report.css">'
    body = (
        '<a href="../tables/islands.tsv">islands</a>'
        '<img src="../plots/V01_circos.svg" alt="V01">'
        '<a href="../plots/data/V01_circos.tsv">data</a>'
        '<script src="./assets/table.js"></script>'
    )
    assert _checks().relative_links(_page(body, head)) == [
        "assets/report.css",
        "../tables/islands.tsv",
        "../plots/V01_circos.svg",
        "../plots/data/V01_circos.tsv",
        "./assets/table.js",
    ]


def test_relative_links_exclude_absolute_protocol_relative_data_and_fragments() -> None:
    body = (
        '<a href="#detection">jump</a>'
        '<a href="https://example.org/doc">doc</a>'
        '<img src="//cdn.example.org/a.png">'
        '<img src="data:image/png;base64,AAAA">'
        '<a href="file:///tmp/x.tsv">x</a>'
        '<a href="mailto:someone@example.org">mail</a>'
        '<a href="javascript:void(0)">js</a>'
        '<a href="">empty</a>'
        '<svg><use xlink:href="#m1"/><path clip-path="url(#c1)"/></svg>'
        '<a href="../tables/regions.tsv">regions</a>'
    )
    assert _checks().relative_links(_page(body)) == ["../tables/regions.tsv"]


def test_relative_links_strip_fragment_and_query_and_decode() -> None:
    body = (
        '<a href="../plots/V03_locus_ISL0001.svg#view">a</a>'
        '<a href="../tables/genes.tsv?v=2">b</a>'
        '<a href="../plots/data/V13%20novel.tsv">c</a>'
    )
    assert _checks().relative_links(_page(body)) == [
        "../plots/V03_locus_ISL0001.svg",
        "../tables/genes.tsv",
        "../plots/data/V13 novel.tsv",
    ]


def test_relative_links_include_css_urls() -> None:
    head = "<style>@import url(assets/print.css); .a { background: url('img/bg.png'); }</style>"
    body = '<div style="background: url(../plots/V02_abundance_landscape.png)">x</div>'
    assert _checks().relative_links(_page(body, head)) == [
        "assets/print.css",
        "img/bg.png",
        "../plots/V02_abundance_landscape.png",
    ]


def test_relative_links_each_once_in_document_order() -> None:
    body = (
        '<a href="../tables/islands.tsv">1</a><img src="../plots/V01_circos.svg">'
        '<a href="../tables/islands.tsv#row3">2</a>'
    )
    assert _checks().relative_links(_page(body)) == [
        "../tables/islands.tsv",
        "../plots/V01_circos.svg",
    ]


def test_relative_links_of_a_page_without_links_is_empty() -> None:
    assert _checks().relative_links(_page("<p>no links</p>")) == []


def test_relative_links_ignore_commented_out_markup() -> None:
    body = '<!-- <a href="../tables/old.tsv">old</a> --><a href="../tables/tus.tsv">tus</a>'
    assert _checks().relative_links(_page(body)) == ["../tables/tus.tsv"]


# ---------------------------------------------------------------- properties

_SEG = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-", min_size=1, max_size=8)
_REL = st.builds(
    lambda up, parts, ext: "../" * up + "/".join(parts) + ext,
    st.integers(0, 2),
    st.lists(_SEG, min_size=1, max_size=3),
    st.sampled_from([".svg", ".png", ".tsv", ".css", ".js"]),
)
_EXT = st.builds(
    lambda pre, host, path: f"{pre}{host}.org/{path}",
    st.sampled_from(["https://", "http://", "//", "HTTPS://"]),
    _SEG,
    _SEG,
)
_FILE = st.builds(lambda p: f"file:///tmp/{p}", _SEG)
_DATA = st.just("data:image/png;base64,iVBORw0KGgo=")
_FRAG = st.builds(lambda s: f"#{s}", _SEG)
_TAG = st.sampled_from(
    [
        ('<img src="{u}">', True),
        ('<script src="{u}"></script>', True),
        ('<link rel="stylesheet" href="{u}">', True),
        ('<iframe src="{u}"></iframe>', True),
        ("<style>.k {{ background: url({u}); }}</style>", True),
        ('<a href="{u}">x</a>', False),
    ]
)
_KIND = st.sampled_from(["rel", "ext", "file", "data", "frag"])


@st.composite
def _documents(draw: st.DrawFn) -> tuple[str, list[str], list[str]]:
    """A page of random references; returns (html, expected relative, expected external)."""
    n = draw(st.integers(0, 12))
    parts: list[str] = []
    rel: list[str] = []
    ext: list[str] = []
    for _ in range(n):
        (tmpl, loads), kind = draw(_TAG), draw(_KIND)
        if kind == "frag" and tmpl.startswith("<style>"):
            kind = "data"  # url(#id) is covered above; keep the generator's CSS simple
        url = draw({"rel": _REL, "ext": _EXT, "file": _FILE, "data": _DATA, "frag": _FRAG}[kind])
        parts.append(tmpl.format(u=url))
        if kind == "rel" and url not in rel:
            rel.append(url)
        if kind in ("ext", "file") and loads and url not in ext:
            ext.append(url)
    return _page("".join(parts)), rel, ext


@settings(max_examples=150, deadline=None)
@given(_documents())
def test_property_detectors_match_the_generated_references(
    doc: tuple[str, list[str], list[str]],
) -> None:
    html, rel, ext = doc
    checks = _checks()
    assert checks.relative_links(html) == rel
    assert checks.external_loads(html) == ext


@settings(max_examples=100, deadline=None)
@given(_documents())
def test_property_relative_links_and_external_loads_are_disjoint(
    doc: tuple[str, list[str], list[str]],
) -> None:
    html, _, _ = doc
    checks = _checks()
    assert not set(checks.relative_links(html)) & set(checks.external_loads(html))


@settings(max_examples=100, deadline=None)
@given(_documents())
def test_property_pure_and_deterministic(doc: tuple[str, list[str], list[str]]) -> None:
    html, _, _ = doc
    checks = _checks()
    assert checks.relative_links(html) == checks.relative_links(html)
    assert checks.external_loads(html) == checks.external_loads(html)
