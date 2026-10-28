"""M9-25 work item 3: the structure snapshot of ``report.html`` and ``compare.html`` on both
frozen fixtures (docs/specs/M9-25-snapshot-baselines.md example 3; GIPPA_SPEC.md §10.3 "A snapshot
of the HTML structure (headings and IDs)").

The ordered ``[tag, id]`` list (``helpers.snapshot_fixture.structure``: every ``<section>``,
``<h2>``, ``<h3>``, plot ``<div>``, id'd ``<div>``, outermost ``<svg>`` and figure stem) of each
page built on the fixture equals ``tests/snapshot/report_structure.json``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from helpers import snapshot_fixture as sf

KEYS = sorted(f"{f}/{p}" for f in sf.FIXTURES for p in ("report.html", "compare.html"))


def _stored() -> dict[str, list[list[str | None]]]:
    if not sf.STRUCTURE_JSON.is_file():
        pytest.fail("tests/snapshot/report_structure.json is not committed")
    data: dict[str, list[list[str | None]]] = json.loads(
        sf.STRUCTURE_JSON.read_text(encoding="utf-8")
    )
    return data


def _page(key: str, tmp: Path) -> str:
    fixture, page = key.split("/")
    why = sf.missing(fixture)
    if why is not None:
        pytest.fail(why)
    build = sf.build_report if page == "report.html" else sf.build_compare
    return build(fixture, tmp).read_text(encoding="utf-8")


def test_structure_json_has_the_four_pages() -> None:
    assert sorted(_stored()) == KEYS


@pytest.mark.parametrize("key", KEYS)
def test_page_structure_matches_the_snapshot(key: str, tmp_path: Path) -> None:
    stored = _stored()
    assert key in stored, f"{key} missing from report_structure.json"
    assert sf.structure(_page(key, tmp_path)) == stored[key]


def test_renamed_section_id_fails_the_snapshot(tmp_path: Path) -> None:
    """Example 3: rename a report section id -> the structure differs."""
    key = f"{sf.RICH}/report.html"
    stored = _stored()[key]
    html = _page(key, tmp_path)
    assert sf.structure(html) == stored
    renamed = html.replace('<section id="detection">', '<section id="detection-renamed">', 1)
    assert renamed != html
    assert sf.structure(renamed) != stored


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_stored_report_sections_are_the_13_fixed_ids(fixture: str) -> None:
    """Index decision 21: the report's sections have fixed IDs, in order."""
    from gippa.report.build import SECTIONS

    stored = _stored()[f"{fixture}/report.html"]
    assert [i for t, i in stored if t == "section"] == [s.id for s in SECTIONS]
    assert len(SECTIONS) == 13


@pytest.mark.parametrize("fixture", sf.FIXTURES)
def test_stored_compare_sections(fixture: str) -> None:
    stored = _stored()[f"{fixture}/compare.html"]
    sections = [i for t, i in stored if t == "section"]
    assert sections == ["runs", "detection", "abundance", "stability", "provenance"]
    assert ["div", "plot:V16"] in stored


def test_rich_report_inlines_a_figure_for_every_plot_but_v16() -> None:
    """On the rich fixture every report plot slot holds at least one inlined figure."""
    stored = _stored()[f"{sf.RICH}/report.html"]
    slots: dict[str, int] = {}
    current = None
    for tag, ident in stored:
        if tag == "div" and str(ident).startswith("plot:"):
            current = str(ident).removeprefix("plot:")
            slots[current] = 0
        elif tag == "figure" and current is not None:
            slots[current] += 1
        elif tag == "section":
            current = None
    from gippa.report.build import SECTIONS

    assert sorted(slots) == sorted(p for s in SECTIONS for p in s.plots)
    assert "V16" not in slots
    assert {p: n for p, n in slots.items() if n == 0} == {}


# ---------------------------------------------------------------- the extractor itself

_HTML = """<html><body><h1>t</h1>
<section id="a"><h2>Alpha  one</h2>
<div class="plot" data-plot="V01"><figure><div id="plotly-V01" class="twin"></div>
<svg id="x"><g id="inner"><svg id="nested"></svg></g><text>no</text></svg>
<figcaption class="links">V01_circos <a href="#">SVG</a></figcaption></figure></div>
<h3 id="h">ignored text</h3><div class="table-wrap"></div></section></body></html>"""


def test_structure_extractor_on_a_small_page() -> None:
    assert sf.structure(_HTML) == [
        ["section", "a"],
        ["h2", "Alpha one"],
        ["div", "plot:V01"],
        ["div", "plotly-V01"],
        ["svg", "x"],
        ["figure", "V01_circos"],
        ["h3", "h"],
    ]
