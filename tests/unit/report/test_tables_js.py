"""Sortable and filterable report tables (M9-14 behaviour 4; examples 4 and 5; degenerate row
"a table with 0 rows"; GIPPA_SPEC.md §10.3 "Every table section shows the top rows, sortable and
filterable"; index decision 24).

Lead decision 1: no headless browser, so the rules are tested through a pure-Python mirror of the
JS rules plus a structural check of the inlined ``tables.js``.

Interface assumed beyond the spec text (M9-14 test author):
- the mirror is ``gippa.report.tables`` (pure, no I/O), the reference the JS follows:
  ``NA = "NA"``; ``is_number(cell) -> bool``; ``numeric_column(cells) -> bool`` (every
  non-``NA`` cell is a number); ``sort_order(cells, descending=False) -> list[int]`` (row indices
  in display order: numeric sort when ``numeric_column``, else string sort; ties keep their order;
  ``NA`` always last in both directions); ``row_matches(cells, query) -> bool``
  (case-insensitive substring over the row; an empty query matches every row);
- ``src/gippa/report/templates/tables.js`` is inlined verbatim (stripped) in one ``<script>``
  without ``src``, under both ``report.interactive`` settings (it is not Plotly);
- each ``<table data-table=K>`` has exactly one ``<input data-filter=K>`` filter box.

Only cases where the JS and the mirror cannot differ are pinned: numbers are plain decimals,
integers or ``1e-05`` style; strings are distinct lowercase ASCII (no locale question); filter
queries lie inside one cell.
"""

from __future__ import annotations

import itertools
import re
from importlib import resources
from pathlib import Path
from types import ModuleType

import pytest
from helpers.interactive_fixtures import page_of, register_twin, run_steps, swap_twins
from helpers.report_fixtures import register_fake, run_out, swap_registry, table_rows
from helpers.viz_out_tree import write_schema_table
from hypothesis import given
from hypothesis import strategies as st


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _mirror() -> ModuleType:
    import gippa.report.tables as tables

    return tables


def _sorted(cells: list[str], descending: bool = False) -> list[str]:
    return [cells[i] for i in _mirror().sort_order(cells, descending=descending)]


# ---------------------------------------------------------------- example 4: sorting


def test_example_4_numeric_ascending_na_last() -> None:
    assert _sorted(["2", "10", "NA", "1.5"]) == ["1.5", "2", "10", "NA"]


def test_numeric_descending_keeps_na_last() -> None:
    assert _sorted(["2", "10", "NA", "1.5"], descending=True) == ["10", "2", "1.5", "NA"]


def test_numeric_column_detection() -> None:
    t = _mirror()
    assert t.NA == "NA"
    assert t.numeric_column(["2", "10", "NA", "1.5", "-3", "1e-05", "0"])
    assert not t.numeric_column(["2", "x", "10"])
    assert not t.numeric_column(["1.2.3", "4"])
    assert not t.numeric_column(["12 kb", "4"])
    for good in ("2", "-3", "0.5", "1e-05", "2.5E+03"):
        assert t.is_number(good), good
    for bad in ("NA", "abc", "1.2.3", "12 kb", "chr1"):
        assert not t.is_number(bad), bad


def test_mixed_column_sorts_as_strings() -> None:
    assert _sorted(["2", "x", "10"]) == ["10", "2", "x"]
    assert _sorted(["2", "NA", "x", "10"]) == ["10", "2", "x", "NA"]


def test_string_sort_both_directions_na_last() -> None:
    assert _sorted(["b", "NA", "a", "c"]) == ["a", "b", "c", "NA"]
    assert _sorted(["b", "NA", "a", "c"], descending=True) == ["c", "b", "a", "NA"]


def test_numeric_is_by_value_not_text() -> None:
    assert _sorted(["1e-05", "0.5", "-3", "10"]) == ["-3", "1e-05", "0.5", "10"]


def test_degenerate_columns() -> None:
    t = _mirror()
    assert t.sort_order([]) == []
    assert _sorted(["NA", "NA"]) == ["NA", "NA"]
    assert t.sort_order(["NA", "NA"]) == [0, 1]
    assert _sorted(["7"]) == ["7"]


def test_ties_keep_their_order() -> None:
    t = _mirror()
    assert t.sort_order(["2", "1", "2", "1"]) == [1, 3, 0, 2]
    assert t.sort_order(["b", "a", "b"], descending=True) == [0, 2, 1]


_NUM = st.one_of(
    st.integers(-10_000, 10_000).map(str),
    st.floats(-1e6, 1e6, allow_nan=False, allow_infinity=False).map(lambda f: f"{f:.6g}"),
)
_COLUMN = st.lists(st.one_of(_NUM, st.just("NA")), max_size=30)
_WORDS = st.lists(st.one_of(st.text("abcdefgh", min_size=1, max_size=5), st.just("NA")),
                  max_size=30)  # fmt: skip


@given(cells=_COLUMN, descending=st.booleans())
def test_property_numeric_sort(cells: list[str], descending: bool) -> None:
    order = _mirror().sort_order(cells, descending=descending)
    assert sorted(order) == list(range(len(cells)))
    shown = [cells[i] for i in order]
    k = sum(c != "NA" for c in cells)
    assert all(c == "NA" for c in shown[k:]) and all(c != "NA" for c in shown[:k])
    values = [float(c) for c in shown[:k]]
    assert values == sorted(values, reverse=descending)


@given(cells=_WORDS, descending=st.booleans())
def test_property_string_sort_stable(cells: list[str], descending: bool) -> None:
    t = _mirror()
    order = t.sort_order(cells, descending=descending)
    assert sorted(order) == list(range(len(cells)))
    shown = [cells[i] for i in order]
    k = sum(c != "NA" for c in cells)
    assert all(c == "NA" for c in shown[k:])
    if not t.numeric_column(cells):
        assert shown[:k] == sorted(shown[:k], reverse=descending)
    for a, b in itertools.pairwise(order):
        if cells[a] == cells[b]:
            assert a < b, "ties must keep their order"


# ---------------------------------------------------------------- example 5: filtering

REGIONS = [
    ["R1", "chr1", "prophage", "12"],
    ["R2", "chr1", "genomic_island", "3"],
    ["R3", "chr2", "Prophage_remnant", "NA"],
    ["R4", "chr2", "integron", "1"],
]


def test_example_5_filter_prophage() -> None:
    t = _mirror()
    visible = [r[0] for r in REGIONS if t.row_matches(r, "prophage")]
    assert visible == ["R1", "R3"]


def test_filter_is_case_insensitive_and_any_cell() -> None:
    t = _mirror()
    assert [r[0] for r in REGIONS if t.row_matches(r, "PROPHAGE")] == ["R1", "R3"]
    assert [r[0] for r in REGIONS if t.row_matches(r, "chr2")] == ["R3", "R4"]
    assert [r[0] for r in REGIONS if t.row_matches(r, "island")] == ["R2"]
    assert [r[0] for r in REGIONS if t.row_matches(r, "nothing here")] == []


def test_empty_filter_shows_every_row() -> None:
    t = _mirror()
    assert all(t.row_matches(r, "") for r in REGIONS)
    assert t.row_matches([], "")


@given(row=st.lists(st.text("abcXYZ_ 1", max_size=6), min_size=1, max_size=5),
       data=st.data())  # fmt: skip
def test_property_any_cell_substring_matches(row: list[str], data: st.DataObject) -> None:
    cell = data.draw(st.sampled_from(row))
    i = data.draw(st.integers(0, len(cell)))
    j = data.draw(st.integers(i, len(cell)))
    q = cell[i:j]
    t = _mirror()
    assert t.row_matches(row, q)
    assert t.row_matches(row, q.upper()) and t.row_matches(row, q.lower())


# ---------------------------------------------------------------- tables.js, structurally


def _js() -> str:
    return (resources.files("gippa.report") / "templates" / "tables.js").read_text(encoding="utf-8")


def test_tables_js_is_small_vanilla_and_offline() -> None:
    js = _js()
    lines = [ln for ln in js.splitlines() if ln.strip()]
    assert 5 <= len(lines) <= 80, f"{len(lines)} non-blank lines (about 40 expected)"
    for banned in ("http:", "https:", "fetch(", "XMLHttpRequest", "import ", "require(",
                   "eval(", "document.write", "jQuery", "$("):  # fmt: skip
        assert banned not in js, banned


def test_tables_js_implements_the_documented_rules() -> None:
    js = _js()
    assert "addEventListener" in js
    assert re.search(r"""["']click["']""", js), "no header click handler"
    assert re.search(r"""["']input["']""", js), "no filter input handler"
    assert re.search(r"""["']NA["']""", js), "NA is not special-cased"
    assert "toLowerCase" in js, "the filter is not case-insensitive"
    assert re.search(r"Number\(|parseFloat\(|isFinite\(|isNaN\(", js), "no numeric parse"
    assert "data-table" in js or "dataset.table" in js, "tables not found by data-table"
    assert "data-filter" in js or "dataset.filter" in js, "filter boxes not found"


def _inlined(html_path: Path) -> list[str]:
    page = page_of(html_path)
    js = _js().strip()
    return [s.text() for s in page.root.by_tag("script")
            if not s.attrs.get("src") and js in s.text()]  # fmt: skip


def test_tables_js_is_inlined_once(fake_registry: ModuleType, tmp_path: Path) -> None:
    html = run_steps(tmp_path, out=run_out(tmp_path, islands=3))
    assert len(_inlined(html)) == 1


def test_tables_js_is_inlined_without_interactive(
    fake_registry: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    swap_twins(monkeypatch)
    register_fake("V17")
    register_twin("V17")
    html = run_steps(tmp_path, "report.interactive=false", out=run_out(tmp_path, islands=3))
    assert len(_inlined(html)) == 1
    assert "plotly.js" not in html.read_text(encoding="utf-8").lower()


# ---------------------------------------------------------------- filter boxes


def test_every_embedded_table_has_one_filter_box(fake_registry: ModuleType, tmp_path: Path) -> None:
    out = run_out(tmp_path, islands=25)
    write_schema_table(out / "tables" / "sensitivity.tsv", "sensitivity")
    page = page_of(run_steps(tmp_path, out=out))
    keys = [t.attrs["data-table"] for t in page.root.by_tag("table") if "data-table" in t.attrs]
    assert {"islands", "sensitivity"} <= set(keys)
    for key in keys:
        boxes = [n for n in page.root.by_tag("input") if n.attrs.get("data-filter") == key]
        assert len(boxes) == 1, (key, len(boxes))
        assert boxes[0].attrs.get("type", "text") in ("text", "search"), boxes[0].attrs


def test_zero_row_table_is_header_only_with_a_filter_box(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    write_schema_table(out / "tables" / "sensitivity.tsv", "sensitivity")
    page = page_of(run_steps(tmp_path, out=out))
    table = page.table("sensitivity")
    assert table is not None
    header, body = table_rows(table)
    assert header and body == []
    assert [n for n in page.root.by_tag("input") if n.attrs.get("data-filter") == "sensitivity"]


def test_tables_still_embed_at_most_20_rows(fake_registry: ModuleType, tmp_path: Path) -> None:
    page = page_of(run_steps(tmp_path, out=run_out(tmp_path, islands=25)))
    table = page.table("islands")
    assert table is not None
    assert len(table_rows(table)[1]) == 20
