"""#221 review round: INF cells sort numerically (ADR-0006), the JS and Python number patterns
are one pattern, ASCII digits only, the size budget is decimal 25 MB, the twin id a JS literal."""

from __future__ import annotations

import re
from pathlib import Path

from giap.report import tables

SRC = Path(__file__).resolve().parents[3] / "src" / "giap"
JS = SRC / "report" / "templates" / "tables.js"


def test_inf_cells_keep_a_column_numeric() -> None:
    cells = ["10", "INF", "2", "NA", "-INF", "1e-05"]
    assert tables.numeric_column(cells)
    assert tables.sort_order(cells) == [4, 5, 2, 0, 1, 3]
    assert tables.sort_order(cells, descending=True) == [1, 0, 2, 5, 4, 3]


def test_unicode_digits_are_not_numbers() -> None:
    assert not tables.is_number("٣")  # ARABIC-INDIC DIGIT THREE: JS \d does not match it
    assert not tables.is_number("inf") and not tables.is_number("Infinity")


def test_js_number_regex_is_the_python_pattern() -> None:
    m = re.search(r"var NUM = /(.*)/;", JS.read_text())
    assert m is not None
    assert m.group(1) == "^" + tables._NUM.pattern + "$"


def test_js_maps_inf_to_infinity() -> None:
    js = JS.read_text()
    assert "Infinity" in js and "-Infinity" in js
    assert "num(cells[a])" in js and "num(cells[b])" in js


def test_report_budget_is_decimal_25_mb() -> None:
    from giap.steps.plots import REPORT_BUDGET

    assert REPORT_BUDGET == 25_000_000


def test_twin_id_is_a_json_literal_in_the_script() -> None:
    block = (SRC / "report" / "templates" / "plotly_block.html.j2").read_text()
    assert "{{ i.twin.id | tojson }}" in block
    assert 'var id = "{{' not in block
