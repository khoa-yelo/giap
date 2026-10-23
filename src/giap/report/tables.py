"""The sort and filter rules of the report tables (M9-14 behaviour 4): the pure reference that
``templates/tables.js`` follows line for line.

A column sorts numerically when every non-``NA`` cell is a plain decimal number (the ``_NUM``
pattern, shared with the JS), else by string (code-point order); ties keep their order and ``NA``
is last in both directions. A row matches a filter when the query is a case-insensitive substring
of any one cell; an empty query matches every row.
"""

from __future__ import annotations

import re
from collections.abc import Sequence

NA = "NA"
_NUM = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")


def is_number(cell: str) -> bool:
    return _NUM.fullmatch(cell) is not None


def numeric_column(cells: Sequence[str]) -> bool:
    return all(is_number(c) for c in cells if c != NA)


def sort_order(cells: Sequence[str], descending: bool = False) -> list[int]:
    """Row indices in display order (stable; ``NA`` last)."""
    keep = [i for i, c in enumerate(cells) if c != NA]
    na = [i for i, c in enumerate(cells) if c == NA]
    if numeric_column(cells):
        keep.sort(key=lambda i: float(cells[i]), reverse=descending)
    else:
        keep.sort(key=lambda i: cells[i], reverse=descending)
    return keep + na


def row_matches(cells: Sequence[str], query: str) -> bool:
    q = query.lower()
    return q == "" or any(q in c.lower() for c in cells)
