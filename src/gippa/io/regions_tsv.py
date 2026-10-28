"""Writer for ``tables/regions_stats.tsv`` (GIPPA_SPEC.md §9.3; docs/specs/M4-18).

Columns ``gippa.context.regions_table.REGIONS_STATS_COLUMNS``; ``start`` converted to 1-based
(``end`` unchanged: half-open 0-based end = closed 1-based end; origin-crossing regions keep
``end > length``, §6.5). Floats with round-trip precision (as ``detection_summary.tsv``), ``NA``
for missing, ``flags`` ``;``-joined (``NA`` when none).
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from gippa.context.regions_table import REGIONS_STATS_COLUMNS
from gippa.io.tsv import write_tsv

__all__ = ["write_regions_stats"]

_INTS = frozenset({"n_genes", "n_eligible", "n_perm", "end"})
_FLOATS = frozenset({"observed", "expected", "o_minus_e", "o_over_e", "hl_shift_log2", "effect",
                     "ci_low", "ci_high", "p_value", "q_value"})  # fmt: skip


def _num(v: object) -> object:
    if v is None or pd.isna(v):
        return None
    f = float(v)  # type: ignore[arg-type]
    if math.isinf(f):
        return f
    return int(f) if f.is_integer() and abs(f) < 2**53 else repr(f)


def _flags(v: object) -> str | None:
    if isinstance(v, str):
        return v or None
    if isinstance(v, list | tuple):
        return ";".join(map(str, v)) or None
    return None


def _cell(col: str, v: object) -> object:
    if col == "flags":
        return _flags(v)
    if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)):
        return None
    if col == "start":
        return int(v) + 1  # type: ignore[call-overload]
    if col in _INTS:
        return int(v)  # type: ignore[call-overload]
    if col in _FLOATS:
        return _num(v)
    return str(v)


def write_regions_stats(path: Path, df: pd.DataFrame) -> None:
    """Write an assembled ``regions_stats`` frame
    (``gippa.context.regions_table.regions_stats``)."""
    cols = list(REGIONS_STATS_COLUMNS)
    rows = [[_cell(c, v) for c, v in zip(cols, r, strict=True)]
            for r in df[cols].itertuples(index=False, name=None)]  # fmt: skip
    write_tsv(path, cols, rows)
