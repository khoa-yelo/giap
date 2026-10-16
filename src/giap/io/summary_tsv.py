"""Writer for ``detection_summary.tsv`` (GIAP_SPEC.md §9.3; docs/specs/M4-16).

Long format with the columns of ``giap.context.summary.SUMMARY_COLUMNS``. Floats are written with
round-trip precision (shortest ``repr``; integral values as integers) so values and CIs read back
exactly; missing values are ``NA``; ``flags`` is ``;``-joined (empty when none).
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from giap.context.summary import SUMMARY_COLUMNS
from giap.io.tsv import write_tsv

__all__ = ["write_detection_summary"]


def _num(v: object) -> object:
    if v is None or pd.isna(v):
        return None
    f = float(v)  # type: ignore[arg-type]
    if math.isinf(f):
        return f
    return int(f) if f.is_integer() and abs(f) < 2**53 else repr(f)


def write_detection_summary(path: Path, df: pd.DataFrame) -> None:
    """Write the assembled detection summary (``detection_summary`` output) to ``path``."""
    rows = []
    for r in df[list(SUMMARY_COLUMNS)].itertuples(index=False):
        n = None if pd.isna(r.n) else int(r.n)
        flags = "" if r.flags is None or pd.isna(r.flags) else str(r.flags)
        rows.append([str(r.section), str(r.group_type), str(r.group), str(r.metric),
                     _num(r.value), _num(r.ci_low), _num(r.ci_high), n, flags])  # fmt: skip
    write_tsv(path, SUMMARY_COLUMNS, rows)
