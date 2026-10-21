"""Writer for ``detection_summary.tsv`` and its Parquet parts (GIAP_SPEC.md §9.3; M4-16).

Long format with the columns of ``giap.context.summary.SUMMARY_COLUMNS``. Floats are written with
round-trip precision (shortest ``repr``; integral values as integers) so values and CIs read back
exactly; missing values are ``NA``; ``flags`` is ``;``-joined (empty when none).

``write_summary_part`` writes one section's rows (the ``tier`` rows of M4-12, M4-20 lead
decision 3) as a byte-reproducible Parquet file that ``context_stats`` reads back and hands to
``giap.context.summary.detection_summary``.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.context.summary import SUMMARY_COLUMNS
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import write_tsv

__all__ = ["SUMMARY_PART_TYPES", "write_detection_summary", "write_summary_part"]

#: Arrow type per ``SUMMARY_COLUMNS`` column of a ``detection_summary`` part file.
SUMMARY_PART_TYPES: tuple[tuple[str, pa.DataType], ...] = tuple(
    zip(SUMMARY_COLUMNS, (pa.string(),) * 4 + (pa.float64(),) * 3 + (pa.int64(), pa.string()),
        strict=True)
)  # fmt: skip


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


def write_summary_part(path: Path, df: pd.DataFrame) -> None:
    """Write the ``detection_summary`` rows ``df`` (``SUMMARY_COLUMNS``) to ``path`` as Parquet."""
    rows = []
    for r in df[list(SUMMARY_COLUMNS)].itertuples(index=False):
        flags = "" if r.flags is None or pd.isna(r.flags) else str(r.flags)
        rows.append([
            str(r.section), str(r.group_type), str(r.group), str(r.metric),
            *(None if pd.isna(v) else float(v) for v in (r.value, r.ci_low, r.ci_high)),
            None if pd.isna(r.n) else int(r.n), flags,
        ])  # fmt: skip
    write_parquet_deterministic(path, SUMMARY_PART_TYPES, rows)
