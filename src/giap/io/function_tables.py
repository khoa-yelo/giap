"""Writers of ``tables/category_summary.tsv``, ``enrichment_detection.tsv`` and
``enrichment_gsea.tsv`` (GIAP_SPEC.md §9.2, §9.3; ADR-0006; docs/specs/M8-06 item 4).

Cells follow ``giap.io.tsv`` (``NA``, ``%.6g``, ``true``/``false``, ``;``-joined ``leading_edge``
and ``flags``; an empty list is ``NA``). Columns are each module's ``COLUMNS``, verbatim; extra
frame columns are dropped. The summary and ORA frames keep the order they are given (already
collection, term, direction); GSEA rows are sorted stably by (collection in ``COLLECTIONS``
order, term). The tables carry no coordinates. The inputs are not modified.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.function.collections import COLLECTIONS
from giap.function.gsea import COLUMNS as GSEA_COLUMNS
from giap.function.ora import COLUMNS as ORA_COLUMNS
from giap.function.summary import COLUMNS as SUMMARY_COLUMNS
from giap.io.tsv import write_tsv

__all__ = ["write_category_summary_tsv", "write_enrichment_detection_tsv",
           "write_enrichment_gsea_tsv"]  # fmt: skip

_ORDER = {c: i for i, c in enumerate(COLLECTIONS)}


def _cell(v: object) -> object:
    if isinstance(v, np.ndarray | list | tuple):
        return [str(x) for x in v]
    return v


def _write(path: Path, df: pd.DataFrame, cols: Sequence[str]) -> None:
    rows = [[_cell(v) for v in r] for r in df[list(cols)].itertuples(index=False, name=None)]
    write_tsv(path, cols, rows)


def write_category_summary_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/category_summary.tsv`` from ``giap.function.summary.category_summary``."""
    _write(path, df, SUMMARY_COLUMNS)


def write_enrichment_detection_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/enrichment_detection.tsv`` from ``giap.function.ora.detection_enrichment``."""
    _write(path, df, ORA_COLUMNS)


def write_enrichment_gsea_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/enrichment_gsea.tsv`` from ``giap.function.gsea.gsea_table``, sorted stably by
    (collection in ``COLLECTIONS`` order, term)."""
    if len(df):
        keys = list(zip(df["collection"].map(_ORDER), df["term"].astype(str), strict=True))
        order = sorted(range(len(df)), key=keys.__getitem__)  # sorted() is stable
        df = df.iloc[order]
    _write(path, df, GSEA_COLUMNS)
