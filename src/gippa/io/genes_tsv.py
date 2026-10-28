"""``tables/genes.tsv`` writer and reader, the gene-table Parquet files and the COG 2020 table
(docs/specs/M4-15-gene-table-assembly.md; GIPPA_SPEC.md §6.5, §6.6, §9.2, §9.3; ADR-0006).

The frame is the internal representation of ``gippa.genetable.assemble`` (0-based half-open, lists
in list columns). The TSV is 1-based closed (``start + 1``, ``end`` unchanged), lists ``;``-joined,
``NA`` for missing, written through M1-04 ``write_tsv``; ``read_genes_tsv`` is the inverse.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from gippa.genetable.assemble import GENES_COLUMNS, GENES_TYPES, ColumnType
from gippa.io.parquet import write_parquet_deterministic
from gippa.io.tsv import read_tsv, write_tsv

__all__ = ["read_cog_def", "read_genes_parquet", "read_genes_tsv", "write_genes_parquet",
           "write_genes_tsv"]  # fmt: skip

_ARROW: dict[ColumnType, Any] = {
    "str": pa.string(), "int": pa.int64(), "float": pa.float64(), "bool": pa.bool_(),
    "list": pa.list_(pa.string()),
}  # fmt: skip


def _na(v: Any) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v))


def _cell(v: Any, kind: ColumnType) -> Any:
    if kind == "list":
        return list(v) if isinstance(v, list | tuple | np.ndarray) else []
    if v is None or v is pd.NA:
        return None
    if kind == "float":
        return None if _na(v) else float(v)
    if kind == "int":
        return None if isinstance(v, float) and np.isnan(v) else int(v)
    if kind == "bool":
        return bool(v)
    return None if isinstance(v, float) and np.isnan(v) else str(v)


def write_genes_tsv(path: Path, df: pd.DataFrame) -> None:
    """Write ``df`` (columns ``GENES_COLUMNS``, any order) as ``genes.tsv`` (1-based closed)."""
    kinds = [GENES_TYPES[c] for c in GENES_COLUMNS]
    frame = df[list(GENES_COLUMNS)]
    i_start = GENES_COLUMNS.index("start")
    rows = []
    for r in frame.itertuples(index=False, name=None):
        cells = [_cell(v, k) for v, k in zip(r, kinds, strict=True)]
        cells[i_start] = cells[i_start] + 1
        rows.append(cells)
    write_tsv(path, GENES_COLUMNS, rows)


def _column(values: list[str | None], kind: ColumnType) -> pd.Series:
    if kind == "list":
        return pd.Series([[] if v is None else v.split(";") for v in values], dtype=object)
    if kind == "int":
        return pd.Series([None if v is None else int(v) for v in values], dtype="Int64")
    if kind == "float":
        return pd.Series([np.nan if v is None else float(v) for v in values], dtype="float64")
    if kind == "bool":
        return pd.Series([v == "true" for v in values], dtype=bool)
    return pd.Series(values, dtype=object)


def read_genes_tsv(path: Path) -> pd.DataFrame:
    """Read ``genes.tsv`` back into the internal representation (``start - 1``)."""
    header, rows = read_tsv(path)
    if tuple(header) != GENES_COLUMNS:
        raise ValueError(f"{path}: header is not the genes.tsv column list")
    data = {c: _column([r[j] for r in rows], GENES_TYPES[c]) for j, c in enumerate(header)}
    data["start"] = data["start"] - 1
    return pd.DataFrame(data, columns=list(GENES_COLUMNS))


def write_genes_parquet(path: Path, df: pd.DataFrame) -> None:
    """Gene-table columns of ``df`` (``GENES_COLUMNS`` order; absent ones skipped) as Parquet."""
    names = [c for c in GENES_COLUMNS if c in df.columns]
    cols = [(c, _ARROW[GENES_TYPES[c]]) for c in names]
    kinds = [GENES_TYPES[c] for c in names]
    rows = [
        [_cell(v, k) for v, k in zip(r, kinds, strict=True)]
        for r in df[names].itertuples(index=False, name=None)
    ]
    write_parquet_deterministic(path, cols, rows)


def read_genes_parquet(path: Path) -> pd.DataFrame:
    """A gene-table Parquet file in the internal representation (list cells as lists, integers
    ``Int64``)."""
    df = pd.read_parquet(path, engine="pyarrow")
    for c in df.columns:
        kind = GENES_TYPES.get(str(c))
        if kind == "list":
            df[c] = pd.Series([list(v) if v is not None else [] for v in df[c]], dtype=object)
        elif kind == "int":
            df[c] = df[c].astype("Int64")
        elif kind == "str":  # missing stays missing (NaN / pd.NA never become "nan")
            df[c] = pd.Series([None if _na(v) else str(v) for v in df[c]], dtype=object)
    return df


def read_cog_def(path: Path) -> pd.DataFrame:
    """NCBI COG 2020 ``cog-20.def.tab`` (no header; COG ID, category letters, name, ...) as a frame
    with columns ``cog_id``, ``cog_category``, ``name``. The file is not always UTF-8."""
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    rows = []
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0].strip():
            rows.append((parts[0].strip(), parts[1].strip(), parts[2] if len(parts) > 2 else ""))
    return pd.DataFrame(rows, columns=["cog_id", "cog_category", "name"])
