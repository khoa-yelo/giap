"""Writers and readers for ``tables/autocorrelation.tsv`` and ``tables/islands.tsv``.

GIAP_SPEC.md §9.2 (table conventions and test columns), §9.3 (the two column dictionaries), §9.5
(the Frictionless resources in ``giap/schemas/tables/``), §6.5 (0-based half-open internally,
1-based closed in the writer) and ADR-0006 (``NA``, ``;`` lists, ``%.6g``); M5-09
(docs/specs/M5-09-context-tables-step.md).

* ``autocorrelation.tsv`` carries no coordinates, so the writer is a pure formatting pass.
* ``islands.tsv`` is 1-based closed: ``start + 1`` with ``end`` unchanged (a half-open 0-based
  end equals a closed 1-based end), and an origin-crossing island keeps ``end > length``.
* An empty frame writes the header line only; a flag outside
  :data:`giap.context.vocab.CONTEXT_FLAGS` raises ``ValueError`` naming it (§9.4 is closed).
* The readers are the exact inverses: ``integer`` fields of the schema come back as pandas
  nullable ``Int64``, ``number`` fields as ``float64``, string fields as ``object`` holding
  ``str``, and the list-valued cells (``flags``, ``genes``) as tuples of ``str``.

The ``flags`` column name collides with ``pandas.Series.flags``, so every access here is
``row["flags"]`` / positional, never ``row.flags``.
"""

from __future__ import annotations

__all__ = [
    "read_autocorrelation",
    "read_islands",
    "write_autocorrelation",
    "write_islands",
]

import math
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.context.autocorrelation import AUTOCORR_COLUMNS
from giap.context.islands import ISLAND_COLUMNS
from giap.context.vocab import require_known_flags
from giap.io.tsv import read_tsv, write_tsv

#: ``integer`` fields of ``autocorrelation.schema.json`` (``lag`` carries NA).
AUTOCORR_INTS: frozenset[str] = frozenset({"lag", "n_perm"})
#: ``number`` fields of ``autocorrelation.schema.json``.
AUTOCORR_FLOATS: frozenset[str] = frozenset(
    {"value", "null_mean", "null_lo", "null_hi", "effect", "ci_low", "ci_high", "p_value",
     "q_value"}
)  # fmt: skip
#: ``integer`` fields of ``islands.schema.json``.
ISLAND_INTS: frozenset[str] = frozenset({"start", "end", "n_tus", "n_genes", "block_len", "n_perm"})
#: ``number`` fields of ``islands.schema.json`` (``stat`` is the scan statistic S, a number).
ISLAND_FLOATS: frozenset[str] = frozenset(
    {"mean_abundance_log2", "detected_frac", "scan_s", "stat", "effect", "ci_low", "ci_high",
     "p_value", "q_value"}
)  # fmt: skip
#: ``;``-joined list cells per table.
AUTOCORR_LISTS: tuple[str, ...] = ("flags",)
ISLAND_LISTS: tuple[str, ...] = ("genes", "flags")


def _is_na(value: object) -> bool:
    """Whether a frame cell is missing (``None``, ``pandas.NA`` or NaN)."""
    if value is None or value is pd.NA:
        return True
    return isinstance(value, float | np.floating) and bool(math.isnan(float(value)))


def _list_cell(value: object) -> tuple[str, ...]:
    """A list-valued cell as a tuple of ``str`` (a missing cell is the empty tuple)."""
    if _is_na(value):
        return ()
    if isinstance(value, str):
        return tuple(part for part in value.split(";") if part)
    if isinstance(value, Sequence | np.ndarray | pd.Series):
        return tuple(str(v) for v in value)
    raise TypeError(f"cannot read {type(value).__name__} {value!r} as a list cell")


def _format_row(
    row: Sequence[object], columns: Sequence[str], *, lists: Sequence[str], table: str
) -> list[object]:
    """One frame row as TSV cells: ``start`` shifted, lists joined, flags checked.

    Everything else is handed to ``giap.io.tsv.format_value`` by value (``%.6g`` floats,
    integers as integers, strings as they are), so a cell that is not of its schema type is
    written as it stands rather than silently coerced.
    """
    cells: list[object] = []
    for column, value in zip(columns, row, strict=True):
        if column == "flags":
            cells.append(list(require_known_flags(_list_cell(value), where=table)) or None)
        elif column in lists:
            cells.append(list(_list_cell(value)) or None)
        elif _is_na(value):
            cells.append(None)
        elif column == "start":
            cells.append(int(value) + 1)  # type: ignore[call-overload]
        else:
            cells.append(value)
    return cells


def _write(
    path: Path, df: pd.DataFrame, columns: tuple[str, ...], *, lists: Sequence[str], table: str
) -> None:
    frame = df[list(columns)]
    rows = [
        _format_row(row, columns, lists=lists, table=table)
        for row in frame.itertuples(index=False, name=None)
    ]
    write_tsv(path, columns, rows)


def write_autocorrelation(path: Path, df: pd.DataFrame) -> None:
    """Write an ``autocorrelation`` frame as ``tables/autocorrelation.tsv`` (§9.3).

    Parameters
    ----------
    path : pathlib.Path
        Destination; a ``.gz`` suffix gzips it deterministically (M1-04).
    df : pandas.DataFrame
        Frame with the columns ``giap.context.autocorrelation.AUTOCORR_COLUMNS`` (any order);
        the rows are written in frame order. ``flags`` holds a sequence of §9.4 flags or NA.

    Raises
    ------
    ValueError
        If a ``flags`` item is outside ``giap.context.vocab.CONTEXT_FLAGS``.
    """
    _write(path, df, AUTOCORR_COLUMNS, lists=AUTOCORR_LISTS, table="autocorrelation")


def write_islands(path: Path, df: pd.DataFrame) -> None:
    """Write an ``islands`` frame as ``tables/islands.tsv``, 1-based closed (§9.3, §6.5).

    Parameters
    ----------
    path : pathlib.Path
        Destination; a ``.gz`` suffix gzips it deterministically (M1-04).
    df : pandas.DataFrame
        Frame with the columns ``giap.context.islands.ISLAND_COLUMNS`` (any order), 0-based
        half-open. ``start`` is written ``start + 1`` and ``end`` unchanged, so an island that
        crosses the origin keeps ``end`` beyond the replicon length.

    Raises
    ------
    ValueError
        If a ``flags`` item is outside ``giap.context.vocab.CONTEXT_FLAGS``.
    """
    _write(path, df, ISLAND_COLUMNS, lists=ISLAND_LISTS, table="islands")


def _parse(cell: str | None, column: str, kind: type[int] | type[float]) -> int | float | None:
    """One numeric cell; ``ValueError`` names the column when the cell is not a number."""
    if cell is None:
        return None
    try:
        return kind(cell)
    except ValueError:
        raise ValueError(
            f"column {column!r} is {kind.__name__} in the table schema, got {cell!r}"
        ) from None


def _column(
    values: list[str | None],
    column: str,
    *,
    ints: frozenset[str],
    floats: frozenset[str],
    lists: Sequence[str],
) -> pd.Series:
    """One read column, typed as the table schema declares it."""
    if column in lists:
        return pd.Series([_list_cell(v) for v in values], dtype=object)
    if column in ints:
        return pd.Series([_parse(v, column, int) for v in values], dtype="Int64")
    if column in floats:
        parsed = [_parse(v, column, float) for v in values]
        return pd.Series([math.nan if v is None else v for v in parsed], dtype="float64")
    return pd.Series([None if v is None else str(v) for v in values], dtype=object)


def _read(
    path: Path,
    columns: tuple[str, ...],
    *,
    ints: frozenset[str],
    floats: frozenset[str],
    lists: Sequence[str],
    table: str,
) -> pd.DataFrame:
    header, rows = read_tsv(path)
    if tuple(header) != columns:
        raise ValueError(f"{path}: {table} header is {header}, expected {list(columns)}")
    data = {
        column: _column([r[i] for r in rows], column, ints=ints, floats=floats, lists=lists)
        for i, column in enumerate(columns)
    }
    out = pd.DataFrame(data, columns=list(columns))
    if "start" in columns:
        out["start"] = out["start"] - 1  # back to 0-based half-open (§6.5)
    return out


def read_autocorrelation(path: Path) -> pd.DataFrame:
    """Read ``tables/autocorrelation.tsv`` back into a frame (inverse of the writer).

    Parameters
    ----------
    path : pathlib.Path
        A table written by :func:`write_autocorrelation`.

    Returns
    -------
    pandas.DataFrame
        ``AUTOCORR_COLUMNS`` with the schema's dtypes: ``Int64`` for ``lag`` and ``n_perm``,
        ``float64`` for the number fields, ``object`` of ``str`` elsewhere and a tuple of
        ``str`` in ``flags``.

    Raises
    ------
    ValueError
        If the header is not ``AUTOCORR_COLUMNS``.
    """
    return _read(path, AUTOCORR_COLUMNS, ints=AUTOCORR_INTS, floats=AUTOCORR_FLOATS,
                 lists=AUTOCORR_LISTS, table="autocorrelation")  # fmt: skip


def read_islands(path: Path) -> pd.DataFrame:
    """Read ``tables/islands.tsv`` back into a 0-based half-open frame (inverse of the writer).

    Parameters
    ----------
    path : pathlib.Path
        A table written by :func:`write_islands`.

    Returns
    -------
    pandas.DataFrame
        ``ISLAND_COLUMNS`` with the schema's dtypes (``Int64`` integers, ``float64`` numbers,
        tuples of ``str`` in ``genes`` and ``flags``) and ``start`` converted back to 0-based.

    Raises
    ------
    ValueError
        If the header is not ``ISLAND_COLUMNS``.
    """
    return _read(path, ISLAND_COLUMNS, ints=ISLAND_INTS, floats=ISLAND_FLOATS,
                 lists=ISLAND_LISTS, table="islands")  # fmt: skip
