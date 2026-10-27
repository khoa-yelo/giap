"""Writers and readers for ``tables/autocorrelation.tsv`` and ``tables/islands.tsv``.

GIPPA_SPEC.md §9.2, §9.3, §9.5 (the shipped Frictionless schemas), §6.5 and ADR-0006; M5-09
(docs/specs/M5-09-context-tables-step.md, docs/methods/context_tables.md).

Cells go to M1-04 ``write_tsv`` **by value** (``%.6g`` floats, integers as integers, ``NA`` for
missing, ``;`` lists), so nothing is silently coerced; only ``start`` is converted, to 1-based
(``end`` unchanged, so an origin-crossing island keeps ``end > length``). An empty frame writes
the header line only, and a flag outside :data:`gippa.context.vocab.CONTEXT_FLAGS` raises
``ValueError`` naming it. The readers are the inverses and type each column from the schema:
``integer`` → nullable ``Int64``, ``number`` → ``float64``, ``string`` → ``object`` of ``str``,
and the list cells (``flags``, ``genes``) → tuples of ``str``.

``pandas.Series.flags`` shadows the ``flags`` *column*, so every access here is by name or
position, never ``row.flags``.
"""

from __future__ import annotations

__all__ = ["read_autocorrelation", "read_islands", "write_autocorrelation", "write_islands"]

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from gippa.context.autocorrelation import AUTOCORR_COLUMNS
from gippa.context.islands import ISLAND_COLUMNS
from gippa.context.vocab import require_known_flags
from gippa.io.tsv import read_tsv, write_tsv


@dataclass(frozen=True, slots=True)
class _Table:
    """A table's name, column order and the field types of its schema."""

    name: str
    columns: tuple[str, ...]
    ints: frozenset[str]
    floats: frozenset[str]
    lists: tuple[str, ...]


#: ``autocorrelation.schema.json`` (no coordinates; ``lag`` is an integer carrying NA).
AUTOCORRELATION = _Table(
    "autocorrelation", AUTOCORR_COLUMNS, frozenset({"lag", "n_perm"}),
    frozenset({"value", "null_mean", "null_lo", "null_hi", "effect", "ci_low", "ci_high",
               "p_value", "q_value"}),
    ("flags",),
)  # fmt: skip
#: ``islands.schema.json`` (1-based closed; ``stat`` is a name, its value is ``scan_s``).
ISLANDS = _Table(
    "islands", ISLAND_COLUMNS,
    frozenset({"start", "end", "n_tus", "n_genes", "block_len", "n_perm"}),
    frozenset({"mean_abundance_log2", "detected_frac", "scan_s", "effect", "ci_low", "ci_high",
               "p_value", "q_value"}),
    ("genes", "flags"),
)  # fmt: skip


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


def _row(row: Sequence[object], table: _Table) -> list[object]:
    """One frame row as TSV cells: ``start`` shifted, lists joined, flags checked."""
    cells: list[object] = []
    for column, value in zip(table.columns, row, strict=True):
        if column == "flags":
            cells.append(list(require_known_flags(_list_cell(value), where=table.name)) or None)
        elif column in table.lists:
            cells.append(list(_list_cell(value)) or None)
        elif _is_na(value):
            cells.append(None)
        elif column == "start":
            cells.append(int(value) + 1)  # type: ignore[call-overload]
        else:
            cells.append(value)
    return cells


def _write(path: Path, df: pd.DataFrame, table: _Table) -> None:
    frame = df[list(table.columns)]
    rows = [_row(r, table) for r in frame.itertuples(index=False, name=None)]
    write_tsv(path, table.columns, rows)


def _column(values: list[str | None], column: str, table: _Table) -> pd.Series:
    """One read column, typed as the table schema declares it."""
    if column in table.lists:
        return pd.Series([_list_cell(v) for v in values], dtype=object)
    if column in table.ints:
        return pd.Series([None if v is None else int(v) for v in values], dtype="Int64")
    if column in table.floats:
        return pd.Series([math.nan if v is None else float(v) for v in values], dtype="float64")
    return pd.Series([None if v is None else str(v) for v in values], dtype=object)


def _read(path: Path, table: _Table) -> pd.DataFrame:
    header, rows = read_tsv(path)
    if tuple(header) != table.columns:
        raise ValueError(f"{path}: {table.name} header is {header}, expected {table.columns}")
    out = pd.DataFrame(
        {c: _column([r[i] for r in rows], c, table) for i, c in enumerate(table.columns)},
        columns=list(table.columns),
    )
    if "start" in table.columns:
        out["start"] = out["start"] - 1  # back to 0-based half-open (§6.5)
    return out


def write_autocorrelation(path: Path, df: pd.DataFrame) -> None:
    """Write an ``AUTOCORR_COLUMNS`` frame as ``tables/autocorrelation.tsv`` (§9.3), rows in
    frame order; a ``flags`` item outside ``CONTEXT_FLAGS`` raises ``ValueError``."""
    _write(path, df, AUTOCORRELATION)


def write_islands(path: Path, df: pd.DataFrame) -> None:
    """Write a 0-based ``ISLAND_COLUMNS`` frame as ``tables/islands.tsv``, 1-based closed
    (``start + 1``, ``end`` unchanged; §9.3, §6.5); an unknown ``flags`` item raises."""
    _write(path, df, ISLANDS)


def read_autocorrelation(path: Path) -> pd.DataFrame:
    """``tables/autocorrelation.tsv`` as an ``AUTOCORR_COLUMNS`` frame with the schema's dtypes
    (inverse of :func:`write_autocorrelation`); a wrong header raises."""
    return _read(path, AUTOCORRELATION)


def read_islands(path: Path) -> pd.DataFrame:
    """``tables/islands.tsv`` as an ``ISLAND_COLUMNS`` frame with the schema's dtypes and
    ``start`` back to 0-based (inverse of :func:`write_islands`)."""
    return _read(path, ISLANDS)
