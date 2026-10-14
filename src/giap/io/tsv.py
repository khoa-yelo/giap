"""TSV conventions (GIAP_SPEC.md §6.6, §9.2, §9.3).

UTF-8, one header row, snake_case column names, ``NA`` for missing, floats ``%.6g``, integers as
integers, booleans ``true``/``false``, lists ``;``-joined, ``\\n`` line endings. A path ending in
``.gz`` is written with a fixed gzip header (``mtime=0``, no file name) so bytes are reproducible.
"""

from __future__ import annotations

import gzip
import io
import math
import numbers
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

NA = "NA"
COLUMN_NAME = re.compile(r"^[a-z][a-z0-9_]*$")
Cell = str | int | float | bool | None | Sequence[str]

_FORBIDDEN = re.compile(r"[\t\n\r]")
_GZIP_MAGIC = b"\x1f\x8b"


def _is_missing(v: object) -> bool:
    # ``pandas.NA`` is detected by type name so this module does not import pandas.
    return v is None or type(v).__name__ == "NAType"


def format_value(v: object) -> str:
    """Format one cell following §6.6 (see module docstring)."""
    if _is_missing(v):
        return NA
    if isinstance(v, str):
        if _FORBIDDEN.search(v):
            raise ValueError(f"tab or newline inside TSV cell: {v!r}")
        return v
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, numbers.Integral):
        return str(int(v))
    if isinstance(v, numbers.Real):
        f = float(v)
        if math.isnan(f):
            return NA
        if f == 0:
            return "0"
        return f"{f:.6g}"
    if isinstance(v, list | tuple):
        return ";".join(format_value(x) for x in v) if v else NA
    raise TypeError(f"cannot format {type(v).__name__} value {v!r} as a TSV cell")


def write_text(path: Path, text: str) -> None:
    """Write UTF-8 ``text``; a ``.gz`` path gets a reproducible gzip header."""
    data = text.encode("utf-8")
    if path.name.endswith(".gz"):
        with (
            path.open("wb") as raw,
            gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz,
        ):
            gz.write(data)
    else:
        path.write_bytes(data)


def open_text(path: Path) -> io.TextIOWrapper:
    """Open ``path`` for UTF-8 reading, decompressing gzip detected by magic bytes."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    # The caller owns the returned handle (used as a context manager).
    raw = gzip.open(path, "rb") if magic == _GZIP_MAGIC else path.open("rb")  # noqa: SIM115
    return io.TextIOWrapper(raw, encoding="utf-8-sig", newline="")


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    """Write a TSV table; column names must be snake_case and every row as wide as the header."""
    for c in columns:
        if not COLUMN_NAME.match(c):
            raise ValueError(f"column name {c!r} is not snake_case ([a-z][a-z0-9_]*)")
    width = len(columns)
    lines = ["\t".join(columns)]
    for i, row in enumerate(rows, start=1):
        if len(row) != width:
            raise ValueError(f"row {i} has {len(row)} cells, header has {width}")
        lines.append("\t".join(format_value(v) for v in row))
    write_text(path, "\n".join(lines) + "\n")


def read_tsv(path: Path) -> tuple[list[str], list[list[str | None]]]:
    """Read a TSV table; ``NA`` and empty cells become ``None``; ragged rows raise ``ValueError``.

    Lines are split on ``\\n`` only, so other Unicode line separators stay inside cells.
    """
    with open_text(path) as fh:
        text = fh.read()
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise ValueError(f"{path}: empty TSV (no header)")
    header = lines[0].rstrip("\r").split("\t")
    rows: list[list[str | None]] = []
    for lineno, line in enumerate(lines[1:], start=2):
        cells = line.rstrip("\r").split("\t")
        if len(cells) != len(header):
            raise ValueError(
                f"{path}: line {lineno} has {len(cells)} fields, header has {len(header)}"
            )
        rows.append([None if c in (NA, "") else c for c in cells])
    return header, rows
