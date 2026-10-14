"""Open text inputs that may be gzip-compressed (GIAP_SPEC.md §5.3: "gzip input is accepted")."""

from __future__ import annotations

import gzip
import io
from pathlib import Path
from typing import IO, TextIO, cast

GZIP_MAGIC = b"\x1f\x8b"


def open_binary(path: Path) -> IO[bytes]:
    """Open ``path`` for reading as bytes; gzip is detected by magic bytes, not by suffix."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == GZIP_MAGIC:
        return cast("IO[bytes]", gzip.open(path, "rb"))
    return path.open("rb")


def open_text(path: Path, *, encoding: str = "utf-8-sig", newline: str | None = None) -> TextIO:
    """Open ``path`` for reading as UTF-8 text (a leading BOM is dropped); gzip by magic bytes.

    ``encoding`` and ``newline`` are passed to :class:`io.TextIOWrapper`.
    """
    return io.TextIOWrapper(open_binary(path), encoding=encoding, newline=newline)


def open_text_write(path: Path, *, compress: bool | None = None) -> TextIO:
    """Open ``path`` for writing UTF-8 text with ``\\n`` newlines.

    Output is gzip-compressed when ``compress`` is true, or (by default) when ``path`` ends with
    ``.gz``; the gzip header is reproducible (``mtime=0``, empty file name).
    """
    if compress is None:
        compress = path.name.endswith(".gz")
    if not compress:
        return path.open("w", encoding="utf-8", newline="\n")
    fh = io.FileIO(path, "wb")
    raw = gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0)
    raw.myfileobj = fh  # GzipFile.close() then also closes the underlying file
    return io.TextIOWrapper(raw, encoding="utf-8", newline="\n")
