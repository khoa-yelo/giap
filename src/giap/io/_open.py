"""Open text inputs that may be gzip-compressed (GIAP_SPEC.md §5.3: "gzip input is accepted")."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO

GZIP_MAGIC = b"\x1f\x8b"


def open_text(path: Path) -> TextIO:
    """Open ``path`` for reading as UTF-8 text; gzip is detected by magic bytes, not by suffix."""
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == GZIP_MAGIC:
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")
