"""Per-line byte decoding shared by the FASTA readers (GIAP_SPEC.md §5.3, §7.0)."""

from __future__ import annotations

import gzip
import zlib
from collections.abc import Iterator
from pathlib import Path

from giap.io._open import open_binary


class FastaFormatError(ValueError):
    """Malformed FASTA input; ``line_no`` is the 1-based line where the problem was found."""

    def __init__(self, message: str, line_no: int) -> None:
        super().__init__(f"line {line_no}: {message}")
        self.line_no = line_no


def read_lines(path: Path) -> Iterator[tuple[int, str]]:
    """Yield ``(line_no, line)`` with line endings (LF, CR, CRLF) removed.

    Lines are decoded one at a time (UTF-8, leading BOM dropped on line 1) so a decoding error
    reports the line that holds the bad bytes; read errors become :class:`FastaFormatError`.
    """
    line_no = 0
    try:
        with open_binary(path) as fh:
            for chunk in fh:
                for raw in chunk.splitlines():
                    line_no += 1
                    yield line_no, raw.decode("utf-8-sig" if line_no == 1 else "utf-8")
    except UnicodeDecodeError as exc:
        raise FastaFormatError(f"cannot decode {path}: {exc}", line_no) from exc
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise FastaFormatError(f"cannot read {path}: {exc}", line_no + 1) from exc
