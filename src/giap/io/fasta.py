"""Genome FASTA reader and writer (GIAP_SPEC.md §5.3, §7.0; ADR-0005 item 1).

The reader reports facts (IDs, sequences, topology tags); validation rules such as duplicate IDs,
invalid letters and "at least one record" are applied by preflight checks (M1-09).
"""

from __future__ import annotations

import gzip
import re
import zlib
from collections import Counter
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from giap.io._open import open_text

IUPAC_DNA: frozenset[str] = frozenset("ACGTNRYSWKMBDHV")
SAFE_ID = re.compile(r"[A-Za-z0-9._-]")

_UNSAFE_ID = re.compile(r"[^A-Za-z0-9._-]")
_CIRCULAR_TAG = re.compile(r"(?:^|\s)circular=true(?=\s|$)|\[topology=circular\]", re.IGNORECASE)
_WHITESPACE = re.compile(r"\s+")


class FastaFormatError(ValueError):
    """Malformed FASTA input; ``line_no`` is the 1-based line where the problem was found."""

    def __init__(self, message: str, line_no: int) -> None:
        super().__init__(f"line {line_no}: {message}")
        self.line_no = line_no


@dataclass(frozen=True, slots=True)
class GenomeRecord:
    id: str
    original_id: str
    header: str
    sequence: str
    circular_tag: bool


def sanitize_id(raw: str) -> str:
    """Map every character outside ``[A-Za-z0-9._-]`` to ``_`` (length-preserving)."""
    return _UNSAFE_ID.sub("_", raw)


def invalid_letters(seq: str) -> dict[str, int]:
    """Count characters of ``seq`` that are not upper-case IUPAC nucleotide codes."""
    return dict(Counter(c for c in seq if c not in IUPAC_DNA))


def _lines(path: Path) -> Iterator[tuple[int, str]]:
    """Yield ``(line_no, line)`` with line endings removed; decoding errors become format errors."""
    line_no = 0
    try:
        with open_text(path) as fh:
            for line_no, line in enumerate(fh, start=1):
                yield line_no, line.rstrip("\r\n")
    except (EOFError, gzip.BadGzipFile, zlib.error, UnicodeDecodeError) as exc:
        raise FastaFormatError(f"cannot read {path}: {exc}", line_no + 1) from exc


def _make_record(header: str, chunks: list[str]) -> GenomeRecord:
    original_id = header.split(maxsplit=1)[0]
    return GenomeRecord(
        id=sanitize_id(original_id),
        original_id=original_id,
        header=header,
        sequence="".join(chunks).upper(),
        circular_tag=_CIRCULAR_TAG.search(header) is not None,
    )


def read_genome_fasta(path: Path) -> list[GenomeRecord]:
    """Read a (possibly gzipped) genome FASTA; records are returned in file order."""
    records: list[GenomeRecord] = []
    header: str | None = None
    chunks: list[str] = []
    for line_no, line in _lines(path):
        if line.startswith(">"):
            if header is not None:
                records.append(_make_record(header, chunks))
            header = line[1:].strip()
            if not header:
                raise FastaFormatError("header has no ID token", line_no)
            chunks = []
        elif not line.strip():
            continue
        elif header is None:
            raise FastaFormatError("sequence data before the first '>' header", line_no)
        else:
            chunks.append(_WHITESPACE.sub("", line))
    if header is not None:
        records.append(_make_record(header, chunks))
    return records


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 80) -> None:
    """Write ``(id, sequence)`` pairs as uncompressed FASTA wrapped at ``width`` columns."""
    if width < 1:
        raise ValueError(f"width must be >= 1, got {width}")
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for rec_id, seq in records:
            fh.write(f">{rec_id}\n")
            for start in range(0, len(seq), width):
                fh.write(seq[start : start + width])
                fh.write("\n")
