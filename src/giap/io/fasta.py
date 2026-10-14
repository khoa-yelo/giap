"""Genome FASTA reader and writer (GIAP_SPEC.md §5.3, §7.0; ADR-0005 item 1).

The reader reports facts (IDs, sequences, topology tags); validation rules such as duplicate IDs,
invalid letters and "at least one record" are applied by preflight checks (M1-09).
"""

from __future__ import annotations

import re
import string
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from giap.io._lines import FastaFormatError, read_lines
from giap.io._open import open_text_write

IUPAC_DNA: frozenset[str] = frozenset("ACGTNRYSWKMBDHV")
SAFE_ID = re.compile(r"[A-Za-z0-9._-]")

_UNSAFE_ID = re.compile(r"[^A-Za-z0-9._-]")
_SAFE_ID_FULL = re.compile(r"[A-Za-z0-9._-]+")
_CIRCULAR_TAG = re.compile(
    r"(?:^|(?<=[\s\[]))circular=true(?=$|[\s;,\]])|\[topology=circular\]", re.IGNORECASE
)
_ASCII_UPPER = str.maketrans(string.ascii_lowercase, string.ascii_uppercase)
_WHITESPACE = re.compile(r"\s+")

__all__ = [
    "IUPAC_DNA",
    "SAFE_ID",
    "FastaFormatError",
    "GenomeRecord",
    "invalid_letters",
    "read_genome_fasta",
    "sanitize_id",
    "write_fasta",
]


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


def _make_record(header: str, chunks: list[str]) -> GenomeRecord:
    original_id = header.split(maxsplit=1)[0]
    return GenomeRecord(
        id=sanitize_id(original_id),
        original_id=original_id,
        header=header,
        sequence="".join(chunks).translate(_ASCII_UPPER),
        circular_tag=_CIRCULAR_TAG.search(header) is not None,
    )


def read_genome_fasta(path: Path) -> list[GenomeRecord]:
    """Read a (possibly gzipped) genome FASTA; records are returned in file order."""
    records: list[GenomeRecord] = []
    header: str | None = None
    chunks: list[str] = []
    for line_no, line in read_lines(path):
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
    """Write ``(id, sequence)`` pairs as uncompressed FASTA wrapped at ``width`` columns.

    IDs must consist only of ``[A-Za-z0-9._-]`` and sequences must not contain line breaks.
    """
    if width < 1:
        raise ValueError(f"width must be >= 1, got {width}")
    with open_text_write(path, compress=False) as fh:
        for rec_id, seq in records:
            if not _SAFE_ID_FULL.fullmatch(rec_id):
                raise ValueError(f"unsafe FASTA ID {rec_id!r}; expected {_SAFE_ID_FULL.pattern}")
            if "\n" in seq or "\r" in seq:
                raise ValueError(f"sequence for {rec_id!r} contains a line break")
            fh.write(f">{rec_id}\n")
            for start in range(0, len(seq), width):
                fh.write(seq[start : start + width])
                fh.write("\n")
