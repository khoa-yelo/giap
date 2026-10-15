"""Protein FASTA reader (GIAP_SPEC.md §5.3 Proteins FASTA, §7.0; ADR-0005 item 1).

The reader reports facts (resolved IDs, sequences, letter counts); validation rules such as
duplicate IDs, warnings for non-standard letters and length-class reporting are applied by
preflight checks (M1-10).
"""

from __future__ import annotations

import re
import string
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from giap.io._lines import FastaFormatError, read_lines
from giap.io._open import open_text_write

STANDARD_AA: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWYU")
COUNTED_AA: frozenset[str] = frozenset("XBZJ")
Route = Literal["peptide", "protein", "dropped"]

_UNIPROT = re.compile(r"^(sp|tr)\|([^|]+)\|")
# ASCII-only upper-casing: str.upper() would turn e.g. "ß" into "SS" (valid residues).
_ASCII_UPPER = str.maketrans(string.ascii_lowercase, string.ascii_uppercase)
_WHITESPACE = re.compile(r"\s+")


class LetterCounts(tuple[tuple[str, int], ...]):
    """Hashable sorted ``(letter, count)`` pairs; ``dict(counts)`` gives the mapping view."""

    __slots__ = ()

    def values(self) -> tuple[int, ...]:
        return tuple(n for _, n in self)


@dataclass(frozen=True, slots=True)
class ProteinRecord:
    """One protein record; ``nonstandard``/``other`` are sorted ``(letter, count)`` pairs."""

    id: str
    members: tuple[str, ...]
    header: str
    sequence: str
    internal_stops: int
    nonstandard: LetterCounts
    other: LetterCounts


def resolve_protein_id(token: str) -> tuple[str, tuple[str, ...]]:
    """Split a protein group on ``;`` and resolve UniProt ``sp|ACC|NAME`` / ``tr|ACC|NAME``.

    Empty and repeated members are dropped; the first occurrence keeps its position.
    """
    members: list[str] = []
    for part in token.split(";"):
        part = part.strip()
        match = _UNIPROT.match(part)
        if match is not None:
            part = match.group(2).strip()
        if part and part not in members:
            members.append(part)
    return ";".join(members), tuple(members)


def route_by_length(length_aa: int, peptide_len: tuple[int, int] = (7, 40)) -> Route:
    """``< min`` → dropped; ``min ≤ L ≤ max`` → peptide; ``> max`` → protein."""
    lo, hi = peptide_len
    if length_aa < lo:
        return "dropped"
    if length_aa <= hi:
        return "peptide"
    return "protein"


def _make_record(header: str, chunks: list[str], line_no: int) -> ProteinRecord:
    tokens = header.split(maxsplit=1)
    group_id, members = resolve_protein_id(tokens[0] if tokens else "")
    if not members:
        raise FastaFormatError("header has no protein ID", line_no)
    sequence = "".join(chunks).translate(_ASCII_UPPER)
    if sequence.endswith("*"):
        sequence = sequence[:-1]
    counts = Counter(sequence)
    stops = counts.pop("*", 0)
    nonstandard = LetterCounts(sorted((c, n) for c, n in counts.items() if c in COUNTED_AA))
    other = LetterCounts(
        sorted((c, n) for c, n in counts.items() if c not in STANDARD_AA and c not in COUNTED_AA)
    )
    return ProteinRecord(
        id=group_id,
        members=members,
        header=header,
        sequence=sequence,
        internal_stops=stops,
        nonstandard=nonstandard,
        other=other,
    )


def read_protein_fasta(path: Path) -> list[ProteinRecord]:
    """Read a (possibly gzipped) protein FASTA; records are returned in file order."""
    records: list[ProteinRecord] = []
    header: str | None = None
    header_line = 0
    chunks: list[str] = []
    for line_no, line in read_lines(path):
        if line.startswith(">"):
            if header is not None:
                records.append(_make_record(header, chunks, header_line))
            header = line[1:].strip()
            header_line = line_no
            chunks = []
        elif not line.strip():
            continue
        elif header is None:
            raise FastaFormatError("sequence data before the first '>' header", line_no)
        else:
            chunks.append(_WHITESPACE.sub("", line))
    if header is not None:
        records.append(_make_record(header, chunks, header_line))
    return records


def write_protein_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> None:
    """Write ``(header, sequence)`` pairs (header without ``>``) wrapped at ``width`` residues."""
    if width < 1:
        raise ValueError(f"width must be >= 1, got {width}")
    lines: list[str] = []
    for header, seq in records:
        if not header or "\n" in header or "\r" in header:
            raise ValueError(f"invalid FASTA header {header!r}")
        if not seq:
            raise ValueError(f"empty sequence for {header!r}")
        lines.append(f">{header}")
        lines.extend(seq[i : i + width] for i in range(0, len(seq), width))
    with open_text_write(path, compress=False) as fh:
        fh.write("".join(f"{line}\n" for line in lines))
