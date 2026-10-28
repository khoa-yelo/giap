"""Shared readers and the per-table translation oracle of the M11-16c simulator tests."""

from __future__ import annotations

import csv
from pathlib import Path

STOPS_4_25 = frozenset({"TAA", "TAG"})
STARTS = frozenset({"ATG", "GTG", "TTG"})
_COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def codons(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]


def translate(nt: str, table: int) -> str:
    """CDS protein under ``table``: initiator read as M, the final stop dropped."""
    from gippasim.codons import code

    t = code(table)
    return "M" + "".join(t[c] for c in codons(nt)[1:-1])


def fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = ""
        elif line:
            out[name] += line.strip()
    return out


def tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def slice_(seq: str, s0: int, e0: int, strand: str) -> str:
    """``[s0, e0)`` of a contig (wrapping past its end), read on ``strand``."""
    n = len(seq)
    nt = "".join(seq[i % n] for i in range(s0, e0))
    return rc(nt) if strand == "-" else nt


def codon_at(seq: str, c0: int, strand: str) -> str:
    return slice_(seq, c0, c0 + 3, strand)
