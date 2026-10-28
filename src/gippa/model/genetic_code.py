"""The NCBI genetic codes GIPPA supports: 4, 11 and 25 (ADR-0067).

The tables are built from the NCBI TCAG strings, not imported from Biopython, so numba kernels
and LUT builders do not depend on it. Only TGA differs between them: a stop under 11, W under 4
and G under 25. Pure: no I/O, results cached.
"""

from __future__ import annotations

import itertools
from collections.abc import Mapping
from functools import cache
from types import MappingProxyType

__all__ = [
    "SUPPORTED_TABLES",
    "codon_table",
    "minus_stop_codons",
    "recoded_codons",
    "stop_codons",
    "table_name",
]

SUPPORTED_TABLES: tuple[int, ...] = (4, 11, 25)

_TCAG: dict[int, str] = {
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    25: "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
_NAMES: dict[int, str] = {
    11: "Bacterial, Archaeal and Plant Plastid",
    4: "Mycoplasma and Spiroplasma",
    25: "Candidate Division SR1 and Gracilibacteria",
}
_CODONS = tuple("".join(c) for c in itertools.product("TCAG", repeat=3))
_COMP = str.maketrans("ACGT", "TGCA")


def _check(table: int) -> None:
    if table not in SUPPORTED_TABLES:
        raise ValueError(f"unsupported translation table {table}; GIPPA supports 4, 11 and 25")


@cache
def codon_table(table: int) -> Mapping[str, str]:
    """The 64 upper-case ACGT codons -> one-letter amino acid, ``*`` for a stop."""
    _check(table)
    return MappingProxyType(dict(zip(_CODONS, _TCAG[table], strict=True)))


@cache
def stop_codons(table: int) -> frozenset[str]:
    """The stop codons of ``table`` on the coding strand."""
    return frozenset(c for c, aa in codon_table(table).items() if aa == "*")


@cache
def minus_stop_codons(table: int) -> frozenset[str]:
    """The reverse complements of :func:`stop_codons`, as read on the plus string."""
    return frozenset(c.translate(_COMP)[::-1] for c in stop_codons(table))


@cache
def recoded_codons(table: int) -> Mapping[str, str]:
    """Bakta's recoded letter -> the stop codon it must sit over: Sec (U) over TGA only where
    TGA is a stop (table 11); Pyl (O) over TAG in every table."""
    stops = stop_codons(table)
    pairs = {"U": "TGA", "O": "TAG"}
    return MappingProxyType({k: v for k, v in pairs.items() if v in stops})


@cache
def table_name(table: int) -> str:
    """The short NCBI name of ``table`` (used in messages and the report)."""
    _check(table)
    return _NAMES[table]
