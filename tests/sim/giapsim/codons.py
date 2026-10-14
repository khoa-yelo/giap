"""Genetic code (table 11) and a bundled E. coli-like codon usage table."""

from __future__ import annotations

BASES = "TCAG"
_AA_11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

#: Table-11 translation of every codon.
CODE_11: dict[str, str] = {
    a + b + c: _AA_11[16 * i + 4 * j + k]
    for i, a in enumerate(BASES)
    for j, b in enumerate(BASES)
    for k, c in enumerate(BASES)
}

STOPS = frozenset(c for c, aa in CODE_11.items() if aa == "*")

#: E. coli K-12-like usage, per thousand codons (sense codons only; stops are drawn separately).
SENSE_USAGE: dict[str, float] = {
    "TTT": 22.1, "TTC": 16.0, "TTA": 14.3, "TTG": 13.0,
    "TCT": 10.4, "TCC": 9.1, "TCA": 8.9, "TCG": 8.5,
    "TAT": 17.5, "TAC": 12.2, "TGT": 5.2, "TGC": 6.1, "TGG": 13.9,
    "CTT": 11.9, "CTC": 10.2, "CTA": 4.2, "CTG": 48.4,
    "CCT": 7.5, "CCC": 5.4, "CCA": 8.6, "CCG": 20.9,
    "CAT": 12.5, "CAC": 9.3, "CAA": 14.6, "CAG": 28.4,
    "CGT": 21.0, "CGC": 22.0, "CGA": 3.8, "CGG": 5.9,
    "ATT": 29.8, "ATC": 23.7, "ATA": 6.8, "ATG": 26.4,
    "ACT": 10.3, "ACC": 22.0, "ACA": 9.3, "ACG": 13.7,
    "AAT": 20.6, "AAC": 24.4, "AAA": 33.6, "AAG": 12.1,
    "AGT": 9.9, "AGC": 15.2, "AGA": 3.6, "AGG": 2.1,
    "GTT": 19.8, "GTC": 14.3, "GTA": 11.6, "GTG": 24.4,
    "GCT": 17.1, "GCC": 24.2, "GCA": 21.2, "GCG": 30.1,
    "GAT": 32.7, "GAC": 19.2, "GAA": 39.1, "GAG": 18.7,
    "GGT": 23.7, "GGC": 25.9, "GGA": 9.9, "GGG": 11.0,
}  # fmt: skip

assert len(SENSE_USAGE) == 61 and not STOPS & SENSE_USAGE.keys()
