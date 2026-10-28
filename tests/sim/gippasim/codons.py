"""Genetic codes 11, 4 and 25 (M11-16) and a bundled E. coli-like codon usage table.

Codes 4 and 25 differ from 11 only at TGA (Trp at 4, Gly at 25). ``sense_usage`` moves
``tga_sense_share`` of that amino acid's per-mille weight onto TGA and scales the other codons
of the family by ``1 - share``, so every amino acid keeps its table-11 weight.
"""

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

#: The tables the simulator supports, TGA's residue in each and the default ``tga_sense_share``
#: (placeholders, spec Q-L1).
TABLES = (4, 11, 25)
_TGA = {4: "W", 25: "G"}
DEFAULT_TGA_SHARE = {4: 0.6, 25: 0.05}


def _check(table: int) -> None:
    if table not in TABLES:
        raise ValueError(f"genetic code {table} is not supported; use one of {TABLES}")


def code(table: int) -> dict[str, str]:
    """Translation of every codon under NCBI ``table`` (``*`` for stops)."""
    _check(table)
    return CODE_11 if table == 11 else {**CODE_11, "TGA": _TGA[table]}


def stops(table: int) -> frozenset[str]:
    """Stop codons of NCBI ``table``."""
    _check(table)
    return STOPS if table == 11 else STOPS - {"TGA"}


def sense_usage(table: int, tga_sense_share: float | None = None) -> dict[str, float]:
    """Per-mille sense-codon usage under ``table``; ``tga_sense_share`` (``None`` = the table's
    default) is TGA's share of its amino acid at 4/25 and is ignored at 11."""
    _check(table)
    if table == 11:
        return SENSE_USAGE
    share = DEFAULT_TGA_SHARE[table] if tga_sense_share is None else tga_sense_share
    if not 0 < share < 1:
        raise ValueError(f"tga_sense_share must be in (0, 1), got {share}")
    aa = _TGA[table]
    fam = [c for c in SENSE_USAGE if CODE_11[c] == aa]
    out = {c: (1 - share) * w if c in fam else w for c, w in SENSE_USAGE.items()}
    out["TGA"] = share * sum(SENSE_USAGE[c] for c in fam)
    return out


_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def translate_cds(codons: list[str], table: int = 11) -> str:
    """Protein of a CDS given as codons (start ... stop): initiator read as M, stop dropped."""
    t = code(table)
    return "M" + "".join(t[c] for c in codons[1:-1])


def codons_of(nt: str) -> list[str]:
    return [nt[i : i + 3] for i in range(0, len(nt) - len(nt) % 3, 3)]
