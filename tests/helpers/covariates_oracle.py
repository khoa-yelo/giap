"""Independent oracle for the M4-10 protein covariates (docs/specs/M4-10-protein-covariates.md).

Written from the spec text, not from ``giap.genetable.covariates``: digestion by regular expression,
Kyte-Doolittle sliding means by cumulative sums (numpy), translation with a literal NCBI table 11,
and the GRAVY/pI/MW reference straight from Biopython ``ProteinAnalysis`` on the cleaned sequence.
"""

from __future__ import annotations

import math
import re

import numpy as np

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte & Doolittle (1982) hydropathy index; any other letter scores 0 (spec).
KD: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}  # fmt: skip

_BASES = "TCAG"
_AA11 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
TABLE11: dict[str, str] = {
    a + b + c: _AA11[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
STOPS = frozenset(c for c, aa in TABLE11.items() if aa == "*")


def digest(aa: str, enzyme: str, missed: int = 0) -> list[str]:
    """Fragments in sequence order (missed = 0), then concatenations of 2..missed+1 fragments."""
    pattern = r"(?<=[KR])" if enzyme == "trypsin_p" else r"(?<=[KR])(?!P)"
    frags = [f for f in re.split(pattern, aa) if f]
    out = list(frags)
    for span in range(2, missed + 2):
        out += ["".join(frags[i : i + span]) for i in range(len(frags) - span + 1)]
    return out


def n_observable(aa: str, enzyme: str, length: tuple[int, int] = (7, 30)) -> int:
    lo, hi = length
    return len({p for p in digest(aa, enzyme) if lo <= len(p) <= hi})


def tm_segments(
    aa: str, window: int = 19, threshold: float = 1.6, cap: int = 10, eps: float = 1e-9
) -> int:
    """Maximal runs of consecutive windows whose Kyte-Doolittle mean is >= threshold."""
    if len(aa) < window:
        return 0
    v = np.array([KD.get(c, 0.0) for c in aa])
    c = np.concatenate([[0.0], np.cumsum(v)])
    means = (c[window:] - c[:-window]) / window
    hit = means >= threshold - eps
    starts = int(hit[0]) + int(np.sum(hit[1:] & ~hit[:-1]))
    return min(cap, starts)


def clean(aa: str) -> str:
    """Strip a leading/trailing ``*`` and drop letters outside the 20 standard amino acids."""
    return "".join(ch for ch in aa.strip("*") if ch in STANDARD_AA)


def biopython_props(aa: str) -> tuple[float, float, float]:
    """(gravy, pi, mw_da) from Biopython on the cleaned sequence; NaN when empty."""
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    s = clean(aa)
    if not s:
        return (math.nan, math.nan, math.nan)
    p = ProteinAnalysis(s)
    return (p.gravy(), p.isoelectric_point(), p.molecular_weight())


def translate_pseudogene(nt: str) -> str:
    """Table 11, whole codons only, first codon read as M, stop before the first stop codon."""
    out = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        codon = nt[i : i + 3]
        if i == 0:
            out.append("M")
            continue
        if codon in STOPS:
            break
        out.append(TABLE11.get(codon, "X"))
    return "".join(out)


# One codon per amino acid for back-translation of test fixtures.
BACK: dict[str, str] = {}
for _codon, _aa in TABLE11.items():
    BACK.setdefault(_aa, _codon)
