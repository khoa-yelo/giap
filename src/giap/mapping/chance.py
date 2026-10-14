"""Chance-hit model E(L, m) and the per-peptide ranking score E_i (GIAP_SPEC.md §7.2.5, G3,
Appendix B.2).

E(L, m) = N * sum_{j <= m} C(L, j) (1 - s)^j s^(L - j), with s = sum f_a^2 and N the number of
six-frame positions. E_i = N * P(<= m mismatches) where position i mismatches with probability
1 - f_{a_i} (Poisson-binomial). E_i is a ranking score: the decoy search calibrates the FDR.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from giap.mapping.sixframe import STOP, Frame, SixFrameIndex, encode_query
from giap.stats.poisson_binomial import poisson_binomial_cdf

_N_RES = STOP  # residue codes 0..18 (I folded into L)


@dataclass(frozen=True)
class ResidueFreqs:
    """Six-frame residue frequencies (19 codes, sum 1, no ``*``/``X``) and position count N."""

    f: np.ndarray
    n_positions: int

    @property
    def s(self) -> float:
        """Per-position match probability of two random residues, sum f_a^2."""
        return float(np.sum(self.f**2))


def _core_residues(frame: Frame, length: int) -> np.ndarray:
    """Residues of ``frame`` that are distinct searchable positions of the replicon.

    A plus-strand codon is kept when its first nucleotide lies in the replicon. A minus-strand codon
    is kept when its lowest coordinate (``first_nt - 2``) does, so codons crossing the origin of a
    circular replicon count once: a stop-free circular replicon of length n gives exactly 2n.
    """
    k = np.arange(len(frame.aa))
    if frame.strand == "+":
        keep = frame.offset + 3 * k < length
    else:
        first_nt = frame.ext_nt - frame.offset - 3 * k - 1
        keep = first_nt - 2 < length
    out: np.ndarray = frame.aa[keep]
    return out


def residue_freqs(index: SixFrameIndex) -> ResidueFreqs:
    """Residue frequencies and N over the target frames, circular-extension duplicates excluded."""
    counts = np.zeros(_N_RES, dtype=np.int64)
    n = 0
    for frame in index.frames:
        core = _core_residues(frame, index.replicons[frame.contig].length)
        n += int(core.size)
        counts += np.bincount(core[core < _N_RES], minlength=_N_RES)
    total = counts.sum()
    f = counts / total if total else np.zeros(_N_RES)
    return ResidueFreqs(f=f.astype(np.float64), n_positions=n)


def _check_m(L: int, m: int) -> None:
    if not 0 <= m <= L:
        raise ValueError(f"mismatch budget m={m} must lie in [0, L={L}]")


def expected_hits_uniform(L: int, m: int, s: float, n_positions: float) -> float:
    """E(L, m): expected chance windows within <= m substitutions of a random L-mer."""
    _check_m(L, m)
    tail = math.fsum(math.comb(L, j) * (1.0 - s) ** j * s ** (L - j) for j in range(m + 1))
    return float(n_positions * tail)


def expected_hits(query: str, m: int, freqs: ResidueFreqs) -> float:
    """E_i: N times the Poisson-binomial P(<= m mismatches), match probability f_{a_i}."""
    codes = encode_query(query)
    _check_m(len(codes), m)
    f_ext = np.zeros(256, dtype=np.float64)
    f_ext[:_N_RES] = freqs.f
    mismatch = 1.0 - f_ext[codes]
    return float(freqs.n_positions * poisson_binomial_cdf(mismatch, m))
