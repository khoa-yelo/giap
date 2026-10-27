"""Genome builders for the frameshift tests (issue M3-08); test tooling, not the simulator."""

from __future__ import annotations

import numpy as np

from gippasim.codons import CODE_11

CODONS: dict[str, list[str]] = {}
for _c, _a in sorted(CODE_11.items()):
    CODONS.setdefault(_a, []).append(_c)
SENSE_RESIDUES = "".join(sorted(a for a in CODONS if a != "*"))  # 20 letters, I included
_RC = str.maketrans("ACGTN", "TGCAN")


def rc(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def background(n: int, seed: int, n_rate: float = 0.0) -> str:
    """Random ACGT of length ``n`` (with ``N`` at rate ``n_rate``)."""
    rng = np.random.default_rng(seed)
    out = rng.choice(list("ACGT"), size=n)
    if n_rate:
        out[rng.random(n) < n_rate] = "N"
    return "".join(out)


def codons_for(peptide: str, rng: np.random.Generator | None = None) -> list[str]:
    """One sense codon per residue (the first one, or a random one with ``rng``)."""
    out = []
    for a in peptide:
        opts = CODONS[a]
        out.append(opts[int(rng.integers(len(opts)))] if rng is not None else opts[0])
    return out


def join_split(prefix: list[str], suffix: list[str], delta: int, extra: str = "A") -> str | None:
    """Nucleotides of a ``delta`` split: ``+1`` inserts ``extra``; ``-1`` reads one base twice.

    For ``-1`` the last base of the last prefix codon must equal the first base of the first
    suffix codon; the codon of either side is re-chosen among synonymous ones when needed.
    Returns ``None`` when no synonymous pair fits.
    """
    if delta == 1:
        return "".join(prefix) + extra + "".join(suffix)
    a, b = prefix[-1], suffix[0]
    for ca in [a, *CODONS[CODE_11[a]]]:
        for cb in [b, *CODONS[CODE_11[b]]]:
            if ca[2] == cb[0]:
                return "".join(prefix[:-1]) + ca + cb[1:] + "".join(suffix[1:])
    return None
