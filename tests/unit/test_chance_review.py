"""Lead decision 1 after statistician review (PR #28; docs/specs/M3-09-chance-hit-model.md).

On a circular replicon N counts the minus-strand codons that cross the origin: for a stop-free,
X-free circular replicon of length n the searchable six-frame positions total exactly 2 * n
(one codon start per nucleotide per strand), not 2 * n - 2.

Stop-freeness by construction: every stop codon (TAA, TAG, TGA) starts with T, and a minus-strand
stop reads TTA, CTA or TCA on the plus strand, all of which contain T. A sequence over {A, C, G}
therefore has no stop in any frame on either strand, across the origin included.
"""

from __future__ import annotations

import pytest

from gippa.mapping.chance import residue_freqs
from gippa.mapping.sixframe import ALPHABET, STOP, SixFrameIndex, build_sixframe
from gippa.model.coords import Replicon

X_CODE = ALPHABET.index("X")
CODONS = ("GCA", "ACG", "CAG", "GAC", "AGC", "CCG", "GGA", "AAC")


def _stop_free_seq(n: int) -> str:
    seq = "".join(CODONS[i % len(CODONS)] for i in range(n // 3 + 1))[:n]
    assert len(seq) == n and "T" not in seq
    return seq


def _circular_index(n: int) -> SixFrameIndex:
    seq = _stop_free_seq(n)
    idx = build_sixframe({"c": seq}, {"c": Replicon("c", n, "circular")})
    for frame in idx.frames:  # extension included: covers the origin-crossing codons
        assert not (frame.aa == STOP).any(), (frame.strand, frame.offset)
        assert not (frame.aa == X_CODE).any(), (frame.strand, frame.offset)
    return idx


@pytest.mark.parametrize("n", [30, 31, 32, 300, 301, 302])
def test_stop_free_circular_n_is_exactly_two_n(n: int) -> None:
    fr = residue_freqs(_circular_index(n))
    assert fr.n_positions == 2 * n


@pytest.mark.parametrize("n", [30, 31, 32, 301])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_each_strand_counts_one_codon_per_nucleotide(n: int, strand: str) -> None:
    """The minus strand keeps its origin-crossing codons (first_nt - 2 < length)."""
    idx = _circular_index(n)
    one = SixFrameIndex(
        replicons=idx.replicons, frames=tuple(f for f in idx.frames if f.strand == strand)
    )
    assert residue_freqs(one).n_positions == n
