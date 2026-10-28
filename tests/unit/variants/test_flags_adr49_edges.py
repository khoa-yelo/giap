"""ADR-0049 A6a edge cases of ``gippa.variants.flags.left_block_site`` and ``left_block_read``
(review of PR #269, items 2a-2d and the two guard nits; docs/adr/ADR-0049-m7-variant-precision.md).

``left_block_site`` reads the reference over the site codon and the ``q + 1`` codons before it
(transcription order) and asks :func:`left_block_read` whether the ALT peptide also aligns as the
left 2->1 block. These tests pin the window geometry where it is easy to break: across the origin
of a circular replicon (both strands, every rotation, genome lengths 0, 1 and 2 mod 3), at the
ends of a linear contig, on the minus strand with ``q >= 1``, with a second substitution in the
peptide, and with alignments that do not belong to the site (other contig, strand or op kind).

Alignments come from ``gippa.variants.align.align_hit`` (0-based half-open hits; a hit across the
origin has ``end > length``), except the op-kind case, which builds its ops by hand.
"""

from __future__ import annotations

from typing import Any

import pytest

from gippa.config.model import GippaConfig
from gippa.model import Replicon
from gippa.variants.align import HitAlignment, ResidueOp, align_hit
from gippa.variants.flags import left_block_read, left_block_site, mass_equivalent

TOL = GippaConfig().variants.mass_equivalent_tol_da
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGT", "TGCA")

# (protein, ALT peptide, protein index of the peptide's first residue, protein index of the site)
NTERM = ("WAGSVATEPIK", "QSVATEPIK", 2, 2)  # q = 0, AG -> Q
INTERNAL = ("DKAAGSVATAEPIKR", "AQSVATAEPIK", 3, 4)  # q = 1, A + (AG -> Q)
LOOSE = ("DKAGTWNHYFR", "AQTWNHYF", 2, 3)  # q = 1, prefix A but ref[k-2] = K: not a block


def _nt(prot: str) -> str:
    return "".join(CODON[a] for a in prot)


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _site(contig: str, strand: str, codon: int, alt_type: str = "substitution") -> dict[str, Any]:
    return {"alt_type": alt_type, "contig": contig, "strand": strand, "codon_start": codon}


def _place(
    case: tuple[str, str, int, int], strand: str, topo: str, rot: int, pad: int
) -> tuple[dict[str, Any], HitAlignment, str, Replicon]:
    """The protein of ``case`` on ``strand`` with ``pad`` filler bases each side, the genome
    rotated left by ``rot`` (circular: coordinate x moves to (x - rot) mod L). Returns the site
    row, the hit's alignment, the genome and the replicon."""
    prot, pep, first, k = case
    nt = _nt(prot)
    seq0 = "C" * pad + (nt if strand == "+" else _rc(nt)) + "C" * pad
    n = len(seq0)
    rot %= n
    seq = seq0[rot:] + seq0[:rot]
    rep = Replicon("c", n, topo)  # type: ignore[arg-type]
    if strand == "+":
        s0, c0 = pad + 3 * first, pad + 3 * k
    else:
        s0 = pad + len(nt) - 3 * (first + len(pep))
        c0 = pad + len(nt) - 3 * (k + 1)
    s = (s0 - rot) % n
    a = align_hit(pep, "c", s, s + 3 * len(pep), strand, 1, None, seq, rep)
    return _site("c", strand, (c0 - rot) % n), a, seq, rep


def _flag(site: dict[str, Any], aligns: list[HitAlignment], seq: str, rep: Replicon) -> bool:
    return left_block_site(site, aligns, {rep.id: seq}, {rep.id: rep}, TOL)


# --------------------------------------------------------------------------- a. circular wrap

PADS = (0, 1, 2)  # genome length 45 + 2 * pad: 0, 2 and 1 mod 3


@pytest.mark.parametrize("pad", PADS)
@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("case", [NTERM, INTERNAL], ids=["q0", "q1"])
def test_block_read_is_flagged_at_every_rotation_of_a_circular_origin(
    case: tuple[str, str, int, int], strand: str, pad: int
) -> None:
    """Every rotation puts the origin somewhere in (or next to) the hit and its window: the site
    codon at 0 with the residue before it at ``L - 3`` (``+``), the window end past ``L``
    (``-``), and hits whose op codons run past ``L``."""
    n = 3 * len(case[0]) + 2 * pad
    for rot in range(n):
        site, a, seq, rep = _place(case, strand, "circular", rot, pad)
        assert _flag(site, [a], seq, rep) is True, (strand, pad, rot, site)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_loose_prefix_is_not_flagged_at_any_rotation(strand: str) -> None:
    n = 3 * len(LOOSE[0])
    for rot in range(n):
        site, a, seq, rep = _place(LOOSE, strand, "circular", rot, 0)
        assert _flag(site, [a], seq, rep) is False, (strand, rot, site)


def test_plus_site_at_codon_zero_reads_the_residue_before_the_origin() -> None:
    """``+``: the site G at codon 0, its A at ``L - 3``; the hit (q = 1) starts at ``L - 3`` and
    crosses the origin, so its substitution op sits at codon ``L``."""
    k = INTERNAL[3]
    rot = 3 * k  # G codon to coordinate 0
    site, a, seq, rep = _place(INTERNAL, "+", "circular", rot, 0)
    assert site["codon_start"] == 0
    sub = [o for o in a.ops if o.kind == "substitution"]
    assert [o.codons[0] for o in sub] == [rep.length], a
    assert seq[rep.length - 3 :] == CODON["A"]
    assert _flag(site, [a], seq, rep) is True


def test_plus_site_at_codon_zero_q0_window_starts_before_the_origin() -> None:
    _, _, _, k = NTERM
    site, a, seq, rep = _place(NTERM, "+", "circular", 3 * k, 0)
    assert site["codon_start"] == 0
    assert _flag(site, [a], seq, rep) is True


def test_minus_site_at_the_last_codon_reads_past_the_origin() -> None:
    """``-``: the site G at ``L - 3``; the residue before it (transcription order) at ``[0, 3)``,
    so the window ``[L - 3, L + 3)`` crosses the origin."""
    prot, _, _, k = NTERM
    n = 3 * len(prot)
    c0 = n - 3 * (k + 1)
    rot = (c0 - (n - 3)) % n  # site codon to L - 3
    site, a, seq, rep = _place(NTERM, "-", "circular", rot, 0)
    assert site["codon_start"] == n - 3
    assert _flag(site, [a], seq, rep) is True


def test_window_longer_than_a_tiny_circular_replicon_is_not_flagged() -> None:
    """A 4-codon circular replicon and a 4-residue hit with its substitution at q = 3: the
    ``q + 2`` = 5-codon window is longer than the replicon. No flag and no ValueError."""
    seq = _nt("WAGS")
    rep = Replicon("t", len(seq), "circular")
    a = align_hit("WAGQ", "t", 0, 12, "+", 1, None, seq, rep)
    assert _flag(_site("t", "+", 9), [a], seq, rep) is False


# --------------------------------------------------------------------------- b. linear edges


def test_plus_site_at_the_first_codon_of_a_linear_contig_is_not_flagged() -> None:
    """No residue before codon 0 on a linear contig: False, no ValueError."""
    prot, pep = "GSVATEPIK", "QSVATEPIK"
    seq = _nt(prot)
    rep = Replicon("c", len(seq), "linear")
    a = align_hit(pep, "c", 0, len(seq), "+", 1, None, seq, rep)
    assert _flag(_site("c", "+", 0), [a], seq, rep) is False


def test_minus_site_at_the_last_codon_of_a_linear_contig_is_not_flagged() -> None:
    prot, pep = "GSVATEPIK", "QSVATEPIK"
    seq = _rc(_nt(prot))
    n = len(seq)
    rep = Replicon("c", n, "linear")
    a = align_hit(pep, "c", 0, n, "-", 1, None, seq, rep)
    assert _flag(_site("c", "-", n - 3), [a], seq, rep) is False


@pytest.mark.parametrize("strand", ["+", "-"])
def test_site_one_codon_in_from_a_linear_edge_is_flagged(strand: str) -> None:
    """The residue before the site is the contig's first (``+``) or last (``-``) codon."""
    site, a, seq, rep = _place(("AGSVATEPIK", "QSVATEPIK", 1, 1), strand, "linear", 0, 0)
    assert _flag(site, [a], seq, rep) is True


# --------------------------------------------------------------------------- c. minus, q >= 1


@pytest.mark.parametrize("pad", [0, 30])
def test_minus_strand_internal_read_is_flagged(pad: int) -> None:
    """q = 1 on ``-``: the window is ``q + 2`` = 3 codons upward from the site codon."""
    site, a, seq, rep = _place(INTERNAL, "-", "linear", 0, pad)
    assert _flag(site, [a], seq, rep) is True


def test_minus_strand_loose_prefix_is_not_flagged() -> None:
    site, a, seq, rep = _place(LOOSE, "-", "linear", 0, 30)
    assert _flag(site, [a], seq, rep) is False


def test_minus_strand_longer_shifted_prefix_is_flagged() -> None:
    # q = 2: ...AAAG S..., peptide AAQSWK, as in the + example of test_flags_adr49
    site, a, seq, rep = _place(("DAAAGSWK", "AAQSWK", 2, 4), "-", "linear", 0, 30)
    assert _flag(site, [a], seq, rep) is True


# --------------------------------------------------------------------------- d. op codon match

TWO_SUBS = ("DWAGSVAGTEPIK", "QSVAWTEPIK", 3, 3)  # G3 -> Q (block-readable), G7 -> W (not)


def _two_subs(strand: str) -> tuple[HitAlignment, str, Replicon, int, int]:
    prot, pep, first, _ = TWO_SUBS
    nt = _nt(prot)
    pad = 30
    seq = "C" * pad + (nt if strand == "+" else _rc(nt)) + "C" * pad
    rep = Replicon("c", len(seq), "linear")
    if strand == "+":
        s = pad + 3 * first
        cod = [pad + 3 * i for i in (3, 7)]
    else:
        s = pad + len(nt) - 3 * (first + len(pep))
        cod = [pad + len(nt) - 3 * (i + 1) for i in (3, 7)]
    a = align_hit(pep, "c", s, s + 3 * len(pep), strand, 2, None, seq, rep)
    return a, seq, rep, cod[0], cod[1]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_second_substitution_does_not_borrow_the_block_read(strand: str) -> None:
    """The peptide reads G>Q at q = 0 (``AG -> Q``, a block) and G>W at q = 4 (``AG -> W`` is no
    equal-mass reading). Only the G>Q site is flagged: each op is read at its own codon."""
    a, seq, rep, c_q, c_w = _two_subs(strand)
    subs = sorted((o.pep[0], o.alt) for o in a.ops if o.kind == "substitution")
    assert subs == [(0, "Q"), (4, "W")], a
    assert not mass_equivalent("AG", "W", TOL)
    assert _flag(_site("c", strand, c_q), [a], seq, rep) is True
    assert _flag(_site("c", strand, c_w), [a], seq, rep) is False


# --------------------------------------------------------------------------- alignment filters


def test_alignment_on_another_contig_is_ignored() -> None:
    """Two contigs with the same sequence: the c2 alignment does not flag the c1 site."""
    site, a, seq, rep = _place(NTERM, "+", "linear", 0, 30)
    rep2 = Replicon("c2", rep.length, "linear")
    a2 = align_hit(a.peptide, "c2", 30 + 3 * NTERM[2], 30 + 3 * NTERM[2] + 3 * len(a.peptide),
                   "+", 1, None, seq, rep2)  # fmt: skip
    genome = {"c": seq, "c2": seq}
    reps = {"c": rep, "c2": rep2}
    assert left_block_site(site, [a], genome, reps, TOL) is True
    assert left_block_site(site, [a2], genome, reps, TOL) is False


def test_alignment_on_the_other_strand_is_ignored() -> None:
    """``+`` reads ...A G T... (``GCT GGC ACC``) as AG -> Q at the site codon; ``-`` at the same
    codon reads (transcription order) G then A, and ``GA`` weighs as Q too. A ``+`` alignment
    must not flag the ``-`` site at that codon."""
    pad = 30
    nt = "GCT" + "GGC" + "ACC" + _nt("SVATEPIK")
    seq = "C" * pad + nt + "C" * pad
    rep = Replicon("c", len(seq), "linear")
    c = pad + 3
    a = align_hit("QTSVATEPIK", "c", c, c + 30, "+", 1, None, seq, rep)
    assert _flag(_site("c", "+", c), [a], seq, rep) is True
    assert mass_equivalent("GA", "Q", TOL)
    assert _flag(_site("c", "-", c), [a], seq, rep) is False


@pytest.mark.parametrize("kind", ["match", "il", "swap", "block", "near_block"])
def test_only_substitution_ops_are_read(kind: str) -> None:
    """Hand-built op of another kind at the site codon, otherwise a perfect AG -> Q read."""
    site, a, seq, rep = _place(NTERM, "+", "linear", 0, 30)
    sub = next(o for o in a.ops if o.kind == "substitution")
    other = ResidueOp(kind, sub.pep, sub.codons, sub.ref, sub.alt)  # type: ignore[arg-type]
    ops = tuple(other if o is sub else o for o in a.ops)
    b = HitAlignment(a.peptide, a.contig, a.strand, a.frame, ops, a.frameshift)
    assert _flag(site, [a], seq, rep) is True
    assert _flag(site, [b], seq, rep) is False


@pytest.mark.parametrize("alt_type", ["swap", "block", "insertion", "deletion", "frameshift"])
def test_only_substitution_sites_are_read(alt_type: str) -> None:
    site, a, seq, rep = _place(NTERM, "+", "linear", 0, 30)
    assert _flag(site, [a], seq, rep) is True
    assert _flag({**site, "alt_type": alt_type}, [a], seq, rep) is False


# --------------------------------------------------------------------------- guards


def test_site_at_reference_index_zero_never_reads_a_block() -> None:
    """``k < 1`` guard: with a one-residue ``ref``, ``ref[k-1:k+1]`` would wrap to ``ref[0:1]``
    and N -> D (deamidation) would pass the mass test."""
    assert mass_equivalent("N", "D", TOL)
    assert left_block_read("D", 0, "N", 0, TOL) is False
    assert left_block_read("DSW", 0, "N", 0, TOL) is False


def test_shifted_prefix_longer_than_the_reference_before_the_pair_is_rejected() -> None:
    """``lo >= 0`` guard: q = 1 with only one residue before the site."""
    assert left_block_read("AQ", 1, "AG", 1, TOL) is False
    assert left_block_read("AAAAQ", 4, "XAG", 2, TOL) is False
