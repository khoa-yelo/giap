"""M7-01 residue alignment of every accepted placement (docs/specs/M7-01-residue-alignment.md;
GIPPA_SPEC.md §7.5.2 step 1, §7.2.4, §7.2.1, §6.5).

Examples 1-10 and every degenerate row of the spec, on hand-built genomes: contig ``c1``,
circular, 1,000 nt, a fixed random background with the coding sequence of each example written
at its locus (first table-11 codon of each residue, ``helpers.fs_genomes.codons_for``).

Hand-worked oracles pin:

* the leftmost placement of co-optimal differences (lead decision 2, test-author reading: the
  lexicographically smallest tuple of the differences' peptide start positions, in peptide
  order, among alignments of equal cost and diff-type multiset);
* I/L: ``il`` ops (REF as the genome encodes it, ALT as the peptide gives it), I folded into L
  for the cost only (a swap ``KI``/``LK`` and a block ``GI``/``AV``);
* the recomputed cost equals ``n_diff`` and a mismatch raises;
* minus-strand and origin-wrap codon coordinates (0-based, lowest coordinate, unwrapped);
* determinism and no mutation of the inputs.

``HitAlignment.frame`` is pinned to the ``peptide_hits.tsv`` convention (M3-02
``SixFrameIndex.frame_label``: ``1 + start % 3`` on ``+``, ``-(1 + ((length - end) % length) %
3)`` on ``-``; the prefix block of a frameshift split).
"""

from __future__ import annotations

import copy
import importlib
import re
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers.align_oracle import (
    codon_at,
    expected_ops,
    locus_aa,
    op_tuples,
    place,
    rc,
)
from helpers.fs_genomes import background, codons_for, join_split

from gippa.mapping.hits import hits_frame
from gippa.mapping.outputs import PEPTIDE_HITS_COLUMNS
from gippa.model import Replicon

L = 1000
REP = Replicon("c1", L, "circular")
LIN = Replicon("c1", L, "linear")
BG = background(L, seed=7)


def align() -> ModuleType:
    """The module under test (imported per test, so each test fails on its own)."""
    return importlib.import_module("gippa.variants.align")


def nt(protein: str) -> str:
    return "".join(codons_for(protein))


def genome(*parts: tuple[int, str, str], base: str = BG) -> str:
    """``base`` with each ``(start, protein, strand)`` written at ``start`` (wrapping)."""
    seq = base
    for start, protein, strand in parts:
        seq = place(seq, start, nt(protein), strand)
    return seq


def m(i: int, codon: int, aa: str) -> tuple[str, tuple[int, int], tuple[int, ...], str, str]:
    return ("match", (i, i + 1), (codon,), aa, aa)


def matches(
    peptide: str, start: int, end: int, strand: str, skip: frozenset[int] | set[int] = frozenset()
) -> list[tuple[str, tuple[int, int], tuple[int, ...], str, str]]:
    return [
        m(i, codon_at(i, start, end, strand), a) for i, a in enumerate(peptide) if i not in skip
    ]


# --------------------------------------------------------------------------- translate_locus


def test_translate_locus_plus_minus_and_no_il_folding() -> None:
    seq = genome((100, "MKILVE", "+"), (400, "WDIRLK", "-"))
    assert align().translate_locus(seq, REP, 100, 118, "+") == "MKILVE"
    assert align().translate_locus(seq, REP, 400, 418, "-") == "WDIRLK"


def test_translate_locus_reads_across_the_origin() -> None:
    seq = genome((995, "MKVLSDEK", "+"))
    assert seq[995:] + seq[:19] == nt("MKVLSDEK")
    assert align().translate_locus(seq, REP, 995, 1019, "+") == "MKVLSDEK"
    seq_m = genome((995, "PEPTIDEK", "-"))
    assert align().translate_locus(seq_m, REP, 995, 1019, "-") == "PEPTIDEK"


def test_translate_locus_ambiguous_codon_is_x() -> None:
    seq = genome((100, "MKVLSDEK", "+"))
    seq = seq[:106] + "N" + seq[107:]
    assert align().translate_locus(seq, REP, 100, 124, "+") == "MKXLSDEK"


def test_translate_locus_linear_beyond_end_raises() -> None:
    seq = genome((976, "MKVLSDEK", "+"))
    assert align().translate_locus(seq, LIN, 976, 1000, "+") == "MKVLSDEK"
    with pytest.raises(ValueError):
        align().translate_locus(seq, LIN, 995, 1019, "+")


# --------------------------------------------------------------------------- examples


def test_example1_exact_plus_with_il() -> None:
    seq = genome((100, "PEPTLDEK", "+"))
    a = align().align_hit("PEPTIDEK", "c1", 100, 124, "+", 0, None, seq, REP)
    want = matches("PEPTIDEK", 100, 124, "+", {4})
    want.insert(4, ("il", (4, 5), (112,), "L", "I"))
    assert op_tuples(a) == want
    assert [op.codons[0] for op in a.ops] == list(range(100, 124, 3))
    assert (a.peptide, a.contig, a.strand, a.frame, a.frameshift) == (
        "PEPTIDEK",
        "c1",
        "+",
        2,
        None,
    )


def test_example2_exact_minus_codons_run_from_end() -> None:
    seq = genome((200, "MKVLSDEK", "-"))
    a = align().align_hit("MKVLSDEK", "c1", 200, 224, "-", 0, None, seq, REP)
    assert a.ops[0].codons == (221,)
    assert a.ops[7].codons == (200,)
    assert op_tuples(a) == matches("MKVLSDEK", 200, 224, "-")
    assert [op.codons[0] for op in a.ops] == [221, 218, 215, 212, 209, 206, 203, 200]
    assert (a.strand, a.frame) == ("-", -3)


def test_example3_substitution() -> None:
    seq = genome((300, "AGSNEFKLLR", "+"))
    a = align().align_hit("AGSDEFKLLR", "c1", 300, 330, "+", 1, None, seq, REP)
    want = matches("AGSDEFKLLR", 300, 330, "+", {3})
    want.insert(3, ("substitution", (3, 4), (309,), "N", "D"))
    assert op_tuples(a) == want


def test_example4_swap() -> None:
    seq = genome((300, "AGSDEFKLLR", "+"))
    a = align().align_hit("AGDSEFKLLR", "c1", 300, 330, "+", 1, None, seq, REP)
    want = [m(0, 300, "A"), m(1, 303, "G"), ("swap", (2, 4), (306, 309), "SD", "DS")]
    want += [m(i, 300 + 3 * i, aa) for i, aa in enumerate("AGDSEFKLLR") if i >= 4]
    assert op_tuples(a) == want


def test_example5_block_two_residues_on_one_codon() -> None:
    seq = genome((300, "ANGEFKLLR", "+"))
    a = align().align_hit("AGGGEFKLLR", "c1", 300, 327, "+", 1, None, seq, REP)
    want = [m(0, 300, "A"), ("block", (1, 3), (303,), "N", "GG")]
    want += [m(i, 300 + 3 * (i - 1), aa) for i, aa in enumerate("AGGGEFKLLR") if i >= 3]
    assert op_tuples(a) == want
    blk = a.ops[1]
    assert blk.codons[0] == 303  # "the first codon of the block" (§7.5.2)


def test_example6_block_one_residue_on_two_codons() -> None:
    seq = genome((300, "AGGGEFKLLRK", "+"))
    a = align().align_hit("ANGEFKLLRK", "c1", 300, 333, "+", 1, None, seq, REP)
    want = [m(0, 300, "A"), ("block", (1, 2), (303, 306), "GG", "N")]
    want += [m(i, 300 + 3 * (i + 1), aa) for i, aa in enumerate("ANGEFKLLRK") if i >= 2]
    assert op_tuples(a) == want


def test_example7_near_block() -> None:
    seq = genome((300, "ADEFKLLRAGK", "+"))
    a = align().align_hit("WEFKLLRAGK", "c1", 300, 333, "+", 1, None, seq, REP)
    want = [("near_block", (0, 1), (300, 303), "AD", "W")]
    want += [m(i, 300 + 3 * (i + 1), aa) for i, aa in enumerate("WEFKLLRAGK") if i >= 1]
    assert op_tuples(a) == want


def fs_genome_plus() -> str:
    """Example 8: ``MKLLAE`` at [300, 318) frame 0, one extra base, ``GSTRK`` at [319, 334)."""
    body = join_split(codons_for("MKLLAE"), codons_for("GSTRK"), 1)
    assert body is not None and len(body) == 34
    return place(BG, 300, body, "+")


def test_example8_frameshift_plus() -> None:
    mod = align()
    seq = fs_genome_plus()
    a = mod.align_hit("MKLLAEGSTRK", "c1", 300, 334, "+", 1, "300-318;319-334", seq, REP)
    want = [m(i, 300 + 3 * i, aa) for i, aa in enumerate("MKLLAE")]
    want += [m(6 + i, 319 + 3 * i, aa) for i, aa in enumerate("GSTRK")]
    assert op_tuples(a) == want
    assert a.frameshift == mod.FrameshiftEvent(
        pep_junction=6, last_prefix_codon=315, first_suffix_codon=319, shift=1
    )
    assert a.frame == 1  # the prefix block's frame


def test_example9_exact_plus_across_the_origin() -> None:
    seq = genome((995, "MKVLSDEK", "+"))
    a = align().align_hit("MKVLSDEK", "c1", 995, 1019, "+", 0, None, seq, REP)
    assert a.ops[1].codons == (998,)  # spans [998, 1001)
    assert a.ops[2].codons == (1001,)  # unwrapped
    assert op_tuples(a) == matches("MKVLSDEK", 995, 1019, "+")


def test_example10_cost_above_n_diff_raises() -> None:
    seq = genome((300, "AGSNEFKLLR", "+"))
    with pytest.raises(ValueError, match="AGSDEFKLLW") as exc:
        align().align_hit("AGSDEFKLLW", "c1", 300, 330, "+", 1, None, seq, REP)
    assert re.search(r"n_diff 1, recomputed 2", str(exc.value))


# --------------------------------------------------------------------------- rule 4: cost and ties


def test_cost_below_n_diff_raises_too() -> None:
    seq = genome((300, "AGSNEFKLLRWK", "+"))
    with pytest.raises(ValueError, match="AGSDEFKLLRWK") as exc:
        align().align_hit("AGSDEFKLLRWK", "c1", 300, 336, "+", 2, None, seq, REP)
    assert re.search(r"n_diff 2, recomputed 1", str(exc.value))


def test_leftmost_placement_substitution_then_swap() -> None:
    """``GLG`` on ``LGL``: {substitution, swap} either as sub@4 + swap[5,7) or swap[4,6) + sub@6;
    leftmost puts the first difference at 4 and the next one at 5."""
    seq = genome((300, "MKDELGLFKPRW", "+"))
    a = align().align_hit("MKDEGLGFKPRW", "c1", 300, 336, "+", 2, None, seq, REP)
    diffs = [o for o in op_tuples(a) if o[0] not in ("match", "il")]
    assert diffs == [
        ("substitution", (4, 5), (312,), "L", "G"),
        ("swap", (5, 7), (315, 318), "GL", "LG"),
    ]
    assert op_tuples(a) == expected_ops("MKDEGLGFKPRW", seq, 300, 336, "+")


def test_leftmost_placement_substitution_then_block() -> None:
    """``AGGKW`` on ``QNKW``: sub A/Q + block GG/N, or block AG/Q + sub G/N; leftmost wins."""
    seq = genome((300, "MKDEQNKWFPRTH", "+"))
    a = align().align_hit("MKDEAGGKWFPRTH", "c1", 300, 339, "+", 2, None, seq, REP)
    diffs = [o for o in op_tuples(a) if o[0] not in ("match", "il")]
    assert diffs == [
        ("substitution", (4, 5), (312,), "Q", "A"),
        ("block", (5, 7), (315,), "N", "GG"),
    ]
    assert op_tuples(a) == expected_ops("MKDEAGGKWFPRTH", seq, 300, 339, "+")


def test_leftmost_placement_on_the_minus_strand_is_in_peptide_order() -> None:
    seq = genome((300, "MKDELGLFKPRW", "-"))
    a = align().align_hit("MKDEGLGFKPRW", "c1", 300, 336, "-", 2, None, seq, REP)
    diffs = [o for o in op_tuples(a) if o[0] not in ("match", "il")]
    assert diffs == [
        ("substitution", (4, 5), (336 - 15,), "L", "G"),
        ("swap", (5, 7), (336 - 18, 336 - 21), "GL", "LG"),
    ]


def test_il_is_free_for_the_cost_and_labelled() -> None:
    """A tolerant row: one substitution plus an I under a genome L (an ``il`` op, not a diff)."""
    seq = genome((300, "AGSNEFKLLR", "+"))
    a = align().align_hit("AGSDEFKILR", "c1", 300, 330, "+", 1, None, seq, REP)
    ops = op_tuples(a)
    assert ops[3] == ("substitution", (3, 4), (309,), "N", "D")
    assert ops[7] == ("il", (7, 8), (321,), "L", "I")
    assert [o[0] for o in ops].count("match") == 8


def test_il_genome_isoleucine_kept_as_ref() -> None:
    seq = genome((100, "PEPTIDEK", "+"))
    a = align().align_hit("PEPTLDEK", "c1", 100, 124, "+", 0, None, seq, REP)
    assert op_tuples(a)[4] == ("il", (4, 5), (112,), "I", "L")


def test_il_swap_of_i_and_l_is_two_il_ops_not_a_swap() -> None:
    seq = genome((100, "PEPTLIEK", "+"))
    a = align().align_hit("PEPTILEK", "c1", 100, 124, "+", 0, None, seq, REP)
    ops = op_tuples(a)
    assert ops[4] == ("il", (4, 5), (112,), "L", "I")
    assert ops[5] == ("il", (5, 6), (115,), "I", "L")


def test_il_folded_inside_a_swap() -> None:
    """``KI`` on genome ``LK`` is a swap once I is folded into L; REF/ALT keep the letters."""
    seq = genome((300, "MKLKEFTLLR", "+"))
    a = align().align_hit("MKKIEFTLLR", "c1", 300, 330, "+", 1, None, seq, REP)
    assert op_tuples(a)[2] == ("swap", (2, 4), (306, 309), "LK", "KI")
    assert len(a.ops) == 9


def test_il_folded_inside_a_block() -> None:
    """``GI`` on genome ``AV``: the block {A,V} <-> {G,L} with I folded into L."""
    seq = genome((300, "MKAVEFKLLR", "+"))
    a = align().align_hit("MKGIEFKLLR", "c1", 300, 330, "+", 1, None, seq, REP)
    assert op_tuples(a)[2] == ("block", (2, 4), (306, 309), "AV", "GI")


def test_minus_strand_tolerant_codons() -> None:
    seq = genome((300, "AGSNEFKLLR", "-"))
    a = align().align_hit("AGSDEFKLLR", "c1", 300, 330, "-", 1, None, seq, REP)
    want = matches("AGSDEFKLLR", 300, 330, "-", {3})
    want.insert(3, ("substitution", (3, 4), (330 - 12,), "N", "D"))
    assert op_tuples(a) == want


def test_minus_strand_block_codons_in_peptide_order() -> None:
    """Rule 5: ``codons`` lists the consumed codons in peptide order, so on ``-`` they descend."""
    seq = genome((300, "AGGGEFKLLRK", "-"))
    a = align().align_hit("ANGEFKLLRK", "c1", 300, 333, "-", 1, None, seq, REP)
    assert op_tuples(a)[1] == ("block", (1, 2), (333 - 6, 333 - 9), "GG", "N")


def test_minus_strand_across_the_origin() -> None:
    seq = genome((995, "MKVLSDEK", "-"))
    a = align().align_hit("MKVLSDEK", "c1", 995, 1019, "-", 0, None, seq, REP)
    assert [op.codons[0] for op in a.ops] == [1016, 1013, 1010, 1007, 1004, 1001, 998, 995]
    assert a.frame == -(1 + ((L - 1019) % L) % 3)


def test_tolerant_block_across_the_origin() -> None:
    """A block whose two codons straddle the origin: codons 998 and 1001 (unwrapped)."""
    seq = genome((995, "AGGGEFKLLRK", "+"))
    a = align().align_hit("ANGEFKLLRK", "c1", 995, 1028, "+", 1, None, seq, REP)
    assert op_tuples(a)[1] == ("block", (1, 2), (998, 1001), "GG", "N")


# --------------------------------------------------------------------------- degenerate rows


def test_exact_row_differing_beyond_il_raises() -> None:
    seq = genome((100, "PEPTVDEK", "+"))
    with pytest.raises(ValueError, match="PEPTIDEK"):
        align().align_hit("PEPTIDEK", "c1", 100, 124, "+", 0, None, seq, REP)


def test_ambiguous_genome_codon_under_a_residue_raises() -> None:
    seq = genome((100, "PEPTIDEK", "+"))
    seq_x = seq[:106] + "N" + seq[107:]
    with pytest.raises(ValueError, match="PEPTIDEK"):
        align().align_hit("PEPTIDEK", "c1", 100, 124, "+", 0, None, seq_x, REP)
    seq2 = genome((300, "AGSNEFKLLR", "+"))
    seq2_x = seq2[:315] + "N" + seq2[316:]  # codon of F (residue 5)
    with pytest.raises(ValueError, match="AGSDEFKLLR"):
        align().align_hit("AGSDEFKLLR", "c1", 300, 330, "+", 1, None, seq2_x, REP)


def test_linear_contig_hit_beyond_end_raises() -> None:
    seq = genome((995, "MKVLSDEK", "+"))
    with pytest.raises(ValueError):
        align().align_hit("MKVLSDEK", "c1", 995, 1019, "+", 0, None, seq, LIN)


# --------------------------------------------------------------------------- frameshift route


def test_frameshift_minus_mirrors_example8() -> None:
    """Example 8 on the reverse-complement contig: prefix = the higher-coordinate block."""
    mod = align()
    seq = rc(fs_genome_plus())
    a = mod.align_hit("MKLLAEGSTRK", "c1", 666, 700, "-", 1, "666-681;682-700", seq, REP)
    want = [m(i, 697 - 3 * i, aa) for i, aa in enumerate("MKLLAE")]
    want += [m(6 + i, 678 - 3 * i, aa) for i, aa in enumerate("GSTRK")]
    assert op_tuples(a) == want
    assert a.frameshift == mod.FrameshiftEvent(
        pep_junction=6, last_prefix_codon=682, first_suffix_codon=678, shift=1
    )
    assert a.frame == -(1 + ((L - 700) % L) % 3)


def test_frameshift_minus_one_overlap() -> None:
    mod = align()
    body = join_split(codons_for("MKLLAE"), codons_for("GSTRK"), -1)
    assert body is not None and len(body) == 32
    seq = place(BG, 300, body, "+")
    a = mod.align_hit("MKLLAEGSTRK", "c1", 300, 332, "+", 1, "300-318;317-332", seq, REP)
    assert [op.codons for op in a.ops] == [(300 + 3 * i,) for i in range(6)] + [
        (317 + 3 * i,) for i in range(5)
    ]
    assert a.frameshift == mod.FrameshiftEvent(6, 315, 317, -1)


def test_frameshift_with_a_substitution_on_the_suffix() -> None:
    mod = align()
    body = join_split(codons_for("MKLLAEG"), codons_for("STRKWV"), 1)
    assert body is not None
    seq = place(BG, 300, body, "+")
    blocks = "300-321;322-340"
    a = mod.align_hit("MKLLAEGSTDKWV", "c1", 300, 340, "+", 2, blocks, seq, REP)
    ops = op_tuples(a)
    assert ops[9] == ("substitution", (9, 10), (322 + 6,), "R", "D")
    assert [o[0] for o in ops].count("match") == 12
    assert a.frameshift == mod.FrameshiftEvent(7, 318, 322, 1)
    with pytest.raises(ValueError, match="MKLLAEGSTDKWV") as exc:
        mod.align_hit("MKLLAEGSTDKWV", "c1", 300, 340, "+", 1, blocks, seq, REP)
    assert re.search(r"n_diff 1, recomputed 2", str(exc.value))


# --------------------------------------------------------------------------- align_hits


def ph_row(peptide: str, start: int, end: int, strand: str, route: str, n_diff: int,
           diff_types: str | None, frame: int) -> dict[str, Any]:  # fmt: skip
    """One ``peptide_hits.tsv`` row, already converted to 0-based half-open (rule 1)."""
    row: dict[str, Any] = dict.fromkeys(PEPTIDE_HITS_COLUMNS)
    row.update(
        peptide=peptide, peptide_len=len(peptide), source="database", contig="c1", start=start,
        end=end, strand=strand, frame=frame, crosses_origin=end > L, route=route, n_diff=n_diff,
        diff_types=diff_types, placement_q=0.0, locus_class="cds_inframe", uniqueness="unique",
    )  # fmt: skip
    return row


def ht_row(peptide: str, start: int, end: int, strand: str, blocks: str | None,
           search: str = "target", n_diff: int = 1) -> dict[str, Any]:  # fmt: skip
    """One ``hits_tolerant.parquet`` row (``gippa.mapping.hits.HIT_COLUMNS``)."""
    return {
        "peptide": peptide, "query": peptide.replace("I", "L"), "search": search,
        "decoy_type": None if search == "target" else "reverse", "contig": "c1", "frame_id": 0,
        "aa_start": 0, "aa_end": 0, "start": start, "end": end, "strand": strand, "frame": 1,
        "route": "tolerant", "n_diff": n_diff, "diff_types": "frameshift", "blocks": blocks,
    }  # fmt: skip


def three_loci() -> tuple[str, pd.DataFrame, pd.DataFrame]:
    """Example 1 at 100, example 8 at 300, example 3 at 500 on one contig."""
    body = join_split(codons_for("MKLLAE"), codons_for("GSTRK"), 1)
    assert body is not None
    seq = place(genome((100, "PEPTLDEK", "+"), (500, "AGSNEFKLLR", "+")), 300, body, "+")
    ph = pd.DataFrame(
        [
            ph_row("PEPTIDEK", 100, 124, "+", "exact", 0, None, 2),
            ph_row("MKLLAEGSTRK", 300, 334, "+", "tolerant", 1, "frameshift", 1),
            ph_row("AGSDEFKLLR", 500, 530, "+", "tolerant", 1, "substitution", 3),
        ],
        columns=list(PEPTIDE_HITS_COLUMNS),
    )
    ht = hits_frame(
        [
            ht_row("MKLLAEGSTRK", 300, 334, "+", "300-317;318-334", search="decoy"),
            ht_row("MKLLAEGSTRK", 300, 334, "+", "300-318;319-334"),
            ht_row("MKLLAEGSTRK", 600, 634, "+", "600-617;618-634"),
        ]
    )
    return seq, ph, ht


def test_align_hits_one_alignment_per_row_in_row_order() -> None:
    mod = align()
    seq, ph, ht = three_loci()
    got = mod.align_hits(ph, ht, {"c1": seq}, {"c1": REP})
    assert len(got) == 3
    assert got[0] == mod.align_hit("PEPTIDEK", "c1", 100, 124, "+", 0, None, seq, REP)
    assert got[1] == mod.align_hit(
        "MKLLAEGSTRK", "c1", 300, 334, "+", 1, "300-318;319-334", seq, REP
    )
    assert got[2] == mod.align_hit("AGSDEFKLLR", "c1", 500, 530, "+", 1, None, seq, REP)
    assert got[1].frameshift is not None and got[1].frameshift.first_suffix_codon == 319
    assert [g.frame for g in got] == [2, 1, 3]


def test_align_hits_empty_table() -> None:
    ph = pd.DataFrame(columns=list(PEPTIDE_HITS_COLUMNS))
    assert align().align_hits(ph, None, {"c1": BG}, {"c1": REP}) == []
    assert align().align_hits(ph, hits_frame(), {"c1": BG}, {"c1": REP}) == []


def test_align_hits_peptide_with_several_loci() -> None:
    seq = genome((100, "AGSNEFKLLR", "+"), (600, "AGSNEFKLLR", "-"))
    ph = pd.DataFrame(
        [
            ph_row("AGSDEFKLLR", 100, 130, "+", "tolerant", 1, "substitution", 2),
            ph_row("AGSDEFKLLR", 600, 630, "-", "tolerant", 1, "substitution", -2),
        ],
        columns=list(PEPTIDE_HITS_COLUMNS),
    )
    got = align().align_hits(ph, None, {"c1": seq}, {"c1": REP})
    assert [(g.peptide, g.strand) for g in got] == [("AGSDEFKLLR", "+"), ("AGSDEFKLLR", "-")]
    assert got[0].ops[3].codons == (109,)
    assert got[1].ops[3].codons == (630 - 12,)


def test_align_hits_frameshift_row_without_work_row_raises() -> None:
    mod = align()
    seq, ph, ht = three_loci()
    with pytest.raises(ValueError, match="MKLLAEGSTRK"):
        mod.align_hits(ph, None, {"c1": seq}, {"c1": REP})
    only_decoy_and_other = ht[(ht["search"] == "decoy") | (ht["start"] == 600)]
    with pytest.raises(ValueError, match="MKLLAEGSTRK"):
        mod.align_hits(ph, only_decoy_and_other, {"c1": seq}, {"c1": REP})


def test_align_hits_frameshift_row_without_blocks_raises() -> None:
    seq, ph, _ = three_loci()
    ht = hits_frame([ht_row("MKLLAEGSTRK", 300, 334, "+", None)])
    with pytest.raises(ValueError, match="MKLLAEGSTRK"):
        align().align_hits(ph, ht, {"c1": seq}, {"c1": REP})


def test_align_hits_exact_row_mismatch_and_cost_mismatch_raise() -> None:
    seq, ph, ht = three_loci()
    bad = ph.copy()
    bad.loc[0, "peptide"] = "PEPTVDEK"  # beyond I/L on an exact row
    with pytest.raises(ValueError, match="PEPTVDEK"):
        align().align_hits(bad, ht, {"c1": seq}, {"c1": REP})
    bad = ph.copy()
    bad.loc[2, "n_diff"] = 2
    with pytest.raises(ValueError, match="AGSDEFKLLR"):
        align().align_hits(bad, ht, {"c1": seq}, {"c1": REP})


def test_align_hits_is_deterministic_and_does_not_mutate_inputs() -> None:
    mod = align()
    seq, ph, ht = three_loci()
    genome_map = {"c1": seq}
    reps = {"c1": REP}
    ph0, ht0, g0, r0 = ph.copy(deep=True), ht.copy(deep=True), dict(genome_map), dict(reps)
    first = mod.align_hits(ph, ht, genome_map, reps)
    second = mod.align_hits(ph, ht, genome_map, reps)
    assert first == second
    assert copy.deepcopy(first) == first
    pd.testing.assert_frame_equal(ph, ph0)
    pd.testing.assert_frame_equal(ht, ht0)
    assert genome_map == g0 and reps == r0


def test_results_are_frozen_value_objects() -> None:
    mod = align()
    seq = genome((300, "AGSNEFKLLR", "+"))
    a = mod.align_hit("AGSDEFKLLR", "c1", 300, 330, "+", 1, None, seq, REP)
    assert isinstance(a, mod.HitAlignment)
    assert isinstance(a.ops, tuple) and all(isinstance(o, mod.ResidueOp) for o in a.ops)
    with pytest.raises(AttributeError):
        a.ops[0].kind = "il"  # type: ignore[misc]
    with pytest.raises(AttributeError):
        a.peptide = "X"  # type: ignore[misc]
    assert a == mod.align_hit("AGSDEFKLLR", "c1", 300, 330, "+", 1, None, seq, REP)


def test_oracle_agrees_with_the_hand_worked_examples() -> None:
    """Guards the test oracle itself (no ``gippa.variants`` needed)."""
    seq = genome((300, "AGSNEFKLLR", "+"))
    assert locus_aa(seq, 300, 330, "+") == "AGSNEFKLLR"
    ops = expected_ops("AGSDEFKLLR", seq, 300, 330, "+")
    assert ops[3] == ("substitution", (3, 4), (309,), "N", "D")
    seq6 = genome((300, "AGGGEFKLLRK", "+"))
    assert expected_ops("ANGEFKLLRK", seq6, 300, 333, "+")[1] == (
        "block", (1, 2), (303, 306), "GG", "N"
    )  # fmt: skip
