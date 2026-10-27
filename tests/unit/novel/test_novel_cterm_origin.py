"""``spans_annotated_stop`` across the origin (M7-16b rule 12; code review of PR #335, must-fix 4).

On a circular replicon the gene's stop codon is compared with peptide intervals stored unwrapped
(``start < n``, ``end`` possibly ``> n``), so the codon offset ``lo - start`` must be taken modulo
``n`` when the codon lies past the origin and the peptide starts before it.

World: the toy ``c1`` of ``test_novel_cterm.py`` (G1 ``+`` ``[0, 57)``, TGA ``[54, 57)``, peptide
``P1`` = ``WQLNEFKWVPAGHTL`` ``[33, 78)`` over the stop), rotated onto a **circular** 126-nt
contig so that toy position p sits at ``(p - shift) % 126``. The minus world is the mirror of the
linear toy, rotated the same way.
"""

from __future__ import annotations

from typing import Any

import pytest
from helpers.novel_evidence_fixture import evidence_table, flags_of
from helpers.novel_oracle import Case, cds, hit_row, mirror_case
from helpers.novel_outputs_fixture import orfs_of
from helpers.readthrough_toys import toy_c1

from gippa.model import Replicon
from gippa.model.coords import Interval

N = 126
P1 = ("WQLNEFKWVPAGHTL", 33, 78)


def _hit(pep: str, s: int, e: int, strand: str) -> dict[str, Any]:
    r = hit_row(pep, "c1", s, e, strand, "cterm_extension", N, gene="G1", route="tolerant")
    r["n_diff"], r["diff_types"] = 1, ["substitution"]
    return r


def _linear() -> Case:
    reps = {"c1": Replicon("c1", N, "linear")}
    return Case(dict(toy_c1()[0]), reps, [_hit(*P1, "+")], genes=[cds("G1", "c1", 0, 57, "+")])


def _rotate(case: Case, shift: int) -> Case:
    """Move everything by ``-shift`` onto a circular contig (intervals stay unwrapped)."""
    seq = case.genome["c1"]
    genome = {"c1": seq[shift:] + seq[:shift]}
    reps = {"c1": Replicon("c1", N, "circular")}

    def mv(s: int, e: int) -> tuple[int, int]:
        a = (s - shift) % N
        return a, a + e - s

    hits = [_hit(h["peptide"], *mv(h["start"], h["end"]), h["strand"]) for h in case.hits]
    genes = [
        cds(g.locus_tag, "c1", *mv(g.interval.start, g.interval.end), g.interval.strand)
        for g in case.genes
    ]
    return Case(genome, reps, hits, genes=genes)


def _spans(case: Case) -> bool:
    (orf,) = orfs_of(case)
    assert orf.cls == "cterm_extension"
    ivs = [Interval("c1", h["start"], h["end"], h["strand"]) for h in case.hits]
    (r,) = evidence_table(case, peptide_intervals={orf.orf_id: ivs}).to_dict("records")
    return "readthrough_candidate" in flags_of(r)


def test_plus_peptide_across_the_origin_spans_the_stop_past_it() -> None:
    """Shift 45: G1 ``[81, 138)`` (unwrapped), stop ``[9, 12)``, P1 ``[114, 159)``. The codon
    offset ``9 - 114`` wraps to 21 (in frame, 3' flank of 24 nt after it)."""
    c = _rotate(_linear(), 45)
    (h,) = c.hits
    assert (h["start"], h["end"]) == (114, 159)
    assert _spans(c)


def test_minus_peptide_across_the_origin_spans_the_stop_past_it() -> None:
    """Mirror (G1 ``-`` ``[69, 126)``, stop ``[69, 72)``, P1 ``[48, 93)``) rotated by 60: G1
    ``[9, 66)``, stop ``[9, 12)``, P1 ``[114, 159)``; the 3' side is below the codon, i.e. the
    part of P1 before the origin."""
    c = _rotate(mirror_case(_linear()), 60)
    (g,) = c.genes
    (h,) = c.hits
    assert (g.interval.start, g.interval.end, g.interval.strand) == (9, 66, "-")
    assert (h["start"], h["end"], h["strand"]) == (114, 159, "-")
    assert _spans(c)


@pytest.mark.parametrize("shift", [20, 40, 45, 54, 60, 100])
@pytest.mark.parametrize("minus", [False, True])
def test_rotations_agree_with_the_linear_call(shift: int, minus: bool) -> None:
    """Rotation invariance: the call does not depend on where the origin falls (it crosses P1
    for shifts 34 to 77 on ``+``)."""
    base = mirror_case(_linear()) if minus else _linear()
    assert _spans(base)
    assert _spans(_rotate(base, shift))
