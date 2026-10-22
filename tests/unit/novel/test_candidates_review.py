"""M7-06 review (#161): antisense ORFs name the opposite-strand gene they overlap (rule 5; C11
``gene`` is NA for antisense), the out-of-frame gene is the smallest same-strand tag, and I/L
twins share one identity for the contaminant ban and the multi-locus count."""

from __future__ import annotations

from helpers.novel_oracle import Case, cds, hit_row, il
from unit.novel.test_candidates import C1, P1, L, ex1, ex1_hits, ex1_seq, ex2, run


def test_antisense_overlapping_gene_is_opposite_strand_gene() -> None:
    seq = ex1_seq().text()
    genes = [cds("AS1", "c1", 60, 300, "-"), cds("SAME", "c1", 60, 300, "+")]
    (o,) = run(Case({"c1": seq}, dict(C1), ex1_hits(seq, "antisense", None), genes=genes))
    assert o.cls == "antisense"
    assert o.overlapping_gene == "AS1"


def test_cds_outframe_gene_is_smallest_same_strand_tag() -> None:
    seq = ex1_seq().text()
    genes = [cds("B2", "c1", 60, 300, "+"), cds("A1", "c1", 61, 301, "+"),
             cds("A0", "c1", 60, 300, "-")]  # fmt: skip
    (o,) = run(Case({"c1": seq}, dict(C1), ex1_hits(seq, "cds_outframe", None), genes=genes))
    assert o.overlapping_gene == "A1"


def test_il_twin_of_a_contaminant_is_excluded() -> None:
    case = ex1()
    case.hits.append(hit_row(il(P1), "c1", 72, 102, "+", "intergenic", L))
    case.contaminant = (P1,)
    (o,) = run(case)
    assert P1 not in o.peptides and il(P1) not in o.peptides


def test_il_twins_at_two_loci_are_multi_locus() -> None:
    case = ex2()
    case.hits.append(hit_row(il(P1), "c1", 151, 181, "+", "intergenic", L))
    for o in run(case):
        assert P1 not in o.peptides and il(P1) not in o.peptides
