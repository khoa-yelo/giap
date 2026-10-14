"""Bakta-derived layers: ribosomal set, 16S and anti-SD, tRNA counts, xrefs, nearest genes (M2-04).

Examples 1-9 of ``docs/specs/M2-04-parse-bakta-step.md`` plus its edge cases; GIAP_SPEC.md §7.1.2
(ribosomal-protein regexes, anti-SD rule, tRNA copy numbers), §9.3 ``annotation_source``.
Lead decision M2-index 7: ``annotation_source`` comes from the Bakta DB type, ``user`` only when
the DB type is unknown. M2-04 Q1 (proposed default): a configured anti-SD other than ``auto`` is a
9-nt DNA string used as is and validated ``^[ACGT]{9}$``.
"""

from __future__ import annotations

import pytest

from giap.model import Interval
from giap.model.coords import Replicon
from giap.model.genes import Gene, NcFeature

ECOLI_TAIL = "GATCACCTCCTTA"


def _gene(
    locus: str,
    start: int,
    end: int,
    strand: str = "+",
    *,
    contig: str = "chr1",
    gene: str | None = None,
    product: str = "",
    xrefs: tuple[str, ...] = (),
    feature_type: str = "cds",
) -> Gene:
    return Gene(
        locus_tag=locus,
        feature_type=feature_type,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon=None,
        gene=gene,
        product=product,
        db_xrefs=xrefs,
    )


def _nc(
    fid: str,
    type_: str,
    start: int = 0,
    end: int = 100,
    *,
    contig: str = "chr1",
    strand: str = "+",
    gene: str | None = None,
    product: str = "",
    anticodon: str | None = None,
    amino_acid: str | None = None,
) -> NcFeature:
    return NcFeature(
        id=fid,
        type=type_,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        gene=gene,
        product=product,
        anticodon=anticodon,
        amino_acid=amino_acid,
    )


# ------------------------------------------------------------------ constants


def test_constants() -> None:
    from giap.annotate import bakta_layers as bl

    assert bl.DEFAULT_ANTI_SD == "ACCTCCTTA"
    assert bl.RIBOSOMAL_PRODUCT.search("30S ribosomal protein S1")
    assert bl.RIBOSOMAL_PRODUCT.search("50S ribosomal protein L2")
    assert not bl.RIBOSOMAL_PRODUCT.search("putative 50S ribosomal protein L2")
    assert bl.RIBOSOMAL_GENE.search("rpmH")
    assert not bl.RIBOSOMAL_GENE.search("rpoB")


# ------------------------------------------------------------------ anti-SD (examples 1-5)


def test_example1_seven_copies() -> None:
    from giap.annotate.bakta_layers import anti_sd

    copies = ["ACGT" * 10 + ECOLI_TAIL] * 7
    assert anti_sd(copies) == ("ACCTCCTTA", ())


def test_example2_no_16s() -> None:
    from giap.annotate.bakta_layers import anti_sd

    assert anti_sd([]) == ("ACCTCCTTA", ("anti_sd_default",))


def test_example3_tail_without_cctcc() -> None:
    from giap.annotate.bakta_layers import anti_sd

    assert anti_sd(["ACGTACGT" + "GATCATTTAAG"]) == ("ACCTCCTTA", ("anti_sd_default",))


def test_example4_majority_wins() -> None:
    from giap.annotate.bakta_layers import anti_sd

    copies = ["GGGG" + "ACCTCCTTA", "TTTT" + "ACCTCCTTT", "CCCC" + "ACCTCCTTA"]
    assert anti_sd(copies) == ("ACCTCCTTA", ())


def test_example5_user_override_ignores_16s() -> None:
    from giap.annotate.bakta_layers import anti_sd

    assert anti_sd(["AAAA" + ECOLI_TAIL], configured="TCCTCCACT") == ("TCCTCCACT", ())
    assert anti_sd([], configured="TCCTCCACT") == ("TCCTCCACT", ())


def test_explicit_auto_is_default() -> None:
    from giap.annotate.bakta_layers import anti_sd

    copies = ["AAAA" + ECOLI_TAIL]
    assert anti_sd(copies, configured="auto") == anti_sd(copies)


@pytest.mark.parametrize("bad", ["ACCTCCTT", "ACCTCCTTAA", "ACCUCCUUA", "acctcctta", "NCCTCCTTA"])
def test_override_must_be_nine_nt_dna(bad: str) -> None:
    from giap.annotate.bakta_layers import anti_sd

    with pytest.raises(ValueError):
        anti_sd([], configured=bad)


def test_tie_first_copy_in_genome_order_wins() -> None:
    from giap.annotate.bakta_layers import anti_sd

    assert anti_sd(["AA" + "ACCTCCTTT", "AA" + "ACCTCCTTA"]) == ("ACCTCCTTT", ())
    assert anti_sd(["AA" + "ACCTCCTTA", "AA" + "ACCTCCTTT"]) == ("ACCTCCTTA", ())


def test_short_16s_skipped() -> None:
    from giap.annotate.bakta_layers import anti_sd

    assert anti_sd(["CCTCC", "CCTCCTT", "AAAA" + "TCCTCCACT"]) == ("TCCTCCACT", ())
    assert anti_sd(["CCTCCTTA"]) == ("ACCTCCTTA", ("anti_sd_default",))


def test_most_common_lacks_cctcc_gives_default() -> None:
    """The rule checks the most common 9-mer only; a minority copy with CCTCC does not rescue it."""
    from giap.annotate.bakta_layers import anti_sd

    copies = ["AA" + "TCATTTAAG", "AA" + "TCATTTAAG", "AA" + "TCCTCCACT"]
    assert anti_sd(copies) == ("ACCTCCTTA", ("anti_sd_default",))


# ------------------------------------------------------------------ ribosomal proteins (example 6)


@pytest.mark.parametrize(
    ("gene", "product", "expected"),
    [
        ("rplB", "50S ribosomal protein L2", True),
        ("rpoB", "DNA-directed RNA polymerase subunit beta", False),
        ("rpsA", "hypothetical protein", True),
        (None, "30S ribosomal protein S12", True),
        (None, "50S ribosomal protein L7/L12", True),
        ("rpmH", "", True),
        ("rpsa", "", False),  # gene regex needs an upper-case letter after the prefix
        (None, "ribosomal protein L7Ae-like", False),
        (None, "putative 30S ribosomal protein", False),
        ("xrpsA", "", False),
        (None, "", False),
    ],
)
def test_example6_is_ribosomal_protein(gene: str | None, product: str, expected: bool) -> None:
    from giap.annotate.bakta_layers import is_ribosomal_protein

    assert is_ribosomal_protein(_gene("X_1", 0, 300, gene=gene, product=product)) is expected


def test_ribosomal_pseudogene_still_ribosomal() -> None:
    """Edge case: a gene can be both ribosomal_protein and pseudogene (M2-17 excludes it later)."""
    from giap.annotate.bakta_layers import is_ribosomal_protein

    g = _gene(
        "X_1", 0, 300, gene="rplB", product="50S ribosomal protein L2", feature_type="pseudogene"
    )
    assert is_ribosomal_protein(g) is True


# ------------------------------------------------------------------ 16S selection


def test_rrna16s_selection_by_product_or_gene() -> None:
    from giap.annotate.bakta_layers import rrna16s

    feats = [
        _nc("r1", "rrna", 0, 1500, product="16S ribosomal RNA", gene="rrs"),
        _nc("r2", "rrna", 2000, 4900, product="23S ribosomal RNA", gene="rrl"),
        _nc("r3", "rrna", 5000, 6500, product="", gene="rrs"),
        _nc("r4", "rrna", 7000, 8500, product="16S ribosomal RNA"),
        _nc("r5", "rrna", 9000, 9120, product="5S ribosomal RNA", gene="rrf"),
        _nc("t1", "trna", 9500, 9576, product="16S-like tRNA", gene="rrs"),
        _nc("n1", "ncrna", 9600, 9700, product="16S"),
    ]
    assert [f.id for f in rrna16s(feats)] == ["r1", "r3", "r4"]


def test_rrna16s_empty() -> None:
    from giap.annotate.bakta_layers import rrna16s

    assert rrna16s([]) == []


# ------------------------------------------------------------------ tRNA counts (example 7)


def test_example7_trna_counts() -> None:
    from giap.annotate.bakta_layers import trna_counts

    feats = [
        _nc("t1", "trna", 0, 76, anticodon="CAG", amino_acid="Leu"),
        _nc("t2", "trna", 100, 176, anticodon="CAG", amino_acid="Leu"),
        _nc("t3", "trna", 200, 276, anticodon="GCC", amino_acid="Gly"),
    ]
    assert trna_counts(feats) == {"CAG": 2, "GCC": 1}


def test_trna_counts_sorted_and_only_trna() -> None:
    from giap.annotate.bakta_layers import trna_counts

    feats = [
        _nc("t1", "trna", 0, 76, anticodon="TTC", amino_acid="Glu"),
        _nc("r1", "rrna", 100, 1600, product="16S ribosomal RNA"),
        _nc("t2", "trna", 2000, 2076, anticodon="AAA", amino_acid="Phe"),
        _nc("t3", "trna", 3000, 3076, anticodon="TTC", amino_acid="Glu"),
        _nc("m1", "tmrna", 4000, 4363),
    ]
    counts = trna_counts(feats)
    assert counts == {"AAA": 1, "TTC": 2}
    assert list(counts) == sorted(counts)


# ------------------------------------------------------------------ xrefs (example 8)


def test_example8_split_xrefs() -> None:
    from giap.annotate.bakta_layers import split_xrefs

    g = _gene(
        "b3702",
        0,
        300,
        xrefs=(
            "COG:COG0593",
            "GO:0005524",
            "EC:3.6.4.12",
            "KEGG:K02313",
            "UniRef:UniRef90_P03004",
        ),
    )
    assert split_xrefs([g]) == [
        ("b3702", "COG", "COG0593"),
        ("b3702", "GO", "GO:0005524"),
        ("b3702", "EC", "3.6.4.12"),
        ("b3702", "KEGG", "K02313"),
        ("b3702", "UniRef", "UniRef90_P03004"),
    ]


def test_split_xrefs_gene_order_and_empty() -> None:
    from giap.annotate.bakta_layers import split_xrefs

    genes = [
        _gene("A_1", 0, 300, xrefs=("RefSeq:WP_000001.1",)),
        _gene("A_2", 400, 700),
        _gene("A_3", 800, 1100, xrefs=("SO:0001217", "UniParc:UPI0000000001")),
    ]
    assert split_xrefs(genes) == [
        ("A_1", "RefSeq", "WP_000001.1"),
        ("A_3", "SO", "0001217"),
        ("A_3", "UniParc", "UPI0000000001"),
    ]
    assert split_xrefs([]) == []


def test_split_xrefs_splits_at_first_colon_only() -> None:
    from giap.annotate.bakta_layers import split_xrefs

    g = _gene("A_1", 0, 300, xrefs=("UniRef:UniRef50:P0A7V0",))
    assert split_xrefs([g]) == [("A_1", "UniRef", "UniRef50:P0A7V0")]


# ------------------------------------------------------------------ annotation_source (example 9)


@pytest.mark.parametrize(
    ("db_type", "expected"),
    [("full", "bakta_full"), ("light", "bakta_light"), ("unknown", "user")],
)
def test_example9_annotation_source(db_type: str, expected: str) -> None:
    from giap.annotate.bakta_layers import annotation_source

    assert annotation_source(db_type) == expected


# ------------------------------------------------------------------ nearest_genes


CIRC = Replicon("chr1", 10_000, "circular", "chromosome")
LIN = Replicon("chr1", 10_000, "linear", "unknown")


def _origin_genes() -> list[Gene]:
    return [
        _gene("A", 9_800, 9_950, "+"),  # 150 bp before the oriC across the origin
        _gene("B", 500, 800, "-"),  # 200 bp after
        _gene("C", 5_000, 5_100, "+"),  # far
        _gene("P", 100, 400, "+", contig="plasmid_1"),  # other contig: never a neighbour
    ]


def test_nearest_genes_circular_is_origin_aware() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    assert nearest_genes(oric, _origin_genes(), CIRC) == [("A", 150), ("B", 200)]


def test_nearest_genes_linear_never_wraps() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    assert nearest_genes(oric, _origin_genes(), LIN) == [("B", 200), ("C", 4_700)]
    assert nearest_genes(oric, _origin_genes(), LIN, n=3) == [
        ("B", 200),
        ("C", 4_700),
        ("A", 9_500),
    ]


def test_nearest_genes_overlap_is_zero() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    genes = [_gene("D", 250, 400, "+"), *_origin_genes()]
    assert nearest_genes(oric, genes, CIRC, n=1) == [("D", 0)]


def test_nearest_genes_origin_crossing_gene_overlaps() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    genes = [_gene("O", 9_900, 10_150, "+"), _gene("B", 500, 800, "-")]
    assert nearest_genes(oric, genes, CIRC) == [("O", 0), ("B", 200)]


def test_nearest_genes_adjacent_is_zero_gap() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    genes = [_gene("E", 300, 600, "+"), _gene("F", 40, 100, "-")]
    assert sorted(nearest_genes(oric, genes, LIN)) == [("E", 0), ("F", 0)]


def test_nearest_genes_fewer_than_n() -> None:
    from giap.annotate.bakta_layers import nearest_genes

    oric = _nc("oriC_1", "oric", 100, 300, strand=".")
    assert nearest_genes(oric, [_gene("B", 500, 800)], CIRC, n=5) == [("B", 200)]
    assert nearest_genes(oric, [], CIRC) == []
