"""M3-14 code-review decisions (PR #50): docs/specs/M3-14-loci-and-assignment.md, "Lead decisions
after code review" items 1-4, and docs/adr/ADR-0013-cterm-extension-locus-class.md.

1. An in-frame, same-strand hit running across a CDS/sORF 3' end (past the annotated stop) is
   ``cterm_extension`` and assigned to that gene, also across the origin; out of frame it stays
   ``cds_outframe``.
2. ``multi_identical`` needs a non-empty ``aa`` on both genes: two unrelated pseudogenes
   (``aa == ""``) are ``shared``.
3. Tie-breaks (smallest locus tag, as in M3-12) and the reviewer's missing boundary cases.
4. An unknown contig raises ``ValueError`` naming it.

Fixtures (model, genomes, row builders) come from tests/unit/test_loci.py.

Conventions pinned here where the decisions are silent:

* a hit that starts upstream of a short gene's start and runs past its end (no stop in between) is
  not ``nterm_extension`` (which only "reaches into" the CDS) but is in frame and extends past the
  annotated stop, so it is ``cterm_extension`` of that gene;
* a hit in frame across one gene's 3' end and out of frame with another same-strand CDS is
  ``cterm_extension`` (``cds_outframe`` only when no in-frame relation exists).

The lead's origin geometry ``[180, 210) +`` over G0 ``[9900, 10200)`` is out of frame after the
shift ((10180 - 9900) % 3 == 1), so it is pinned as ``cds_outframe``; ``[179, 209) +`` is the
in-frame origin case.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from test_loci import (  # type: ignore[import-not-found]
    G0,
    G1,
    G2,
    GENOME,
    MODEL,
    ORIGIN_MODEL,
    REPS,
    S1,
    acc,
    annotate,
    background,
    cds,
    classify,
    fill,
    flags_of,
    hit_row,
    iv,
    na,
    peps,
    put,
)

from giap.model.coords import Replicon
from giap.model.genes import GeneModel

# --- 1. cterm_extension ---------------------------------------------------------------------


def test_cterm_extension_is_a_locus_class() -> None:
    from giap.mapping.loci import LOCUS_CLASSES

    assert "cterm_extension" in LOCUS_CLASSES


@pytest.mark.parametrize(
    ("hit", "gm", "expected"),
    [
        # the reviewer's four geometries (in frame across the 3' end / annotated stop)
        pytest.param(iv(370, 415, "+"), MODEL, ("cterm_extension", G1), id="g1-plus"),
        pytest.param(iv(2060, 2090, "+"), MODEL, ("cterm_extension", S1), id="sorf-plus"),
        pytest.param(iv(985, 1015, "-"), MODEL, ("cterm_extension", G2), id="g2-minus"),
        # (10179 - 9900) % 3 == 0: in frame after the +L shift, over G0's end at 10200 (= 200)
        pytest.param(iv(179, 209, "+"), ORIGIN_MODEL, ("cterm_extension", G0), id="g0-origin"),
        # boundaries: 3 nt past the end (the stop codon itself) on both strands
        pytest.param(iv(373, 403, "+"), MODEL, ("cterm_extension", G1), id="g1-3nt-past"),
        pytest.param(iv(997, 1027, "-"), MODEL, ("cterm_extension", G2), id="g2-3nt-past"),
        pytest.param(iv(173, 203, "+"), ORIGIN_MODEL, ("cterm_extension", G0), id="g0-3nt-past"),
        # ending exactly at the CDS end is inside: cds_inframe
        pytest.param(iv(370, 400, "+"), MODEL, ("cds_inframe", G1), id="g1-ends-at-end"),
        pytest.param(iv(1000, 1030, "-"), MODEL, ("cds_inframe", G2), id="g2-ends-at-end"),
        pytest.param(iv(170, 200, "+"), ORIGIN_MODEL, ("cds_inframe", G0), id="g0-ends-at-end"),
        # out of frame across the 3' end: still cds_outframe
        pytest.param(iv(371, 416, "+"), MODEL, ("cds_outframe", None), id="g1-outframe"),
        pytest.param(iv(374, 404, "+"), MODEL, ("cds_outframe", None), id="g1-outframe-3nt"),
        pytest.param(iv(2061, 2091, "+"), MODEL, ("cds_outframe", None), id="sorf-outframe"),
        pytest.param(iv(986, 1016, "-"), MODEL, ("cds_outframe", None), id="g2-outframe"),
        pytest.param(iv(180, 210, "+"), ORIGIN_MODEL, ("cds_outframe", None), id="g0-outframe"),
        # opposite strand across the end: antisense
        pytest.param(iv(370, 415, "-"), MODEL, ("antisense", None), id="g1-antisense"),
    ],
)
def test_cterm_extension_geometries(
    hit: Any, gm: GeneModel, expected: tuple[str, str | None]
) -> None:
    assert tuple(classify(hit, gm)) == expected


def test_cterm_extension_beats_outframe_overlap() -> None:
    # H1 + [100, 400); H2 + [302, 700) is in another frame and overlaps the hit [370, 415)
    gm = GeneModel(REPS, (cds("H1", "chr1", 100, 400, "+"), cds("H2", "chr1", 302, 700, "+")), ())
    assert tuple(classify(iv(370, 415, "+"), gm)) == ("cterm_extension", "H1")
    # out of frame with both
    assert tuple(classify(iv(372, 417, "+"), gm)) == ("cds_outframe", None)


@settings(max_examples=150, deadline=None)
@given(a=st.integers(1, 90), b=st.integers(1, 150), shift=st.integers(0, 2))
def test_property_g1_3prime_overlap(a: int, b: int, shift: int) -> None:
    """Hits [400 - 3a + shift, 400 + 3b + shift): in frame -> cterm G1, else cds_outframe."""
    h = iv(400 - 3 * a + shift, 400 + 3 * b + shift, "+")
    want = ("cterm_extension", G1) if shift == 0 else ("cds_outframe", None)
    assert tuple(classify(h)) == want


def test_cterm_extension_assigned_in_annotate() -> None:
    rows = [
        hit_row("CTP", 370, 415, "+"),
        hit_row("CTM", 985, 1015, "-"),
        hit_row("CTS", 2060, 2090, "+"),
        hit_row("CTO", 371, 416, "+"),
    ]
    out = annotate(rows)
    for k, gene in (
        (("CTP", 370, 415, "+"), G1),
        (("CTM", 985, 1015, "-"), G2),
        (("CTS", 2060, 2090, "+"), S1),
    ):
        assert out[k]["locus_class"] == "cterm_extension"
        assert out[k]["gene"] == gene
        assert out[k]["uniqueness"] == "unique"
        assert "nterm_extension_evidence" not in flags_of(out[k]["flags"])
    o = out[("CTO", 371, 416, "+")]
    assert o["locus_class"] == "cds_outframe" and na(o["gene"])


def test_cterm_extension_origin_in_annotate() -> None:
    out = annotate([hit_row("CTG0", 179, 209, "+")], gm=ORIGIN_MODEL)
    r = out[("CTG0", 179, 209, "+")]
    assert (r["locus_class"], r["gene"]) == ("cterm_extension", G0)
    assert "origin_spanning" not in flags_of(r["flags"])


def test_cterm_and_inframe_of_one_gene_unique_but_multi() -> None:
    rows = [hit_row("CTU", 130, 160, "+"), hit_row("CTU", 370, 400, "+"),
            hit_row("CTU", 373, 403, "+")]  # fmt: skip
    out = annotate(rows)
    assert {r["gene"] for r in out.values()} == {G1}
    assert {r["uniqueness"] for r in out.values()} == {"unique"}
    assert all("multi" in flags_of(r["flags"]) for r in out.values())


# --- 2. pseudogenes with empty aa are shared -------------------------------------------------


def _pseudo_model(aa_a: str, aa_b: str) -> GeneModel:
    return GeneModel(
        REPS,
        (
            cds("PA", "chr1", 3000, 3600, "+", kind="pseudogene", aa=aa_a),
            cds("PB", "chr1", 4500, 4800, "-", kind="pseudogene", aa=aa_b),
        ),
        (),
    )


def test_two_unrelated_pseudogenes_empty_aa_are_shared() -> None:
    rows = [hit_row("PEPX", 3100, 3130, "+"), hit_row("PEPX", 4600, 4630, "-")]
    out = annotate(rows, gm=_pseudo_model("", ""))
    assert {r["gene"] for r in out.values()} == {"PA", "PB"}
    for r in out.values():
        assert r["locus_class"] == "pseudogene"
        assert r["uniqueness"] == "shared"
        assert "multi_identical" not in flags_of(r["flags"])


def test_cds_with_empty_aa_not_identical() -> None:
    gm = GeneModel(
        REPS,
        (cds("EA", "chr1", 100, 400, "+", aa=""), cds("EB", "chr1", 1000, 1300, "-", aa="")),
        (),
    )
    rows = [hit_row("PEPE", 130, 160, "+"), hit_row("PEPE", 1270, 1300, "-")]
    for r in annotate(rows, gm=gm).values():
        assert r["uniqueness"] == "shared"


def test_pseudogenes_with_identical_nonempty_aa_stay_multi_identical() -> None:
    same = "MKTAYIAKQRQISFVKSHFSRQ"
    rows = [hit_row("PEPY", 3100, 3130, "+"), hit_row("PEPY", 4600, 4630, "-")]
    for r in annotate(rows, gm=_pseudo_model(same, same)).values():
        assert r["uniqueness"] == "multi_identical"
        assert "multi_identical" in flags_of(r["flags"])


def test_one_empty_aa_one_not_is_shared() -> None:
    rows = [hit_row("PEPZ", 3100, 3130, "+"), hit_row("PEPZ", 4600, 4630, "-")]
    for r in annotate(rows, gm=_pseudo_model("", "MKT")).values():
        assert r["uniqueness"] == "shared"


# --- 3. tie-breaks -------------------------------------------------------------------------


def _both_orders(*genes: Any) -> list[GeneModel]:
    return [GeneModel(REPS, tuple(genes), ()), GeneModel(REPS, tuple(reversed(genes)), ())]


def test_tiebreak_nested_same_frame_cds() -> None:
    # ZZ + [190, 400) nested in AA + [100, 700), same frame ((190 - 100) % 3 == 0)
    for gm in _both_orders(cds("ZZ", "chr1", 190, 400, "+"), cds("AA", "chr1", 100, 700, "+")):
        assert tuple(classify(iv(220, 250, "+"), gm)) == ("cds_inframe", "AA")
    for gm in _both_orders(cds("AA", "chr1", 190, 400, "+"), cds("ZZ", "chr1", 100, 700, "+")):
        assert tuple(classify(iv(220, 250, "+"), gm)) == ("cds_inframe", "AA")


def test_tiebreak_nested_same_frame_minus() -> None:
    for gm in _both_orders(cds("ZM", "chr1", 1000, 1300, "-"), cds("AM", "chr1", 700, 1300, "-")):
        assert tuple(classify(iv(1270, 1300, "-"), gm)) == ("cds_inframe", "AM")


def test_tiebreak_two_overlapping_pseudogenes() -> None:
    for gm in _both_orders(
        cds("PZ", "chr1", 3000, 3600, "+", kind="pseudogene"),
        cds("PA", "chr1", 3050, 3700, "+", kind="pseudogene"),
    ):
        assert tuple(classify(iv(3100, 3130, "+"), gm)) == ("pseudogene", "PA")


def test_tiebreak_two_genes_sharing_a_start_nterm() -> None:
    # GENOME is stop-free on [40, 100): [70, 100) + is upstream of both
    for gm in _both_orders(cds("YS", "chr1", 100, 400, "+"), cds("XS", "chr1", 100, 700, "+")):
        assert tuple(classify(iv(70, 100, "+"), gm)) == ("nterm_extension", "XS")
    for gm in _both_orders(cds("XS", "chr1", 100, 400, "+"), cds("YS", "chr1", 100, 700, "+")):
        assert tuple(classify(iv(70, 100, "+"), gm)) == ("nterm_extension", "XS")


def test_tiebreak_two_genes_sharing_a_start_nterm_minus() -> None:
    # GENOME is stop-free on [1300, 1360) on the minus strand
    for gm in _both_orders(cds("YM", "chr1", 1000, 1300, "-"), cds("XM", "chr1", 700, 1300, "-")):
        assert tuple(classify(iv(1300, 1330, "-"), gm)) == ("nterm_extension", "XM")


def test_tiebreak_two_genes_sharing_an_end_cterm() -> None:
    for gm in _both_orders(cds("YE", "chr1", 100, 400, "+"), cds("XE", "chr1", 205, 400, "+")):
        assert tuple(classify(iv(370, 415, "+"), gm)) == ("cterm_extension", "XE")


# --- 3. missing boundary cases -------------------------------------------------------------


def _with(genome: Mapping[str, str], pos: int, codon: str) -> dict[str, str]:
    s = list(genome["chr1"])
    put(s, pos, codon)
    return {**genome, "chr1": "".join(s)}


def test_stop_immediately_before_start_plus() -> None:
    gm = GeneModel(REPS, (cds(G1, "chr1", 100, 400, "+"),), ())
    h = iv(67, 97, "+")
    assert tuple(classify(h, gm, GENOME)) == ("nterm_extension", G1)
    # TAA at [97, 100), in frame and adjacent to the start codon
    assert tuple(classify(h, gm, _with(GENOME, 97, "TAA"))) == ("intergenic", None)
    assert tuple(classify(iv(40, 70, "+"), gm, _with(GENOME, 97, "TAA"))) == ("intergenic", None)


def test_stop_immediately_before_start_minus() -> None:
    gm = GeneModel(REPS, (cds(G2, "chr1", 1000, 1300, "-"),), ())
    h = iv(1303, 1333, "-")
    assert tuple(classify(h, gm, GENOME)) == ("nterm_extension", G2)
    # TTA on plus = TAA on minus at [1300, 1303), in frame and adjacent to the start codon
    g = _with(GENOME, 1300, "TTA")
    assert tuple(classify(h, gm, g)) == ("intergenic", None)
    assert tuple(classify(iv(1330, 1360, "-"), gm, g)) == ("intergenic", None)


def test_upstream_hit_past_end_of_short_gene_plus() -> None:
    gm = GeneModel(REPS, (cds("SG", "chr1", 100, 130, "+"),), ())
    assert tuple(classify(iv(70, 130, "+"), gm, GENOME)) == ("nterm_extension", "SG")
    assert tuple(classify(iv(70, 127, "+"), gm, GENOME)) == ("nterm_extension", "SG")
    assert tuple(classify(iv(70, 133, "+"), gm, GENOME)) == ("cterm_extension", "SG")


def test_upstream_hit_past_end_of_short_gene_minus() -> None:
    gm = GeneModel(REPS, (cds("SM", "chr1", 1270, 1300, "-"),), ())
    assert tuple(classify(iv(1270, 1330, "-"), gm, GENOME)) == ("nterm_extension", "SM")
    assert tuple(classify(iv(1273, 1330, "-"), gm, GENOME)) == ("nterm_extension", "SM")
    assert tuple(classify(iv(1267, 1330, "-"), gm, GENOME)) == ("cterm_extension", "SM")


@pytest.mark.parametrize("split_at", [2998, 2999], ids=["2+1", "1+2"])
def test_stop_codon_split_across_origin(split_at: int) -> None:
    n = 3000
    gstart = 31 if split_at == 2998 else 32  # in frame: (split_at - (gstart + n)) % 3 == 0
    rep = {"c5": Replicon("c5", n, "circular", "plasmid")}
    gm = GeneModel(rep, (cds("GS", "c5", gstart, gstart + 300, "+"),), ())
    s = background(n)
    fill(s, 2900, n + gstart)
    clean = {"c5": "".join(s)}
    hs = 2950 + (gstart - 31)  # in frame with the gene
    h = iv(hs, hs + 30, "+", "c5")
    assert (split_at - (gstart + n)) % 3 == 0 and (hs - (gstart + n)) % 3 == 0
    assert tuple(classify(h, gm, clean)) == ("nterm_extension", "GS")
    put(s, split_at, "TAA")
    stopped = "".join(s)
    assert (stopped[split_at:] + stopped[:3])[:3] == "TAA"  # really split across the origin
    assert tuple(classify(h, gm, {"c5": stopped})) == ("intergenic", None)


def test_origin_spanning_end_equal_to_length() -> None:
    out = annotate([hit_row("ENDL", 9970, 10000, "+"), hit_row("ENDP", 9973, 10003, "+")])
    assert "origin_spanning" not in flags_of(out[("ENDL", 9970, 10000, "+")]["flags"])
    assert "origin_spanning" in flags_of(out[("ENDP", 9973, 10003, "+")]["flags"])


# --- 4. unknown contig ---------------------------------------------------------------------


def test_classify_unknown_contig_raises_value_error() -> None:
    with pytest.raises(ValueError, match="chrX"):
        classify(iv(100, 130, "+", "chrX"))


def test_annotate_unknown_contig_raises_value_error() -> None:
    from giap.mapping.loci import annotate_peptide_hits

    rows = [hit_row("PEPU", 130, 160, "+"), hit_row("PEPV", 100, 130, "+", contig="plasmidQ")]
    with pytest.raises(ValueError, match="plasmidQ"):
        annotate_peptide_hits(acc(rows), MODEL, GENOME, peps(rows), "related")
