"""M7-03 VCF records (docs/specs/M7-03-variant-tables-and-vcf.md; GIPPA_SPEC.md §7.5.5, G14,
Appendix B.7, §6.5; ADR-0005 item 1).

Examples 1-11 and the degenerate rows of the spec, on the fixture genome of
``helpers.variant_tables_fixture``. Inputs are built with the public M7-02 ``pileup``/``sites``
and passed through untouched, so these tests do not depend on how the pileup carries its ALT
detail (lead decision 2). Records are written with the merged M1-05 ``write_vcf`` (with
``sequences=genome``, so every REF is checked against the genome).

Interface reading (test author; see the report to the lead):

* ``vcf_records(sites, pile, genome, replicons) -> (records, table)``: ``records`` are
  ``gippa.io.vcf.VcfRecord``\\ s (any order); ``table`` holds one row per ``sites`` row, same
  ``variant_id``\\ s in the same order, with ``nt_changes``, ``alt_codon_vcf`` and ``flags``.
* M7-05 flags reach ``vcf_records`` as a ``flags`` column of tuples on ``sites``; a missing
  column means no flags. A ``class_error_rate`` column (float or ``unestimated``) feeds
  ``CLASS_ERR``; a missing column means ``unestimated``.
* A two-codon event (swap, 2<->2 block) gets ``POS`` = the lowest 1-based coordinate of its two
  codons (on ``-`` the site's ``codon_start`` is the peptide-order first codon, the higher one),
  ``AA_REF`` = the two REF residues and ``AA_ALT`` = the two ALT residues.
* A two-codon event whose 6-nt REF crosses the origin is unrepresentable like a crossing codon.
* The ``sample = auto`` degenerate row is enforced by ``write_vcf`` (the only VCF writer).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers.variant_tables_fixture import (
    CODONS64,
    GENOME,
    GM_,
    GP,
    REPS,
    RESIDUES21,
    SAMPLE,
    SEQ_C2,
    SEQ_CHR,
    SHA,
    VERSION,
    aa_at,
    al,
    body_lines,
    flag_set,
    gene_codon_nt,
    gene_peptide,
    info_dict,
    is_na,
    make_hit,
    match_op,
    minimal_codon,
    peptides_at,
    pl,
    raw_hit,
    rc,
    records,
    ref_peptides,
    site_row,
    sub,
    translate,
    two,
    vcf_io,
    vr,
    write,
)

from gippa.model import Replicon

# --------------------------------------------------------------------------- minimal_alt_codon


def test_example2_j_over_i_and_l_codons() -> None:
    assert minimal_codon("GTT", "J") == ("ATT", 1)
    assert vr().minimal_alt_codon("GTT", "J") == ("ATT", 1)


def test_example3_tie_is_lexicographic() -> None:
    assert minimal_codon("AAA", "N") == ("AAC", 1)
    assert vr().minimal_alt_codon("AAA", "N") == ("AAC", 1)


def test_example6_and_8_oracle_values() -> None:
    """The spec's worked values, re-derived by the brute-force oracle."""
    assert minimal_codon("AAT", "D") == ("GAT", 1)
    assert minimal_codon("AAT", "E") == ("GAA", 2)
    assert minimal_codon("TCT", "D") == ("GAT", 2)
    assert minimal_codon("GAT", "S") == ("AGT", 2)
    assert minimal_codon("GCT", "S") == ("TCT", 1)


def test_minimal_alt_codon_exhaustive_table() -> None:
    """§7.5.6: 64 codons x 21 residues (20 + ``*``) plus ``J`` against the brute-force oracle."""
    f = vr().minimal_alt_codon
    assert len(CODONS64) == 64 and len(RESIDUES21) == 21
    bad = [
        (c, aa, f(c, aa), minimal_codon(c, aa))
        for c in CODONS64
        for aa in (*RESIDUES21, "J")
        if tuple(f(c, aa)) != minimal_codon(c, aa)
    ]
    assert not bad, bad[:10]


def test_minimal_alt_codon_of_the_ref_residue_is_the_ref_codon() -> None:
    f = vr().minimal_alt_codon
    for c in CODONS64:
        assert tuple(f(c, translate(c))) == (c, 0)


def test_minimal_alt_codon_stop() -> None:
    """``*`` cannot come from a peptide but is defined (degenerate row 3)."""
    f = vr().minimal_alt_codon
    assert tuple(f("TGG", "*")) == ("TAG", 1)  # TAG, TGA both 1 change; lexicographic
    assert minimal_codon("TGG", "*") == ("TAG", 1)


# --------------------------------------------------------------------------- INFO / FORMAT


def _decl(fields: Any) -> list[tuple[str, str, str]]:
    return [(f.id, f.number, f.type) for f in fields]


def test_info_fields_rule_6() -> None:
    assert _decl(vr().INFO_FIELDS) == [
        ("GENE", "1", "String"),
        ("AA_POS", "1", "Integer"),
        ("AA_REF", "1", "String"),
        ("AA_ALT", "A", "String"),
        ("CALL", "1", "String"),
        ("CONF", "A", "String"),
        ("NT_CHANGES", "A", "Integer"),
        ("ALT_TYPE", "A", "String"),
        ("FLAGS", ".", "String"),
        ("CLASS_ERR", "A", "String"),
    ]
    assert all(f.description.strip() for f in vr().INFO_FIELDS)


def test_format_fields_rule_7() -> None:
    assert _decl(vr().FORMAT_FIELDS) == [
        ("GT", "1", "String"),
        ("AD", "R", "Integer"),
        ("DP", "1", "Integer"),
    ]
    assert all(f.description.strip() for f in vr().FORMAT_FIELDS)


# --------------------------------------------------------------------------- examples


def _one(recs: list[Any]) -> Any:
    assert len(recs) == 1, recs
    return recs[0]


def test_records_are_vcf_records() -> None:
    recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    assert all(isinstance(r, vcf_io().VcfRecord) for r in recs)


def test_example1_appendix_b7_record(tmp_path: Path) -> None:
    aligns = peptides_at(GP, 44, sub(GP, 44, "S"), 3)
    recs, table, _, _ = records(aligns, {"flags": [("mass_equivalent", "snv")]})
    r = _one(recs)
    assert (r.contig, r.pos0, r.ref, tuple(r.alts)) == ("chr", 1000, "GCT", ("TCT",))
    assert (r.gt, tuple(r.ad), r.dp, r.filter) == ("1", (0, 3), 3, "PASS")
    write(tmp_path / "v.vcf", recs)
    (line,) = body_lines(tmp_path / "v.vcf")
    assert line == [
        "chr", "1001", ".", "GCT", "TCT", ".", "PASS",
        "GENE=NISS_00012;AA_POS=45;AA_REF=A;AA_ALT=S;CALL=ALT;CONF=high;NT_CHANGES=1;"
        "ALT_TYPE=substitution;FLAGS=mass_equivalent|snv;CLASS_ERR=unestimated",
        "GT:AD:DP", "1:0,3:3",
    ]  # fmt: skip
    (row,) = site_row(table, "chr", 1000, "+")
    assert (int(row["nt_changes"]), row["alt_codon_vcf"]) == (1, "TCT")
    assert flag_set(row["flags"]) == {"mass_equivalent", "snv"}


def test_example4_minus_strand_codon_conversion(tmp_path: Path) -> None:
    """§7.5.6 "minus-strand codon conversion": gene strand GCT -> TCT; forward AGC -> AGA."""
    assert SEQ_CHR[1567:1570] == "AGC" and gene_codon_nt(GM_, 10) == "GCT"
    recs, table, _, _ = records(peptides_at(GM_, 10, sub(GM_, 10, "S"), 1))
    r = _one(recs)
    assert (r.contig, r.pos0, r.ref, tuple(r.alts)) == ("chr", 1567, "AGC", ("AGA",))
    (row,) = site_row(table, "chr", 1567, "-")
    assert (row["alt_codon_vcf"], int(row["nt_changes"])) == ("AGA", 1)
    write(tmp_path / "v.vcf", recs)
    (line,) = body_lines(tmp_path / "v.vcf")
    assert line[1] == "1568"
    info = info_dict(line[7])
    assert (info["GENE"], info["AA_POS"], info["AA_REF"], info["AA_ALT"]) == (
        "NISS_00020", "11", "A", "S",
    )  # fmt: skip
    assert info["CONF"] == "candidate" and line[6] == "PASS"


def test_example5_mixed_gt_missing(tmp_path: Path) -> None:
    aligns = ref_peptides(GP, 50, 1) + peptides_at(GP, 50, sub(GP, 50, "D"), 2)
    recs, *_ = records(aligns)
    r = _one(recs)
    assert (r.gt, tuple(r.ad), r.dp, tuple(r.alts)) == (".", (1, 2), 3, ("GAT",))
    write(tmp_path / "v.vcf", recs)
    (line,) = body_lines(tmp_path / "v.vcf")
    assert info_dict(line[7])["CALL"] == "mixed" and line[9] == ".:1,2:3"


def test_example6_multiallelic(tmp_path: Path) -> None:
    aligns = peptides_at(GP, 50, sub(GP, 50, "D"), 2) + peptides_at(GP, 50, sub(GP, 50, "E"), 1)
    recs, table, _, _ = records(aligns)
    r = _one(recs)
    assert (r.pos0, r.ref, tuple(r.alts)) == (1018, "AAT", ("GAT", "GAA"))
    assert (r.gt, tuple(r.ad), r.dp) == ("1", (0, 2, 1), 3)
    write(tmp_path / "v.vcf", recs)
    (line,) = body_lines(tmp_path / "v.vcf")
    info = info_dict(line[7])
    assert info["AA_ALT"] == "D,E" and info["CONF"] == "high,candidate"
    assert info["NT_CHANGES"] == "1,2" and info["ALT_TYPE"] == "substitution,substitution"
    assert info["CLASS_ERR"] == "unestimated,unestimated" and info["CALL"] == "ALT"
    assert line[4] == "GAT,GAA" and line[9] == "1:0,2,1:3"
    by_alt = {r["alt_aa"]: r for r in site_row(table, "chr", 1018, "+")}
    assert (by_alt["D"]["alt_codon_vcf"], int(by_alt["D"]["nt_changes"])) == ("GAT", 1)
    assert (by_alt["E"]["alt_codon_vcf"], int(by_alt["E"]["nt_changes"])) == ("GAA", 2)


def test_alleles_ordered_by_n_alt_descending() -> None:
    aligns = peptides_at(GP, 50, sub(GP, 50, "D"), 1) + peptides_at(GP, 50, sub(GP, 50, "E"), 2)
    r = _one(records(aligns)[0])
    assert tuple(r.alts) == ("GAA", "GAT") and r.gt == "1" and tuple(r.ad) == (0, 2, 1)


def test_example7_gt_tie_is_missing() -> None:
    aligns = peptides_at(GP, 50, sub(GP, 50, "D"), 1) + peptides_at(GP, 50, sub(GP, 50, "E"), 1)
    r = _one(records(aligns)[0])
    assert tuple(r.alts) == ("GAT", "GAA") and r.gt == "." and tuple(r.ad) == (0, 1, 1)


def test_example8_swap_six_nt_record(tmp_path: Path) -> None:
    assert SEQ_CHR[1048:1054] == "TCTGAT"
    recs, table, _, _ = records(peptides_at(GP, 60, two(GP, 60, "swap", "DS"), 1))
    r = _one(recs)
    assert (r.pos0, r.ref, tuple(r.alts)) == (1048, "TCTGAT", ("GATAGT",))
    write(tmp_path / "v.vcf", recs)
    (line,) = body_lines(tmp_path / "v.vcf")
    info = info_dict(line[7])
    assert (info["NT_CHANGES"], info["ALT_TYPE"], info["AA_ALT"]) == ("4", "swap", "DS")
    assert info["AA_REF"] == "SD"
    (row,) = site_row(table, "chr", 1048, "+")
    assert (row["alt_codon_vcf"], int(row["nt_changes"])) == ("GATAGT", 4)


def test_swap_on_minus_strand_starts_at_the_lower_codon() -> None:
    """Gene strand TCT GAT (codons 1537, 1534) -> GAT AGT; forward REF/ALT reverse-complemented."""
    assert rc(SEQ_CHR[1534:1540]) == "TCTGAT"
    recs, table, _, _ = records(peptides_at(GM_, 20, two(GM_, 20, "swap", "DS"), 1))
    r = _one(recs)
    assert (r.pos0, r.ref, tuple(r.alts)) == (1534, "ATCAGA", (rc("GATAGT"),))
    (row,) = site_row(table, "chr", 1537, "-")
    assert row["alt_codon_vcf"] == "ACTATC" and int(row["nt_changes"]) == 4


def test_two_by_two_block_enters_the_vcf() -> None:
    """GT -> AS at GGT ACT: G->A GCT (1), T->S AGT/TCT (1 each, tie -> AGT)."""
    assert minimal_codon("GGT", "A") == ("GCT", 1) and minimal_codon("ACT", "S") == ("AGT", 1)
    recs, table, _, _ = records(peptides_at(GP, 70, two(GP, 70, "block", "AS"), 2))
    r = _one(recs)
    assert (r.pos0, r.ref, tuple(r.alts), r.gt) == (1078, "GGTACT", ("GCTAGT",), "1")
    assert r.info["ALT_TYPE"] in ("block", ("block",), ["block"])
    (row,) = site_row(table, "chr", 1078, "+")
    assert (row["alt_codon_vcf"], int(row["nt_changes"])) == ("GCTAGT", 2)


def test_example9_one_to_two_block_not_in_vcf() -> None:
    op = ("block", (GP.codon(50),), "N", "GG")
    recs, table, _, _ = records(peptides_at(GP, 50, op, 1))
    assert recs == []
    (row,) = site_row(table, "chr", 1018, "+", "block")
    assert is_na(row["alt_codon_vcf"])


def test_two_to_one_block_not_in_vcf() -> None:
    op = ("block", (GP.codon(6), GP.codon(7)), aa_at("chr", GP.codon(6), "+")
          + aa_at("chr", GP.codon(7), "+"), "N")  # fmt: skip
    recs, table, _, _ = records(peptides_at(GP, 6, op, 1))
    assert recs == []
    (row,) = site_row(table, "chr", GP.codon(6), "+", "block")
    assert is_na(row["alt_codon_vcf"])


def test_frameshift_not_in_vcf() -> None:
    ev = al().FrameshiftEvent(4, GP.codon(91), GP.codon(92) + 1, 1)
    ops = [match_op("chr", GP.codon(i), "+") for i in range(88, 92)]
    recs, table, _, sites = records([raw_hit("chr", "+", ops, frameshift=ev)])
    assert (sites["alt_type"] == "frameshift").sum() == 1
    assert recs == []
    (row,) = [r for _, r in table.iterrows() if r["alt_type"] == "frameshift"]
    assert is_na(row["alt_codon_vcf"])


def _c2_sub(codon: int) -> list[Any]:
    """One peptide on ``c2`` ``+``: matches, a substitution at ``codon``, a match after it."""
    ref = aa_at("c2", codon, "+")
    alt = "W" if ref != "W" else "C"
    ops = [match_op("c2", codon - 3, "+"), ("substitution", (codon,), ref, alt),
           match_op("c2", codon + 3, "+")]  # fmt: skip
    return [raw_hit("c2", "+", ops)]


@pytest.mark.parametrize("codon", [998, 999])
def test_example10_origin_crossing_codon_flagged_not_written(codon: int) -> None:
    recs, table, _, _ = records(_c2_sub(codon))
    assert recs == []
    (row,) = site_row(table, "c2", codon, "+")
    assert is_na(row["alt_codon_vcf"])
    assert "vcf_unrepresentable" in flag_set(row["flags"])


def test_origin_crossing_keeps_m7_05_flags() -> None:
    _recs, table, _, _ = records(_c2_sub(999), {"flags": [("snv",)]})
    (row,) = site_row(table, "c2", 999, "+")
    assert flag_set(row["flags"]) == {"snv", "vcf_unrepresentable"}


def test_codon_ending_at_the_contig_end_is_representable(tmp_path: Path) -> None:
    recs, table, _, _ = records(_c2_sub(997))
    r = _one(recs)
    assert (r.contig, r.pos0, r.ref) == ("c2", 997, SEQ_C2[997:1000])
    (row,) = site_row(table, "c2", 997, "+")
    assert "vcf_unrepresentable" not in flag_set(row["flags"])
    write(tmp_path / "v.vcf", recs)
    assert body_lines(tmp_path / "v.vcf")[0][1] == "998"


def test_swap_whose_six_nt_cross_the_origin_is_unrepresentable() -> None:
    r1, r2 = aa_at("c2", 996, "+"), aa_at("c2", 999, "+")
    assert (r1, r2) == ("S", "D")
    ops = [match_op("c2", 993, "+"), ("swap", (996, 999), r1 + r2, r2 + r1),
           match_op("c2", 1002, "+")]  # fmt: skip
    recs, table, _, _ = records([raw_hit("c2", "+", ops)])
    assert recs == []
    (row,) = site_row(table, "c2", 996, "+", "swap")
    assert is_na(row["alt_codon_vcf"]) and "vcf_unrepresentable" in flag_set(row["flags"])


def test_representable_sites_are_not_flagged() -> None:
    _, table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 2))
    assert all("vcf_unrepresentable" not in flag_set(v) for v in table["flags"])


# --------------------------------------------------------------------------- rules 5-8


def test_substitution_and_swap_at_one_codon_are_two_records(tmp_path: Path) -> None:
    aligns = peptides_at(GP, 60, sub(GP, 60, "A"), 1) + peptides_at(
        GP, 60, two(GP, 60, "swap", "DS"), 1, shift=3
    )
    recs, *_ = records(aligns)
    assert [(r.pos0, len(r.ref)) for r in recs] == [(1048, 3), (1048, 6)]
    write(tmp_path / "v.vcf", recs)
    assert [ln[1] for ln in body_lines(tmp_path / "v.vcf")] == ["1049", "1049"]


def test_dp_counts_il_peptides() -> None:
    """Decision 4: DP = n_ref + n_il + sum n_alt (L codon CTT: REF, IL, two F)."""
    c = GP.codon(80)
    assert SEQ_CHR[c : c + 3] == "CTT"
    il = ("il", (c,), "L", "I")
    aligns = [
        make_hit(al(), "chr", "+", gene_peptide(GP, 79, 82, {}), length=2400),
        make_hit(al(), "chr", "+", gene_peptide(GP, 78, 81, {80: il}), length=2400),
        *peptides_at(GP, 80, sub(GP, 80, "F"), 2),
    ]
    recs, _, _pile, _ = records(aligns)
    r = _one(recs)
    assert (r.gt, tuple(r.ad), r.dp) == (".", (1, 2), 4)


def test_filter_is_pass_for_candidates() -> None:
    recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1))
    r = _one(recs)
    assert r.filter == "PASS"
    assert r.info["CONF"] in ("candidate", ("candidate",), ["candidate"])


def test_j_alt(tmp_path: Path) -> None:
    recs, _table, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "I"), 1))
    r = _one(recs)
    want, n = minimal_codon("GCT", "J")
    assert tuple(r.alts) == (want,)
    write(tmp_path / "v.vcf", recs)
    info = info_dict(body_lines(tmp_path / "v.vcf")[0][7])
    assert info["AA_ALT"] == "J" and info["NT_CHANGES"] == str(n)


def test_gene_and_aa_pos_omitted_off_gene(tmp_path: Path) -> None:
    ref = aa_at("chr", 200, "+")
    alt = "W" if ref != "W" else "C"
    ops = [match_op("chr", 197, "+"), ("substitution", (200,), ref, alt),
           match_op("chr", 203, "+")]  # fmt: skip
    recs, *_ = records([raw_hit("chr", "+", ops)])
    write(tmp_path / "v.vcf", recs)
    info = info_dict(body_lines(tmp_path / "v.vcf")[0][7])
    assert "GENE" not in info and "AA_POS" not in info
    assert info["AA_REF"] == ref and info["AA_ALT"] == alt


def test_flags_omitted_when_no_allele_has_one(tmp_path: Path) -> None:
    for extra in (None, {"flags": [()]}):
        recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), extra)
        write(tmp_path / "v.vcf", recs)
        assert "FLAGS" not in info_dict(body_lines(tmp_path / "v.vcf")[0][7])


def test_class_err_formats_a_rate(tmp_path: Path) -> None:
    recs, *_ = records(
        peptides_at(GP, 44, sub(GP, 44, "S"), 1), {"class_error_rate": [0.0123456789]}
    )
    write(tmp_path / "v.vcf", recs)
    assert info_dict(body_lines(tmp_path / "v.vcf")[0][7])["CLASS_ERR"] == "0.0123457"


def test_table_keeps_every_site_in_order() -> None:
    aligns = (
        peptides_at(GP, 44, sub(GP, 44, "S"), 2)
        + peptides_at(GP, 50, ("block", (GP.codon(50),), "N", "GG"), 1)
        + peptides_at(GM_, 10, sub(GM_, 10, "S"), 1)
        + _c2_sub(999)
    )
    _, table, _, sites = records(aligns)
    assert table["variant_id"].tolist() == sites["variant_id"].tolist()
    for col in ("nt_changes", "alt_codon_vcf", "flags"):
        assert col in table.columns


# --------------------------------------------------------------------------- degenerate cases


def test_no_sites_gives_no_records_and_a_valid_vcf(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    recs, table, _, _ = records([])
    assert recs == [] and len(table) == 0
    write(tmp_path / "v.vcf", recs)
    assert body_lines(tmp_path / "v.vcf") == []
    with pysam.VariantFile(str(tmp_path / "v.vcf")) as vf:
        assert list(vf) == [] and list(vf.header.samples) == [SAMPLE]


def test_ref_only_pileup_gives_no_records() -> None:
    recs, table, _, _ = records(ref_peptides(GP, 44, 2))
    assert recs == [] and len(table) == 0


def test_sample_auto_is_rejected(tmp_path: Path) -> None:
    recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1))
    with pytest.raises(ValueError, match="auto"):
        write(tmp_path / "v.vcf", recs, sample="auto")


def test_unsafe_contig_id_raises_in_the_writer(tmp_path: Path) -> None:
    v = vr()
    with pytest.raises(ValueError):
        vcf_io().write_vcf(
            tmp_path / "v.vcf", [], [Replicon("chr 1", 10, "linear")], gippa_version=VERSION,
            genome_sha256=SHA, sample=SAMPLE, info_fields=v.INFO_FIELDS,
            format_fields=v.FORMAT_FIELDS,
        )  # fmt: skip


# --------------------------------------------------------------------------- example 11


def _mixed_bag() -> list[Any]:
    return (
        peptides_at(GP, 44, sub(GP, 44, "S"), 3)
        + peptides_at(GP, 50, sub(GP, 50, "D"), 2)
        + peptides_at(GP, 50, sub(GP, 50, "E"), 1)
        + peptides_at(GP, 60, two(GP, 60, "swap", "DS"), 1)
        + peptides_at(GP, 70, two(GP, 70, "block", "AS"), 2)
        + peptides_at(GM_, 10, sub(GM_, 10, "S"), 1)
        + ref_peptides(GM_, 10, 1)
        + peptides_at(GM_, 20, two(GM_, 20, "swap", "DS"), 1)
        + _c2_sub(997)
        + _c2_sub(999)
    )


def test_example11_pysam_round_trip(tmp_path: Path, capfd: pytest.CaptureFixture[str]) -> None:
    pysam = pytest.importorskip("pysam")
    recs, _table, _, _ = records(_mixed_bag())
    path = tmp_path / "variants.vcf"
    write(path, recs)
    capfd.readouterr()
    with pysam.VariantFile(str(path)) as vf:
        hdr = vf.header
        assert list(hdr.samples) == [SAMPLE]
        assert {c: hdr.contigs[c].length for c in hdr.contigs} == {"chr": 2400, "c2": 1000}
        got = list(vf)
    err = capfd.readouterr().err
    assert err == "", err
    assert len(got) == len(recs) == 7
    by_pos = {(g.chrom, g.pos, len(g.ref)): g for g in got}
    g = by_pos[("chr", 1019, 3)]
    assert g.alts == ("GAT", "GAA")
    assert g.info["AA_ALT"] == ("D", "E") and g.info["CONF"] == ("high", "candidate")
    assert g.info["NT_CHANGES"] == (1, 2) and g.info["CALL"] == "ALT"
    assert g.info["GENE"] == "NISS_00012" and g.info["AA_POS"] == 51
    s = g.samples[SAMPLE]
    assert s["GT"] == (1,) and s["AD"] == (0, 2, 1) and s["DP"] == 3
    m = by_pos[("chr", 1568, 3)]
    assert m.ref == "AGC" and m.alts == ("AGA",) and m.samples[SAMPLE]["GT"] == (None,)
    assert m.info["CALL"] == "mixed" and m.samples[SAMPLE]["AD"] == (1, 1)
    w = by_pos[("chr", 1049, 6)]
    assert w.ref == "TCTGAT" and w.alts == ("GATAGT",) and w.info["NT_CHANGES"] == (4,)
    assert ("c2", 1000, 3) not in by_pos and ("c2", 998, 3) in by_pos
    assert all(list(x.filter) == ["PASS"] for x in got)
    assert all(x.info["CLASS_ERR"] == tuple("unestimated" for _ in x.alts) for x in got)
    text = path.read_text(encoding="utf-8")
    assert text.startswith("##fileformat=VCFv4.3\n")
    assert f"##source=gippa-{VERSION}\n" in text and f"##reference=sha256:{SHA}\n" in text


def test_every_written_ref_matches_the_genome(tmp_path: Path) -> None:
    """``write_vcf(..., sequences=genome)`` raises on a REF mismatch; also re-check here."""
    recs, *_ = records(_mixed_bag())
    write(tmp_path / "v.vcf", recs)
    for r in recs:
        seq = GENOME[r.contig]
        assert seq[r.pos0 : r.pos0 + len(r.ref)] == r.ref
        assert r.pos0 + len(r.ref) <= REPS[r.contig].length


def test_pileup_module_is_the_m7_02_one() -> None:
    """Guard: these tests build their inputs through the public M7-02 functions."""
    assert hasattr(pl(), "pileup") and hasattr(pl(), "sites")
