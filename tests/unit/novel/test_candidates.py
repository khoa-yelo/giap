"""Novel-ORF candidates, ORF reconstruction, grouping, IDs, second-round windows, OSTIR rows and the
pass-1 minimal call (issue M7-06; GIAP_SPEC.md §7.6.1, §7.6.2, §7.6.4 "Neither", §7.1.4, §6.5).

Examples 1-10, the degenerate cases and the acceptance criteria of
docs/specs/M7-06-novel-candidates-and-orfs.md, with the lead decisions of M7-index.md (M7-06
Q1-Q6 accepted: placed peptides are checked against the contaminant library by exact I = L
substring and C9 ``is_contaminant``; multi-locus candidate peptides are excluded; an ORF takes the
highest-priority M3-14 class; the ID is ``NOVEL_<contig>_<start1>_<p|m>`` from the lowest 1-based
coordinate of the longest interval; ``novel_candidate_id`` stays NA in pass 1).

Conventions these tests pin where the spec is silent (see also tests/helpers/novel_oracle.py):

* ``candidate_rows`` returns the retained C11 rows (every input column kept, possibly more);
  ``reconstruct`` takes its output, so the q-value join of rule 9 is internal to the module;
* ``Orf.stop``/``starts``/``start_longest`` are lowest 0-based codon coordinates wrapped into
  ``[0, L)``; ``starts`` are ordered 5'->3' on the ORF strand (descending on ``-``);
* ``Orf.frame`` equals the C11 ``frame`` of its peptides (checked only off the origin);
* ``best_q`` is ``None`` (not NaN) when no peptide has a q-value;
* ``promoter_windows`` returns windows unique by ``(contig, strand, end3p)``; ``Window.owner`` is
  the ``orf_id``; ``Window.k`` is relative to the start the window was built for;
* the minimal call is ``minimal_call(orfs, cfg) -> DataFrame`` with columns ``orf_id``,
  ``n_peptides``, ``called``, ``class_fdr``, ``flags`` (one row per ORF, input order; flags
  ``;``-joined alphabetically: ``class_fdr_unestimated`` always, ``no_start_codon`` without a
  start, ``partial`` when either end is open).
"""

from __future__ import annotations

import copy
import dataclasses

import pandas as pd
import pytest
from helpers.novel_oracle import (
    Case,
    Seq,
    cds,
    hit_row,
    il,
    mirror_case,
    mirror_codon,
    oracle_orfs,
    revcomp,
    translate,
    wrapped,
)

from giap.annotate.ostir import OstirRow
from giap.annotate.promoters import Window, gene_windows
from giap.config.model import GiapConfig
from giap.io.protein_fasta import LetterCounts, ProteinRecord
from giap.model.coords import Replicon
from giap.novel.candidates import (
    Orf,
    candidate_rows,
    minimal_call,
    ostir_rows,
    promoter_windows,
    reconstruct,
)

L = 1000
C1 = {"c1": Replicon("c1", L, "circular", "chromosome")}
LIN = {"lin": Replicon("lin", 300, "linear", "plasmid")}
ANTI_SD = "ACCTCCTTA"

# Example 1's coding codons from 72 (frame 0): P1 = [72, 102), P2 = [90, 111).
EX1_CODONS = [
    "GCT",
    "ATT",
    "GAA",
    "AAA",
    "CGT",
    "GAT",
    "TCT",
    "CCG",
    "CAT",
    "TGG",
    "TTC",
    "AAC",
    "CAG",
]
P1 = translate("".join(EX1_CODONS[:10]))  # AIEKRDSPHW (contains I: folded ALEKRDSPHW)
P2 = translate("".join(EX1_CODONS[6:]))  # SPHWFNQ
P3_CODONS = ["GAC", "TAC", "AAG", "CGC", "GTT", "GAG", "CAC", "TTT", "ACC", "ATC"]
P3 = translate("".join(P3_CODONS))  # DYKRVEHFTI


def cfg(**kw: object) -> GiapConfig:
    return GiapConfig(seed=1, **kw)  # type: ignore[arg-type]


def ex1_seq(*, atg: bool = True, gtg: bool = True, m_first: bool = False) -> Seq:
    s = Seq(L).fill(33, 120).put(30, "TAA").put(120, "TGA")
    if atg:
        s.put(42, "ATG")
    if gtg:
        s.put(60, "GTG")
    codons = list(EX1_CODONS)
    if m_first:
        codons[0] = "ATG"
    return s.put(72, "".join(codons))


def ex1_hits(seq: str, cls: str = "intergenic", gene: str | None = None) -> list[dict]:
    return [
        hit_row(translate(seq[72:102]), "c1", 72, 102, "+", cls, L, gene=gene),
        hit_row(translate(seq[90:111]), "c1", 90, 111, "+", cls, L, gene=gene),
    ]


def ex1(**kw: bool) -> Case:
    seq = ex1_seq(**kw).text()
    return Case({"c1": seq}, dict(C1), ex1_hits(seq))


def ex2() -> Case:
    s = ex1_seq().fill(133, 190).put(130, "TAA").put(139, "ATG").put(151, "".join(P3_CODONS))
    seq = s.put(190, "TAG").text()
    hits = [*ex1_hits(seq), hit_row(P3, "c1", 151, 181, "+", "intergenic", L)]
    return Case({"c1": seq}, dict(C1), hits)


def ex5(*, atg: bool = True) -> Case:
    s = Seq(300).fill(0, 60).put(6, "".join(EX1_CODONS[:10])).put(60, "TAA")
    if atg:
        s.put(3, "ATG")
    seq = s.text()
    return Case({"lin": seq}, dict(LIN), [hit_row(P1, "lin", 6, 36, "+", "intergenic", 300)])


def open3_case() -> Case:
    """Linear 300 nt: stop TAA [198, 201), ATG [210, 213), peptides [240, 270), [270, 288);
    no stop downstream before the contig end (last whole codon [297, 300))."""
    tail = ["TTC", "AAC", "CAG", "GAC", "TAC", "AAG"]
    s = Seq(300).put(198, "TAA").fill(201, 300).put(210, "ATG")
    seq = s.put(240, "".join(EX1_CODONS[:10])).put(270, "".join(tail)).text()
    hits = [
        hit_row(translate(seq[240:270]), "lin", 240, 270, "+", "intergenic", 300),
        hit_row(translate(seq[270:288]), "lin", 270, 288, "+", "intergenic", 300),
    ]
    return Case({"lin": seq}, dict(LIN), hits)


def ex6() -> Case:
    s = Seq(L).fill(243, 897).put(240, "TAA").put(249, "TTG").put(300, "ATG")
    seq = s.put(255, "".join(EX1_CODONS[:10])).put(897, "TAA").text()
    g = cds("G", "c1", 300, 900, "+")
    hit = hit_row(translate(seq[255:285]), "c1", 255, 285, "+", "nterm_extension", L, gene="G")
    return Case({"c1": seq}, dict(C1), [hit], genes=[g])


def wrap_case() -> Case:
    """+ ORF across the origin of c1: stop TAA [960, 963), ATG [972, 975), peptides [990, 1011)
    (crossing, unwrapped) and [14, 35) (= unwrapped [1014, 1035)), stop TGA at [50, 53)."""
    a = ["GCT", "ATT", "GAA", "AAA", "CGT", "GAT", "TCT"]
    b = ["CCG", "CAT", "TGG", "TTC", "AAC", "CAG", "GAC"]
    s = Seq(L).fill(963, 1050).put(960, "TAA").put(972, "ATG")
    seq = s.put(990, "".join(a)).put(1014, "".join(b)).put(1050, "TGA").text()
    hits = [
        hit_row(translate("".join(a)), "c1", 990, 1011, "+", "intergenic", L),
        hit_row(translate("".join(b)), "c1", 14, 35, "+", "intergenic", L),
    ]
    return Case({"c1": seq}, dict(C1), hits)


def run(case: Case, config: GiapConfig | None = None, contaminants=None) -> list[Orf]:
    rows = candidate_rows(case.hits_df(), case.peptides_df(), contaminants)
    return reconstruct(rows, case.gm, case.genome, config or cfg())


def by_stop(orfs: list[Orf]) -> dict[int, Orf]:
    return {o.stop: o for o in orfs}


def keratin(seq: str) -> ProteinRecord:
    return ProteinRecord(
        id="CON_KRT1_HUMAN",
        members=("CON_KRT1_HUMAN",),
        header="CON_KRT1_HUMAN Keratin, type II cytoskeletal 1",
        sequence=seq,
        internal_stops=0,
        nonstandard=LetterCounts(),
        other=LetterCounts(),
    )


# ---------------------------------------------------------------- the Orf record


def test_orf_is_a_frozen_dataclass_with_the_spec_fields() -> None:
    names = [f.name for f in dataclasses.fields(Orf)]
    assert names == [
        "orf_id", "cls", "contig", "strand", "frame", "stop", "starts", "start_longest",
        "end_5p_open", "end_3p_open", "peptides", "best_q", "overlapping_gene",
    ]  # fmt: skip
    (orf,) = run(ex1())
    with pytest.raises(dataclasses.FrozenInstanceError):
        orf.stop = 0  # type: ignore[misc]


# ---------------------------------------------------------------- examples 1-10


def test_example1_one_orf_two_starts_id_and_length() -> None:
    (orf,) = run(ex1())
    assert (orf.contig, orf.strand, orf.stop) == ("c1", "+", 120)
    assert orf.starts == (42, 60)
    assert orf.start_longest == 42
    assert orf.orf_id == "NOVEL_c1_43_p"
    assert (orf.stop - orf.start_longest) // 3 == 26  # length from 42, stop excluded
    assert orf.cls == "intergenic"
    assert orf.frame == 1
    assert orf.peptides == tuple(sorted({il(P1), il(P2)}))
    assert orf.peptides == ("ALEKRDSPHW", "SPHWFNQ")
    assert (orf.end_5p_open, orf.end_3p_open) == (False, False)
    assert orf.overlapping_gene is None
    call = minimal_call([orf], cfg())
    assert list(call.columns) == ["orf_id", "n_peptides", "called", "class_fdr", "flags"]
    assert call.to_dict("records") == [
        {"orf_id": "NOVEL_c1_43_p", "n_peptides": 2, "called": True,
         "class_fdr": "unestimated", "flags": "class_fdr_unestimated"}
    ]  # fmt: skip


def test_example2_other_frame_is_a_second_orf() -> None:
    orfs = run(ex2())
    assert sorted((o.contig, o.strand, o.stop) for o in orfs) == [
        ("c1", "+", 120),
        ("c1", "+", 190),
    ]
    o3 = by_stop(orfs)[190]
    assert o3.peptides == (il(P3),)
    assert o3.frame == 2
    assert o3.starts == (139,)
    assert o3.orf_id == "NOVEL_c1_140_p"
    assert by_stop(orfs)[120].peptides == ("ALEKRDSPHW", "SPHWFNQ")


def test_example3_peptide_also_cds_inframe_elsewhere_is_excluded() -> None:
    case = ex1()
    case.genes.append(cds("X", "c1", 498, 700, "+"))
    case.hits.append(hit_row(P2, "c1", 501, 522, "+", "cds_inframe", L, gene="X"))
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P1}
    (orf,) = reconstruct(rows, case.gm, case.genome, cfg())
    assert orf.peptides == (il(P1),)
    call = minimal_call([orf], cfg())
    assert call["n_peptides"].tolist() == [1]
    assert call["called"].tolist() == [False]
    assert call["class_fdr"].tolist() == ["unestimated"]


def test_example4_minus_strand_mirror() -> None:
    case = mirror_case(ex1())
    assert case.hits[0]["start"] == 898 and case.hits[0]["end"] == 928
    (orf,) = run(case)
    assert (orf.contig, orf.strand, orf.stop) == ("c1", "-", 877)
    assert orf.starts == (955, 937)  # 5'->3' on minus: ATG [955, 958), GTG [937, 940)
    assert orf.start_longest == 955
    assert orf.orf_id == "NOVEL_c1_878_m"  # the interval [877, 958): lowest base 878 (1-based)
    assert orf.frame == -1
    assert orf.peptides == ("ALEKRDSPHW", "SPHWFNQ")
    assert (orf.start_longest + 3 - orf.stop) // 3 - 1 == 26


def test_example5_linear_start_without_upstream_stop_is_open() -> None:
    (orf,) = run(ex5())
    assert orf.end_5p_open is True
    assert orf.end_3p_open is False
    assert orf.stop == 60
    assert orf.starts == (3,)  # codons 0, 3, 6 searched; only [3, 6) is a start
    assert orf.orf_id == "NOVEL_lin_4_p"
    (row,) = minimal_call([orf], cfg()).to_dict("records")
    assert row["flags"] == "class_fdr_unestimated;partial"
    assert not row["called"]  # one peptide


def test_example5_codon_at_contig_start_is_eligible() -> None:
    case = ex5(atg=False)
    seq = Seq(300, case.genome["lin"]).put(0, "GTG").text()
    case.genome["lin"] = seq
    (orf,) = run(case)
    assert orf.starts == (0,)
    assert orf.orf_id == "NOVEL_lin_1_p"


def test_example6_nterm_extension_uses_the_gene_stop() -> None:
    (orf,) = run(ex6())
    assert orf.cls == "nterm_extension"
    assert orf.stop == 897  # G's own stop codon [897, 900)
    assert orf.starts == (249,)  # G's ATG at 300 is downstream of the peptide: not a candidate
    assert orf.start_longest == 249
    assert orf.overlapping_gene == "G"
    assert orf.orf_id == "NOVEL_c1_250_p"


def test_example7_no_start_codon() -> None:
    (orf,) = run(ex1(atg=False, gtg=False))
    assert orf.starts == ()
    assert orf.start_longest is None
    assert orf.stop == 120
    assert orf.orf_id == "NOVEL_c1_73_p"  # interval from the first peptide
    (row,) = minimal_call([orf], cfg()).to_dict("records")
    assert row["flags"] == "class_fdr_unestimated;no_start_codon"
    assert row["n_peptides"] == 2


def test_example8_promoter_windows_union() -> None:
    case = ex1()
    orfs = run(case)
    wins = promoter_windows(orfs, case.gm, case.genome, cfg())
    assert all(isinstance(w, Window) for w in wins)
    keys = [(w.contig, w.strand, w.end3p) for w in wins]
    assert len(keys) == len(set(keys)) == 239  # deduplicated by (contig, strand, end3p)
    expected = {(s + k - 1) % L + 1 for s in (42, 60) for k in range(-200, 21)}
    assert {k[2] for k in keys} == expected
    assert {w.strand for w in wins} == {"+"}
    seq = case.genome["c1"]
    for w in wins:
        assert len(w.seq) == 40
        assert w.seq == wrapped(seq, w.end3p - 40, w.end3p)
        assert w.owner == "NOVEL_c1_43_p"


def test_example9_ostir_row_spans_both_starts() -> None:
    case = ex1()
    orfs = run(case)
    (row,) = ostir_rows(orfs, case.gm, case.genome, ANTI_SD, cfg())
    assert isinstance(row, OstirRow)
    seq = case.genome["c1"]
    assert row.id == "NOVEL_c1_43_p"
    assert row.seq == seq[992:] + seq[:110]  # genome[-8, 110) wrapped
    assert len(row.seq) == 118
    assert (row.start, row.end) == (51, 69)
    assert row.seq[row.start - 1 : row.start + 2] == "ATG"
    assert row.seq[row.end - 1 : row.end + 2] == "GTG"
    assert row.anti_sd == ANTI_SD


def test_example9_minus_strand_ostir_row() -> None:
    case = mirror_case(ex1())
    orfs = run(case)
    (row,) = ostir_rows(orfs, case.gm, case.genome, ANTI_SD, cfg())
    plus = ex1().genome["c1"]
    assert row.id == "NOVEL_c1_878_m"
    assert row.seq == plus[992:] + plus[:110]  # the reverse complement of the mirrored window
    assert (row.start, row.end) == (51, 69)


def test_example10_contaminant_substring_with_i_equals_l_is_excluded() -> None:
    case = ex1()
    lib = [keratin("MSRQFSSRSGYRSGGGFSSGSAGL" + "ALEKRDSPHW" + "GGSGGSY")]  # P1 has I here
    rows = candidate_rows(case.hits_df(), case.peptides_df(), lib)
    assert set(rows["peptide"]) == {P2}
    (orf,) = reconstruct(rows, case.gm, case.genome, cfg())
    assert orf.peptides == (P2,)
    assert minimal_call([orf], cfg())["called"].tolist() == [False]


# ---------------------------------------------------------------- exclusion rules (rule 1)


def test_exclusion_c9_is_contaminant() -> None:
    """Rule 1: C9 marks the peptide a contaminant (engine ID with a contaminant prefix)."""
    case = ex1()
    case.contaminant = (P1,)
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P2}


def test_exclusion_rna_overlap_row() -> None:
    """Rule 1: a peptide with an ``rna_overlap`` row never supports a novel ORF."""
    case = ex1()
    case.hits.append(hit_row(P1, "c1", 600, 630, "+", "rna_overlap", L))
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P2}


def test_exclusion_multi_candidate_locus() -> None:
    """Decision 2 (Q2): a peptide with two candidate loci is excluded from every ORF."""
    case = ex1()
    case.hits.append(hit_row(P2, "c1", 600, 621, "-", "antisense", L, uniqueness="shared"))
    case.hits[1]["uniqueness"] = "shared"
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P1}


def test_exclusion_non_candidate_classes_only() -> None:
    """Rows of ``cds_inframe``/``rna_overlap`` are never candidate rows themselves."""
    case = ex1()
    case.hits.append(hit_row("QQQQQQQK", "c1", 600, 624, "+", "cds_inframe", L, gene="X"))
    case.hits.append(hit_row("WWWWWWWK", "c1", 700, 724, "+", "rna_overlap", L))
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P1, P2}
    assert set(rows["locus_class"]) == {"intergenic"}


def test_contaminant_superstring_is_not_excluded() -> None:
    """Only peptides that are substrings of a library protein are excluded, not the reverse."""
    case = ex1()
    lib = [keratin(P2[:4])]
    rows = candidate_rows(case.hits_df(), case.peptides_df(), lib)
    assert set(rows["peptide"]) == {P1, P2}


def test_no_library_still_applies_is_contaminant() -> None:
    case = ex1()
    case.contaminant = (P2,)
    assert set(candidate_rows(case.hits_df(), case.peptides_df(), [])["peptide"]) == {P1}
    assert set(candidate_rows(case.hits_df(), case.peptides_df(), None)["peptide"]) == {P1}


@pytest.mark.parametrize(
    "cls", ["nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
            "pseudogene", "sorf"],
)  # fmt: skip
def test_every_candidate_class_is_retained(cls: str) -> None:
    case = ex1()
    for h in case.hits:
        h["locus_class"] = cls
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert set(rows["peptide"]) == {P1, P2}


def test_candidate_rows_keep_the_input_columns() -> None:
    case = ex1()
    hits = case.hits_df()
    rows = candidate_rows(hits, case.peptides_df(), None)
    assert set(hits.columns) <= set(rows.columns)
    got = rows[["peptide", "contig", "start", "end", "strand"]].sort_values("peptide")
    exp = hits[["peptide", "contig", "start", "end", "strand"]].sort_values("peptide")
    pd.testing.assert_frame_equal(got.reset_index(drop=True), exp.reset_index(drop=True))
    assert rows["novel_candidate_id"].isna().all()  # Q6: NA in pass 1


# ---------------------------------------------------------------- degenerate cases


def test_no_candidate_rows() -> None:
    case = ex1()
    for h in case.hits:
        h["locus_class"] = "cds_inframe"
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    assert len(rows) == 0
    orfs = reconstruct(rows, case.gm, case.genome, cfg())
    assert orfs == []
    assert promoter_windows(orfs, case.gm, case.genome, cfg()) == []
    assert ostir_rows(orfs, case.gm, case.genome, ANTI_SD, cfg()) == []
    assert len(minimal_call(orfs, cfg())) == 0


def test_empty_peptide_hits() -> None:
    case = Case(ex1().genome, dict(C1), [])
    assert run(case) == []


def test_in_frame_stop_between_two_peptides_gives_two_orfs() -> None:
    s = ex1_seq().put(102, "TAA")
    seq = s.text()
    hits = [
        hit_row(translate(seq[72:102]), "c1", 72, 102, "+", "intergenic", L),
        hit_row(translate(seq[105:120]), "c1", 105, 120, "+", "intergenic", L),
    ]
    orfs = run(Case({"c1": seq}, dict(C1), hits))
    assert sorted(o.stop for o in orfs) == [102, 120]
    a, b = by_stop(orfs)[102], by_stop(orfs)[120]
    assert a.starts == (42, 60) and a.orf_id == "NOVEL_c1_43_p"
    assert b.starts == () and b.orf_id == "NOVEL_c1_106_p"  # (102, 105] holds no start codon


def test_pseudogene_peptides_in_two_frames_give_two_orfs() -> None:
    s = Seq(L).fill(400, 700).put(402, "TAG").put(409, "TGA")
    s.put(420, "".join(EX1_CODONS[:10])).put(481, "".join(P3_CODONS))
    seq = s.put(600, "TAA").put(652, "TGA").text()
    ps = cds("PS1", "c1", 400, 700, "+", kind="pseudogene")
    hits = [
        hit_row(translate(seq[420:450]), "c1", 420, 450, "+", "pseudogene", L, gene="PS1"),
        hit_row(translate(seq[481:511]), "c1", 481, 511, "+", "pseudogene", L, gene="PS1"),
    ]
    case = Case({"c1": seq}, dict(C1), hits, genes=[ps])
    orfs = run(case)
    assert len(orfs) == 2
    assert len({o.stop for o in orfs}) == 2
    assert {o.cls for o in orfs} == {"pseudogene"}
    assert {o.overlapping_gene for o in orfs} == {"PS1"}
    expected = {(o.stop, o.starts) for o in oracle_orfs(hits, case.genome, case.replicons)}
    assert {(o.stop, o.starts) for o in orfs} == expected


def test_linear_contig_end_without_stop_is_3p_open() -> None:
    (orf,) = run(open3_case())
    assert orf.end_3p_open is True
    assert orf.end_5p_open is False
    assert orf.starts == (210,)
    assert orf.orf_id == "NOVEL_lin_211_p"
    assert len(orf.peptides) == 2  # both peptides of the open frame are one ORF
    (row,) = minimal_call([orf], cfg()).to_dict("records")
    assert row["called"]
    assert row["flags"] == "class_fdr_unestimated;partial"


def test_linear_minus_strand_open_ends() -> None:
    (orf,) = run(mirror_case(ex5()))
    assert orf.strand == "-"
    assert orf.end_5p_open is True and orf.end_3p_open is False
    assert orf.stop == 237  # mirror of the TAA [60, 63)
    assert orf.starts == (294,)  # CAT [294, 297) = ATG on minus
    assert orf.orf_id == "NOVEL_lin_238_m"
    (orf3,) = run(mirror_case(open3_case()))
    assert orf3.end_3p_open is True and orf3.end_5p_open is False
    assert orf3.starts == (mirror_codon(210, 300),)


def test_orf_crossing_the_origin_plus() -> None:
    (orf,) = run(wrap_case())
    assert (orf.strand, orf.stop) == ("+", 50)  # unwrapped [1050, 1053) wrapped
    assert orf.starts == (972,)
    assert orf.start_longest == 972
    assert orf.orf_id == "NOVEL_c1_973_p"
    assert len(orf.peptides) == 2  # both sides of the origin in one ORF


def test_orf_crossing_the_origin_minus() -> None:
    (orf,) = run(mirror_case(wrap_case()))
    assert (orf.strand, orf.stop) == ("-", 947)
    assert orf.starts == (25,)
    assert orf.orf_id == "NOVEL_c1_948_m"  # interval [947, 1028): lowest coordinate near L
    assert len(orf.peptides) == 2


def test_circular_frame_without_any_stop_raises() -> None:
    seq = "GCA" * 100
    rep = {"ring": Replicon("ring", 300, "circular")}
    hits = [hit_row(translate(seq[30:60]), "ring", 30, 60, "+", "intergenic", 300)]
    with pytest.raises(ValueError) as err:
        run(Case({"ring": seq}, rep, hits))
    assert "ring" in str(err.value)
    assert "frame" in str(err.value).lower()


def test_no_start_codon_gives_no_windows_and_no_ostir_row() -> None:
    case = ex1(atg=False, gtg=False)
    orfs = run(case)
    assert promoter_windows(orfs, case.gm, case.genome, cfg()) == []
    assert ostir_rows(orfs, case.gm, case.genome, ANTI_SD, cfg()) == []


def test_start_near_linear_contig_start_drops_windows_and_clips_ostir() -> None:
    case = ex5()
    orfs = run(case)
    assert promoter_windows(orfs, case.gm, case.genome, cfg()) == []  # every window runs off
    (row,) = ostir_rows(orfs, case.gm, case.genome, ANTI_SD, cfg())
    assert row.seq == case.genome["lin"][0:53]
    assert (row.start, row.end) == (4, 4)


def test_first_peptide_codon_is_an_eligible_start() -> None:
    (orf,) = run(ex1(m_first=True))
    assert orf.starts == (42, 60, 72)
    assert orf.start_longest == 42


# ---------------------------------------------------------------- class and overlapping gene


def test_class_is_the_highest_priority_among_the_peptides() -> None:
    case = ex1()
    case.genes.append(cds("C2", "c1", 103, 400, "+"))
    case.hits[1]["locus_class"] = "cds_outframe"  # P2 overlaps C2 in another frame
    (orf,) = run(case)
    assert orf.cls == "cds_outframe"  # cds_outframe > intergenic
    assert orf.overlapping_gene == "C2"


@pytest.mark.parametrize(
    ("cls", "kind", "gene"),
    [("sorf", "sorf", "S1"), ("pseudogene", "pseudogene", "PS1"), ("intergenic", None, None)],
)
def test_overlapping_gene_by_class(cls: str, kind: str | None, gene: str | None) -> None:
    seq = ex1_seq().text()
    genes = [] if kind is None else [cds(gene, "c1", 42, 123, "+", kind=kind)]  # type: ignore[arg-type]
    case = Case({"c1": seq}, dict(C1), ex1_hits(seq, cls, gene), genes=genes)
    (orf,) = run(case)
    assert orf.cls == cls
    assert orf.overlapping_gene == gene


# ---------------------------------------------------------------- peptides, q and the call


def test_il_variants_count_once() -> None:
    case = ex1()
    case.hits.append(hit_row(il(P1), "c1", 72, 102, "+", "intergenic", L, route="tolerant"))
    case.q = {P1: 0.003, il(P1): 0.0001, P2: 0.002}
    (orf,) = run(case)
    assert orf.peptides == ("ALEKRDSPHW", "SPHWFNQ")
    assert orf.best_q == pytest.approx(0.0001)
    assert minimal_call([orf], cfg())["n_peptides"].tolist() == [2]


def test_il_variants_alone_are_one_distinct_peptide() -> None:
    case = ex1()
    del case.hits[1]
    case.hits.append(hit_row(il(P1), "c1", 72, 102, "+", "intergenic", L, route="tolerant"))
    (orf,) = run(case)
    assert orf.peptides == ("ALEKRDSPHW",)
    call = minimal_call([orf], cfg())
    assert call["n_peptides"].tolist() == [1]
    assert call["called"].tolist() == [False]


def test_best_q_is_the_minimum_input_q() -> None:
    case = ex1()
    case.q = {P1: 0.004, P2: 0.0007}
    (orf,) = run(case)
    assert orf.best_q == pytest.approx(0.0007)


def test_best_q_is_none_for_de_novo_input() -> None:
    case = ex1()
    case.q = {P1: None, P2: None}
    (orf,) = run(case)
    assert orf.best_q is None


def test_min_peptides_is_read_from_config() -> None:
    (orf,) = run(ex1())
    strict = cfg(novel_orfs={"min_peptides_without_decoys": 3})
    assert minimal_call([orf], strict)["called"].tolist() == [False]
    lax = cfg(novel_orfs={"min_peptides_without_decoys": 1})
    (o3,) = [o for o in run(ex2()) if o.stop == 190]
    assert minimal_call([o3], lax)["called"].tolist() == [True]


def test_minimal_call_keeps_input_order() -> None:
    orfs = run(ex2())
    call = minimal_call(list(reversed(orfs)), cfg())
    assert call["orf_id"].tolist() == [o.orf_id for o in reversed(orfs)]
    assert set(call["class_fdr"]) == {"unestimated"}


# ---------------------------------------------------------------- config is honoured


def test_start_codons_from_config() -> None:
    (orf,) = run(ex1(), cfg(annotation={"start_codons": ["ATG"]}))
    assert orf.starts == (42,)
    (orf,) = run(ex1(), cfg(annotation={"start_codons": ["GTG"]}))
    assert orf.starts == (60,)
    assert orf.orf_id == "NOVEL_c1_61_p"


def test_window_range_from_config() -> None:
    case = ex6()
    c = cfg(annotation={"promoters": {"window_3p_range": [-10, 5]}})
    wins = promoter_windows(run(case, c), case.gm, case.genome, c)
    assert sorted(w.end3p for w in wins) == [249 + k for k in range(-10, 6)]
    assert sorted(w.k for w in wins) == list(range(-10, 6))


def test_ostir_flanks_from_config() -> None:
    case = ex1()
    c = cfg(annotation={"ostir": {"upstream_nt": 20, "downstream_nt": 10}})
    (row,) = ostir_rows(run(case, c), case.gm, case.genome, ANTI_SD, c)
    assert row.seq == case.genome["c1"][22:70]
    assert (row.start, row.end) == (21, 39)


# ---------------------------------------------------------------- one generalised window function


@pytest.mark.parametrize("mirror", [False, True])
def test_windows_equal_gene_windows_at_an_annotated_start(mirror: bool) -> None:
    case = mirror_case(ex6()) if mirror else ex6()
    (orf,) = run(case)
    (s,) = orf.starts
    g = case.gm.by_locus("G")
    assert (g.interval.start if not mirror else g.interval.end - 3) != s  # G's start differs
    gene = cds("ANN", "c1", s, 900, "+") if not mirror else cds("ANN", "c1", 100, s + 3, "-")
    expected = gene_windows(gene, case.replicons["c1"], case.genome["c1"], (-200, 20))
    got = promoter_windows([orf], case.gm, case.genome, cfg())

    def key(w: Window) -> tuple:
        return (w.contig, w.strand, w.end3p, w.seq, w.k)

    assert sorted(map(key, got)) == sorted(map(key, expected))
    assert {w.owner for w in got} == {orf.orf_id}


def test_minus_strand_window_sequences_are_reverse_complements() -> None:
    plus, minus = ex1(), mirror_case(ex1())
    wp = promoter_windows(run(plus), plus.gm, plus.genome, cfg())
    wm = promoter_windows(run(minus), minus.gm, minus.genome, cfg())
    assert sorted(w.seq for w in wp) == sorted(w.seq for w in wm)
    assert {w.strand for w in wm} == {"-"}
    mseq = minus.genome["c1"]
    for w in wm:  # 3'-most base at end3p (1-based) is the window's lowest base on minus
        assert w.seq == revcomp(wrapped(mseq, w.end3p - 1, w.end3p + 39))


# ---------------------------------------------------------------- determinism and purity


def test_inputs_are_not_mutated_and_results_are_deterministic() -> None:
    case = ex2()
    hits, peps = case.hits_df(), case.peptides_df()
    hits0, peps0 = hits.copy(deep=True), peps.copy(deep=True)
    genome0 = copy.deepcopy(case.genome)
    gm = case.gm
    rows = candidate_rows(hits, peps, None)
    rows0 = rows.copy(deep=True)
    orfs = reconstruct(rows, gm, case.genome, cfg())
    wins = promoter_windows(orfs, gm, case.genome, cfg())
    ost = ostir_rows(orfs, gm, case.genome, ANTI_SD, cfg())
    pd.testing.assert_frame_equal(hits, hits0)
    pd.testing.assert_frame_equal(peps, peps0)
    pd.testing.assert_frame_equal(rows, rows0)
    assert case.genome == genome0
    assert candidate_rows(hits, peps, None).equals(rows)
    assert reconstruct(rows, gm, case.genome, cfg()) == orfs
    assert promoter_windows(orfs, gm, case.genome, cfg()) == wins
    assert ostir_rows(orfs, gm, case.genome, ANTI_SD, cfg()) == ost


def test_result_does_not_depend_on_row_order() -> None:
    case = ex2()
    ref = run(case)
    hits = case.hits_df().iloc[::-1].reset_index(drop=True)
    peps = case.peptides_df().iloc[::-1].reset_index(drop=True)
    got = reconstruct(candidate_rows(hits, peps, None), case.gm, case.genome, cfg())
    assert got == ref


def test_orf_ids_are_unique() -> None:
    orfs = run(ex2())
    assert len({o.orf_id for o in orfs}) == len(orfs)
