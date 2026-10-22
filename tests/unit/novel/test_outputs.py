"""Novel-ORF outputs: ``novel_orfs.tsv``, ``novel_orfs.gff3``, ``novel_orfs.faa`` and
``annotation_patch.gff3`` (issue M7-07; GIAP_SPEC.md §7.6.6, §9.1, §9.3, §6.5; ADR-0012
decision 3, ADR-0013 decision 3).

Examples 1-6 and 8 and every degenerate row of docs/specs/M7-07-novel-orf-outputs.md with its lead
decisions (minus-strand ``start_*``/``end`` in transcription direction; no OSTIR score -> NA
``start_best_ostir``/``ostir_tir``, ``*_best`` from the longest start, flag ``ostir_missing``;
FASTA suffixes ``|longest``/``|best_ostir``, merged when equal; no ``cterm_extension`` patch row).
Example 7 (``giap makedb --patch``) is tests/integration/test_patch_makedb_roundtrip.py.

Conventions pinned by the test author (see tests/helpers/novel_outputs_fixture.py):

* ``orf_table(orfs, scored, genome=..., replicons=...)``; ``scored`` is keyed by ``orf_id`` (the
  M7-06 ``minimal_call`` columns plus optional ``start_best_ostir``, a 0-based codon of
  ``Orf.starts``, and any other table column); the TSV is written by ``write_novel_orfs``;
* TSV cells follow giap.io.tsv: ``NA``, ``true``/``false``, ``%.6g``; ``peptides`` and ``flags``
  are ``;``-joined (flags in the schema vocabulary order, M7-09 decision 4);
* GFF3 attributes use the column names; the first three keys are ``ID``, ``class``, ``called``;
  ``start_longest`` is written like the TSV (1-based, transcription direction);
* FASTA records follow the table rows, ``|longest`` before ``|best_ostir``; every record starts
  with ``M``; ``length_aa_longest``/``length_aa_best`` equal the lengths of the proteins;
* without a start codon, ``length_aa_longest`` counts the codons from the first peptide to the
  stop (the M7-06 interval; the GFF3 feature needs it on ``-``, where the ID gives the stop side).
"""

from __future__ import annotations

import copy
import dataclasses
import random
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.novel_outputs_fixture import (
    NOVEL_ORFS_COLUMNS,
    PATCH_PRODUCT,
    Case,
    attr_keys,
    attrs_of,
    check_gff3,
    ex1,
    ex1_ex6,
    ex6,
    fasta_records,
    feature_lines,
    mirror_case,
    novel_io,
    open3_case,
    orfs_of,
    outputs,
    plus_minus,
    protein_from,
    proteins,
    rotate_case,
    scored,
    table,
    tsv_rows,
    two_contigs,
    wrap_case,
    write_all,
)

from giap.io.gff3 import GffFeature, read_gff3, write_gff3
from giap.io.protein_fasta import ProteinRecord

ID1, ID6, ID1M = "NOVEL_c1_43_p", "NOVEL_c1_250_p", "NOVEL_c1_878_m"
GFF_ATTRS = {
    "ID", "class", "called", "rank", "n_peptides", "start_longest", "start_codon",
    "promoter_supported", "terminator_downstream", "terminator_conf", "terminator_distance_nt",
    "Note", "partial",
}  # fmt: skip
NEW_ORF_CLASSES = ("intergenic", "cds_outframe", "antisense")


def _ex1_table(best: int | None = 60, **kw: dict[str, Any]) -> tuple[Case, list[Any], Any]:
    case = ex1()
    orfs = orfs_of(case)
    tir = {ID1: 1234.5} if best is not None else {}
    sc = scored(orfs, best={ID1: best}, ostir_tir=tir, **kw)
    return case, orfs, table(case, orfs, sc)


def _written(tmp_path: Path, case: Case, t: pd.DataFrame) -> dict[str, Path]:
    return write_all(tmp_path / "tables", case, t)


def _row(path: Path, orf_id: str) -> dict[str, str | None]:
    return next(r for r in tsv_rows(path) if r["orf_id"] == orf_id)


def _features(path: Path) -> dict[str, list[str]]:
    return {attrs_of(c[8])["ID"]: c for c in feature_lines(path)}


# ---------------------------------------------------------------- columns


def test_columns_constant_is_section_9_3_plus_adr_0012() -> None:
    cols = tuple(outputs().NOVEL_ORFS_COLUMNS)
    assert cols == NOVEL_ORFS_COLUMNS
    i = cols.index("promoter_score")
    assert cols[i + 1 : i + 5] == (
        "promoter_supported", "terminator_downstream", "terminator_conf", "terminator_distance_nt",
    )  # fmt: skip


def test_table_has_every_column_in_order_never_omitted() -> None:
    _, orfs, t = _ex1_table()
    assert tuple(t.columns) == NOVEL_ORFS_COLUMNS
    assert len(t) == len(orfs) == 1


# ---------------------------------------------------------------- examples 1-6


def test_example1_tsv_row(tmp_path: Path) -> None:
    case, _, t = _ex1_table()
    paths = _written(tmp_path, case, t)
    (row,) = tsv_rows(paths["novel_orfs.tsv"])
    assert row["orf_id"] == ID1
    assert row["class"] == "intergenic"
    assert (row["contig"], row["strand"], row["frame"]) == ("c1", "+", "1")
    assert (row["start_longest"], row["start_best_ostir"], row["end"]) == ("43", "61", "123")
    assert (row["length_aa_longest"], row["length_aa_best"]) == ("26", "20")
    assert row["start_codon_best"] == "GTG"
    assert row["n_peptides"] == "2"
    assert row["peptides"] == "ALEKRDSPHW;SPHWFNQ"
    assert row["best_peptide_q"] == "0.001"
    assert row["ostir_tir"] == "1234.5"
    assert row["class_fdr"] == "unestimated"
    assert row["called"] == "true"
    assert row["overlapping_gene"] is None
    assert row["flags"] == "class_fdr_unestimated"


def test_example2_fasta_two_start_variants(tmp_path: Path) -> None:
    case, _, t = _ex1_table()
    seq = case.genome["c1"]
    longest, best = protein_from(seq, 42, 120, "+"), protein_from(seq, 60, 120, "+")
    assert (len(longest), len(best)) == (26, 20)
    assert proteins(case, t) == {f"{ID1}|longest": longest, f"{ID1}|best_ostir": best}
    paths = _written(tmp_path, case, t)
    recs = fasta_records(paths["novel_orfs.faa"])
    assert recs == [(f"{ID1}|longest", longest), (f"{ID1}|best_ostir", best)]
    assert best[0] == "M" and seq[60:63] == "GTG"  # the GTG initiator is read as M


def test_protein_records_are_protein_records() -> None:
    case, _, t = _ex1_table()
    recs = outputs().protein_records(t, case.genome, case.replicons)
    assert all(isinstance(r, ProteinRecord) for r in recs)
    assert [r.header for r in recs] == [f"{ID1}|longest", f"{ID1}|best_ostir"]
    assert all("*" not in r.sequence and r.internal_stops == 0 for r in recs)


def test_example3_gff3_line(tmp_path: Path) -> None:
    case, _, t = _ex1_table()
    feats = outputs().gff_features(t)
    assert all(isinstance(f, GffFeature) for f in feats)
    paths = _written(tmp_path, case, t)
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[:8] == ["c1", "GIAP", "CDS", "61", "123", ".", "+", "0"]
    assert attr_keys(cols[8])[:3] == ["ID", "class", "called"]
    assert set(attr_keys(cols[8])) <= GFF_ATTRS
    a = attrs_of(cols[8])
    assert (a["ID"], a["class"], a["called"]) == (ID1, "intergenic", "true")
    assert (a["n_peptides"], a["start_longest"], a["start_codon"]) == ("2", "43", "GTG")
    assert a["Note"] == "intergenic"


def test_example4_patch_line(tmp_path: Path) -> None:
    case, _, t = _ex1_table()
    paths = _written(tmp_path, case, t)
    assert feature_lines(paths["annotation_patch.gff3"]) == [
        ["c1", "GIAP", "CDS", "61", "123", ".", "+", "0",
         f"ID={ID1};product={PATCH_PRODUCT}"],
    ]  # fmt: skip


def test_example5_nterm_extension_patch(tmp_path: Path) -> None:
    case = ex6()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, called={ID6: True}))
    paths = _written(tmp_path, case, t)
    assert feature_lines(paths["annotation_patch.gff3"]) == [
        ["c1", "GIAP", "CDS", "250", "900", ".", "+", "0", "ID=G;locus_tag=G;Note=nterm_extension"],
    ]  # fmt: skip
    row = _row(paths["novel_orfs.tsv"], ID6)
    assert (row["class"], row["overlapping_gene"]) == ("nterm_extension", "G")
    assert (row["start_longest"], row["end"], row["start_codon_best"]) == ("250", "900", "TTG")
    assert row["length_aa_longest"] == str((897 - 249) // 3)


def test_example5_mirror_nterm_patch_uses_the_gene_end(tmp_path: Path) -> None:
    case = mirror_case(ex6())  # G - [100, 700); TTG [748, 751)
    orfs = orfs_of(case)
    (orf,) = orfs
    t = table(case, orfs, scored(orfs, called={orf.orf_id: True}))
    paths = _written(tmp_path, case, t)
    assert feature_lines(paths["annotation_patch.gff3"]) == [
        ["c1", "GIAP", "CDS", "101", "751", ".", "-", "0", "ID=G;locus_tag=G;Note=nterm_extension"],
    ]  # fmt: skip
    row = _row(paths["novel_orfs.tsv"], orf.orf_id)
    assert (row["start_longest"], row["end"]) == ("751", "101")


def test_example6_minus_strand_transcription_direction(tmp_path: Path) -> None:
    case = mirror_case(ex1())
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, best={ID1M: None}))
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], ID1M)
    assert (row["strand"], row["frame"]) == ("-", "-1")
    assert (row["start_longest"], row["end"]) == ("958", "878")  # start > end on -
    assert row["start_best_ostir"] is None
    assert (row["length_aa_longest"], row["length_aa_best"]) == ("26", "26")
    assert row["start_codon_best"] == "ATG"
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[:8] == ["c1", "GIAP", "CDS", "878", "958", ".", "-", "0"]
    assert attrs_of(cols[8])["start_longest"] == "958"


def test_example6_minus_strand_best_ostir_start(tmp_path: Path) -> None:
    case = mirror_case(ex1())
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, best={ID1M: 937}))  # the mirrored GTG [937, 940)
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], ID1M)
    assert (row["start_longest"], row["start_best_ostir"], row["end"]) == ("958", "940", "878")
    assert (row["length_aa_best"], row["start_codon_best"]) == ("20", "GTG")
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:7] == ["878", "940", ".", "-"]
    plus = ex1()
    assert proteins(case, t) == {
        f"{ID1M}|longest": protein_from(plus.genome["c1"], 42, 120, "+"),
        f"{ID1M}|best_ostir": protein_from(plus.genome["c1"], 60, 120, "+"),
    }


# ---------------------------------------------------------------- example 8: GFF3 round trip


@pytest.mark.parametrize("name", ["novel_orfs.gff3", "annotation_patch.gff3"])
def test_example8_gff3_round_trip(tmp_path: Path, name: str) -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, best={ID1: 60}, called={ID6: True}))
    feats = outputs().gff_features(t) if name == "novel_orfs.gff3" else (
        outputs().patch_features(t, case.gm)
    )  # fmt: skip
    path = tmp_path / name
    write_gff3(path, feats, list(case.replicons.values()))
    back = read_gff3(path, case.replicons)

    def norm(f: GffFeature) -> tuple[object, ...]:
        attrs = {k: v if isinstance(v, str) else tuple(v) for k, v in f.attributes.items()}
        return (f.interval, f.source, f.type, f.score, f.phase, attrs)

    assert sorted(map(norm, back), key=repr) == sorted(map(norm, feats), key=repr)
    assert len(back) == 2


def test_both_gff3_files_pass_the_validity_check(tmp_path: Path) -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, best={ID1: 60}, called={ID6: True}))
    paths = _written(tmp_path, case, t)
    assert len(check_gff3(paths["novel_orfs.gff3"], case.replicons)) == 2
    assert len(check_gff3(paths["annotation_patch.gff3"], case.replicons)) == 2


# ---------------------------------------------------------------- rule 2: no OSTIR score


def test_no_ostir_score_falls_back_to_the_longest_start(tmp_path: Path) -> None:
    case, _, t = _ex1_table(best=None)
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], ID1)
    assert row["start_best_ostir"] is None and row["ostir_tir"] is None
    assert (row["length_aa_best"], row["start_codon_best"]) == ("26", "ATG")
    assert row["flags"] is not None
    assert set(row["flags"].split(";")) == {"class_fdr_unestimated", "ostir_missing"}
    assert row["flags"].split(";") == sorted(row["flags"].split(";"))
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:5] == ["43", "123"]  # longest start -> stop
    assert attrs_of(cols[8])["start_codon"] == "ATG"
    seq = case.genome["c1"]
    assert fasta_records(paths["novel_orfs.faa"]) == [
        (f"{ID1}|longest", protein_from(seq, 42, 120, "+"))
    ]
    assert feature_lines(paths["annotation_patch.gff3"])[0][3:5] == ["43", "123"]


def test_scored_none_writes_na_best_columns(tmp_path: Path) -> None:
    case = ex1()
    orfs = orfs_of(case)
    t = table(case, orfs, None)
    assert tuple(t.columns) == NOVEL_ORFS_COLUMNS and len(t) == 1
    novel_io().write_novel_orfs(tmp_path / "novel_orfs.tsv", t)
    (row,) = tsv_rows(tmp_path / "novel_orfs.tsv")
    assert row["start_best_ostir"] is None and row["ostir_tir"] is None
    assert row["start_longest"] == "43" and row["length_aa_best"] == "26"
    assert "ostir_missing" in (row["flags"] or "").split(";")


def test_best_equal_to_longest_gives_one_merged_record(tmp_path: Path) -> None:
    case, _, t = _ex1_table(best=42)
    paths = _written(tmp_path, case, t)
    seq = case.genome["c1"]
    assert fasta_records(paths["novel_orfs.faa"]) == [
        (f"{ID1}|longest,best_ostir", protein_from(seq, 42, 120, "+"))
    ]
    row = _row(paths["novel_orfs.tsv"], ID1)
    assert row["start_best_ostir"] == row["start_longest"] == "43"
    assert "ostir_missing" not in (row["flags"] or "")


# ---------------------------------------------------------------- rule 1 and 3: cells


def test_pass1_columns_are_na(tmp_path: Path) -> None:
    case, _, t = _ex1_table()
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], ID1)
    for col in ("promoter_score", "promoter_supported", "terminator_downstream",
                "terminator_conf", "terminator_distance_nt", "rbs_supported", "rank"):  # fmt: skip
        assert row[col] is None, col
    text = paths["novel_orfs.tsv"].read_text(encoding="utf-8")
    assert len(text.splitlines()[1].split("\t")) == len(NOVEL_ORFS_COLUMNS)
    assert "\tNA\t" in text


def test_scored_values_reach_tsv_and_gff3(tmp_path: Path) -> None:
    case, _, t = _ex1_table(
        promoter_score={ID1: 0.87},
        promoter_supported={ID1: True},
        terminator_downstream={ID1: False},
        terminator_conf={ID1: 0.25},
        terminator_distance_nt={ID1: 35},
        rbs_supported={ID1: True},
        rank={ID1: 1},
        class_fdr={ID1: 0.0123456789},
    )
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], ID1)
    assert (row["promoter_score"], row["promoter_supported"]) == ("0.87", "true")
    assert (row["terminator_downstream"], row["terminator_conf"]) == ("false", "0.25")
    assert (row["terminator_distance_nt"], row["rbs_supported"], row["rank"]) == ("35", "true", "1")
    assert row["class_fdr"] == "0.0123457"  # %.6g
    a = attrs_of(feature_lines(paths["novel_orfs.gff3"])[0][8])
    assert a["promoter_supported"] == "true"
    assert a["terminator_downstream"] == "false"
    assert a["terminator_conf"] == "0.25"
    assert a["terminator_distance_nt"] == "35"
    assert a["rank"] == "1"


def test_class_fdr_is_unestimated_literal_or_float(tmp_path: Path) -> None:
    case = two_contigs()
    orfs = orfs_of(case)
    ids = [o.orf_id for o in orfs]
    t = table(case, orfs, scored(orfs, class_fdr={ids[0]: 0.05, ids[1]: "unestimated",
                                                   ids[2]: 1e-7}))  # fmt: skip
    novel_io().write_novel_orfs(tmp_path / "novel_orfs.tsv", t)
    got = {r["orf_id"]: r["class_fdr"] for r in tsv_rows(tmp_path / "novel_orfs.tsv")}
    assert got == {ids[0]: "0.05", ids[1]: "unestimated", ids[2]: "1e-07"}


# ---------------------------------------------------------------- rule 4: order


def test_rows_follow_fasta_contig_order_then_coordinates_not_rank(tmp_path: Path) -> None:
    case = two_contigs()  # contig order zz, c1 (not alphabetical)
    orfs = orfs_of(case)
    assert [o.contig for o in orfs] == ["c1", "c1", "zz"]  # M7-06 sorts by name
    sc = scored(orfs, rank={"NOVEL_zz_43_p": 3, "NOVEL_c1_43_p": 2, "NOVEL_c1_140_p": 1})
    t = table(case, list(reversed(orfs)), sc.iloc[::-1].reset_index(drop=True))
    want = ["NOVEL_zz_43_p", "NOVEL_c1_43_p", "NOVEL_c1_140_p"]
    assert list(t["orf_id"]) == want
    paths = _written(tmp_path, case, t)
    assert [r["orf_id"] for r in tsv_rows(paths["novel_orfs.tsv"])] == want
    assert [r["rank"] for r in tsv_rows(paths["novel_orfs.tsv"])] == ["3", "2", "1"]
    heads = [h.split("|")[0] for h, _ in fasta_records(paths["novel_orfs.faa"])]
    assert list(dict.fromkeys(heads)) == want


@pytest.mark.parametrize("shift", [0, 100, 150, 900])
def test_mixed_strands_are_ordered_by_lowest_coordinate(shift: int) -> None:
    case = rotate_case(plus_minus(), shift)
    orfs = orfs_of(case)
    t = table(case, orfs, None)
    lowest = {o.orf_id: int(o.orf_id.split("_")[2]) for o in orfs}
    assert sorted(lowest) == sorted(t["orf_id"])
    assert list(t["orf_id"]) == sorted(lowest, key=lambda i: lowest[i])
    assert {o.strand for o in orfs} == {"+", "-"}


# ---------------------------------------------------------------- degenerate cases


def test_no_orfs(tmp_path: Path) -> None:
    case = ex1()
    out = outputs()
    t = table(case, [], None)
    assert tuple(t.columns) == NOVEL_ORFS_COLUMNS and len(t) == 0
    assert out.gff_features(t) == []
    assert out.patch_features(t, case.gm) == []
    assert out.protein_records(t, case.genome, case.replicons) == []
    paths = _written(tmp_path, case, t)
    lines = paths["novel_orfs.tsv"].read_text(encoding="utf-8").splitlines()
    assert lines == ["\t".join(NOVEL_ORFS_COLUMNS)]
    for name in ("novel_orfs.gff3", "annotation_patch.gff3"):
        text = paths[name].read_text(encoding="utf-8")
        assert all(ln.startswith("##") or "\tregion\t" in ln for ln in text.splitlines())
        assert feature_lines(paths[name]) == []
        check_gff3(paths[name], case.replicons)
    assert paths["novel_orfs.faa"].read_bytes() == b""


def test_no_called_orf(tmp_path: Path) -> None:
    case = two_contigs()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, called={o.orf_id: False for o in orfs}))
    paths = _written(tmp_path, case, t)
    assert feature_lines(paths["annotation_patch.gff3"]) == []
    assert len(tsv_rows(paths["novel_orfs.tsv"])) == 3
    assert {r["called"] for r in tsv_rows(paths["novel_orfs.tsv"])} == {"false"}
    assert len(feature_lines(paths["novel_orfs.gff3"])) == 3
    assert {attrs_of(c[8])["called"] for c in feature_lines(paths["novel_orfs.gff3"])} == {"false"}
    heads = {h.split("|")[0] for h, _ in fasta_records(paths["novel_orfs.faa"])}
    assert heads == {o.orf_id for o in orfs}


def test_orf_crossing_the_origin(tmp_path: Path) -> None:
    case = wrap_case()
    orfs = orfs_of(case)
    (orf,) = orfs
    t = table(case, orfs, scored(orfs, called={orf.orf_id: True}))
    paths = _written(tmp_path, case, t)
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:7] == ["973", "1053", ".", "+"]  # end > length (M1-05)
    assert feature_lines(paths["annotation_patch.gff3"])[0][3:5] == ["973", "1053"]
    row = _row(paths["novel_orfs.tsv"], orf.orf_id)
    assert (row["start_longest"], row["end"], row["length_aa_longest"]) == ("973", "53", "26")
    want = protein_from(case.genome["c1"], 972, 1050, "+")
    assert len(want) == 26
    assert fasta_records(paths["novel_orfs.faa"]) == [(f"{orf.orf_id}|longest", want)]
    check_gff3(paths["novel_orfs.gff3"], case.replicons)


def test_orf_open_at_a_linear_end(tmp_path: Path) -> None:
    case = open3_case()
    orfs = orfs_of(case)
    (orf,) = orfs
    assert orf.end_3p_open
    t = table(case, orfs, scored(orfs))
    paths = _written(tmp_path, case, t)
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:7] == ["211", "300", ".", "+"]  # to the last whole codon
    assert attrs_of(cols[8])["partial"] == "true"
    seq = case.genome["lin"]
    want = protein_from(seq, 210, 300, "+")
    assert len(want) == 30 and "*" not in want
    assert fasta_records(paths["novel_orfs.faa"]) == [(f"{orf.orf_id}|longest", want)]
    row = _row(paths["novel_orfs.tsv"], orf.orf_id)
    assert (row["start_longest"], row["end"], row["length_aa_longest"]) == ("211", "300", "30")
    assert "partial" in (row["flags"] or "").split(";")
    assert "\tregion\t" not in paths["novel_orfs.gff3"].read_text(encoding="utf-8")


def test_orf_without_a_start_codon(tmp_path: Path) -> None:
    case = ex1(atg=False, gtg=False)
    orfs = orfs_of(case)
    (orf,) = orfs
    assert orf.orf_id == "NOVEL_c1_73_p" and not orf.starts
    t = table(case, orfs, scored(orfs, called={orf.orf_id: True}))
    paths = _written(tmp_path, case, t)
    row = _row(paths["novel_orfs.tsv"], orf.orf_id)
    assert row["start_longest"] is None and row["start_best_ostir"] is None
    assert row["start_codon_best"] is None
    assert row["end"] == "123" and row["called"] == "true"
    assert row["length_aa_longest"] == "16"  # the interval from the first peptide (72) to the stop
    assert "no_start_codon" in (row["flags"] or "").split(";")
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:7] == ["73", "123", ".", "+"]  # from the first peptide
    assert attrs_of(cols[8])["start_codon"] == "NA"
    assert outputs().protein_records(t, case.genome, case.replicons) == []
    assert paths["novel_orfs.faa"].read_bytes() == b""
    assert feature_lines(paths["annotation_patch.gff3"]) == []  # never patched


def test_minus_orf_without_a_start_codon(tmp_path: Path) -> None:
    case = mirror_case(ex1(atg=False, gtg=False))  # first peptide [898, 928), stop [877, 880)
    orfs = orfs_of(case)
    (orf,) = orfs
    assert orf.orf_id == ID1M and not orf.starts
    t = table(case, orfs, scored(orfs))
    paths = _written(tmp_path, case, t)
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert cols[3:7] == ["878", "928", ".", "-"]
    row = _row(paths["novel_orfs.tsv"], ID1M)
    assert (row["start_longest"], row["end"], row["length_aa_longest"]) == (None, "878", "16")
    assert paths["novel_orfs.faa"].read_bytes() == b""


@pytest.mark.parametrize("cls", ["cterm_extension", "pseudogene", "sorf"])
def test_classes_without_patch_rows(tmp_path: Path, cls: str) -> None:
    case = ex6()
    (orf,) = orfs_of(case)
    o = dataclasses.replace(orf, cls=cls)
    t = table(case, [o], scored([o], called={o.orf_id: True}))
    assert outputs().patch_features(t, case.gm) == []
    paths = _written(tmp_path, case, t)
    assert _row(paths["novel_orfs.tsv"], o.orf_id)["class"] == cls
    (cols,) = feature_lines(paths["novel_orfs.gff3"])
    assert (attrs_of(cols[8])["class"], attrs_of(cols[8])["Note"]) == (cls, cls)


@pytest.mark.parametrize("cls", NEW_ORF_CLASSES)
def test_new_orf_classes_are_patched_as_cds(tmp_path: Path, cls: str) -> None:
    case, (orf,), _ = _ex1_table()
    o = dataclasses.replace(orf, cls=cls)
    t = table(case, [o], scored([o], best={ID1: 60}))
    (f,) = outputs().patch_features(t, case.gm)
    assert (f.type, f.source, f.phase) == ("CDS", "GIAP", 0)
    assert (f.interval.start, f.interval.end, f.interval.strand) == (60, 123, "+")
    assert dict(f.attributes) == {"ID": ID1, "product": PATCH_PRODUCT}


def test_uncalled_orfs_are_not_patched() -> None:
    case, orfs, _ = _ex1_table()
    t = table(case, orfs, scored(orfs, best={ID1: 60}, called={ID1: False}))
    assert outputs().patch_features(t, case.gm) == []
    assert len(outputs().gff_features(t)) == 1


def test_two_orfs_with_one_id_are_rejected() -> None:
    case, (orf,), _ = _ex1_table()
    twin = dataclasses.replace(orf, stop=orf.stop + 3)
    with pytest.raises((AssertionError, ValueError)):
        t = table(case, [orf, twin], scored([orf, twin]))
        outputs().patch_features(t, case.gm)


# ---------------------------------------------------------------- rule 5/6 invariants


def test_every_gff_feature_is_a_giap_cds_with_phase_0_and_note() -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs))
    for f in outputs().gff_features(t):
        assert (f.source, f.type, f.phase, f.score) == ("GIAP", "CDS", 0, None)
        assert f.attributes["Note"] == f.attributes["class"]
        assert f.interval.strand in ("+", "-")


def test_fasta_lengths_match_the_table(tmp_path: Path) -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, best={ID1: 60}))
    paths = _written(tmp_path, case, t)
    recs = dict(fasta_records(paths["novel_orfs.faa"]))
    for row in tsv_rows(paths["novel_orfs.tsv"]):
        oid = row["orf_id"]
        longest = recs.get(f"{oid}|longest") or recs[f"{oid}|longest,best_ostir"]
        assert str(len(longest)) == row["length_aa_longest"]
        best = recs.get(f"{oid}|best_ostir", longest)
        assert str(len(best)) == row["length_aa_best"]
        assert longest.endswith(best[1:])
    assert all(s[0] == "M" and "*" not in s for s in recs.values())


# ---------------------------------------------------------------- TE7 and purity


def _build(tmp: Path, shuffle: int | None = None) -> dict[str, bytes]:
    case = ex1_ex6()
    orfs = orfs_of(case)
    sc = scored(orfs, best={ID1: 60}, called={ID6: True}, ostir_tir={ID1: 1234.5})
    if shuffle is not None:
        random.Random(shuffle).shuffle(orfs)
        sc = sc.sample(frac=1.0, random_state=shuffle).reset_index(drop=True)
    paths = write_all(tmp, case, table(case, orfs, sc))
    return {name: p.read_bytes() for name, p in paths.items()}


def test_te7_four_files_byte_identical_across_reruns(tmp_path: Path) -> None:
    a, b = _build(tmp_path / "a"), _build(tmp_path / "b")
    assert a == b
    assert all(a[name] for name in a if name != "novel_orfs.faa") and a["novel_orfs.faa"]


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_outputs_do_not_depend_on_input_order(tmp_path: Path, seed: int) -> None:
    assert _build(tmp_path / "a") == _build(tmp_path / f"s{seed}", shuffle=seed)


def test_inputs_are_not_mutated() -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    sc = scored(orfs, best={ID1: 60}, called={ID6: True})
    orfs_before, sc_before = list(orfs), sc.copy(deep=True)
    genome_before, reps_before = copy.deepcopy(case.genome), dict(case.replicons)
    t = table(case, orfs, sc)
    assert orfs == orfs_before
    pd.testing.assert_frame_equal(sc, sc_before)
    t_before = t.copy(deep=True)
    out = outputs()
    gm = case.gm
    out.gff_features(t)
    out.patch_features(t, gm)
    out.protein_records(t, case.genome, case.replicons)
    pd.testing.assert_frame_equal(t, t_before)
    assert case.genome == genome_before and dict(case.replicons) == reps_before
