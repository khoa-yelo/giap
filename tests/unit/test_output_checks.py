"""In-pipeline output checks ``check_tsv`` and ``check_vcf`` (M7-15 item 2; review #278 B2, lead
decision 11; GIPPA_SPEC.md §6.3, §9.5).

One hand-made file per branch. The cases kill the review's surviving mutants: POS >= 1, the
undeclared contig, required-NA and the enum check. The decision-11 nits are pinned too:
``check_tsv`` honours the schema ``pattern`` (full match), and ``check_vcf`` checks
POS <= the ``##contig`` length.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from gippa.io.output_checks import check_tsv, check_vcf

# ---------------------------------------------------------------- TSV (real schemas)

VARIANT_ROW: dict[str, str] = {
    "variant_id": "VAR_1",
    "contig": "chr1",
    "codon_start": "10",
    "strand": "+",
    "gene": "NA",
    "aa_pos": "4",
    "ref_codon": "ATG",
    "ref_aa": "M",
    "alt_aa": "I",
    "alt_type": "substitution",
    "call": "ALT",
    "confidence": "high",
    "n_ref": "0",
    "n_alt": "2",
    "nt_changes": "1",
    "alt_codon_vcf": "ATA",
    "class_error_rate": "0.01",
    "class_error_n": "5",
    "source": "database",
    "flags": "snv",
}

ORF_ROW: dict[str, str] = {
    "orf_id": "NOVEL_chr1_100_p",
    "class": "intergenic",
    "contig": "chr1",
    "strand": "+",
    "frame": "1",
    "start_longest": "10",
    "start_best_ostir": "10",
    "end": "100",
    "length_aa_longest": "30",
    "length_aa_best": "30",
    "start_codon_best": "ATG",
    "n_peptides": "2",
    "peptides": "PEPTIDEK;LLAMAR",
    "best_peptide_q": "0.001",
    "ostir_tir": "120.5",
    "promoter_score": "NA",
    "promoter_supported": "true",
    "terminator_downstream": "false",
    "terminator_conf": "NA",
    "terminator_distance_nt": "NA",
    "rbs_supported": "true",
    "class_fdr": "0.01",
    "called": "true",
    "rank": "1",
    "overlapping_gene": "NA",
    "flags": "rbs_supported",
}


def _tsv(path: Path, row: dict[str, str], **change: str) -> Path:
    cells = {**row, **change}
    text = "\t".join(cells) + "\n" + "\t".join(cells.values()) + "\n"
    if path.suffix == ".gz":
        path.write_bytes(gzip.compress(text.encode(), mtime=0))
    else:
        path.write_text(text, encoding="utf-8")
    return path


def _variants(tmp_path: Path, **change: str) -> Path:
    return _tsv(tmp_path / "variants.tsv", VARIANT_ROW, **change)


def _orfs(tmp_path: Path, **change: str) -> Path:
    return _tsv(tmp_path / "novel_orfs.tsv", ORF_ROW, **change)


def _fails(key: str, path: Path, stem: str, *needles: str) -> str:
    with pytest.raises(ValueError, match=r"^output ") as info:
        check_tsv(key, path, stem)
    msg = str(info.value)
    assert msg.startswith(f"output {key} ({path})"), msg
    for n in needles:
        assert n in msg, msg
    return msg


def test_valid_rows_pass(tmp_path: Path) -> None:
    check_tsv("variants", _variants(tmp_path), "variants")
    check_tsv("novel_orfs", _orfs(tmp_path), "novel_orfs")


def test_gzipped_tsv_is_read(tmp_path: Path) -> None:
    check_tsv("variants", _tsv(tmp_path / "v.tsv.gz", VARIANT_ROW), "variants")


def test_na_in_an_optional_column_passes(tmp_path: Path) -> None:
    check_tsv("variants", _variants(tmp_path, gene="NA", aa_pos="", flags="NA"), "variants")


def test_na_in_a_required_column_fails(tmp_path: Path) -> None:
    """Kills the 'required-NA check removed' mutant."""
    _fails("variants", _variants(tmp_path, contig="NA"), "variants", "line 2", "contig", "NA")


def test_empty_cell_in_a_required_column_fails(tmp_path: Path) -> None:
    _fails("variants", _variants(tmp_path, call=""), "variants", "call")


def test_required_column_with_field_missing_values_accepts_na(tmp_path: Path) -> None:
    # ref_aa declares its own missingValues [] (an empty residue is a deletion), so NA is not
    # a required-NA failure there.
    check_tsv("variants", _variants(tmp_path, ref_aa="NA"), "variants")


def test_value_outside_the_enum_fails(tmp_path: Path) -> None:
    """Kills the 'enum check removed' mutant (strand is a plain string otherwise)."""
    _fails("variants", _variants(tmp_path, strand="x"), "variants", "strand", "'x'")


def test_second_enum_column_fails(tmp_path: Path) -> None:
    _fails("variants", _variants(tmp_path, call="REF"), "variants", "call", "'REF'")


def test_value_below_the_minimum_fails(tmp_path: Path) -> None:
    _fails("variants", _variants(tmp_path, codon_start="0"), "variants", "codon_start", "minimum")


def test_minimum_holds_at_the_bound(tmp_path: Path) -> None:
    check_tsv("variants", _variants(tmp_path, codon_start="1", n_ref="0"), "variants")


def test_non_integer_fails(tmp_path: Path) -> None:
    _fails("variants", _variants(tmp_path, n_alt="1.5"), "variants", "n_alt", "integer")


def test_non_number_fails(tmp_path: Path) -> None:
    _fails("novel_orfs", _orfs(tmp_path, ostir_tir="high"), "novel_orfs", "ostir_tir", "number")


def test_bad_boolean_fails(tmp_path: Path) -> None:
    _fails("novel_orfs", _orfs(tmp_path, called="yes"), "novel_orfs", "called", "boolean")


def test_missing_schema_column_fails(tmp_path: Path) -> None:
    row = {k: v for k, v in VARIANT_ROW.items() if k != "confidence"}
    path = _tsv(tmp_path / "variants.tsv", row)
    _fails("variants", path, "variants", "missing", "confidence")


def test_ragged_row_fails_naming_the_output(tmp_path: Path) -> None:
    path = tmp_path / "variants.tsv"
    path.write_text("\t".join(VARIANT_ROW) + "\nVAR_1\tchr1\n", encoding="utf-8")
    _fails("variants", path, "variants")


def test_error_names_the_offending_line(tmp_path: Path) -> None:
    path = tmp_path / "variants.tsv"
    good = "\t".join(VARIANT_ROW.values())
    bad = "\t".join({**VARIANT_ROW, "strand": "?"}.values())
    path.write_text("\t".join(VARIANT_ROW) + f"\n{good}\n{good}\n{bad}\n", encoding="utf-8")
    _fails("variants", path, "variants", "line 4", "strand")


# -- schema pattern (decision 11 nit; review N3): fails until check_tsv honours `pattern`


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("variant_id", "VAR_0"),  # VAR_[1-9][0-9]*
        ("variant_id", "VAR_12x"),  # full match, not a prefix match
        ("ref_codon", "ATX"),  # [ACGTN]+
        ("alt_codon_vcf", "atg"),
        ("source", "database;"),
        ("flags", "snv;bogus"),
    ],
)
def test_value_not_matching_the_schema_pattern_fails(
    tmp_path: Path, column: str, value: str
) -> None:
    _fails("variants", _variants(tmp_path, **{column: value}), "variants", column, repr(value))


def test_orf_id_pattern_is_checked(tmp_path: Path) -> None:
    _fails("novel_orfs", _orfs(tmp_path, orf_id="ORF_1"), "novel_orfs", "orf_id")


def test_values_matching_the_pattern_pass(tmp_path: Path) -> None:
    path = _variants(
        tmp_path,
        variant_id="VAR_120",
        ref_codon="ATN",
        source="database;denovo",
        flags="snv;strain_difference",
        class_error_rate="unestimated",
    )
    check_tsv("variants", path, "variants")


# ---------------------------------------------------------------- VCF

HEADER = [
    "##fileformat=VCFv4.3",
    "##source=GIPPA",
    "##contig=<ID=chr1,length=100>",
    "##contig=<ID=plasmid_2,length=50>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample",
]


def _record(chrom: str = "chr1", pos: str = "10") -> str:
    return "\t".join([chrom, pos, "VAR_1", "A", "G", ".", "PASS", "AA=M", "GT", "1"])


def _vcf(tmp_path: Path, lines: list[str]) -> Path:
    path = tmp_path / "variants.vcf"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _vcf_fails(path: Path, *needles: str) -> str:
    with pytest.raises(ValueError, match=r"^output variants_vcf \(") as info:
        check_vcf("variants_vcf", path)
    msg = str(info.value)
    assert str(path) in msg, msg
    for n in needles:
        assert n in msg, msg
    return msg


def test_valid_vcf_passes(tmp_path: Path) -> None:
    check_vcf("variants_vcf", _vcf(tmp_path, [*HEADER, _record(), _record("plasmid_2", "50")]))


def test_header_only_vcf_passes(tmp_path: Path) -> None:
    check_vcf("variants_vcf", _vcf(tmp_path, HEADER))


def test_gzipped_vcf_is_read(tmp_path: Path) -> None:
    path = tmp_path / "variants.vcf.gz"
    path.write_bytes(gzip.compress(("\n".join([*HEADER, _record()]) + "\n").encode(), mtime=0))
    check_vcf("variants_vcf", path)


def test_non_vcf_first_line_fails(tmp_path: Path) -> None:
    _vcf_fails(_vcf(tmp_path, ["##source=GIPPA", *HEADER[2:], _record()]), "fileformat")


def test_empty_file_fails(tmp_path: Path) -> None:
    path = tmp_path / "variants.vcf"
    path.write_text("", encoding="utf-8")
    _vcf_fails(path, "fileformat")


def test_record_before_the_chrom_line_fails(tmp_path: Path) -> None:
    _vcf_fails(_vcf(tmp_path, [*HEADER[:4], _record(), HEADER[4]]), "line 5", "#CHROM")


def test_missing_chrom_line_fails(tmp_path: Path) -> None:
    _vcf_fails(_vcf(tmp_path, HEADER[:4]), "#CHROM")


def test_too_few_columns_fails(tmp_path: Path) -> None:
    short = "\t".join(_record().split("\t")[:7])
    _vcf_fails(_vcf(tmp_path, [*HEADER, short]), "line 6", "7 columns")


def test_header_with_fewer_than_eight_columns_fails(tmp_path: Path) -> None:
    head = "\t".join(HEADER[4].split("\t")[:7])
    rec = "\t".join(_record().split("\t")[:7])
    _vcf_fails(_vcf(tmp_path, [*HEADER[:4], head, rec]), "line 6")


def test_undeclared_contig_fails(tmp_path: Path) -> None:
    """Kills the 'undeclared-contig check removed' mutant: the record is otherwise valid."""
    _vcf_fails(_vcf(tmp_path, [*HEADER, _record("chr9", "10")]), "line 6", "'chr9'")


def test_pos_zero_fails(tmp_path: Path) -> None:
    """Kills the 'POS >= 1 check removed' mutant: the record has the full column count."""
    _vcf_fails(_vcf(tmp_path, [*HEADER, _record("chr1", "0")]), "line 6", "POS", "'0'")


@pytest.mark.parametrize("pos", ["-3", "1.5", "x", ""])
def test_non_positive_integer_pos_fails(tmp_path: Path, pos: str) -> None:
    _vcf_fails(_vcf(tmp_path, [*HEADER, _record("chr1", pos)]), "POS")


# -- POS <= contig length (decision 11 nit; review N3): fails until check_vcf reads `length=`


def test_pos_at_the_contig_length_passes(tmp_path: Path) -> None:
    check_vcf("variants_vcf", _vcf(tmp_path, [*HEADER, _record("chr1", "100")]))


@pytest.mark.parametrize(("chrom", "pos"), [("chr1", "101"), ("plasmid_2", "51")])
def test_pos_past_the_contig_length_fails(tmp_path: Path, chrom: str, pos: str) -> None:
    _vcf_fails(_vcf(tmp_path, [*HEADER, _record(chrom, pos)]), "line 6", pos)
