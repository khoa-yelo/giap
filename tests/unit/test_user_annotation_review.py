"""M2-11 review decisions (PR #64): multi-row CDS gaps, required contig order, surviving mutants.

docs/specs/M2-11-user-annotation-id-map.md, "Lead decisions after code review (PR #64, changes
requested)", items 1-3:

1. Multi-row GFF3 CDS sharing one ``ID`` merge only when consecutive parts (after unwrapping on
   circular contigs) overlap or leave a gap of at most 3 nt; a larger gap raises ``ValueError``
   naming the ID. On a linear contig, rows split at the origin raise ``ValueError``.
2. ``build_id_map`` requires the contig order (FASTA order, §6.5); no alphabetical fallback.
3. One test per surviving mutant of the review.

Convention pinned here: the contig order is the third parameter ``contig_order``, given as the
sequence of contig IDs in input FASTA order (the replicon/contig list), passed by keyword.
"""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from gippa.annotate.user_annotation import (
    UserCds,
    build_id_map,
    detect_format,
    read_user_cds,
    write_bakta_regions,
)
from gippa.model import Gene, Interval, Replicon

CHR = Replicon("chr1", 100_000, "circular", "chromosome")
PL = Replicon("plasmid_1", 8_000, "circular", "plasmid")
LIN = Replicon("linear_1", 5_000, "linear", "unknown")
MINI = Replicon("mini_1", 1_000, "circular", "plasmid")
REPS = {r.id: r for r in (CHR, PL, LIN, MINI)}
IDENT = {r: r for r in REPS}
ORDER = list(REPS)


def gene(tag: str, contig: str, start: int, end: int, strand: str, ftype: str = "cds") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type=ftype,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def row(seqid: str, s1: int, e1: int, strand: str, attrs: str) -> str:
    return "\t".join([seqid, "test", "CDS", str(s1), str(e1), ".", strand, "0", attrs])


def gff(path: Path, *rows: str) -> Path:
    path.write_text("\n".join(["##gff-version 3", *rows]) + "\n")
    return path


def read(path: Path) -> list[UserCds]:
    return read_user_cds(path, "gff3", REPS, IDENT)


def gbk(path: Path, *, locus: str, version: str, length: int, topology: str) -> Path:
    path.write_text(
        "\n".join(
            [
                f"LOCUS       {locus:<16}{length:>12} bp    DNA     {topology:<8} BCT 14-OCT-2026",
                "DEFINITION  test record.",
                f"ACCESSION   {version}",
                f"VERSION     {version}",
                "FEATURES             Location/Qualifiers",
                f"     source          1..{length}",
                '                     /mol_type="genomic DNA"',
                "     CDS             101..400",
                '                     /locus_tag="G1"',
                "ORIGIN",
                "//",
            ]
        )
        + "\n"
    )
    return path


# ============================================================ item 1: multi-row CDS gaps
def test_shared_id_far_apart_raises_naming_the_id(tmp_path: Path) -> None:
    """Two distinct genes that reuse ``ID=cds`` must not merge into one 5.2 kb span."""
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 400, "+", "ID=cds"),
        row("chr1", 5001, 5300, "+", "ID=cds"),
    )
    with pytest.raises(ValueError, match=r"\bcds\b"):
        read(p)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_programmed_frameshift_overlap_one_nt_merges(tmp_path: Path, strand: str) -> None:
    """prfB-like: two parts overlapping by 1 nt (75 nt + 1026 nt) are one CDS."""
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 175, strand, "ID=cds-prfB;locus_tag=PRFB"),
        row("chr1", 175, 1200, strand, "ID=cds-prfB;locus_tag=PRFB"),
    )
    assert read(p) == [UserCds("PRFB", Interval("chr1", 100, 1200, strand), False)]  # type: ignore[arg-type]


@pytest.mark.parametrize("gap", [0, 1, 2, 3])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_parts_with_gap_up_to_3nt_merge(tmp_path: Path, gap: int, strand: str) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 400, strand, "ID=c1"),
        row("chr1", 401 + gap, 700, strand, "ID=c1"),
    )
    assert read(p) == [UserCds("c1", Interval("chr1", 100, 700, strand), False)]  # type: ignore[arg-type]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_parts_with_gap_4nt_raise(tmp_path: Path, strand: str) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 400, strand, "ID=c1"),
        row("chr1", 405, 700, strand, "ID=c1"),
    )
    with pytest.raises(ValueError, match=r"\bc1\b"):
        read(p)


def test_gap_rule_applies_to_any_consecutive_pair(tmp_path: Path) -> None:
    """Three parts: the first pair touches, the second pair is 4 nt apart."""
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 400, "+", "ID=c1"),
        row("chr1", 401, 700, "+", "ID=c1"),
        row("chr1", 705, 900, "+", "ID=c1"),
    )
    with pytest.raises(ValueError, match=r"\bc1\b"):
        read(p)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_split_rows_on_linear_contig_raise(tmp_path: Path, strand: str) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("linear_1", 4901, 5000, strand, "ID=c1"),
        row("linear_1", 1, 150, strand, "ID=c1"),
    )
    with pytest.raises(ValueError):
        read(p)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_origin_split_rows_on_circular_contig_cross_origin(tmp_path: Path, strand: str) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 1, 150, strand, "ID=c1"),
        row("chr1", 99901, 100000, strand, "ID=c1"),
    )
    assert read(p) == [UserCds("c1", Interval("chr1", 99900, 100150, strand), False)]  # type: ignore[arg-type]


def test_origin_split_gap_counted_across_the_origin(tmp_path: Path) -> None:
    """After unwrapping, a 3-nt gap at the origin merges and a 4-nt gap raises."""
    ok = gff(
        tmp_path / "ok.gff3",
        row("chr1", 99901, 99997, "+", "ID=c1"),
        row("chr1", 1, 150, "+", "ID=c1"),
    )
    assert read(ok) == [UserCds("c1", Interval("chr1", 99900, 100150, "+"), False)]
    bad = gff(
        tmp_path / "bad.gff3",
        row("chr1", 99901, 99996, "+", "ID=c1"),
        row("chr1", 1, 150, "+", "ID=c1"),
    )
    with pytest.raises(ValueError, match=r"\bc1\b"):
        read(bad)


# ============================================================ item 2: required contig order
def test_build_id_map_without_contig_order_is_type_error() -> None:
    cds = [UserCds("U1", Interval("chr1", 100, 400, "+"), False)]
    with pytest.raises(TypeError):
        build_id_map(cds, [gene("G1", "chr1", 100, 400, "+")])  # type: ignore[call-arg]


def test_build_id_map_contig_order_parameter_has_no_default() -> None:
    param = inspect.signature(build_id_map).parameters["contig_order"]
    assert param.default is inspect.Parameter.empty


def test_build_id_map_follows_fasta_order_not_alphabetical() -> None:
    cds = [
        UserCds("A_CHR", Interval("chr1", 100, 400, "+"), False),
        UserCds("B_PL", Interval("plasmid_1", 100, 400, "+"), False),
        UserCds("C_LIN", Interval("linear_1", 100, 400, "+"), False),
    ]
    rows = build_id_map(cds, [], contig_order=["plasmid_1", "linear_1", "chr1"])
    assert [r.source_id for r in rows] == ["B_PL", "C_LIN", "A_CHR"]
    rows = build_id_map(cds[::-1], [], contig_order=["chr1", "plasmid_1", "linear_1"])
    assert [r.source_id for r in rows] == ["A_CHR", "B_PL", "C_LIN"]


# ============================================================ item 3: surviving mutants
def test_circular_guard_linear_parts_touching_both_ends_not_unwrapped(tmp_path: Path) -> None:
    """Mutant: drop ``rep.topology == "circular"``. On a linear contig parts are never unwrapped,
    even when the gap at the contig ends (0 nt) is smaller than the gap between parts (2 nt)."""
    p = gff(
        tmp_path / "u.gff3",
        row("linear_1", 1, 2000, "+", "ID=c1"),
        row("linear_1", 2003, 5000, "+", "ID=c1"),
    )
    assert read(p) == [UserCds("c1", Interval("linear_1", 0, 5000, "+"), False)]


@pytest.mark.parametrize("other", ["sorf", "pseudogene"])
def test_id_map_prefers_cds_over_other_gene_types(other: str) -> None:
    """Mutant: index genes by locus tag only. A sORF/pseudogene with an earlier tag on the same
    coordinates must not win over the CDS."""
    cds = [UserCds("U1", Interval("chr1", 100, 400, "+"), False)]
    genes = [gene("ZZZ_CDS", "chr1", 100, 400, "+"), gene("AAA_X", "chr1", 100, 400, "+", other)]
    for g in (genes, genes[::-1]):
        (r,) = build_id_map(cds, g, contig_order=ORDER)
        assert (r.locus_tag, r.match) == ("ZZZ_CDS", "identical")


@pytest.mark.parametrize("which", [0, 1])
def test_pseudo_flag_from_any_row(tmp_path: Path, which: int) -> None:
    """Mutant: take the pseudo flag from the first row only."""
    attrs = ["ID=c1", "ID=c1"]
    attrs[which] += ";pseudo=true"
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 175, "+", attrs[0]),
        row("chr1", 175, 1200, "+", attrs[1]),
    )
    (c,) = read(p)
    assert c.pseudo is True


def test_cds_rows_on_both_strands_raise(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", 101, 175, "+", "ID=c1"),
        row("chr1", 175, 1200, "-", "ID=c1"),
    )
    with pytest.raises(ValueError, match=r"\bc1\b"):
        read(p)


def test_gbff_locus_name_preferred_over_version_id(tmp_path: Path) -> None:
    """LOCUS name is a known contig, VERSION id is not: the LOCUS name is used."""
    p = gbk(tmp_path / "a.gbff", locus="chr1", version="NZ_XYZ01.1", length=100_000,
            topology="circular")  # fmt: skip
    (c,) = read_user_cds(p, "gbff", REPS, IDENT)
    assert c.interval == Interval("chr1", 100, 400, "+")


def test_gbff_falls_back_to_version_id(tmp_path: Path) -> None:
    """LOCUS name is unknown (e.g. truncated to 16 chars), VERSION id is the contig."""
    p = gbk(tmp_path / "a.gbff", locus="unknown_locus", version="plasmid_1", length=8_000,
            topology="circular")  # fmt: skip
    (c,) = read_user_cds(p, "gbff", REPS, IDENT)
    assert c.interval == Interval("plasmid_1", 100, 400, "+")


def test_gbff_version_id_through_contig_map(tmp_path: Path) -> None:
    p = gbk(tmp_path / "a.gbff", locus="unknown_locus", version="NZ_AB1.1", length=8_000,
            topology="circular")  # fmt: skip
    (c,) = read_user_cds(p, "gbff", REPS, {**IDENT, "NZ_AB1.1": "plasmid_1"})
    assert c.interval.contig == "plasmid_1"


def test_nested_part_does_not_shorten_reach(tmp_path: Path) -> None:
    """Mutant: ``reach = part.end`` instead of ``max``. A part nested in the first part must not
    open a fake 300-nt gap (larger than the 20-nt gap at the origin of a 1 kb circle)."""
    p = gff(
        tmp_path / "u.gff3",
        row("mini_1", 11, 400, "+", "ID=c1"),
        row("mini_1", 51, 100, "+", "ID=c1"),
        row("mini_1", 401, 990, "+", "ID=c1"),
    )
    assert read(p) == [UserCds("c1", Interval("mini_1", 10, 990, "+"), False)]


@pytest.mark.parametrize(
    ("body", "fmt"),
    [
        ("##gff-version 3\n", "gff3"),
        ("LOCUS       chr1\n", "gbff"),
        ('{"genome": {}}\n', "bakta_json"),
    ],
)
def test_detect_format_strips_utf8_bom(tmp_path: Path, body: str, fmt: str) -> None:
    p = tmp_path / "ann.txt"
    p.write_bytes(b"\xef\xbb\xbf" + body.encode())
    assert detect_format(p) == fmt


@pytest.mark.parametrize(
    "iv",
    [
        Interval("linear_1", 4900, 5100, "+"),  # linear contigs never wrap
        Interval("plasmid_1", 8000, 8300, "+"),  # starts beyond the replicon
        Interval("plasmid_1", 100, 8200, "-"),  # longer than the circular replicon
    ],
)
def test_regions_writer_checks_intervals(tmp_path: Path, iv: Interval) -> None:
    with pytest.raises(ValueError):
        write_bakta_regions(tmp_path / "r.gff3", [UserCds("X", iv, False)], list(REPS.values()))
