"""``--annotation`` GFF3/GBFF: Bakta ``--regions`` input and ``annotation/id_map.tsv`` (M2-11).

docs/specs/M2-11-user-annotation-id-map.md, examples 1-8, edge cases and the simulator round trip;
GIAP_SPEC.md §5.3 ("Existing annotation"), §7.1.1 (Bakta ``--regions`` takes CDS only), §7.1.8
(``annotation/id_map.tsv``), §6.5 (0-based half-open internally, 1-based closed in writers,
origin-crossing intervals stored unwrapped with ``end > length``).

Conventions pinned here (ambiguities resolved by the test author, spec text first):

* ``id_map.tsv`` header is exactly the 7 columns of the spec (lead Q1 default), 1-based closed,
  ``NA`` for a missing locus tag; an origin-crossing row keeps ``end > length``.
* GFF3 ``source_id`` fallback ``<contig>_<start1>_<strand>`` is only checked for its
  ``<contig>_<start1>_`` prefix and for being strand-specific (the strand token is not pinned).
* Bakta's own pseudogene marker (``pseudogene=<type>`` in GFF3, ``/pseudogene="<type>"`` in GBFF;
  INSDC: ``/pseudogene`` implies ``/pseudo``) counts as pseudo, like ``pseudo=true`` and ``/pseudo``
  (real K-12 recording, ``KOIAIB_00745``).
* ``detect_format`` raises ``ValueError`` on content that is none of JSON, GFF3 or GBFF.
* Real Bakta 1.12.2 full-DB outputs (tests/data/user_annotation/PROVENANCE.md) are realistic
  "other annotator" inputs whose locus tags and coordinates are known.
"""

from __future__ import annotations

import csv
import gzip
import random
import re
from collections.abc import Mapping
from pathlib import Path

import pytest

from giap.annotate.user_annotation import (
    IdMapRow,
    UserCds,
    build_id_map,
    detect_format,
    read_user_cds,
    write_bakta_regions,
    write_id_map,
)
from giap.io.gff3 import read_gff3
from giap.model import Gene, Interval, Replicon

DATA = Path(__file__).resolve().parents[1] / "data" / "user_annotation"

CHR = Replicon("chr1", 100_000, "circular", "chromosome")
PL = Replicon("plasmid_1", 8_000, "circular", "plasmid")
LIN = Replicon("linear_1", 5_000, "linear", "unknown")
REPS = {r.id: r for r in (CHR, PL, LIN)}
IDENT = {r: r for r in REPS}
ORDER = list(REPS)  # input FASTA order (build_id_map requires it, review decision 2)
ID_MAP_HEADER = ["source_id", "locus_tag", "contig", "start", "end", "strand", "match"]


# --------------------------------------------------------------- helpers
def gene(tag: str, contig: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def gff(path: Path, *rows: str, pragmas: tuple[str, ...] = ()) -> Path:
    """GFF3 file from tab-separated feature rows (``|`` in a row is NOT a column separator)."""
    path.write_text("\n".join(["##gff-version 3", *pragmas, *rows]) + "\n")
    return path


def row(seqid: str, ftype: str, s1: int, e1: int, strand: str, attrs: str) -> str:
    phase = "0" if ftype == "CDS" else "."
    return "\t".join([seqid, "test", ftype, str(s1), str(e1), ".", strand, phase, attrs])


def gbk(
    path: Path,
    features: list[tuple[str, list[str]]],
    *,
    name: str = "chr1",
    length: int = 100_000,
    topology: str = "circular",
) -> Path:
    """Minimal GenBank record: (location, qualifier lines) CDS features, empty ORIGIN."""
    lines = [
        f"LOCUS       {name:<16}{length:>12} bp    DNA     {topology:<8} BCT 14-OCT-2026",
        "DEFINITION  test record.",
        f"ACCESSION   {name}",
        f"VERSION     {name}",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{length}",
        '                     /mol_type="genomic DNA"',
    ]
    for loc, quals in features:
        lines.append(f"     CDS             {loc}")
        lines.extend(f"                     {q}" for q in quals)
    lines += ["ORIGIN", "//"]
    path.write_text("\n".join(lines) + "\n")
    return path


def tsv_rows(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text().splitlines()]


def feature_lines(path: Path) -> list[list[str]]:
    return [ln.split("\t") for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]


def bakta_tsv_genes(path: Path) -> list[Gene]:
    """CDS rows of a real Bakta TSV as Genes (1-based closed -> 0-based half-open)."""
    out = []
    with path.open() as fh:
        for cells in csv.reader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            if cells[1] == "cds":
                out.append(gene(cells[5], cells[0], int(cells[2]) - 1, int(cells[3]), cells[4]))
    return out


SIM_REPS = {
    "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
    "plasmid_1": Replicon("plasmid_1", 8_000, "circular", "plasmid"),
}
K12 = Replicon("U00096.3", 4_641_652, "linear", "chromosome")


# --------------------------------------------------------------- example 1
def ex1_gff(tmp_path: Path) -> Path:
    return gff(
        tmp_path / "user.gff3",
        row("chr1", "gene", 101, 400, "+", "ID=gene-1;locus_tag=ABC_0001"),
        row("chr1", "mRNA", 101, 400, "+", "ID=rna-1;Parent=gene-1"),
        row("chr1", "CDS", 101, 400, "+", "ID=cds-1;Parent=rna-1;locus_tag=ABC_0001"),
        pragmas=("##sequence-region chr1 1 100000",),
    )


def test_example1_gff3_reads_only_cds(tmp_path: Path) -> None:
    got = read_user_cds(ex1_gff(tmp_path), "gff3", REPS, IDENT)
    assert got == [UserCds("ABC_0001", Interval("chr1", 100, 400, "+"), False)]


def test_example1_regions_file(tmp_path: Path) -> None:
    cds = read_user_cds(ex1_gff(tmp_path), "gff3", REPS, IDENT)
    out = tmp_path / "regions.gff3"
    write_bakta_regions(out, cds, [CHR])
    lines = out.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert "##sequence-region chr1 1 100000" in lines
    feats = feature_lines(out)
    assert len(feats) == 1
    f = feats[0]
    assert len(f) == 9
    assert (f[0], f[2:]) == ("chr1", ["CDS", "101", "400", ".", "+", "0", "ID=ABC_0001"])


def test_regions_file_has_sequence_region_per_replicon_and_cds_only(tmp_path: Path) -> None:
    cds = [
        UserCds("B", Interval("plasmid_1", 10, 310, "-"), False),
        UserCds("A", Interval("chr1", 500, 800, "+"), False),
        UserCds("C", Interval("chr1", 100, 400, "-"), True),
    ]
    out = tmp_path / "regions.gff3"
    write_bakta_regions(out, cds, [CHR, PL])
    text = out.read_text()
    assert "##sequence-region chr1 1 100000" in text
    assert "##sequence-region plasmid_1 1 8000" in text
    feats = feature_lines(out)
    assert {f[2] for f in feats} == {"CDS"}
    # sorted per §6.5: contig in replicon order, then start
    assert [(f[0], f[3]) for f in feats] == [("chr1", "101"), ("chr1", "501"), ("plasmid_1", "11")]
    ids = [re.match(r"ID=([^;]+)", f[8]).group(1) for f in feats]  # type: ignore[union-attr]
    assert ids == ["C", "A", "B"]


def test_regions_file_reads_back_with_m1_05_reader(tmp_path: Path) -> None:
    cds = [
        UserCds("A", Interval("chr1", 100, 400, "+"), False),
        UserCds("B", Interval("plasmid_1", 7000, 7600, "-"), False),
        UserCds("L", Interval("linear_1", 0, 300, "+"), False),
    ]
    out = tmp_path / "regions.gff3"
    write_bakta_regions(out, cds, [CHR, PL, LIN])
    back = [f for f in read_gff3(out, REPS) if f.type == "CDS"]
    assert sorted((str(f.attributes["ID"]), f.interval) for f in back) == [
        (c.source_id, c.interval) for c in cds
    ]


# --------------------------------------------------------------- examples 2-3
def test_example2_gbff_complement(tmp_path: Path) -> None:
    p = gbk(tmp_path / "a.gbff", [("complement(1001..1300)", ['/locus_tag="ABC_0002"'])])
    got = read_user_cds(p, "gbff", REPS, IDENT)
    assert got == [UserCds("ABC_0002", Interval("chr1", 1000, 1300, "-"), False)]


def test_example3_gbff_join_across_origin(tmp_path: Path) -> None:
    p = gbk(tmp_path / "a.gbff", [("join(99901..100000,1..150)", ['/locus_tag="ABC_0003"'])])
    (cds,) = read_user_cds(p, "gbff", REPS, IDENT)
    assert cds.interval == Interval("chr1", 99900, 100150, "+")
    assert cds.interval.crosses_origin(CHR)


def test_gbff_complement_join_across_origin(tmp_path: Path) -> None:
    p = gbk(
        tmp_path / "a.gbff",
        [("complement(join(99901..100000,1..150))", ['/locus_tag="ABC_0004"'])],
    )
    (cds,) = read_user_cds(p, "gbff", REPS, IDENT)
    assert cds.interval == Interval("chr1", 99900, 100150, "-")


def test_origin_crossing_cds_written_once_to_regions(tmp_path: Path) -> None:
    out = tmp_path / "regions.gff3"
    write_bakta_regions(out, [UserCds("X", Interval("chr1", 99900, 100150, "+"), False)], [CHR])
    feats = feature_lines(out)
    assert len(feats) == 1
    assert feats[0][3] == "99901"
    assert feats[0][8].startswith("ID=X")


# --------------------------------------------------------------- examples 4-5
def test_example4_identical_row_and_tsv(tmp_path: Path) -> None:
    cds = [UserCds("ABC_0001", Interval("chr1", 100, 400, "+"), False)]
    rows = build_id_map(cds, [gene("GSIM_00005", "chr1", 100, 400, "+")], contig_order=ORDER)
    assert rows == [IdMapRow("ABC_0001", "GSIM_00005", "chr1", 100, 400, "+", "identical")]
    out = tmp_path / "id_map.tsv"
    write_id_map(out, rows)
    assert tsv_rows(out) == [
        ID_MAP_HEADER,
        ["ABC_0001", "GSIM_00005", "chr1", "101", "400", "+", "identical"],
    ]


def test_example5_start_differs_unmatched(tmp_path: Path) -> None:
    cds = [UserCds("ABC_0001", Interval("chr1", 97, 400, "+"), False)]
    rows = build_id_map(cds, [gene("GSIM_00005", "chr1", 100, 400, "+")], contig_order=ORDER)
    assert rows == [IdMapRow("ABC_0001", None, "chr1", 97, 400, "+", "unmatched")]
    out = tmp_path / "id_map.tsv"
    write_id_map(out, rows)
    assert tsv_rows(out)[1] == ["ABC_0001", "NA", "chr1", "98", "400", "+", "unmatched"]


@pytest.mark.parametrize(
    "iv",
    [
        Interval("chr1", 100, 403, "+"),  # stop differs
        Interval("chr1", 100, 400, "-"),  # strand differs
        Interval("plasmid_1", 100, 400, "+"),  # contig differs
    ],
)
def test_only_identical_coordinates_match(iv: Interval) -> None:
    rows = build_id_map(
        [UserCds("U1", iv, False)], [gene("G1", "chr1", 100, 400, "+")], contig_order=ORDER
    )
    assert [(r.locus_tag, r.match) for r in rows] == [(None, "unmatched")]
    assert (rows[0].contig, rows[0].start, rows[0].end, rows[0].strand) == (
        iv.contig,
        iv.start,
        iv.end,
        iv.strand,
    )


def test_one_row_per_user_cds_sorted_and_unused_genes_ignored(tmp_path: Path) -> None:
    cds = [
        UserCds("U3", Interval("chr1", 5000, 5300, "-"), False),
        UserCds("U1", Interval("chr1", 100, 400, "+"), False),
        UserCds("U2", Interval("chr1", 1000, 1300, "+"), False),
    ]
    genes = [
        gene("G9", "chr1", 9000, 9300, "+"),
        gene("G3", "chr1", 5000, 5300, "-"),
        gene("G1", "chr1", 100, 400, "+"),
    ]
    rows = build_id_map(cds, genes, contig_order=ORDER)
    assert [(r.source_id, r.locus_tag, r.match) for r in rows] == [
        ("U1", "G1", "identical"),
        ("U2", None, "unmatched"),
        ("U3", "G3", "identical"),
    ]


def test_origin_crossing_identical_and_tsv_keeps_unwrapped_end(tmp_path: Path) -> None:
    iv = Interval("chr1", 99900, 100150, "-")
    rows = build_id_map(
        [UserCds("OX", iv, False)],
        [gene("G_ORI", "chr1", 99900, 100150, "-")],
        contig_order=ORDER,
    )
    assert rows == [IdMapRow("OX", "G_ORI", "chr1", 99900, 100150, "-", "identical")]
    out = tmp_path / "id_map.tsv"
    write_id_map(out, rows)
    assert tsv_rows(out)[1] == ["OX", "G_ORI", "chr1", "99901", "100150", "-", "identical"]


# --------------------------------------------------------------- examples 6-7
def test_example6_contig_map_translates_gff3_seqid(tmp_path: Path) -> None:
    reps = {"chr_1_x": Replicon("chr_1_x", 100_000, "circular", "chromosome")}
    p = gff(tmp_path / "u.gff3", row("chr 1|x", "CDS", 101, 400, "+", "ID=c1;locus_tag=ABC_1"))
    got = read_user_cds(p, "gff3", reps, {"chr 1|x": "chr_1_x"})
    assert got == [UserCds("ABC_1", Interval("chr_1_x", 100, 400, "+"), False)]


def test_contig_map_translates_gbff_locus_name(tmp_path: Path) -> None:
    reps = {"chr_1": Replicon("chr_1", 100_000, "circular", "chromosome")}
    p = gbk(tmp_path / "a.gbff", [("1001..1300", ['/locus_tag="T1"'])], name="NZ_CP0001")
    got = read_user_cds(p, "gbff", reps, {"NZ_CP0001": "chr_1"})
    assert got == [UserCds("T1", Interval("chr_1", 1000, 1300, "+"), False)]


def test_example7_unknown_gff3_contig(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "CDS", 101, 400, "+", "ID=a"),
        row("plasmidZ", "CDS", 101, 400, "+", "ID=b"),
    )
    with pytest.raises(ValueError, match=re.escape("annotation contig plasmidZ not in genome")):
        read_user_cds(p, "gff3", REPS, IDENT)


def test_unknown_gbff_contig(tmp_path: Path) -> None:
    p = gbk(tmp_path / "a.gbff", [("1..300", ['/locus_tag="T1"'])], name="plasmidZ")
    with pytest.raises(ValueError, match=re.escape("annotation contig plasmidZ not in genome")):
        read_user_cds(p, "gbff", REPS, IDENT)


# --------------------------------------------------------------- example 8
def test_example8_detect_gzipped_gbff(tmp_path: Path) -> None:
    p = gbk(tmp_path / "a.gbff", [("1..300", ['/locus_tag="T1"'])])
    gz = tmp_path / "a.gbff.gz"
    gz.write_bytes(gzip.compress(p.read_bytes()))
    assert detect_format(gz) == "gbff"


def test_example8_detect_real_gzipped_bakta_gbff() -> None:
    assert detect_format(DATA / "sim100k_bakta.gbff.gz") == "gbff"


@pytest.mark.parametrize("compress", [False, True])
@pytest.mark.parametrize(
    ("text", "fmt"),
    [
        ('{"sequences": [], "features": []}\n', "bakta_json"),
        ("##gff-version 3\n##sequence-region chr1 1 100\n", "gff3"),
        ("LOCUS       chr1      100 bp    DNA     circular BCT 14-OCT-2026\n//\n", "gbff"),
    ],
)
def test_detect_format_by_content_not_name(
    tmp_path: Path, text: str, fmt: str, compress: bool
) -> None:
    p = tmp_path / "annotation.dat"  # extension says nothing
    raw = text.encode()
    p.write_bytes(gzip.compress(raw) if compress else raw)
    assert detect_format(p) == fmt


def test_detect_format_real_files() -> None:
    assert detect_format(DATA / "sim100k_bakta.gff3") == "gff3"
    assert detect_format(DATA / "k12_excerpt.gbff") == "gbff"


def test_detect_format_unknown_content_raises(tmp_path: Path) -> None:
    p = tmp_path / "x.txt"
    p.write_text(">chr1\nACGT\n")
    with pytest.raises(ValueError):
        detect_format(p)


# --------------------------------------------------------------- source_id rules
@pytest.mark.parametrize(
    ("attrs", "expected"),
    [
        ("ID=cds-1;Name=dnaA;locus_tag=ABC_0001", "ABC_0001"),
        ("ID=cds-1;Name=dnaA", "cds-1"),
        ("Name=dnaA", "dnaA"),
    ],
)
def test_gff3_source_id_precedence(tmp_path: Path, attrs: str, expected: str) -> None:
    p = gff(tmp_path / "u.gff3", row("chr1", "CDS", 101, 400, "+", attrs))
    (cds,) = read_user_cds(p, "gff3", REPS, IDENT)
    assert cds.source_id == expected


def test_gff3_source_id_fallback_contig_start(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "CDS", 101, 400, "+", "."),
        row("chr1", "CDS", 101, 400, "-", "product=x"),
    )
    got = read_user_cds(p, "gff3", REPS, IDENT)
    assert len(got) == 2
    assert all(c.source_id.startswith("chr1_101_") for c in got)
    assert got[0].source_id != got[1].source_id
    assert {c.interval.strand for c in got} == {"+", "-"}


@pytest.mark.parametrize(
    ("quals", "expected"),
    [
        (['/gene="dnaA"', '/locus_tag="ABC_0001"', '/protein_id="WP_1.1"'], "ABC_0001"),
        (['/gene="dnaA"', '/protein_id="WP_1.1"'], "WP_1.1"),
        (['/gene="dnaA"'], "dnaA"),
    ],
)
def test_gbff_source_id_precedence(tmp_path: Path, quals: list[str], expected: str) -> None:
    p = gbk(tmp_path / "a.gbff", [("101..400", quals)])
    (cds,) = read_user_cds(p, "gbff", REPS, IDENT)
    assert cds.source_id == expected


# --------------------------------------------------------------- GFF3 merging
@pytest.mark.parametrize("strand", ["+", "-"])
def test_gff3_origin_split_rows_merge(tmp_path: Path, strand: str) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "CDS", 99901, 100000, strand, "ID=cds-9;locus_tag=ABC_9"),
        row("chr1", "CDS", 1, 150, strand, "ID=cds-9;locus_tag=ABC_9"),
    )
    got = read_user_cds(p, "gff3", REPS, IDENT)
    assert got == [UserCds("ABC_9", Interval("chr1", 99900, 100150, strand), False)]  # type: ignore[arg-type]


def test_gff3_unwrapped_origin_row(tmp_path: Path) -> None:
    p = gff(tmp_path / "u.gff3", row("chr1", "CDS", 99901, 100150, "+", "ID=cds-9"))
    got = read_user_cds(p, "gff3", REPS, IDENT)
    assert got == [UserCds("cds-9", Interval("chr1", 99900, 100150, "+"), False)]


def test_gff3_multi_row_cds_sharing_id_merges(tmp_path: Path) -> None:
    """A programmed frameshift (prfB-like) written as two CDS rows with one ID."""
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "CDS", 2001, 2075, "+", "ID=cds-prfB;locus_tag=ABC_PRFB"),
        row("chr1", "CDS", 2075, 3100, "+", "ID=cds-prfB;locus_tag=ABC_PRFB"),
        row("chr1", "CDS", 5001, 5300, "+", "ID=cds-other"),
    )
    got = read_user_cds(p, "gff3", REPS, IDENT)
    assert UserCds("ABC_PRFB", Interval("chr1", 2000, 3100, "+"), False) in got
    assert len(got) == 2


def test_gff3_ignores_non_cds_types(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "region", 1, 100000, ".", "ID=chr1;Is_circular=true"),
        row("chr1", "gene", 101, 400, "+", "ID=g1;locus_tag=A1"),
        row("chr1", "tRNA", 501, 576, "-", "ID=t1;locus_tag=A2"),
        row("chr1", "rRNA", 1001, 2500, "+", "ID=r1;locus_tag=A3"),
        row("chr1", "pseudogene", 3001, 3300, "+", "ID=p1;locus_tag=A4"),
        row("chr1", "CDS", 101, 400, "+", "ID=c1;Parent=g1;locus_tag=A1"),
    )
    assert [c.source_id for c in read_user_cds(p, "gff3", REPS, IDENT)] == ["A1"]


def test_gff3_with_fasta_section(tmp_path: Path) -> None:
    p = tmp_path / "u.gff3"
    p.write_text(
        "##gff-version 3\n"
        + row("chr1", "CDS", 101, 400, "+", "ID=c1")
        + "\n##FASTA\n>chr1\nACGT\n"
    )
    assert [c.source_id for c in read_user_cds(p, "gff3", REPS, IDENT)] == ["c1"]


def test_gff3_gzipped_input(tmp_path: Path) -> None:
    plain = ex1_gff(tmp_path)
    gz = tmp_path / "user.gff3.gz"
    gz.write_bytes(gzip.compress(plain.read_bytes()))
    assert read_user_cds(gz, "gff3", REPS, IDENT) == read_user_cds(plain, "gff3", REPS, IDENT)


# --------------------------------------------------------------- pseudo
def test_gff3_pseudo_true_flagged_and_written(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "CDS", 101, 400, "+", "ID=c1;locus_tag=P1;pseudo=true"),
        row("chr1", "CDS", 1001, 1300, "+", "ID=c2;locus_tag=P2"),
    )
    got = read_user_cds(p, "gff3", REPS, IDENT)
    assert [(c.source_id, c.pseudo) for c in got] == [("P1", True), ("P2", False)]
    out = tmp_path / "regions.gff3"
    write_bakta_regions(out, got, [CHR])
    ids = [re.match(r"ID=([^;]+)", f[8]).group(1) for f in feature_lines(out)]  # type: ignore[union-attr]
    assert ids == ["P1", "P2"]


def test_gbff_pseudo_flagged(tmp_path: Path) -> None:
    p = gbk(
        tmp_path / "a.gbff",
        [
            ("101..400", ['/locus_tag="P1"', "/pseudo"]),
            ("1001..1300", ['/locus_tag="P2"']),
        ],
    )
    got = read_user_cds(p, "gbff", REPS, IDENT)
    assert [(c.source_id, c.pseudo) for c in got] == [("P1", True), ("P2", False)]


def test_real_k12_bakta_gbff_pseudogene_and_partial_location() -> None:
    got = read_user_cds(DATA / "k12_excerpt.gbff", "gbff", {K12.id: K12}, {K12.id: K12.id})
    assert got == [
        UserCds("KOIAIB_00744", Interval("U00096.3", 715411, 715906, "+"), False),
        UserCds("KOIAIB_00745", Interval("U00096.3", 716387, 716597, "-"), True),
        UserCds("KOIAIB_00746", Interval("U00096.3", 716945, 718265, "-"), False),
    ]


def test_real_k12_bakta_gff3_pseudogene() -> None:
    got = read_user_cds(DATA / "k12_excerpt.gff3", "gff3", {K12.id: K12}, {K12.id: K12.id})
    assert got == [
        UserCds("KOIAIB_00744", Interval("U00096.3", 715411, 715906, "+"), False),
        UserCds("KOIAIB_00745", Interval("U00096.3", 716387, 716597, "-"), True),
        UserCds("KOIAIB_00746", Interval("U00096.3", 716945, 718265, "-"), False),
    ]


# --------------------------------------------------------------- edge cases
def test_empty_gff3(tmp_path: Path) -> None:
    p = gff(
        tmp_path / "u.gff3",
        row("chr1", "gene", 101, 400, "+", "ID=g1"),
        pragmas=("##sequence-region chr1 1 100000",),
    )
    cds = read_user_cds(p, "gff3", REPS, IDENT)
    assert cds == []
    regions = tmp_path / "regions.gff3"
    write_bakta_regions(regions, cds, [CHR, PL])
    lines = regions.read_text().splitlines()
    assert lines[0] == "##gff-version 3"
    assert all(ln.startswith("#") for ln in lines)
    idm = tmp_path / "id_map.tsv"
    write_id_map(idm, build_id_map(cds, [gene("G1", "chr1", 100, 400, "+")], contig_order=ORDER))
    assert idm.read_text() == "\t".join(ID_MAP_HEADER) + "\n"


def test_duplicate_coordinates_both_map_to_same_locus() -> None:
    iv = Interval("chr1", 100, 400, "+")
    cds = [UserCds("DUP_B", iv, False), UserCds("DUP_A", iv, False)]
    rows = build_id_map(cds, [gene("G1", "chr1", 100, 400, "+")], contig_order=ORDER)
    assert len(rows) == 2
    assert {r.source_id for r in rows} == {"DUP_A", "DUP_B"}
    assert {(r.locus_tag, r.match) for r in rows} == {("G1", "identical")}


def test_gff3_linear_contig_beyond_length_raises(tmp_path: Path) -> None:
    p = gff(tmp_path / "u.gff3", row("linear_1", "CDS", 4901, 5150, "+", "ID=c1"))
    with pytest.raises(ValueError):
        read_user_cds(p, "gff3", REPS, IDENT)


def test_gbff_join_on_linear_record_raises(tmp_path: Path) -> None:
    p = gbk(
        tmp_path / "a.gbff",
        [("join(4901..5000,1..150)", ['/locus_tag="L1"'])],
        name="linear_1",
        length=5_000,
        topology="linear",
    )
    with pytest.raises(ValueError):
        read_user_cds(p, "gbff", REPS, IDENT)


def test_gbff_multi_record(tmp_path: Path) -> None:
    a = gbk(tmp_path / "a.gbff", [("101..400", ['/locus_tag="C1"'])])
    b = gbk(
        tmp_path / "b.gbff",
        [("complement(11..310)", ['/locus_tag="P1"'])],
        name="plasmid_1",
        length=8_000,
    )
    both = tmp_path / "both.gbff"
    both.write_text(a.read_text() + b.read_text())
    got = read_user_cds(both, "gbff", REPS, IDENT)
    assert sorted(got, key=lambda c: c.source_id) == [
        UserCds("C1", Interval("chr1", 100, 400, "+"), False),
        UserCds("P1", Interval("plasmid_1", 10, 310, "-"), False),
    ]


def test_gbff_non_cds_features_ignored(tmp_path: Path) -> None:
    p = tmp_path / "a.gbff"
    text = gbk(p, [("101..400", ['/locus_tag="C1"'])]).read_text()
    extra = (
        "     gene            101..400\n"
        '                     /locus_tag="C1"\n'
        "     tRNA            complement(501..576)\n"
        '                     /locus_tag="T1"\n'
    )
    p.write_text(text.replace("     CDS ", extra + "     CDS ", 1))
    assert [c.source_id for c in read_user_cds(p, "gbff", REPS, IDENT)] == ["C1"]


def test_dataclasses_frozen() -> None:
    c = UserCds("A", Interval("chr1", 0, 3, "+"), False)
    r = IdMapRow("A", None, "chr1", 0, 3, "+", "unmatched")
    with pytest.raises(AttributeError):
        c.source_id = "B"  # type: ignore[misc]
    with pytest.raises(AttributeError):
        r.match = "identical"  # type: ignore[misc]


# --------------------------------------------------------------- real Bakta outputs
def test_real_bakta_gff3_maps_identically_to_its_locus_tags() -> None:
    cds = read_user_cds(DATA / "sim100k_bakta.gff3", "gff3", SIM_REPS, {k: k for k in SIM_REPS})
    genes = bakta_tsv_genes(DATA / "sim100k_bakta.tsv")
    assert len(cds) == len(genes) == 86
    assert not any(c.pseudo for c in cds)
    rows = build_id_map(cds, genes, contig_order=list(SIM_REPS))
    assert len(rows) == 86
    assert all(r.match == "identical" and r.locus_tag == r.source_id for r in rows)
    assert [r.locus_tag for r in rows] == [g.locus_tag for g in genes]  # Bakta TSV is §6.5 order


def test_real_bakta_gbff_equals_gff3() -> None:
    ident = {k: k for k in SIM_REPS}
    from_gbff = read_user_cds(DATA / "sim100k_bakta.gbff.gz", "gbff", SIM_REPS, ident)
    from_gff3 = read_user_cds(DATA / "sim100k_bakta.gff3", "gff3", SIM_REPS, ident)
    key = lambda c: (c.interval.contig, c.interval.start, c.source_id)  # noqa: E731
    assert sorted(from_gbff, key=key) == sorted(from_gff3, key=key)


def test_real_bakta_gff3_through_regions_and_back(tmp_path: Path) -> None:
    ident = {k: k for k in SIM_REPS}
    cds = read_user_cds(DATA / "sim100k_bakta.gff3", "gff3", SIM_REPS, ident)
    regions = tmp_path / "regions.gff3"
    write_bakta_regions(regions, cds, list(SIM_REPS.values()))
    again = read_user_cds(regions, "gff3", SIM_REPS, ident)
    key = lambda c: (c.interval.contig, c.interval.start, c.source_id)  # noqa: E731
    assert sorted(again, key=key) == sorted(cds, key=key)
    assert detect_format(regions) == "gff3"


# --------------------------------------------------------------- determinism
def test_writers_are_byte_identical_and_order_independent(tmp_path: Path) -> None:
    ident = {k: k for k in SIM_REPS}
    cds = read_user_cds(DATA / "sim100k_bakta.gff3", "gff3", SIM_REPS, ident)
    genes = bakta_tsv_genes(DATA / "sim100k_bakta.tsv")
    shuffled = list(cds)
    random.Random(7).shuffle(shuffled)
    shuffled_genes = list(genes)
    random.Random(8).shuffle(shuffled_genes)
    outs: list[tuple[bytes, bytes]] = []
    for i, (c, g) in enumerate([(cds, genes), (cds, genes), (shuffled, shuffled_genes)]):
        r, m = tmp_path / f"r{i}.gff3", tmp_path / f"m{i}.tsv"
        write_bakta_regions(r, c, list(SIM_REPS.values()))
        write_id_map(m, build_id_map(c, g, contig_order=list(SIM_REPS)))
        outs.append((r.read_bytes(), m.read_bytes()))
    assert outs[0] == outs[1] == outs[2]


def test_read_user_cds_deterministic() -> None:
    ident = {k: k for k in SIM_REPS}
    a = read_user_cds(DATA / "sim100k_bakta.gbff.gz", "gbff", SIM_REPS, ident)
    b = read_user_cds(DATA / "sim100k_bakta.gbff.gz", "gbff", SIM_REPS, ident)
    assert a == b


# --------------------------------------------------------------- simulator round trip
def _sim_round_trip(cfg: object, seed: int, tmp: Path) -> None:
    from giap.annotate.bakta import parse_bakta_json
    from giapsim import simulate

    sim = simulate(cfg, seed, tmp)  # type: ignore[arg-type]
    ann = parse_bakta_json(sim.bakta_json)
    reps: Mapping[str, Replicon] = ann.model.replicons
    assert sim.annotation_gff3 is not None
    cds = read_user_cds(sim.annotation_gff3, "gff3", reps, {k: k for k in reps})
    assert cds
    rows = build_id_map(cds, ann.model.genes, contig_order=list(reps))
    assert len(rows) == len(cds)
    assert all(r.match == "identical" for r in rows)
    assert all(r.locus_tag == r.source_id for r in rows)  # giapsim: ID = locus_tag


def test_simulator_round_trip_all_identical(tmp_path: Path) -> None:
    from giapsim import SimConfig

    cfg = SimConfig(
        genome_length=40_000,
        n_plasmids=1,
        plasmid_length=(6_000, 6_000),
        linear_contig_length=5_000,
        origin_crossing_gene=True,
        n_pseudogenes=2,
        n_sorfs=2,
    )
    _sim_round_trip(cfg, 11, tmp_path)


def test_simulator_round_trip_sim100k(tmp_path: Path) -> None:
    from giapsim import SIM100K, SIM100K_SEED

    _sim_round_trip(SIM100K, SIM100K_SEED, tmp_path)
