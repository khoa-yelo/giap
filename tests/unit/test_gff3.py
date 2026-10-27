"""GFF3 writer and minimal reader (issue M1-05; GIPPA_SPEC.md §6.5, §9.5, decision 16).

GFF3 uses 1-based closed coordinates (``start + 1``, ``end``). Origin-crossing features on circular
replicons keep ``end > length``; circular replicons get a region line with ``Is_circular=true``.
"""

from __future__ import annotations

import gzip

import pytest

from gippa.io.gff3 import GffFeature, escape_attr, read_gff3, write_gff3
from gippa.model import Interval, Replicon

C = Replicon("chr", 2000, "circular")
L = Replicon("chr", 2000, "linear")
P = Replicon("pl1", 500, "circular", "plasmid")
LIN = Replicon("lin", 300, "linear")

REGION_C = "chr\tgippa\tregion\t1\t2000\t.\t.\t.\tID=chr;Is_circular=true"


def _lines(path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    return text[:-1].split("\n")


def _feature_lines(path) -> list[str]:
    return [ln for ln in _lines(path) if ln and not ln.startswith("#")]


def _g1() -> GffFeature:
    return GffFeature(
        Interval("chr", 99, 300, "+"),
        "gippa",
        "CDS",
        phase=0,
        attributes={"ID": "g1", "Name": "dnaA"},
    )


# --- Example 1: 1-based closed columns ----------------------------------------------------------


def test_example_01_cds_line(tmp_path) -> None:
    p = tmp_path / "a.gff3"
    write_gff3(p, [_g1()], [C])
    assert "chr\tgippa\tCDS\t100\t300\t.\t+\t0\tID=g1;Name=dnaA" in _feature_lines(p)


def test_example_01_linear_full_file(tmp_path) -> None:
    p = tmp_path / "a.gff3"
    write_gff3(p, [_g1()], [L])
    assert p.read_text() == (
        "##gff-version 3\n"
        "##sequence-region chr 1 2000\n"
        "chr\tgippa\tCDS\t100\t300\t.\t+\t0\tID=g1;Name=dnaA\n"
    )


# --- Example 2: origin-crossing feature on a circular replicon ----------------------------------


def test_example_02_origin_crossing_full_file(tmp_path) -> None:
    f = GffFeature(
        Interval("chr", 1950, 2050, "-"), "gippa", "CDS", phase=0, attributes={"ID": "g2"}
    )
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [C])
    assert p.read_text() == (
        "##gff-version 3\n"
        "##sequence-region chr 1 2000\n"
        f"{REGION_C}\n"
        "chr\tgippa\tCDS\t1951\t2050\t.\t-\t0\tID=g2\n"
    )


def test_first_line_is_gff_version(tmp_path) -> None:
    p = tmp_path / "a.gff3"
    write_gff3(p, [_g1()], [C])
    assert _lines(p)[0] == "##gff-version 3"


def test_region_lines_only_for_circular_and_before_features(tmp_path) -> None:
    feats = [
        GffFeature(Interval("lin", 10, 20, "+"), "gippa", "gene", attributes={"ID": "a"}),
        GffFeature(Interval("pl1", 490, 510, "+"), "gippa", "gene", attributes={"ID": "b"}),
        _g1(),
    ]
    p = tmp_path / "a.gff3"
    write_gff3(p, feats, [C, LIN, P])
    lines = _lines(p)
    sr = [ln for ln in lines if ln.startswith("##sequence-region")]
    assert sr == [
        "##sequence-region chr 1 2000",
        "##sequence-region lin 1 300",
        "##sequence-region pl1 1 500",
    ]
    regions = [ln for ln in lines if "\tregion\t" in ln]
    assert regions == [REGION_C, "pl1\tgippa\tregion\t1\t500\t.\t.\t.\tID=pl1;Is_circular=true"]
    # each region line follows its own sequence-region line; all precede the first feature
    first_feature = min(i for i, ln in enumerate(lines) if "\tgene\t" in ln or "\tCDS\t" in ln)
    for rid, region in (("chr", regions[0]), ("pl1", regions[1])):
        assert lines.index(f"##sequence-region {rid} 1 {'2000' if rid == 'chr' else '500'}") < (
            lines.index(region)
        )
        assert lines.index(region) < first_feature
    assert "##sequence-region pl1 1 500" not in lines[first_feature:]


# --- Example 3: attribute escaping ---------------------------------------------------------------


def test_example_03_escape_attr() -> None:
    assert escape_attr("tRNA-Leu; anticodon=CAG") == "tRNA-Leu%3B anticodon%3DCAG"


@pytest.mark.parametrize(
    ("raw", "escaped"),
    [
        (";", "%3B"),
        ("=", "%3D"),
        ("&", "%26"),
        (",", "%2C"),
        ("\t", "%09"),
        ("\n", "%0A"),
        ("\r", "%0D"),
        ("%", "%25"),
        ("\x01", "%01"),
        ("\x7f", "%7F"),
        ("a b-c_d.e:f/g", "a b-c_d.e:f/g"),
    ],
)
def test_escape_attr_characters(raw: str, escaped: str) -> None:
    assert escape_attr(raw) == escaped


def test_escape_attr_percent_is_not_double_decoded() -> None:
    assert escape_attr("50%3B") == "50%253B"


def test_example_03_written_and_read_back(tmp_path) -> None:
    product = "tRNA-Leu; anticodon=CAG"
    f = GffFeature(
        Interval("chr", 9, 90, "+"), "gippa", "tRNA", attributes={"ID": "t1", "product": product}
    )
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p) == [
        "chr\tgippa\ttRNA\t10\t90\t.\t+\t.\tID=t1;product=tRNA-Leu%3B anticodon%3DCAG"
    ]
    (back,) = read_gff3(p, {"chr": L})
    assert back.attributes["product"] == product


def test_list_values_joined_by_comma_after_escaping(tmp_path) -> None:
    f = GffFeature(
        Interval("chr", 9, 90, "+"),
        "gippa",
        "CDS",
        phase=0,
        attributes={"ID": "c1", "Dbxref": ["EC:1.2", "a,b"]},
    )
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p)[0].endswith("\tID=c1;Dbxref=EC:1.2,a%2Cb")
    (back,) = read_gff3(p, {"chr": L})
    assert list(back.attributes["Dbxref"]) == ["EC:1.2", "a,b"]


def test_attribute_insertion_order_kept(tmp_path) -> None:
    attrs = {"ID": "c1", "zeta": "1", "alpha": "2", "Name": "x"}
    f = GffFeature(Interval("chr", 9, 90, "+"), "gippa", "CDS", phase=0, attributes=attrs)
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p)[0].endswith("\tID=c1;zeta=1;alpha=2;Name=x")
    (back,) = read_gff3(p, {"chr": L})
    assert list(back.attributes) == ["ID", "zeta", "alpha", "Name"]


# --- Columns: score, phase, strand ---------------------------------------------------------------


@pytest.mark.parametrize(
    ("score", "col"), [(None, "."), (12.5, "12.5"), (0.123456789, "0.123457"), (1e-10, "1e-10")]
)
def test_score_column(tmp_path, score, col) -> None:
    f = GffFeature(
        Interval("chr", 0, 9, "+"), "gippa", "CDS", score=score, phase=0, attributes={"ID": "x"}
    )
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p)[0].split("\t")[5] == col


def test_cds_phase_none_raises_on_write(tmp_path) -> None:
    # ADR-0005 item 1 (after the M1-05 review): a CDS must carry a phase (GFF3 3.1.26).
    f = GffFeature(Interval("chr", 0, 9, "+"), "gippa", "CDS", attributes={"ID": "x"})
    with pytest.raises(ValueError):
        write_gff3(tmp_path / "a.gff3", [f], [L])


@pytest.mark.parametrize("phase", [0, 1, 2])
def test_phase_written(tmp_path, phase) -> None:
    f = GffFeature(Interval("chr", 0, 9, "-"), "gippa", "CDS", phase=phase, attributes={"ID": "x"})
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p)[0].split("\t")[6:8] == ["-", str(phase)]


def test_unstranded_feature(tmp_path) -> None:
    f = GffFeature(Interval("chr", 0, 9, "."), "gippa", "repeat_region", attributes={"ID": "r"})
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [L])
    assert _feature_lines(p) == ["chr\tgippa\trepeat_region\t1\t9\t.\t.\t.\tID=r"]


# --- Sorting (§6.5: contig order, start, end, strand, ID) ----------------------------------------


def test_features_sorted_by_contig_order_then_start(tmp_path) -> None:
    def g(contig: str, s: int, e: int, strand: str, id_: str) -> GffFeature:
        return GffFeature(Interval(contig, s, e, strand), "gippa", "gene", attributes={"ID": id_})

    feats = [
        g("chr", 500, 600, "+", "d"),
        g("pl1", 10, 20, "+", "e"),
        g("chr", 100, 300, "-", "c"),
        g("chr", 100, 300, "+", "b"),
        g("chr", 100, 200, "+", "a"),
    ]
    p = tmp_path / "a.gff3"
    write_gff3(p, feats, [P, C])  # pl1 is first in the given order
    ids = [ln.rsplit("ID=", 1)[1] for ln in _feature_lines(p) if "\tgene\t" in ln]
    assert ids == ["e", "a", "b", "c", "d"]


# --- Validation (strict, ADR-0006 style) ----------------------------------------------------------


def test_writer_unknown_contig_raises(tmp_path) -> None:
    f = GffFeature(Interval("nope", 0, 9, "+"), "gippa", "CDS", phase=0, attributes={"ID": "x"})
    with pytest.raises(ValueError):
        write_gff3(tmp_path / "a.gff3", [f], [C])


def test_writer_linear_overrun_raises(tmp_path) -> None:
    f = GffFeature(Interval("chr", 1950, 2050, "+"), "gippa", "CDS", phase=0, attributes={"ID": "x"})
    with pytest.raises(ValueError):
        write_gff3(tmp_path / "a.gff3", [f], [L])


# --- Edge cases: no features; FASTA section ------------------------------------------------------


def test_no_features_circular(tmp_path) -> None:
    p = tmp_path / "a.gff3"
    write_gff3(p, [], [C])
    assert p.read_text() == f"##gff-version 3\n##sequence-region chr 1 2000\n{REGION_C}\n"
    assert read_gff3(p, {"chr": C}) == []


def test_no_features_linear(tmp_path) -> None:
    p = tmp_path / "a.gff3"
    write_gff3(p, [], [LIN])
    assert p.read_text() == "##gff-version 3\n##sequence-region lin 1 300\n"


def test_fasta_section_at_end(tmp_path) -> None:
    seq = "ACGT" * 75
    p = tmp_path / "a.gff3"
    feat = GffFeature(Interval("lin", 0, 9, "+"), "gippa", "CDS", phase=0, attributes={"ID": "x"})
    write_gff3(p, [feat], [LIN], fasta={"lin": seq})
    lines = _lines(p)
    i = lines.index("##FASTA")
    assert all(not ln.startswith("lin\t") for ln in lines[i:])
    assert lines[i + 1] == ">lin"
    assert "".join(lines[i + 2 :]) == seq
    assert read_gff3(p, {"lin": LIN}) == [feat]


# --- Example 4: round trip -----------------------------------------------------------------------


def test_example_04_round_trip(tmp_path) -> None:
    feats = [
        _g1(),
        GffFeature(
            Interval("chr", 400, 700, "-"),
            "gippa",
            "gene",
            score=3.5,
            attributes={"ID": "g3", "Note": "a;b=c,d&e%f"},
        ),
        GffFeature(
            Interval("chr", 1950, 2050, "-"), "gippa", "CDS", phase=2, attributes={"ID": "g2"}
        ),
    ]
    p = tmp_path / "a.gff3"
    write_gff3(p, feats, [C])
    assert read_gff3(p, {"chr": C}) == feats


def test_two_writes_identical_bytes(tmp_path) -> None:
    feats = [
        GffFeature(
            Interval("chr", 1950, 2050, "-"), "gippa", "CDS", phase=0, attributes={"ID": "g2"}
        ),
        _g1(),
    ]
    a, b = tmp_path / "a.gff3", tmp_path / "b.gff3"
    write_gff3(a, feats, [C])
    write_gff3(b, list(reversed(feats)), [C])
    assert a.read_bytes() == b.read_bytes()


def test_gzip_output_reproducible(tmp_path) -> None:
    a, b, plain = tmp_path / "a.gff3.gz", tmp_path / "b.gff3.gz", tmp_path / "p.gff3"
    for path in (a, b, plain):
        write_gff3(path, [_g1()], [C])
    assert a.read_bytes() == b.read_bytes()
    assert gzip.decompress(a.read_bytes()) == plain.read_bytes()
    assert read_gff3(a, {"chr": C}) == [_g1()]


# --- Example 5 and reader validation -------------------------------------------------------------


def _write(tmp_path, body: str):
    p = tmp_path / "in.gff3"
    p.write_text("##gff-version 3\n" + body, encoding="utf-8")
    return p


def test_example_05_end_beyond_length_on_linear_raises(tmp_path) -> None:
    p = _write(tmp_path, "chr\t.\tCDS\t1951\t2050\t.\t-\t0\tID=x\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": L})


def test_end_beyond_length_on_circular_accepted(tmp_path) -> None:
    p = _write(tmp_path, "chr\t.\tCDS\t1951\t2050\t.\t-\t0\tID=x\n")
    (f,) = read_gff3(p, {"chr": C})
    assert f.interval == Interval("chr", 1950, 2050, "-")
    assert (f.source, f.type, f.score, f.phase) == (".", "CDS", None, 0)


def test_start_greater_than_end_raises_even_on_circular(tmp_path) -> None:
    p = _write(tmp_path, "chr\t.\tCDS\t1951\t50\t.\t+\t0\tID=x\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


def test_unknown_seqid_raises(tmp_path) -> None:
    p = _write(tmp_path, "chrX\t.\tCDS\t1\t9\t.\t+\t0\tID=x\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


@pytest.mark.parametrize(
    "line",
    [
        "chr\t.\tCDS\t1\t9\t.\t+\t0",  # 8 columns
        "chr\t.\tCDS\t0\t9\t.\t+\t0\tID=x",  # start < 1
        "chr\t.\tCDS\ta\t9\t.\t+\t0\tID=x",  # non-integer start
        "chr\t.\tCDS\t1\t9\t.\t+\t3\tID=x",  # phase out of range
        "chr\t.\tCDS\t1\t9\t.\t*\t0\tID=x",  # bad strand
        "chr\t.\tCDS\t1\t9\tabc\t+\t0\tID=x",  # bad score
    ],
)
def test_malformed_lines_raise(tmp_path, line: str) -> None:
    p = _write(tmp_path, line + "\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


def test_reader_third_party_file(tmp_path) -> None:
    p = _write(
        tmp_path,
        "##sequence-region chr 1 2000\n"
        "# a free comment\n"
        "\n"
        "##species https://example.org\n"
        "chr\tProdigal\tCDS\t100\t300\t12.5\t+\t0\tID=cds1;Note=a%3Bb%2Cc;Dbxref=x:1\n"
        "###\n"
        "chr\tProdigal\tCDS\t1951\t2050\t.\t-\t.\tID=cds2\n"
        "##FASTA\n"
        ">chr\n"
        "ACGT\n",
    )
    got = read_gff3(p, {"chr": C})
    assert got == [
        GffFeature(
            Interval("chr", 99, 300, "+"),
            "Prodigal",
            "CDS",
            score=12.5,
            phase=0,
            attributes={"ID": "cds1", "Note": "a;b,c", "Dbxref": "x:1"},
        ),
        GffFeature(Interval("chr", 1950, 2050, "-"), "Prodigal", "CDS", attributes={"ID": "cds2"}),
    ]


def test_reader_stops_at_fasta(tmp_path) -> None:
    p = _write(
        tmp_path,
        "chr\t.\tgene\t1\t9\t.\t+\t.\tID=a\n##FASTA\n>chr\nchr\t.\tgene\t1\t9\t.\t+\t.\tID=b\n",
    )
    assert [f.attributes["ID"] for f in read_gff3(p, {"chr": C})] == ["a"]
