"""VCF 4.3 writer (issue M1-05; GIAP_SPEC.md §6.5, §7.5.5, §9.5, B.7, G14; ADR-0005 item 1).

POS is 1-based within ``[1, length]``; REF/ALT are forward-strand codons supplied by the caller
(minus-strand conversion is the caller's job, G14). A codon crossing the origin raises.
"""

from __future__ import annotations

import gzip
import importlib.util

import pytest

from giap.io.vcf import VcfHeaderField, VcfRecord, write_vcf
from giap.model import Replicon

C = Replicon("chr", 2000, "circular")
L = Replicon("lin", 900, "linear")
SHA = "0123456789abcdef" * 4

INFO_FIELDS = [
    VcfHeaderField("GENE", "1", "String", "Locus tag"),
    VcfHeaderField("AA_POS", "1", "Integer", "Residue position in the protein"),
    VcfHeaderField("AA_REF", "1", "String", "Reference residue"),
    VcfHeaderField("AA_ALT", "A", "String", "Alternative residue(s)"),
    VcfHeaderField("CALL", "1", "String", "Call"),
    VcfHeaderField("CONF", "1", "String", "Confidence"),
    VcfHeaderField("NT_CHANGES", "A", "Integer", "Nucleotide changes"),
    VcfHeaderField("ALT_TYPE", "1", "String", "ALT type"),
    VcfHeaderField("FLAGS", "1", "String", "Flags joined by |"),
    VcfHeaderField("CLASS_ERR", "1", "String", "Class error rate"),
]
FORMAT_FIELDS = [
    VcfHeaderField("GT", "1", "String", "Genotype"),
    VcfHeaderField("AD", "R", "Integer", "Distinct peptides for REF and each ALT"),
    VcfHeaderField("DP", "1", "Integer", "Total distinct peptides"),
]

B7_INFO = {
    "GENE": "NISS_00012",
    "AA_POS": 45,
    "AA_REF": "A",
    "AA_ALT": "S",
    "CALL": "ALT",
    "CONF": "high",
    "NT_CHANGES": 1,
    "ALT_TYPE": "substitution",
    "FLAGS": ("mass_equivalent", "snv"),
    "CLASS_ERR": "unestimated",
}
B7_LINE = (
    "chr\t1001\t.\tGCT\tTCT\t.\tPASS\t"
    "GENE=NISS_00012;AA_POS=45;AA_REF=A;AA_ALT=S;CALL=ALT;CONF=high;NT_CHANGES=1;"
    "ALT_TYPE=substitution;FLAGS=mass_equivalent|snv;CLASS_ERR=unestimated\t"
    "GT:AD:DP\t1:0,3:3"
)
CHROM_LINE = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpilot_A3"


def _b7(**kw) -> VcfRecord:
    base = dict(
        contig="chr", pos0=1000, ref="GCT", alts=("TCT",), info=B7_INFO, gt="1", ad=(0, 3), dp=3
    )
    base.update(kw)
    return VcfRecord(**base)


def _write(path, records, replicons=(C,), **kw):
    args = dict(
        giap_version="0.1.0",
        genome_sha256=SHA,
        sample="pilot_A3",
        info_fields=INFO_FIELDS,
        format_fields=FORMAT_FIELDS,
    )
    args.update(kw)
    write_vcf(path, records, list(replicons), **args)
    return path


def _lines(path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    return text[:-1].split("\n")


def _records(path) -> list[str]:
    return [ln for ln in _lines(path) if not ln.startswith("#")]


# --- Header ---------------------------------------------------------------------------------------


def test_header_exact(tmp_path) -> None:
    p = _write(
        tmp_path / "v.vcf",
        [],
        replicons=(C, L),
        info_fields=INFO_FIELDS[:2],
        format_fields=FORMAT_FIELDS,
    )
    assert _lines(p) == [
        "##fileformat=VCFv4.3",
        "##source=giap-0.1.0",
        f"##reference=sha256:{SHA}",
        "##contig=<ID=chr,length=2000>",
        "##contig=<ID=lin,length=900>",
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Locus tag">',
        '##INFO=<ID=AA_POS,Number=1,Type=Integer,Description="Residue position in the protein">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Distinct peptides for REF and each '
        'ALT">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total distinct peptides">',
        CHROM_LINE,
    ]


@pytest.mark.parametrize("bad", ["", "abc", SHA.upper(), SHA[:-1], "g" * 64, "sha256:" + SHA])
def test_bad_genome_sha256_raises(tmp_path, bad: str) -> None:
    with pytest.raises(ValueError):
        _write(tmp_path / "v.vcf", [], genome_sha256=bad)


# --- Example 6: Appendix B.7 record byte for byte -----------------------------------------------


def test_example_06_appendix_b7_record(tmp_path) -> None:
    p = _write(tmp_path / "v.vcf", [_b7()])
    lines = _lines(p)
    assert lines[-2] == CHROM_LINE
    assert lines[-1] == B7_LINE


def test_example_06_full_file_with_b7_record(tmp_path) -> None:
    p = _write(tmp_path / "v.vcf", [_b7()], info_fields=[], format_fields=[])
    assert p.read_text() == (
        "##fileformat=VCFv4.3\n"
        "##source=giap-0.1.0\n"
        f"##reference=sha256:{SHA}\n"
        "##contig=<ID=chr,length=2000>\n"
        f"{CHROM_LINE}\n"
        f"{B7_LINE}\n"
    )


# --- Example 7 and origin handling (ADR-0005 item 1) --------------------------------------------


def test_example_07_unwrapped_pos_wraps(tmp_path) -> None:
    p = _write(tmp_path / "v.vcf", [_b7(pos0=2003)])
    assert _records(p)[0].split("\t")[:2] == ["chr", "4"]


def test_codon_ending_exactly_at_contig_end_ok(tmp_path) -> None:
    p = _write(tmp_path / "v.vcf", [_b7(pos0=1997)])
    assert _records(p)[0].split("\t")[1] == "1998"


@pytest.mark.parametrize("pos0", [1998, 1999, 3999])
def test_codon_crossing_origin_raises(tmp_path, pos0: int) -> None:
    with pytest.raises(ValueError):
        _write(tmp_path / "v.vcf", [_b7(pos0=pos0)])


@pytest.mark.parametrize("pos0", [898, 900, 5000])
def test_linear_overrun_raises(tmp_path, pos0: int) -> None:
    with pytest.raises(ValueError):
        _write(tmp_path / "v.vcf", [_b7(contig="lin", pos0=pos0)], replicons=(C, L))


def test_unknown_contig_raises(tmp_path) -> None:
    with pytest.raises(ValueError):
        _write(tmp_path / "v.vcf", [_b7(contig="nope")])


def test_ad_length_must_match_alleles(tmp_path) -> None:
    with pytest.raises(ValueError):
        _write(tmp_path / "v.vcf", [_b7(alts=("TCT", "ACT"), ad=(0, 3))])


# --- Sorting --------------------------------------------------------------------------------------


def test_records_sorted_by_contig_order_then_pos(tmp_path) -> None:
    recs = [
        _b7(pos0=1000),
        _b7(contig="lin", pos0=10),
        _b7(pos0=2003),  # wraps to POS 4
        _b7(pos0=30),
    ]
    p = _write(tmp_path / "v.vcf", recs, replicons=(L, C))
    assert [tuple(r.split("\t")[:2]) for r in _records(p)] == [
        ("lin", "11"),
        ("chr", "4"),
        ("chr", "31"),
        ("chr", "1001"),
    ]


def test_two_writes_identical_bytes(tmp_path) -> None:
    recs = [_b7(pos0=1000), _b7(pos0=30), _b7(contig="lin", pos0=10)]
    a = _write(tmp_path / "a.vcf", recs, replicons=(C, L))
    b = _write(tmp_path / "b.vcf", list(reversed(recs)), replicons=(C, L))
    assert a.read_bytes() == b.read_bytes()


def test_gzip_output_reproducible(tmp_path) -> None:
    a = _write(tmp_path / "a.vcf.gz", [_b7()])
    b = _write(tmp_path / "b.vcf.gz", [_b7()])
    plain = _write(tmp_path / "p.vcf", [_b7()])
    assert a.read_bytes() == b.read_bytes()
    assert gzip.decompress(a.read_bytes()) == plain.read_bytes()


# --- Records: minus strand, multi-allelic, mixed, columns ---------------------------------------


def test_minus_strand_codon_written_as_given_forward_strand(tmp_path) -> None:
    # Gene on the minus strand: residue codon GCT (Ala) is AGC on the forward strand; S via TCT
    # is AGA on the forward strand. The writer writes the forward-strand bases verbatim.
    info = dict(B7_INFO, GENE="NISS_00099")
    p = _write(tmp_path / "v.vcf", [_b7(pos0=500, ref="AGC", alts=("AGA",), info=info)])
    assert _records(p)[0].split("\t")[1:5] == ["501", ".", "AGC", "AGA"]


def test_multi_allelic_and_mixed(tmp_path) -> None:
    info = dict(B7_INFO, AA_ALT=("S", "T"), CALL="mixed", NT_CHANGES=(1, 1))
    rec = _b7(alts=("TCT", "ACT"), info=info, gt=".", ad=(2, 2, 1), dp=5)
    cols = _records(_write(tmp_path / "v.vcf", [rec]))[0].split("\t")
    assert cols[4] == "TCT,ACT"
    assert "AA_ALT=S,T" in cols[7].split(";")
    assert "NT_CHANGES=1,1" in cols[7].split(";")
    assert cols[8:] == ["GT:AD:DP", ".:2,2,1:5"]


def test_id_qual_filter_columns(tmp_path) -> None:
    rec = _b7(id="v1", qual="30", filter="low_conf")
    cols = _records(_write(tmp_path / "v.vcf", [rec]))[0].split("\t")
    assert cols[2] == "v1"
    assert cols[5:7] == ["30", "low_conf"]


def test_info_float_six_significant_digits(tmp_path) -> None:
    info = dict(B7_INFO, CLASS_ERR=0.0123456789)
    cols = _records(_write(tmp_path / "v.vcf", [_b7(info=info)]))[0].split("\t")
    assert cols[7].endswith(";CLASS_ERR=0.0123457")


def test_info_insertion_order_kept(tmp_path) -> None:
    info = {"CALL": "ALT", "GENE": "g", "AA_POS": 3}
    cols = _records(_write(tmp_path / "v.vcf", [_b7(info=info)]))[0].split("\t")
    assert cols[7] == "CALL=ALT;GENE=g;AA_POS=3"


@pytest.mark.parametrize(
    ("raw", "escaped"),
    [
        ("a:b", "a%3Ab"),
        ("a;b", "a%3Bb"),
        ("a=b", "a%3Db"),
        ("50%", "50%25"),
        ("a,b", "a%2Cb"),
        ("a b", "a%20b"),
        ("a\tb", "a%09b"),
        ("a\nb", "a%0Ab"),
        ("a\rb", "a%0Db"),
    ],
)
def test_info_value_percent_encoding(tmp_path, raw: str, escaped: str) -> None:
    info = {"GENE": raw}
    cols = _records(_write(tmp_path / "v.vcf", [_b7(info=info)]))[0].split("\t")
    assert cols[7] == f"GENE={escaped}"


def test_flags_items_escaped_then_joined_by_pipe(tmp_path) -> None:
    info = {"FLAGS": ("a b", "c;d")}
    cols = _records(_write(tmp_path / "v.vcf", [_b7(info=info)]))[0].split("\t")
    assert cols[7] == "FLAGS=a%20b|c%3Bd"


# --- Edge: no records ----------------------------------------------------------------------------


def test_no_records_header_only(tmp_path) -> None:
    p = _write(tmp_path / "v.vcf", [])
    lines = _lines(p)
    assert lines[0] == "##fileformat=VCFv4.3"
    assert lines[-1] == CHROM_LINE
    assert _records(p) == []


# --- Example 8: pysam parses the output (§9.5) --------------------------------------------------

needs_pysam = pytest.mark.skipif(
    importlib.util.find_spec("pysam") is None, reason="pysam (dev extra) not installed"
)


@needs_pysam
def test_example_08_pysam_parses(tmp_path) -> None:
    import pysam

    info = dict(B7_INFO, AA_ALT=("S", "T"), CALL="mixed", NT_CHANGES=(1, 1))
    recs = [
        _b7(),
        _b7(pos0=2003),
        _b7(pos0=30, alts=("TCT", "ACT"), info=info, gt=".", ad=(2, 2, 1), dp=5),
        _b7(contig="lin", pos0=10),
    ]
    p = _write(tmp_path / "v.vcf", recs, replicons=(C, L))
    with pysam.VariantFile(str(p)) as vf:
        header = str(vf.header)
        assert "##contig=<ID=chr,length=2000>" in header
        assert f"##reference=sha256:{SHA}" in header
        assert list(vf.header.samples) == ["pilot_A3"]
        got = list(vf)
    assert [(r.contig, r.pos) for r in got] == [("chr", 4), ("chr", 31), ("chr", 1001), ("lin", 11)]
    b7 = got[2]
    assert b7.pos == 1001
    assert (b7.ref, b7.alts) == ("GCT", ("TCT",))
    assert b7.info["GENE"] == "NISS_00012"
    assert b7.info["FLAGS"] == "mass_equivalent|snv"
    assert b7.samples["pilot_A3"]["AD"] == (0, 3)
    assert b7.samples["pilot_A3"]["DP"] == 3
    assert got[1].alts == ("TCT", "ACT")
    assert got[1].samples["pilot_A3"]["AD"] == (2, 2, 1)


@needs_pysam
def test_example_08_pysam_parses_empty_vcf(tmp_path) -> None:
    import pysam

    p = _write(tmp_path / "v.vcf", [])
    with pysam.VariantFile(str(p)) as vf:
        assert list(vf) == []
        assert list(vf.header.samples) == ["pilot_A3"]
