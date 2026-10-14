"""Lead decisions after the M1-05 review (ADR-0005 item 1, "After the M1-05 review").

GFF3: CDS must carry a phase (GFF3 3.1.26); reader accepts strand ``?``, rejects duplicate
attribute keys, invalid percent-escapes and non-plain-digit coordinates; seqids must be safe IDs.
VCF: INFO/FORMAT keys declared; NaN -> ``.``; bool -> Flag; empty lists and ``|`` in FLAGS items
are errors; REF/ALT uppercase ``[ACGTN]+``, ALT != REF, same length; GT within the ALT count;
no control characters in header text; safe contig IDs; optional REF check against ``sequences``.
Every VCF written here must parse with pysam without htslib warnings.
"""

from __future__ import annotations

import importlib.util
import math
import warnings

import pytest

from giap.io.gff3 import GffFeature, read_gff3, write_gff3
from giap.io.vcf import VcfHeaderField, VcfRecord, write_vcf
from giap.model import Interval, Replicon

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
    VcfHeaderField("SCORE", "1", "Float", "A float score"),
    VcfHeaderField("MASS_EQ", "0", "Flag", "Mass-equivalent substitution"),
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

needs_pysam = pytest.mark.skipif(
    importlib.util.find_spec("pysam") is None, reason="pysam (dev extra) not installed"
)


def _b7(**kw) -> VcfRecord:
    base = dict(
        contig="chr", pos0=1000, ref="GCT", alts=("TCT",), info=B7_INFO, gt="1", ad=(0, 3), dp=3
    )
    base.update(kw)
    return VcfRecord(**base)


def _vcf(path, records, replicons=(C,), **kw):
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


def _records(path) -> list[str]:
    return [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]


def _info_col(path) -> str:
    return _records(path)[0].split("\t")[7]


def _pysam_clean(path, capfd) -> list:
    """Open and iterate ``path`` with pysam; fail on Python warnings or htslib stderr output."""
    pysam = pytest.importorskip("pysam")
    capfd.readouterr()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pysam.VariantFile(str(path)) as vf:
            recs = list(vf)
    err = capfd.readouterr().err
    assert not caught, [str(w.message) for w in caught]
    assert err == "", err
    return recs


def _genome(pos0: int = 1000, codon: str = "GCT", length: int = 2000) -> dict[str, str]:
    seq = ["A"] * length
    seq[pos0 : pos0 + len(codon)] = codon
    return {"chr": "".join(seq)}


# =================================================================================================
# GFF3
# =================================================================================================


def _gff_text(tmp_path, body: str):
    p = tmp_path / "in.gff3"
    p.write_text("##gff-version 3\n" + body, encoding="utf-8")
    return p


def test_gff3_cds_without_phase_raises_on_write(tmp_path) -> None:
    f = GffFeature(Interval("chr", 0, 9, "+"), "giap", "CDS", attributes={"ID": "x"})
    with pytest.raises(ValueError):
        write_gff3(tmp_path / "a.gff3", [f], [C])


def test_gff3_non_cds_without_phase_still_written_dot(tmp_path) -> None:
    f = GffFeature(Interval("chr", 0, 9, "+"), "giap", "gene", attributes={"ID": "x"})
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [C])
    assert "chr\tgiap\tgene\t1\t9\t.\t+\t.\tID=x" in p.read_text().splitlines()


def test_gff3_reader_reads_strand_question_mark_as_unknown(tmp_path) -> None:
    p = _gff_text(tmp_path, "chr\t.\tgene\t1\t9\t.\t?\t.\tID=x\n")
    (f,) = read_gff3(p, {"chr": C})
    assert f.interval == Interval("chr", 0, 9, ".")
    assert f.attributes == {"ID": "x"}


def test_gff3_reader_lenient_cds_without_phase(tmp_path) -> None:
    # Only the writer requires a CDS phase; third-party files with phase "." are read.
    p = _gff_text(tmp_path, "chr\tProdigal\tCDS\t1\t9\t.\t+\t.\tID=x\n")
    (f,) = read_gff3(p, {"chr": C})
    assert (f.type, f.phase) == ("CDS", None)


@pytest.mark.parametrize("attrs", ["ID=a;ID=b", "ID=a;Note=x;Note=y", "ID=a;Note=x;Note=x"])
def test_gff3_reader_duplicate_attribute_key_raises(tmp_path, attrs: str) -> None:
    p = _gff_text(tmp_path, f"chr\t.\tgene\t1\t9\t.\t+\t.\t{attrs}\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


@pytest.mark.parametrize(
    "attrs",
    [
        "ID=a;Note=%FF",  # not valid UTF-8
        "ID=a;Note=x%C3",  # truncated UTF-8 sequence
        "ID=a;Note=%G1",  # not hex
        "ID=a;Note=50%",  # bare % at end
        "ID=a;Note=%4",  # incomplete escape
        "ID=a;No%ZZte=x",  # bad escape in the key
    ],
)
def test_gff3_reader_invalid_percent_escape_raises(tmp_path, attrs: str) -> None:
    p = _gff_text(tmp_path, f"chr\t.\tgene\t1\t9\t.\t+\t.\t{attrs}\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


def test_gff3_reader_valid_utf8_percent_escape_decoded(tmp_path) -> None:
    p = _gff_text(tmp_path, "chr\t.\tgene\t1\t9\t.\t+\t.\tID=a;Note=%C3%A9%3B\n")
    (f,) = read_gff3(p, {"chr": C})
    assert f.attributes["Note"] == "é;"


@pytest.mark.parametrize("bad", ["+5", " 5", "5 ", "1_000", chr(0x665), "5.0", "0x5", "-0"])
@pytest.mark.parametrize("column", ["start", "end"])
def test_gff3_reader_rejects_non_plain_digit_coordinates(tmp_path, bad: str, column: str) -> None:
    start, end = (bad, "1200") if column == "start" else ("1", bad)
    p = _gff_text(tmp_path, f"chr\t.\tgene\t{start}\t{end}\t.\t+\t.\tID=x\n")
    with pytest.raises(ValueError):
        read_gff3(p, {"chr": C})


@pytest.mark.parametrize("seqid", ["chr 1", "chr;1", "chr\t1", "chr|1", "chr,1", "chr=1", ""])
def test_gff3_unsafe_seqid_raises_on_write(tmp_path, seqid: str) -> None:
    with pytest.raises(ValueError):
        rep = Replicon(seqid, 100, "linear")
        f = GffFeature(Interval(seqid, 0, 9, "+"), "giap", "gene", attributes={"ID": "x"})
        write_gff3(tmp_path / "a.gff3", [f], [rep])


@pytest.mark.parametrize("seqid", ["chr 1", "chr;1", "chr,1"])
def test_gff3_unsafe_seqid_raises_on_write_without_features(tmp_path, seqid: str) -> None:
    # The ##sequence-region pragma alone must not carry an unsafe ID either.
    with pytest.raises(ValueError):
        write_gff3(tmp_path / "a.gff3", [], [Replicon(seqid, 100, "circular")])


ATTR_CASES = [
    {"Name": "été"},  # Latin-1 range
    {"product": "β-lactamase"},  # Greek
    {"Note": "中文;=,&%"},  # CJK plus reserved characters
    {"Note": "emoji \U0001f9ec"},  # astral plane
    {"Note": "50%C3%A9"},  # literal percent text that looks like an escape
    {"Dbxref": ("GO:0001", "EC:1.2.3.4")},  # list value
    {"Alias": ("é", "a,b", "c;d", "e=f", "%25")},  # list with non-ASCII and reserved
    {"Note": "tab\there", "Alias": ("line\nbreak", "\x7f")},  # control characters
    {"Name": "Ångström", "Alias": ("ΔG", "→")},
]


@pytest.mark.parametrize("extra", ATTR_CASES)
@pytest.mark.parametrize("ftype", ["gene", "CDS"])
def test_gff3_attribute_round_trip_non_ascii_and_lists(tmp_path, extra, ftype: str) -> None:
    attrs = {"ID": "f1", **extra}
    phase = 0 if ftype == "CDS" else None
    f = GffFeature(Interval("chr", 1950, 2050, "-"), "giap", ftype, phase=phase, attributes=attrs)
    p = tmp_path / "a.gff3"
    write_gff3(p, [f], [C])
    raw = p.read_bytes()
    raw.decode("utf-8")  # written as UTF-8
    (g,) = read_gff3(p, {"chr": C})
    assert g == f
    assert list(g.attributes) == list(attrs)
    for key, value in attrs.items():
        assert type(g.attributes[key]) is type(value)


# =================================================================================================
# VCF
# =================================================================================================


def test_vcf_b7_record_still_byte_exact(tmp_path) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7()])
    assert _records(p) == [B7_LINE]


@needs_pysam
def test_vcf_b7_record_pysam_clean(tmp_path, capfd) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7()])
    (rec,) = _pysam_clean(p, capfd)
    assert rec.info["FLAGS"] == "mass_equivalent|snv"


def test_vcf_undeclared_info_key_raises(tmp_path) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g", "UNDECLARED": "x"})])


def test_vcf_b7_with_no_info_header_raises(tmp_path) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7()], info_fields=[])


@pytest.mark.parametrize("missing", ["GT", "AD", "DP"])
def test_vcf_format_field_must_be_declared(tmp_path, missing: str) -> None:
    fmt = [f for f in FORMAT_FIELDS if f.id != missing]
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7()], format_fields=fmt)


def test_vcf_no_format_header_raises(tmp_path) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7()], format_fields=[])


def test_vcf_nan_info_written_as_dot(tmp_path, capfd) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g", "SCORE": math.nan})])
    assert _info_col(p) == "GENE=g;SCORE=."
    if importlib.util.find_spec("pysam") is not None:
        (rec,) = _pysam_clean(p, capfd)
        assert rec.info.get("SCORE") is None


def test_vcf_flag_true_written_as_bare_key(tmp_path, capfd) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g", "MASS_EQ": True, "AA_POS": 3})])
    assert _info_col(p) == "GENE=g;MASS_EQ;AA_POS=3"
    if importlib.util.find_spec("pysam") is not None:
        (rec,) = _pysam_clean(p, capfd)
        assert rec.info["MASS_EQ"] is True


def test_vcf_flag_false_omitted(tmp_path, capfd) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g", "MASS_EQ": False})])
    assert _info_col(p) == "GENE=g"
    if importlib.util.find_spec("pysam") is not None:
        (rec,) = _pysam_clean(p, capfd)
        assert "MASS_EQ" not in rec.info


def test_vcf_only_false_flag_gives_dot_info(tmp_path) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7(info={"MASS_EQ": False})])
    assert _info_col(p) == "."


@pytest.mark.parametrize("key", ["AA_ALT", "NT_CHANGES", "FLAGS"])
@pytest.mark.parametrize("empty", [(), []])
def test_vcf_empty_list_info_raises(tmp_path, key: str, empty) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g", key: empty})])


@pytest.mark.parametrize("items", [("a|b",), ("ok", "x|"), ("|",)])
def test_vcf_pipe_inside_flags_item_raises(tmp_path, items) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7(info={"FLAGS": items})])


@pytest.mark.parametrize(
    ("ref", "alts"),
    [
        ("gct", ("TCT",)),  # lowercase REF
        ("GCT", ("tct",)),  # lowercase ALT
        ("GCU", ("TCT",)),  # RNA base
        ("GCT", ("TCU",)),
        ("GCR", ("TCT",)),  # IUPAC ambiguity code
        ("GCT", ("*",)),
        ("GCT", ("<DEL>",)),
        ("GCT", ("",)),
        ("GCT", ("GCT",)),  # ALT == REF
        ("GCT", ("TCT", "GCT")),  # one ALT == REF
        ("GCT", ("TC",)),  # length differs
        ("GCT", ("TCTA",)),
        ("GCT", ("TCT", "AC")),
        ("G CT", ("TCT",)),
    ],
)
def test_vcf_bad_ref_alt_raises(tmp_path, ref: str, alts) -> None:
    ad = (0,) * (1 + len(alts))
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7(ref=ref, alts=alts, ad=ad, gt=".", info={"GENE": "g"})])


@needs_pysam
@pytest.mark.parametrize(("ref", "alt"), [("GCT", "TCT"), ("N", "A"), ("ACGTN", "TTTTT")])
def test_vcf_good_ref_alt_pysam_clean(tmp_path, capfd, ref: str, alt: str) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7(ref=ref, alts=(alt,), info={"GENE": "g"})])
    (rec,) = _pysam_clean(p, capfd)
    assert (rec.ref, rec.alts) == (ref, (alt,))


@pytest.mark.parametrize(
    ("gt", "alts"),
    [
        ("2", ("TCT",)),
        ("3", ("TCT", "ACT")),
        ("0", ("TCT",)),
        ("-1", ("TCT",)),
        ("1/1", ("TCT",)),
        ("a", ("TCT",)),
        ("", ("TCT",)),
        ("01", ("TCT",)),
    ],
)
def test_vcf_gt_out_of_range_raises(tmp_path, gt: str, alts) -> None:
    rec = _b7(alts=alts, ad=(0,) * (1 + len(alts)), gt=gt, info={"GENE": "g"})
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [rec])


@pytest.mark.parametrize("gt", [".", "1", "2"])
def test_vcf_gt_in_range_ok(tmp_path, capfd, gt: str) -> None:
    rec = _b7(alts=("TCT", "ACT"), ad=(1, 2, 3), dp=6, gt=gt, info={"GENE": "g"})
    p = _vcf(tmp_path / "v.vcf", [rec])
    assert _records(p)[0].split("\t")[9] == f"{gt}:1,2,3:6"
    if importlib.util.find_spec("pysam") is not None:
        _pysam_clean(p, capfd)


@pytest.mark.parametrize("bad", ["0.1\n", "0.1\r", "0.1\t", "0.1\x00", "0.1\x1b", "0.1\x7f"])
def test_vcf_control_chars_in_giap_version_raise(tmp_path, bad: str) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [], giap_version=bad)


@pytest.mark.parametrize("bad", ["a\nb", "a\rb", "a\tb", "a\x00b", "a\x1bb", "a\x7fb"])
@pytest.mark.parametrize("kind", ["info", "format"])
def test_vcf_control_chars_in_description_raise(tmp_path, bad: str, kind: str) -> None:
    field = VcfHeaderField("GENE" if kind == "info" else "GT", "1", "String", bad)
    if kind == "info":
        kw = {"info_fields": [field]}
    else:
        kw = {"format_fields": [field, *FORMAT_FIELDS[1:]]}
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [], **kw)


@needs_pysam
def test_vcf_description_with_quote_and_backslash_pysam_clean(tmp_path, capfd) -> None:
    field = VcfHeaderField("GENE", "1", "String", 'Locus "tag" \\ é')
    p = _vcf(tmp_path / "v.vcf", [_b7(info={"GENE": "g"})], info_fields=[field])
    _pysam_clean(p, capfd)


@pytest.mark.parametrize("cid", ["chr 1", "chr,1", "chr=1", "chr;1", "chr<1>", "chr\t1", ""])
def test_vcf_unsafe_contig_id_raises(tmp_path, cid: str) -> None:
    with pytest.raises(ValueError):
        rep = Replicon(cid, 2000, "circular")
        _vcf(tmp_path / "v.vcf", [_b7(contig=cid)], replicons=(rep,))


@pytest.mark.parametrize("cid", ["chr 1", "chr,1", "chr<1>"])
def test_vcf_unsafe_contig_id_raises_without_records(tmp_path, cid: str) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [], replicons=(Replicon(cid, 2000, "circular"),))


# --- Optional REF check against the genome ------------------------------------------------------


def test_vcf_b7_passes_with_matching_sequence(tmp_path, capfd) -> None:
    p = _vcf(tmp_path / "v.vcf", [_b7()], sequences=_genome())
    assert _records(p) == [B7_LINE]
    if importlib.util.find_spec("pysam") is not None:
        _pysam_clean(p, capfd)


def test_vcf_ref_mismatch_with_sequence_raises(tmp_path) -> None:
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7()], sequences=_genome(codon="GCA"))


def test_vcf_ref_checked_at_wrapped_pos(tmp_path) -> None:
    # pos0 2003 on a 2000-bp circular contig is POS 4 (0-based 3).
    _vcf(tmp_path / "ok.vcf", [_b7(pos0=2003)], sequences=_genome(pos0=3))
    with pytest.raises(ValueError):
        _vcf(tmp_path / "bad.vcf", [_b7(pos0=2003)], sequences=_genome(pos0=1000))


def test_vcf_ref_not_checked_without_sequences(tmp_path) -> None:
    # Without sequences the writer cannot check REF; the B.7 record is written as given.
    p = _vcf(tmp_path / "v.vcf", [_b7()])
    assert _records(p) == [B7_LINE]


def test_vcf_ref_check_off_by_one_detected(tmp_path) -> None:
    # Codon placed one base later (0-based 1001) must not match POS 1001 (0-based 1000).
    with pytest.raises(ValueError):
        _vcf(tmp_path / "v.vcf", [_b7()], sequences=_genome(pos0=1001))
