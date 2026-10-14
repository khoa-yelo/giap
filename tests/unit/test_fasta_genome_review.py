"""Code-review findings on the genome FASTA reader/writer (M1-02 review; GIAP_SPEC.md §5.3)."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import pytest

from giap.io.fasta import FastaFormatError, invalid_letters, read_genome_fasta, write_fasta


def _write_bytes(tmp_path: Path, data: bytes, name: str = "g.fa") -> Path:
    p = tmp_path / name
    p.write_bytes(data)
    return p


def _read(tmp_path: Path, text: str) -> list:
    return read_genome_fasta(_write_bytes(tmp_path, text.encode("utf-8")))


# --- 1. Non-ASCII letters are not upper-cased into valid bases -----------------------------


@pytest.mark.parametrize(("raw", "bad"), [("AC\u00df", "\u00df"), ("AC\u017f", "\u017f")])
def test_non_ascii_letter_not_uppercased_into_bases(tmp_path: Path, raw: str, bad: str) -> None:
    (rec,) = _read(tmp_path, f">chr1\n{raw}\n")
    assert rec.sequence not in {"ACSS", "ACS"}
    assert rec.sequence.startswith("AC")
    assert invalid_letters(rec.sequence), "non-ASCII letter must be reported as invalid"
    assert bad in invalid_letters(rec.sequence)


# --- 2. UTF-8 BOM at file start is ignored -------------------------------------------------


@pytest.mark.parametrize("compress", [False, True])
def test_utf8_bom_ignored(tmp_path: Path, compress: bool) -> None:
    data = "﻿>chr1 desc\nACGT\n".encode()
    if compress:
        data = gzip.compress(data)
    recs = read_genome_fasta(_write_bytes(tmp_path, data))
    assert [r.id for r in recs] == ["chr1"]
    assert recs[0].original_id == "chr1"
    assert recs[0].sequence == "ACGT"


# --- 3. Decode error reports no wrong line number ------------------------------------------


def test_decode_error_line_number_not_wrong(tmp_path: Path) -> None:
    lines = [b">chr1"] + [b"ACGTACGTAC"] * 9999 + [b"ACG\xffT"] + [b"ACGT"] * 10
    path = _write_bytes(tmp_path, b"\n".join(lines) + b"\n")  # bad byte on line 10001
    with pytest.raises(FastaFormatError) as ei:
        read_genome_fasta(path)
    msg = str(ei.value)
    assert "line 8193" not in msg
    for n in re.findall(r"line (\d+)", msg):
        assert int(n) == 10001


# --- 4. write_fasta rejects unsafe IDs and multi-line sequences ----------------------------


@pytest.mark.parametrize("bad_id", ["a b", "x\ny", ""])
def test_write_fasta_rejects_unsafe_id(tmp_path: Path, bad_id: str) -> None:
    with pytest.raises(ValueError):
        write_fasta(tmp_path / "o.fa", [(bad_id, "ACGT")])


def test_write_fasta_rejects_newline_in_sequence(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        write_fasta(tmp_path / "o.fa", [("chr1", "AC\nGT")], width=80)


# --- 5. Mutation gaps -----------------------------------------------------------------------


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("chr1 xcircular=true", False),
        ("chr1 circular=trueish", False),
        ("chr1 topology=circularish", False),
        ("chr1 [topology=circularish]", False),
        ("chr1 desc circular=true", True),
        ("chr1 circular=true more words", True),
        ("chr1 desc [topology=circular]", True),
        ("chr1 [topology=circular] more words", True),
        ("chr1 circular=true;", True),
        ("chr1 circular=true, more", True),
        ("chr1 [circular=true]", True),
    ],
)
def test_circular_tag_boundaries(tmp_path: Path, header: str, expected: bool) -> None:
    (rec,) = _read(tmp_path, f">{header}\nACGT\n")
    assert rec.circular_tag is expected


def test_corrupt_deflate_stream_is_format_error(tmp_path: Path) -> None:
    payload = "".join(f">chr{i}\n{'ACGT' * 500}\n" for i in range(50)).encode()
    data = bytearray(gzip.compress(payload))
    mid = len(data) // 2
    for i in range(mid, mid + 64):
        data[i] ^= 0xFF
    path = _write_bytes(tmp_path, bytes(data), "g.fa.gz")
    with pytest.raises(FastaFormatError):
        read_genome_fasta(path)


def test_header_whitespace_stripped(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">chr1 desc   \nACGT\n")
    assert rec.header == "chr1 desc"
    assert rec.id == "chr1"


def test_leading_whitespace_after_gt(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">  chr1 desc\nACGT\n")
    assert rec.id == "chr1"
    assert rec.original_id == "chr1"


def test_multi_member_gzip_reads_all_records(tmp_path: Path) -> None:
    data = gzip.compress(b">chr1\nACGT\n") + gzip.compress(b">chr2\nGGCC\n")
    recs = read_genome_fasta(_write_bytes(tmp_path, data, "g.fa.gz"))
    assert [(r.id, r.sequence) for r in recs] == [("chr1", "ACGT"), ("chr2", "GGCC")]


def test_gt_inside_sequence_line_is_invalid_letter(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">chr1\nAC>GT\n")
    assert len(recs) == 1
    assert recs[0].sequence == "AC>GT"
    assert invalid_letters(recs[0].sequence) == {">": 1}
