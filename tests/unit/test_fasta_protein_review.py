"""Review pins for the protein FASTA reader (issue M1-03; GIPPA_SPEC.md §5.3, §7.0)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from gippa.io.fasta import FastaFormatError
from gippa.io.protein_fasta import ProteinRecord, read_protein_fasta


def _write_bytes(tmp_path: Path, data: bytes, name: str = "p.faa") -> Path:
    p = tmp_path / name
    p.write_bytes(data)
    return p


def _read(tmp_path: Path, text: str) -> list[ProteinRecord]:
    return read_protein_fasta(_write_bytes(tmp_path, text.encode("utf-8")))


# --- 1. error line numbers ------------------------------------------------------------------


def test_invalid_utf8_small_file_line_no(tmp_path: Path) -> None:
    data = b">P1\nMKV\nMKV\nMK\xffV\n>P2\nMKV\n"
    p = _write_bytes(tmp_path, data)
    with pytest.raises(FastaFormatError) as info:
        read_protein_fasta(p)
    assert info.value.line_no == 4


def test_invalid_utf8_large_file_line_no(tmp_path: Path) -> None:
    lines = [b">P1"] + [b"MKVLAAGIVG" * 6] * 3000  # lines 1..3001
    lines.append(b"MK\xffV")  # line 3002
    lines.append(b"MKV")
    p = _write_bytes(tmp_path, b"\n".join(lines) + b"\n")
    with pytest.raises(FastaFormatError) as info:
        read_protein_fasta(p)
    assert info.value.line_no == 3002


def test_invalid_utf8_gzip_line_no(tmp_path: Path) -> None:
    data = b">P1\nMKV\nMKV\nMK\xffV\n"
    p = _write_bytes(tmp_path, gzip.compress(data), name="p.faa.gz")
    with pytest.raises(FastaFormatError) as info:
        read_protein_fasta(p)
    assert info.value.line_no == 4


@pytest.mark.parametrize(
    ("text", "line_no"),
    [
        (">\nMKV\n", 1),
        (">P1\nMKV\n\n>   \nMKV\n", 4),
        (">P1\nMKV\n>;;\nMKV\n", 3),
        (">P1\nMKV\n>P2\nMKV\n>;\nMKV\n", 5),
    ],
)
def test_no_protein_id_line_no(tmp_path: Path, text: str, line_no: int) -> None:
    with pytest.raises(FastaFormatError, match="no protein ID") as info:
        _read(tmp_path, text)
    assert info.value.line_no == line_no


# --- 2. non-ASCII letters -------------------------------------------------------------------


def test_non_ascii_kept_and_counted_as_other(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">P1\nMKß\n")
    assert rec.sequence == "MKß"
    assert dict(rec.other) == {"ß": 1}
    assert dict(rec.nonstandard) == {}


# --- 3. whitespace --------------------------------------------------------------------------


def test_internal_whitespace_removed(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">P1\nMK V\tL A\nAG \t G\n")
    assert rec.sequence == "MKVLAAGG"
    assert dict(rec.other) == {}


def test_header_leading_whitespace_resolves_id(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">  sp|P1|X  desc\nMKV\n")
    assert rec.id == "P1"
    assert rec.members == ("P1",)


# --- 4. BOM, line endings, truncated gzip ---------------------------------------------------

_LF = ">sp|P1|X desc\nMKVL\nAAG\n>Q2\nMKV*\n"


def _summary(records: list[ProteinRecord]) -> list[tuple[str, tuple[str, ...], str, str]]:
    return [(r.id, r.members, r.header, r.sequence) for r in records]


def test_bom_ignored(tmp_path: Path) -> None:
    expected = _summary(_read(tmp_path, _LF))
    p = _write_bytes(tmp_path, b"\xef\xbb\xbf" + _LF.encode(), name="bom.faa")
    got = read_protein_fasta(p)
    assert _summary(got) == expected
    assert got[0].id == "P1"


@pytest.mark.parametrize("eol", ["\r\n", "\r"])
def test_line_endings_equivalent(tmp_path: Path, eol: str) -> None:
    expected = _summary(_read(tmp_path, _LF))
    p = _write_bytes(tmp_path, _LF.replace("\n", eol).encode(), name="eol.faa")
    assert _summary(read_protein_fasta(p)) == expected


def test_truncated_gzip_is_format_error(tmp_path: Path) -> None:
    body = ("".join(f">P{i}\n{'MKVLAAGIVG' * 8}\n" for i in range(500))).encode()
    blob = gzip.compress(body)
    p = _write_bytes(tmp_path, blob[: len(blob) // 2], name="t.faa.gz")
    with pytest.raises(FastaFormatError):
        read_protein_fasta(p)


# --- 5. lead decisions ----------------------------------------------------------------------


def test_record_is_hashable_and_counts_compare_as_dicts(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">P1\nMKXXBß*Q\n")
    hash(rec)  # must not raise
    assert dict(rec.nonstandard) == {"X": 2, "B": 1}
    assert dict(rec.other) == {"ß": 1}
    assert rec.internal_stops == 1
    (again,) = _read(tmp_path, ">P1\nMKXXBß*Q\n")
    assert rec == again
    assert hash(rec) == hash(again)
    assert len({rec, again}) == 1


def test_duplicate_group_members_deduplicated(tmp_path: Path) -> None:
    (rec,) = _read(tmp_path, ">sp|P1|X;sp|P1|Y;Q2\nMKV\n")
    assert rec.members == ("P1", "Q2")
    assert rec.id == "P1;Q2"
