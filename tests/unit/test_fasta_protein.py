"""Protein FASTA reader (issue M1-03; GIAP_SPEC.md §5.3, §7.0; ADR-0005 item 1)."""

from __future__ import annotations

import dataclasses
import gzip
from pathlib import Path

import pytest

from giap.io.fasta import FastaFormatError
from giap.io.protein_fasta import (
    COUNTED_AA,
    STANDARD_AA,
    ProteinRecord,
    read_protein_fasta,
    resolve_protein_id,
    route_by_length,
)


def _write(tmp_path: Path, text: str, name: str = "p.faa") -> Path:
    p = tmp_path / name
    p.write_bytes(text.encode("utf-8"))
    return p


def _read(tmp_path: Path, text: str) -> list[ProteinRecord]:
    return read_protein_fasta(_write(tmp_path, text))


# --- Constants and record type ------------------------------------------------------------


def test_standard_aa_is_20_plus_u() -> None:
    assert frozenset("ACDEFGHIKLMNPQRSTVWYU") == STANDARD_AA
    assert len(STANDARD_AA) == 21


def test_counted_aa_is_xbzj() -> None:
    assert frozenset("XBZJ") == COUNTED_AA
    assert not (STANDARD_AA & COUNTED_AA)
    assert "O" not in STANDARD_AA | COUNTED_AA


def test_protein_record_is_frozen(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMKV\n")
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.id = "b"  # type: ignore[misc]


# --- Examples 1-9 (one test each) ---------------------------------------------------------


def test_example1_swissprot_header(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">sp|P0A7V8|RS4_ECOLI 30S ribosomal protein S4\nMARYLGPKLK\n")
    assert r.id == "P0A7V8"
    assert r.members == ("P0A7V8",)
    assert r.header == "sp|P0A7V8|RS4_ECOLI 30S ribosomal protein S4"
    assert r.sequence == "MARYLGPKLK"


def test_example2_trembl_header(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">tr|A0A0K|A0A0K_ECOLI\nMKV\n")
    assert r.id == "A0A0K"
    assert r.members == ("A0A0K",)


def test_example3_plain_group(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">NISS_00012;NISS_00013 group\nMKV\n")
    assert r.id == "NISS_00012;NISS_00013"
    assert r.members == ("NISS_00012", "NISS_00013")


def test_example4_uniprot_group_resolves_each_member(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">sp|P1|X_ECOLI;sp|P2|Y_ECOLI\nMKV\n")
    assert r.id == "P1;P2"
    assert r.members == ("P1", "P2")


def test_example5_trailing_stop_removed(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMKV*\n")
    assert r.sequence == "MKV"
    assert r.internal_stops == 0


def test_example6_internal_stop_and_x(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMK*VLX\n")
    assert r.internal_stops == 1
    assert dict(r.nonstandard) == {"X": 1}
    assert dict(r.other) == {}
    assert r.sequence == "MK*VLX"


def test_example7_lower_case_upper_cased_u_is_standard(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nmkub\n")
    assert r.sequence == "MKUB"
    assert dict(r.nonstandard) == {"B": 1}
    assert dict(r.other) == {}


def test_example8_route_by_length_default_boundaries() -> None:
    assert [route_by_length(n) for n in (6, 7, 40, 41)] == [
        "dropped",
        "peptide",
        "peptide",
        "protein",
    ]


def test_example9_gzip_matches_plain(tmp_path: Path) -> None:
    text = ">sp|P1|X_ECOLI first\nMKV*\n>b\nmk*vlx\n"
    plain = _write(tmp_path, text, "p.faa")
    gz = tmp_path / "p.faa.gz"
    gz.write_bytes(gzip.compress(text.encode("utf-8")))
    got_plain = read_protein_fasta(plain)
    assert len(got_plain) == 2
    assert read_protein_fasta(gz) == got_plain


# --- resolve_protein_id -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("token", "expected"),
    [
        ("sp|P0A7V8|RS4_ECOLI", ("P0A7V8", ("P0A7V8",))),
        ("tr|A0A0K|A0A0K_ECOLI", ("A0A0K", ("A0A0K",))),
        ("NISS_00012;NISS_00013", ("NISS_00012;NISS_00013", ("NISS_00012", "NISS_00013"))),
        ("sp|P1|X_ECOLI;sp|P2|Y_ECOLI", ("P1;P2", ("P1", "P2"))),
        ("sp|P1|X;NISS_1", ("P1;NISS_1", ("P1", "NISS_1"))),
        ("plain", ("plain", ("plain",))),
        ("sp|ACC", ("sp|ACC", ("sp|ACC",))),
        ("A;;B;", ("A;B", ("A", "B"))),
        (" A ; B ", ("A;B", ("A", "B"))),
    ],
)
def test_resolve_protein_id(token: str, expected: tuple[str, tuple[str, ...]]) -> None:
    assert resolve_protein_id(token) == expected


def test_sp_without_third_field_keeps_raw_id(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">sp|ACC desc\nMKV\n")
    assert r.id == "sp|ACC"
    assert r.members == ("sp|ACC",)


def test_id_is_first_header_token_only(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">prot1 some;thing sp|Q|R\nMKV\n")
    assert r.id == "prot1"
    assert r.members == ("prot1",)
    assert r.header == "prot1 some;thing sp|Q|R"


# --- Sequences and alphabet ---------------------------------------------------------------


def test_multiline_sequence_joined(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMKV\nLLA\n\nGG*\n")
    assert r.sequence == "MKVLLAGG"
    assert r.internal_stops == 0


def test_only_one_trailing_stop_removed(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMKV**\n")
    assert r.sequence == "MKV*"
    assert r.internal_stops == 1


def test_only_stop_gives_empty_sequence_and_dropped_route(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\n*\n")
    assert r.sequence == ""
    assert r.internal_stops == 0
    assert route_by_length(len(r.sequence)) == "dropped"


def test_routing_length_excludes_trailing_stop(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\n" + "M" * 40 + "*\n")
    assert len(r.sequence) == 40
    assert route_by_length(len(r.sequence)) == "peptide"


def test_all_counted_letters(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nMXXBZJjx\n")
    assert dict(r.nonstandard) == {
        "X": 3,
        "B": 1,
        "Z": 1,
        "J": 2,
    }  # j is upper-cased (spec: sequence is upper case)
    assert dict(r.other) == {}


def test_pyrrolysine_digits_and_dash_counted_as_other_record_kept(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a\nMKoO1-2-.\n>b\nMKV\n")
    assert [r.id for r in recs] == ["a", "b"]
    a = recs[0]
    assert dict(a.other) == {"O": 2, "1": 1, "2": 1, "-": 2, ".": 1}
    assert dict(a.nonstandard) == {}
    assert a.internal_stops == 0


def test_standard_sequence_has_no_counts(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">a\nACDEFGHIKLMNPQRSTVWYU\n")
    assert dict(r.nonstandard) == {}
    assert dict(r.other) == {}
    assert r.internal_stops == 0


def test_short_records_are_returned_not_dropped_by_reader(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">short\nMKV\n>pep\nMKVLLAG\n>prot\n" + "M" * 41 + "\n")
    assert [r.id for r in recs] == ["short", "pep", "prot"]
    assert [route_by_length(len(r.sequence)) for r in recs] == ["dropped", "peptide", "protein"]


# --- Duplicates, empty file and format errors ---------------------------------------------


def test_duplicate_ids_returned_as_is(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a\nMKV\n>a\nLLL\n")
    assert [(r.id, r.sequence) for r in recs] == [("a", "MKV"), ("a", "LLL")]


def test_duplicate_sequences_returned_as_is(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a\nMKV\n>b\nMKV\n")
    assert [(r.id, r.sequence) for r in recs] == [("a", "MKV"), ("b", "MKV")]


def test_empty_file_returns_empty_list(tmp_path: Path) -> None:
    assert _read(tmp_path, "") == []


@pytest.mark.parametrize("header", [">", ">   ", ">;", ">; ;", ">; desc"])
def test_header_without_members_raises(tmp_path: Path, header: str) -> None:
    with pytest.raises(FastaFormatError):
        _read(tmp_path, f">ok\nMKV\n{header}\nMKV\n")


def test_sequence_before_header_raises(tmp_path: Path) -> None:
    with pytest.raises(FastaFormatError):
        _read(tmp_path, "MKV\n>a\nMKV\n")


# --- route_by_length boundaries -----------------------------------------------------------


@pytest.mark.parametrize(
    ("length", "expected"),
    [(0, "dropped"), (6, "dropped"), (7, "peptide"), (40, "peptide"), (41, "protein")],
)
def test_route_by_length_default(length: int, expected: str) -> None:
    assert route_by_length(length) == expected


@pytest.mark.parametrize(
    ("length", "expected"),
    [
        (6, "dropped"),
        (7, "dropped"),
        (8, "peptide"),
        (30, "peptide"),
        (31, "protein"),
        (40, "protein"),
        (41, "protein"),
    ],
)
def test_route_by_length_custom_bounds(length: int, expected: str) -> None:
    assert route_by_length(length, peptide_len=(8, 30)) == expected
