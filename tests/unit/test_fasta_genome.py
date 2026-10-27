"""Genome FASTA reader and writer (issue M1-02; GIPPA_SPEC.md §5.3, §7.0; ADR-0005 item 1)."""

from __future__ import annotations

import dataclasses
import gzip
import io
from pathlib import Path

import pytest

import gippasim
from gippa.io._open import open_text
from gippa.io.fasta import (
    IUPAC_DNA,
    SAFE_ID,
    FastaFormatError,
    GenomeRecord,
    invalid_letters,
    read_genome_fasta,
    sanitize_id,
    write_fasta,
)

EX1 = ">chr1 Escherichia coli circular=true\nACGT\nacgn\n"


def _write(tmp_path: Path, text: str, name: str = "g.fa") -> Path:
    p = tmp_path / name
    p.write_bytes(text.encode("utf-8"))
    return p


def _read(tmp_path: Path, text: str) -> list[GenomeRecord]:
    return read_genome_fasta(_write(tmp_path, text))


# --- Examples 1-9 (one test each) ---------------------------------------------------------


def test_example1_basic_record_with_circular_tag(tmp_path: Path) -> None:
    recs = _read(tmp_path, EX1)
    assert len(recs) == 1
    r = recs[0]
    assert r.id == "chr1"
    assert r.original_id == "chr1"
    assert r.header == "chr1 Escherichia coli circular=true"
    assert r.sequence == "ACGTACGN"
    assert r.circular_tag is True


def test_example2_topology_bracket_tag(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">p1 [topology=circular] [gcode=11]\nGG\n")
    assert r.id == "p1"
    assert r.circular_tag is True
    assert r.sequence == "GG"


def test_example3_no_tag_is_not_circular(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">ctg_7 len=5000\nAT\n")
    assert r.id == "ctg_7"
    assert r.circular_tag is False


def test_example4_unsafe_id_characters_mapped(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">NZ_CP0|chr:1 x\nA\n")
    assert r.id == "NZ_CP0_chr_1"
    assert r.original_id == "NZ_CP0|chr:1"
    assert r.header == "NZ_CP0|chr:1 x"


def test_example5_gzip_detected_by_magic_bytes(tmp_path: Path) -> None:
    plain = read_genome_fasta(_write(tmp_path, EX1, "plain.fa"))
    gz = tmp_path / "g.fa"  # no .gz suffix
    gz.write_bytes(gzip.compress(EX1.encode("utf-8")))
    assert gz.read_bytes()[:2] == b"\x1f\x8b"
    assert read_genome_fasta(gz) == plain


def test_example6_text_before_first_header_is_error(tmp_path: Path) -> None:
    with pytest.raises(FastaFormatError) as exc:
        _read(tmp_path, "ACGT\n>c\nA\n")
    assert exc.value.line_no == 1


def test_example7_invalid_letters_reported_not_rejected(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">c\nACGTXU-\n")
    assert r.sequence == "ACGTXU-"
    assert invalid_letters("ACGTXU-") == {"X": 1, "U": 1, "-": 1}


def test_example8_empty_file_gives_no_records(tmp_path: Path) -> None:
    assert _read(tmp_path, "") == []


def test_example9_write_fasta_wraps_at_80(tmp_path: Path) -> None:
    p = tmp_path / "out.fa"
    write_fasta(p, [("c1", "A" * 85)])
    assert p.read_bytes() == b">c1\n" + b"A" * 80 + b"\n" + b"A" * 5 + b"\n"


# --- Reader rules and edge cases ----------------------------------------------------------


def test_fasta_format_error_is_value_error() -> None:
    assert issubclass(FastaFormatError, ValueError)


def test_genome_record_is_frozen_slots_dataclass() -> None:
    r = GenomeRecord(id="a", original_id="a", header="a", sequence="A", circular_tag=False)
    assert dataclasses.is_dataclass(r)
    assert not hasattr(r, "__dict__")
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.id = "b"  # type: ignore[misc]


def test_error_line_no_is_one_based_after_blank_lines(tmp_path: Path) -> None:
    with pytest.raises(FastaFormatError) as exc:
        _read(tmp_path, "\n\nNNNN\n>c\nA\n")
    assert exc.value.line_no == 3


def test_blank_lines_before_header_skipped(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, "\n  \n>c\nAC\n\nGT\n\n")
    assert r.sequence == "ACGT"


def test_crlf_line_endings_accepted(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">c1 circular=true\r\nAC\r\nGT\r\n")
    assert r.id == "c1"
    assert r.header == "c1 circular=true"
    assert r.sequence == "ACGT"
    assert r.circular_tag is True


def test_sequence_whitespace_removed_and_upper_cased(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">c\nac gt\tnn \nrY\n")
    assert r.sequence == "ACGTNNRY"


def test_multiple_records_keep_file_order(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">z\nA\n>a\nC\n>m\nG\n")
    assert [r.id for r in recs] == ["z", "a", "m"]
    assert [r.sequence for r in recs] == ["A", "C", "G"]


def test_last_record_without_trailing_newline(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a\nAC\n>b\nGT")
    assert [(r.id, r.sequence) for r in recs] == [("a", "AC"), ("b", "GT")]


@pytest.mark.parametrize("header", [">\n", "> \n", ">\t\n"])
def test_empty_first_token_is_error(tmp_path: Path, header: str) -> None:
    with pytest.raises(FastaFormatError) as exc:
        _read(tmp_path, ">ok\nA\n" + header + "A\n")
    assert exc.value.line_no == 3


def test_header_without_sequence_gives_empty_sequence(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a\n>b\nAC\n>c\n")
    assert [(r.id, r.sequence) for r in recs] == [("a", ""), ("b", "AC"), ("c", "")]


def test_ids_colliding_after_mapping_are_both_returned(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">a|b\nA\n>a:b\nC\n")
    assert [r.id for r in recs] == ["a_b", "a_b"]
    assert [r.original_id for r in recs] == ["a|b", "a:b"]


def test_exact_duplicate_ids_not_rejected(tmp_path: Path) -> None:
    recs = _read(tmp_path, ">x\nA\n>x\nC\n")
    assert [r.id for r in recs] == ["x", "x"]


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("c circular=true", True),
        ("c Circular=True", True),
        ("c CIRCULAR=TRUE", True),
        ("c [topology=circular]", True),
        ("c [Topology=Circular]", True),
        ("c [TOPOLOGY=CIRCULAR] [gcode=11]", True),
        ("c circular=false", False),
        ("c [topology=linear]", False),
        ("c", False),
        ("c plasmid complete", False),
    ],
)
def test_topology_tag_detection(tmp_path: Path, header: str, expected: bool) -> None:
    (r,) = _read(tmp_path, f">{header}\nACGT\n")
    assert r.circular_tag is expected


def test_tag_in_sequence_part_is_ignored(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">c\nACGT\ncircular=true\n")
    assert r.circular_tag is False


def test_trailing_star_and_digits_kept_as_invalid_letters(tmp_path: Path) -> None:
    (r,) = _read(tmp_path, ">c\nACGT*\n12\n")
    assert r.sequence == "ACGT*12"
    assert invalid_letters(r.sequence) == {"*": 1, "1": 1, "2": 1}


@pytest.mark.parametrize("cut", [5, 20, -4])
def test_truncated_gzip_is_format_error(tmp_path: Path, cut: int) -> None:
    data = gzip.compress((EX1 * 50).encode("utf-8"))
    p = tmp_path / "trunc.fa.gz"
    p.write_bytes(data[:cut])
    with pytest.raises(FastaFormatError):
        read_genome_fasta(p)


# --- IUPAC_DNA, invalid_letters, sanitize_id, SAFE_ID -------------------------------------


def test_iupac_dna_alphabet_exact() -> None:
    assert isinstance(IUPAC_DNA, frozenset)
    assert set("ACGTNRYSWKMBDHV") == IUPAC_DNA


@pytest.mark.parametrize("ch", ["U", "-", ".", "*", "X", "a", "n"])
def test_not_in_iupac_dna(ch: str) -> None:
    assert ch not in IUPAC_DNA


def test_invalid_letters_empty_for_full_iupac() -> None:
    assert invalid_letters("ACGTNRYSWKMBDHV" * 3) == {}
    assert invalid_letters("") == {}


def test_invalid_letters_counts_u_dash_dot() -> None:
    # ADR-0005 item 1: U and gap characters are errors.
    assert invalid_letters("UUA-.-") == {"U": 2, "-": 2, ".": 1}


@pytest.mark.parametrize(
    ("raw", "safe"),
    [
        ("chr1", "chr1"),
        ("NZ_CP0|chr:1", "NZ_CP0_chr_1"),
        ("a.b-c_D9", "a.b-c_D9"),
        ("a/b\\c", "a_b_c"),
        ("x;y,z=1", "x_y_z_1"),
        ("contig#1(2)", "contig_1_2_"),
        ("é", "_"),
    ],
)
def test_sanitize_id(raw: str, safe: str) -> None:
    assert sanitize_id(raw) == safe


def test_safe_id_pattern() -> None:
    for ch in "aZ09._-":
        assert SAFE_ID.fullmatch(ch)
    for ch in "|: /\\#é":
        assert not SAFE_ID.fullmatch(ch)


# --- open_text ----------------------------------------------------------------------------


def test_open_text_plain_and_gzip_by_magic(tmp_path: Path) -> None:
    text = "héllo\nworld\n"
    plain = tmp_path / "a.gz"  # misleading suffix, plain content
    plain.write_bytes(text.encode("utf-8"))
    gz = tmp_path / "b.txt"  # misleading suffix, gzip content
    gz.write_bytes(gzip.compress(text.encode("utf-8")))
    for p in (plain, gz):
        with open_text(p) as fh:
            assert isinstance(fh, io.TextIOBase)
            assert fh.read() == text


# --- write_fasta --------------------------------------------------------------------------


def test_write_fasta_exact_multiple_of_width(tmp_path: Path) -> None:
    p = tmp_path / "o.fa"
    write_fasta(p, [("a", "C" * 160), ("b", "G")])
    assert p.read_text() == ">a\n" + "C" * 80 + "\n" + "C" * 80 + "\n>b\nG\n"


def test_write_fasta_custom_width_and_lf_only(tmp_path: Path) -> None:
    p = tmp_path / "o.fa"
    write_fasta(p, [("a", "ACGTA")], width=2)
    data = p.read_bytes()
    assert data == b">a\nAC\nGT\nA\n"
    assert b"\r" not in data


def test_write_fasta_uncompressed_even_with_gz_suffix(tmp_path: Path) -> None:
    p = tmp_path / "o.fa.gz"
    write_fasta(p, [("a", "ACGT")])
    assert p.read_bytes() == b">a\nACGT\n"


def test_write_fasta_accepts_generator(tmp_path: Path) -> None:
    p = tmp_path / "o.fa"
    write_fasta(p, ((f"c{i}", "A" * i) for i in range(1, 4)))
    assert p.read_text() == ">c1\nA\n>c2\nAA\n>c3\nAAA\n"


def test_round_trip_write_then_read(tmp_path: Path) -> None:
    pairs = [("chr1", "ACGT" * 50 + "N"), ("p.1", "RYSWKMBDHV"), ("x-2", "G" * 80)]
    p = tmp_path / "rt.fa"
    write_fasta(p, pairs)
    recs = read_genome_fasta(p)
    assert [(r.id, r.sequence) for r in recs] == pairs
    assert all(r.original_id == r.id and r.header == r.id for r in recs)
    assert not any(r.circular_tag for r in recs)


# --- Sim genome ---------------------------------------------------------------------------


def test_sim_genome_reads_as_one_circular_record(tmp_path: Path) -> None:
    cfg = gippasim.SimConfig(genome_length=20_000)
    gippasim.simulate(cfg, seed=3, outdir=tmp_path)
    recs = read_genome_fasta(tmp_path / "genome.fasta")
    assert len(recs) == 1
    (r,) = recs
    assert r.id == cfg.contig_id
    assert r.circular_tag is True
    assert len(r.sequence) == cfg.genome_length
    assert invalid_letters(r.sequence) == {}
