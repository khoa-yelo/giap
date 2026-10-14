"""Hypothesis properties for the genome FASTA writer/reader (issue M1-02; GIAP_SPEC.md §5.3)."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.io.fasta import SAFE_ID, invalid_letters, read_genome_fasta, sanitize_id, write_fasta

FAST = settings(max_examples=150, deadline=None)
SAFE_CHARS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-"
IDS = st.text(alphabet=SAFE_CHARS, min_size=1, max_size=20)
SEQS = st.text(alphabet="ACGTNRYSWKMBDHV", min_size=0, max_size=400)
RECORDS = st.lists(st.tuples(IDS, SEQS), min_size=0, max_size=8)
WIDTHS = st.integers(min_value=1, max_value=120)


@FAST
@given(records=RECORDS, width=WIDTHS)
def test_write_read_round_trip(records: list[tuple[str, str]], width: int) -> None:
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "g.fa"
        write_fasta(p, records, width=width)
        got = read_genome_fasta(p)
        assert [(r.id, r.sequence) for r in got] == records
        assert all(r.original_id == r.id for r in got)
        assert all(invalid_letters(r.sequence) == {} for r in got)
        lines = p.read_text().split("\n")
        assert lines[-1] == ""
        assert all(len(line) <= width for line in lines if not line.startswith(">"))


@FAST
@given(records=RECORDS)
def test_gzip_round_trip_matches_plain(records: list[tuple[str, str]]) -> None:
    with tempfile.TemporaryDirectory() as d:
        plain = Path(d) / "g.fa"
        write_fasta(plain, records)
        gz = Path(d) / "g2.fa"
        gz.write_bytes(gzip.compress(plain.read_bytes()))
        assert read_genome_fasta(gz) == read_genome_fasta(plain)


@FAST
@given(seq=SEQS)
def test_lower_case_input_is_upper_cased(seq: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "g.fa"
        p.write_text(">c\n" + seq.lower() + "\n")
        (r,) = read_genome_fasta(p)
        assert r.sequence == seq


@FAST
@given(raw=st.text(min_size=1, max_size=30).filter(lambda s: not any(c.isspace() for c in s)))
def test_sanitize_id_is_safe_length_preserving_and_idempotent(raw: str) -> None:
    safe = sanitize_id(raw)
    assert len(safe) == len(raw)
    assert all(SAFE_ID.fullmatch(c) for c in safe)
    assert sanitize_id(safe) == safe
    assert all(s == r for s, r in zip(safe, raw, strict=True) if SAFE_ID.fullmatch(r))
