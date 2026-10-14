"""Hypothesis properties for TSV, BED and proBed writers (issue M1-04; GIAP_SPEC.md §6.5, §6.6)."""

from __future__ import annotations

import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.io.bed import BedRecord, bed_lines
from giap.io.probed import PROBED_EXTRA, ProBedRecord, write_probed
from giap.io.tsv import format_value, read_tsv, write_tsv
from giap.model import Interval, Replicon

FAST = settings(max_examples=150, deadline=None)

TEXT = st.text(
    alphabet=st.characters(blacklist_categories=("Cs", "Cc"), blacklist_characters="\t\n\r"),
    max_size=12,
)
CELL = st.one_of(
    st.none(),
    st.booleans(),
    st.integers(min_value=-(2**63), max_value=2**63),
    st.floats(allow_nan=True, allow_infinity=True),
    TEXT,
    # ADR-0006: list items may not be empty or contain the ";" separator
    st.lists(TEXT.filter(lambda s: s != "" and ";" not in s), max_size=4),
)
COLUMN = st.from_regex(r"[a-z][a-z0-9_]{0,10}", fullmatch=True)


@st.composite
def table(draw: st.DrawFn) -> tuple[list[str], list[list[object]]]:
    cols = draw(st.lists(COLUMN, min_size=1, max_size=5, unique=True))
    rows = draw(st.lists(st.lists(CELL, min_size=len(cols), max_size=len(cols)), max_size=8))
    return cols, rows


@FAST
@given(table(), st.booleans())
def test_tsv_write_read_round_trip(tab: tuple[list[str], list[list[object]]], gz: bool) -> None:
    cols, rows = tab
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / ("t.tsv.gz" if gz else "t.tsv")
        write_tsv(p, cols, rows)
        header, got = read_tsv(p)
    expected = [[None if (s := format_value(c)) in ("NA", "") else s for c in row] for row in rows]
    assert header == cols
    assert got == expected


@FAST
@given(table())
def test_tsv_gzip_reproducible(tab: tuple[list[str], list[list[object]]]) -> None:
    cols, rows = tab
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / "a.tsv.gz", Path(d) / "b.tsv.gz"
        write_tsv(a, cols, rows)
        write_tsv(b, cols, rows)
        assert a.read_bytes() == b.read_bytes()


@st.composite
def circular_interval(draw: st.DrawFn) -> tuple[Replicon, Interval]:
    length = draw(st.integers(min_value=1, max_value=10_000))
    rep = Replicon("c1", length, "circular")
    start = draw(st.integers(min_value=0, max_value=length - 1))
    end = draw(st.integers(min_value=start + 1, max_value=start + length))
    return rep, Interval("c1", start, end, draw(st.sampled_from(["+", "-", "."])))


@FAST
@given(circular_interval())
def test_bed_parts_cover_interval(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    lines = bed_lines(BedRecord(iv, "g"), rep)
    assert len(lines) == (2 if iv.end > rep.length else 1)
    total = 0
    for line in lines:
        f = line.split("\t")
        assert len(f) == 6
        s, e = int(f[1]), int(f[2])
        assert 0 <= s < e <= rep.length
        assert f[5] == iv.strand
        total += e - s
    assert total == iv.length_bp
    if len(lines) == 2:
        assert {x.split("\t")[3] for x in lines} == {"g_part1", "g_part2"}


EXTRA_VALUE = st.one_of(
    st.none(),
    st.integers(min_value=0, max_value=10),
    st.floats(min_value=0, max_value=1e6, allow_nan=False),
    st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=20),
)


@settings(max_examples=20, deadline=None)
@given(
    st.lists(
        st.tuples(
            circular_interval(),
            st.dictionaries(st.sampled_from(PROBED_EXTRA), EXTRA_VALUE, max_size=6),
        ),
        min_size=50,
        max_size=50,
    )
)
def test_probed_lines_have_25_fields(cases: list) -> None:
    rep = Replicon("c1", 10_000, "circular")
    recs = []
    for i, ((_r, iv), extra) in enumerate(cases):
        recs.append(ProBedRecord(BedRecord(iv, f"p{i}"), extra))
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x.probed"
        write_probed(p, recs, {"c1": rep}, "v1")
        lines = [x for x in p.read_text().splitlines() if not x.startswith("#")]
    assert len(lines) >= 50
    for line in lines:
        f = line.split("\t")
        assert len(f) == 25
        assert 0 <= int(f[1]) < int(f[2]) <= rep.length
