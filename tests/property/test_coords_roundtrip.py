"""Hypothesis properties for coordinates (issue M1-01; GIAP_SPEC.md §6.5, §11.2)."""

from __future__ import annotations

from typing import get_args

from giap.model.coords import Strand
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.model import (
    Interval,
    Replicon,
    check_interval,
    from_one_based_closed,
    split_at_origin,
    to_one_based_closed,
    wrap_position,
)

FAST = settings(max_examples=200, deadline=None)
STRANDS = st.sampled_from(["+", "-", "."])


@st.composite
def replicon_and_interval(draw: st.DrawFn) -> tuple[Replicon, Interval]:
    """A valid interval on a random replicon (both topologies, all strands)."""
    length = draw(st.integers(min_value=1, max_value=10_000))
    topology = draw(st.sampled_from(["circular", "linear"]))
    rep = Replicon("c1", length, topology)
    start = draw(st.integers(min_value=0, max_value=length - 1))
    max_end = start + length if topology == "circular" else length
    end = draw(st.integers(min_value=start + 1, max_value=max_end))
    return rep, Interval("c1", start, end, draw(STRANDS))


@st.composite
def circular_and_interval(draw: st.DrawFn) -> tuple[Replicon, Interval]:
    length = draw(st.integers(min_value=1, max_value=10_000))
    rep = Replicon("c1", length, "circular")
    start = draw(st.integers(min_value=0, max_value=length - 1))
    end = draw(st.integers(min_value=start + 1, max_value=start + length))
    return rep, Interval("c1", start, end, draw(STRANDS))


def test_strand_literal_values() -> None:
    assert set(get_args(Strand)) == {"+", "-", "."}


@FAST
@given(replicon_and_interval())
def test_generated_intervals_are_valid(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    check_interval(iv, rep)


@FAST
@given(replicon_and_interval())
def test_one_based_round_trip(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    start1, end1 = to_one_based_closed(iv)
    assert (start1, end1) == (iv.start + 1, iv.end)
    assert from_one_based_closed(iv.contig, start1, end1, iv.strand, rep) == iv


@FAST
@given(replicon_and_interval())
def test_from_one_based_then_back_round_trip(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    start1, end1 = iv.start + 1, iv.end
    out = from_one_based_closed(iv.contig, start1, end1, iv.strand, rep)
    assert to_one_based_closed(out) == (start1, end1)


@FAST
@given(circular_and_interval())
def test_wrapped_one_based_input_unwraps(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    assume(iv.crosses_origin(rep) and iv.end - rep.length < iv.start + 1)  # start1 > end1
    wrapped_end1 = iv.end - rep.length
    assert from_one_based_closed(iv.contig, iv.start + 1, wrapped_end1, iv.strand, rep) == iv


@FAST
@given(replicon_and_interval())
def test_split_preserves_length_and_bounds(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    parts = split_at_origin(iv, rep)
    assert sum(p.length_bp for p in parts) == iv.length_bp
    for p in parts:
        assert 0 <= p.start < p.end <= rep.length
        assert p.contig == iv.contig
        assert p.strand == iv.strand
        check_interval(p, rep)


@FAST
@given(replicon_and_interval())
def test_split_count_matches_crosses_origin(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    parts = split_at_origin(iv, rep)
    assert iv.crosses_origin(rep) == (iv.end > rep.length)
    if rep.topology == "linear":
        assert not iv.crosses_origin(rep)
    if iv.crosses_origin(rep):
        assert parts == (
            Interval(iv.contig, iv.start, rep.length, iv.strand),
            Interval(iv.contig, 0, iv.end - rep.length, iv.strand),
        )
    else:
        assert parts == (iv,)


@FAST
@given(circular_and_interval())
def test_split_then_merge_recovers_interval(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    parts = split_at_origin(iv, rep)
    merged_end = parts[0].end + (parts[1].end if len(parts) == 2 else 0)
    assert Interval(iv.contig, parts[0].start, merged_end, iv.strand) == iv


@FAST
@given(circular_and_interval())
def test_split_covers_same_positions(case: tuple[Replicon, Interval]) -> None:
    rep, iv = case
    expected = [wrap_position(p, rep) for p in range(iv.start, iv.end)]
    got = [p for part in split_at_origin(iv, rep) for p in range(part.start, part.end)]
    assert got == expected
    assert len(set(got)) == len(got)


@FAST
@given(st.integers(min_value=1, max_value=10_000), st.data())
def test_wrap_position_circular(length: int, data: st.DataObject) -> None:
    rep = Replicon("c1", length, "circular")
    pos = data.draw(st.integers(min_value=0, max_value=2 * length - 1))
    w = wrap_position(pos, rep)
    assert 0 <= w < length
    assert w == pos % length


@FAST
@given(st.integers(min_value=1, max_value=10_000), st.data())
def test_wrap_position_linear_identity_inside(length: int, data: st.DataObject) -> None:
    rep = Replicon("c1", length, "linear")
    pos = data.draw(st.integers(min_value=0, max_value=length - 1))
    assert wrap_position(pos, rep) == pos
