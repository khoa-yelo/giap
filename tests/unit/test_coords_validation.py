"""Review findings pinned for gippa.model.coords input validation (M1-01)."""

from __future__ import annotations

import pytest

from gippa.model.coords import (
    Interval,
    Replicon,
    from_one_based_closed,
    interval_sort_key,
    wrap_position,
)

CIRC = Replicon("c1", 1000, "circular", "chromosome")
LIN = Replicon("l1", 1000, "linear", "unknown")


# 1. from_one_based_closed rejects non-positive 1-based coordinates.
@pytest.mark.parametrize(
    ("contig", "start1", "end1", "rep"),
    [
        ("c1", 5, 0, CIRC),
        ("c1", 5, -3, CIRC),
        ("c1", 0, 10, CIRC),
        ("c1", -2, 10, CIRC),
        ("l1", 0, 10, LIN),
        ("l1", 5, 0, LIN),
    ],
)
def test_from_one_based_closed_rejects_nonpositive(contig, start1, end1, rep):
    with pytest.raises(ValueError):
        from_one_based_closed(contig, start1, end1, "+", rep)


# 2. wrap_position rejects negative positions on both topologies.
@pytest.mark.parametrize("rep", [CIRC, LIN], ids=["circular", "linear"])
def test_wrap_position_rejects_negative(rep):
    with pytest.raises(ValueError):
        wrap_position(-1, rep)


# 3. Interval validates field values at construction.
@pytest.mark.parametrize("strand", ["x", "", "++", "1", None])
def test_interval_rejects_bad_strand(strand):
    with pytest.raises(ValueError):
        Interval("c1", 0, 10, strand)


@pytest.mark.parametrize(
    ("start", "end"),
    [(0.5, 10), (0, 10.0), ("3", 10), (0, "10"), (True, 10), (0, True), (False, 5)],
)
def test_interval_rejects_non_int_coords(start, end):
    with pytest.raises((TypeError, ValueError)):
        Interval("c1", start, end, "+")


# 4. interval_sort_key names an unknown contig.
def test_interval_sort_key_unknown_contig_names_it():
    iv = Interval("mystery_contig", 0, 10, "+")
    with pytest.raises(ValueError, match="mystery_contig"):
        interval_sort_key(iv, {"c1": 0, "c2": 1})
