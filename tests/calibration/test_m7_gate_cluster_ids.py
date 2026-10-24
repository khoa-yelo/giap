"""The makedb ID of a clustered-null entry (M7-13b lead decision 7;
``helpers.m7_gate_novel_records.makedb_id`` / ``segment_span``). Unmarked: runs on every commit.

A segment carries ``nt = (nt0, t, n)`` in strand-local coordinates (``nt0`` the first base of the
stop-free stretch, ``t`` the first base of its closing stop codon, ``n`` the contig length) and
``stop``, the 1-based lowest plus-strand base of the stop codon (``_segments``). The ID is
``GIAPORF_<cls>_<contig>_<a>_<b>_<p|m>`` with ``a = nt0 + 1, b = t`` on ``+`` and
``a = n - t + 1, b = n - nt0`` on ``-``: the stretch without its stop codon, 1-based closed, in
plus-strand coordinates (giapsim's makedb convention; statistician probe check_ids.py).
"""

from __future__ import annotations

from typing import Any

from helpers import m7_gate_novel_records as N
from hypothesis import given, settings
from hypothesis import strategies as st


def seg(strand: str, nt0: int, t: int, n: int, contig: str = "c1") -> dict[str, Any]:
    stop0 = t if strand == "+" else n - (t + 3)  # as in _segments
    return {"contig": contig, "strand": strand, "stop": stop0 + 1, "nt": (nt0, t, n)}


def test_plus_strand_id_example() -> None:
    s = seg("+", 10, 40, 1000)
    assert N.makedb_id("antisense", s) == "GIAPORF_antisense_c1_11_40_p"
    assert N.segment_span(s) == (10, 43)
    assert s["stop"] == 41


def test_minus_strand_id_example() -> None:
    s = seg("-", 10, 40, 1000)
    assert N.makedb_id("antisense", s) == "GIAPORF_antisense_c1_961_990_m"
    assert N.segment_span(s) == (957, 990)
    assert s["stop"] == 958


def test_class_and_contig_enter_the_id() -> None:
    assert N.makedb_id("intergenic", seg("+", 0, 30, 99, "chr_2")) == (
        "GIAPORF_intergenic_chr_2_1_30_p"
    )


@st.composite
def segments(draw: Any) -> dict[str, Any]:
    n = draw(st.integers(30, 10**6))
    codons = draw(st.integers(1, (n - 3) // 3))
    nt0 = draw(st.integers(0, n - 3 - 3 * codons))
    return seg(draw(st.sampled_from("+-")), nt0, nt0 + 3 * codons, n)


def _ab(id_: str) -> tuple[int, int, str]:
    *_, a, b, st_ = id_.split("_")
    return int(a), int(b), st_


@settings(max_examples=300, deadline=None)
@given(segments())
def test_id_is_the_stretch_before_the_stop_codon(s: dict[str, Any]) -> None:
    nt0, t, n = s["nt"]
    a, b, tag = _ab(N.makedb_id("antisense", s))
    lo, hi = N.segment_span(s)
    assert 1 <= a <= b <= n and b - a + 1 == t - nt0 and (b - a + 1) % 3 == 0
    assert hi - lo == t - nt0 + 3 and 0 <= lo < hi <= n
    if s["strand"] == "+":
        assert tag == "p" and (lo, hi) == (a - 1, b + 3) and s["stop"] == b + 1
    else:
        assert tag == "m" and (lo, hi) == (a - 4, b) and s["stop"] == a - 3


@settings(max_examples=200, deadline=None)
@given(segments())
def test_minus_strand_id_mirrors_the_plus_strand_coordinates(s: dict[str, Any]) -> None:
    nt0, t, n = s["nt"]
    a, b, _ = _ab(N.makedb_id("antisense", s))
    a2, b2, _ = _ab(N.makedb_id("antisense", seg("-" if s["strand"] == "+" else "+", nt0, t, n)))
    assert (a2, b2) == (n - b + 1, n - a + 1)
