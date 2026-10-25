"""``readthrough.consumed_sites`` (M7-16b rule 15; spec examples 11-16 and the "Edge cases (16b,
16c labels)"; docs/specs/M7-16-stop-codon-readthrough.md).

A ``peptide_hits`` row consumes the site ``S`` when it is a target row without ``blocks``, on
``S``'s contig and strand, and some ``x`` in ``{S.codon_start, S.codon_start + n}`` (the second only
on a circular replicon of length ``n``) has ``start <= x``, ``x + 3 <= end`` and
``(x - start) mod 3 == 0`` on the row's unwrapped 0-based ``[start, end)``. The function returns,
per row, the consumed sites sorted by ``codon_start`` (``()`` for the rest).

Interface pinned here (spec rule 15)::

    consumed_sites(rows, sites, lengths, circular) -> list[tuple[StopSite, ...]]

Toy contig ``c1`` (linear, 126 nt) and gene ``G1`` = CDS ``+`` ``[0, 57)``: its terminal TGA is
``[54, 57)`` (cell p = 18). The TAA at ``[90, 93)`` is not an annotated stop.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd
import pytest

from giap.mapping.readthrough import StopSite

G1_SITE = StopSite("c1", "+", 54, ("G1",), "TGA", "terminal")
C2_SITE = StopSite("c2", "-", 120, ("G2",), "TAA", "terminal")  # example 8
LIN = {"c1": 126, "c2": 300}
LIN_TOPO = {"c1": False, "c2": False}


def _consumed(
    rows: pd.DataFrame,
    sites: Sequence[StopSite],
    lengths: Mapping[str, int] = LIN,
    circular: Mapping[str, bool] = LIN_TOPO,
) -> list[tuple[StopSite, ...]]:
    from giap.mapping.readthrough import consumed_sites

    return consumed_sites(rows, sites, lengths, circular)


def row(
    contig: str,
    start: int,
    end: int,
    strand: str = "+",
    *,
    search: str = "target",
    route: str = "tolerant",
    blocks: str | None = None,
    diff_types: str | None = "substitution",
    peptide: str = "PEPTIDEKAAA",
) -> dict[str, Any]:
    """One hit row (``giap.mapping.hits`` columns; ``frame_id``/``aa_*`` unused by the rule)."""
    return {
        "peptide": peptide,
        "query": peptide,
        "search": search,
        "decoy_type": "reverse" if search == "decoy" else None,
        "contig": contig,
        "frame_id": 0,
        "aa_start": 0,
        "aa_end": (end - start) // 3,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": 1,
        "route": route,
        "n_diff": 0 if route == "exact" else 1,
        "diff_types": None if route == "exact" else diff_types,
        "blocks": blocks,
    }


def frame(*rows: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(list(rows))


# ---------------------------------------------------------------- examples 11-16


def test_example_11_window_across_the_stop_consumes_g1() -> None:
    """``WQLNEFKWVPAGHTL`` at ``[33, 78) +``: (54 - 33) mod 3 = 0 and [54, 57) inside."""
    assert _consumed(frame(row("c1", 33, 78)), [G1_SITE]) == [(G1_SITE,)]


def test_example_12_window_ending_on_the_stop_consumes_g1() -> None:
    """``GRWQLNEFKW`` at ``[27, 57) +`` (Q3): the codon ends the window, ``x + 3 = end``."""
    assert _consumed(frame(row("c1", 27, 57)), [G1_SITE]) == [(G1_SITE,)]


def test_example_13_unannotated_stop_and_exact_rows_consume_nothing() -> None:
    """A window over the TAA at p = 30 (``[69, 111)``) and an exact interior row."""
    rows = frame(row("c1", 69, 111), row("c1", 0, 30, route="exact"))
    assert _consumed(rows, [G1_SITE]) == [(), ()]


def test_example_14_minus_strand_in_frame_and_out_of_frame() -> None:
    """Site ``(c2, -, 120)``: ``[111, 141) -`` consumes it ((120 - 111) mod 3 = 0); the same
    length one base over, ``[112, 142) -`` (another frame), does not."""
    rows = frame(row("c2", 111, 141, "-"), row("c2", 112, 142, "-"))
    assert _consumed(rows, [C2_SITE]) == [(C2_SITE,), ()]


def test_example_15_origin_copy_on_a_circular_replicon() -> None:
    """Circular length 300, ``+`` site at 0; the unwrapped row ``[282, 312) +`` holds
    ``x = 0 + 300`` with (300 - 282) mod 3 = 0. On a linear contig there is no second copy."""
    site = StopSite("o", "+", 0, ("O1",), "TAA", "terminal")
    rows = frame(row("o", 282, 312))
    assert _consumed(rows, [site], {"o": 300}, {"o": True}) == [(site,)]
    assert _consumed(rows, [site], {"o": 300}, {"o": False}) == [()]


def test_example_16_frameshift_rows_are_never_consumers() -> None:
    """A frameshift row (``blocks`` set) whose unwrapped span covers the codon in frame."""
    fs = row("c1", 33, 78, blocks="33-60;59-78", diff_types="frameshift")
    assert _consumed(frame(fs), [G1_SITE]) == [()]


# ---------------------------------------------------------------- edge cases


def test_decoy_rows_are_never_consumers() -> None:
    """Decoy rows carry the mirrored *target* interval; they never reach ``peptide_hits``."""
    assert _consumed(frame(row("c1", 33, 78, search="decoy")), [G1_SITE]) == [()]


def test_strand_and_contig_must_match() -> None:
    rows = frame(row("c1", 33, 78, "-"), row("c2", 33, 78, "+"))
    assert _consumed(rows, [G1_SITE]) == [(), ()]


@pytest.mark.parametrize(
    ("start", "end", "hit"),
    [
        (54, 84, True),  # the codon starts the window
        (24, 54, False),  # ends just before the codon
        (57, 87, False),  # starts just after it
        (55, 85, False),  # overlaps it out of frame
        (30, 56, False),  # cuts the codon (never a real row; containment needs x + 3 <= end)
        (21, 57, True),
    ],
)
def test_containment_and_phase_on_a_linear_contig(start: int, end: int, hit: bool) -> None:
    got = _consumed(frame(row("c1", start, end)), [G1_SITE])
    assert got == [(G1_SITE,) if hit else ()]


def test_origin_crossing_codon_needs_the_unwrapped_end() -> None:
    """Circular length 300, ``+`` stop ``[299, 302)`` (``codon_start`` 299): the row must reach
    302 unwrapped; a row ending at 300 does not hold the codon."""
    site = StopSite("o", "+", 299, ("O2",), "TGA", "terminal")
    rows = frame(row("o", 281, 302), row("o", 270, 300))
    assert _consumed(rows, [site], {"o": 300}, {"o": True}) == [(site,), ()]


def test_minus_strand_origin_copy() -> None:
    """Circular length 300, ``-`` stop ``[0, 3)``: the unwrapped row ``[291, 321) -`` holds the
    copy ``x = 300``; ``[0, 30) -`` holds ``x = 0``."""
    site = StopSite("o", "-", 0, ("O3",), "TAA", "terminal")
    rows = frame(row("o", 291, 321, "-"), row("o", 0, 30, "-"), row("o", 292, 322, "-"))
    assert _consumed(rows, [site], {"o": 300}, {"o": True}) == [(site,), (site,), ()]


def test_two_cells_in_one_row_sorted_by_codon_start() -> None:
    """Adjacent short genes (k = 2): both sites, ascending, whatever the input order."""
    a = StopSite("c1", "+", 15, ("A",), "TAA", "terminal")
    b = StopSite("c1", "+", 30, ("B",), "TGA", "terminal")
    assert _consumed(frame(row("c1", 0, 45)), [b, a]) == [(a, b)]


def test_empty_sites_and_empty_rows() -> None:
    rows = frame(row("c1", 33, 78), row("c1", 0, 30, route="exact"))
    assert _consumed(rows, []) == [(), ()]
    empty = frame(row("c1", 33, 78)).iloc[0:0]
    assert _consumed(empty, [G1_SITE]) == []


def test_one_entry_per_row_in_row_order() -> None:
    rows = frame(
        row("c1", 0, 30, route="exact"),
        row("c1", 33, 78),
        row("c2", 111, 141, "-"),
        row("c1", 27, 57),
    )
    got = _consumed(rows, [C2_SITE, G1_SITE])
    assert isinstance(got, list) and all(isinstance(t, tuple) for t in got)
    assert got == [(), (G1_SITE,), (C2_SITE,), (G1_SITE,)]


def test_recoded_sites_are_consumed_like_terminal_ones() -> None:
    """The rule is geometric: the site's ``kind`` decides the flag (rule 16/19), not the
    consumption."""
    rec = StopSite("c1", "+", 54, ("G1",), "TGA", "recoded")
    assert _consumed(frame(row("c1", 33, 78)), [rec]) == [(rec,)]
