"""Read-through re-alignment: the review-335 cases (M7-16b rules 10-11; code review of PR #335,
must-fix 4 and 5; statistical review item 4).

* ``_residues`` keys an accepted stop by the **wrapped** codon start: a hit stored unwrapped
  across the origin (codon ``>= n``) aligns with ``stop_sites`` keyed by ``codon % n``.
* Over a read-through cell (code 22) the DP uses only move 0, a substitution of cost 1 (rule 10):
  no swap, 1-to-2 or 2-to-1 block shape may involve the cell, even where the residues next to it
  would allow one. In particular, peptide ``...EW VPAG...`` against genome ``...EF * VPAG...``
  (one residue for the two cells ``F *``, a block shape over the cell) has no alignment.
* Error-rate classes (D4): a REF ``*`` site with a standard ALT, and U over TGA / O over TAG
  (``nt_changes`` 0), all fall in ``site_class`` ``(substitution, False)``.

World: residues ``MSDKPLTEAGRWQLNEF`` + ``*`` (cell p = 17, codon ``[51, 54)``) + ``VPAGHTLDMRK``
+ TAA on a 90-nt contig, G1 = CDS ``+`` ``[0, 54)``. The minus world is its reverse complement
(G1 ``-`` ``[36, 90)``, stop ``[36, 39)``).
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import back_translate, gene, model, rc

from gippa.config import GippaConfig
from gippa.model import Replicon
from gippa.variants.align import align_hit
from gippa.variants.error_rates import site_class

BODY, FLANK = "MSDKPLTEAGRWQLNEF", "VPAGHTLDMRK"
P = len(BODY)  # 17, the cell
N = 3 * (len(BODY) + 1 + len(FLANK) + 1)  # 90
WIN = (11, 25)  # aa window WQLNEF*VPAGHTL


def _nt(stop: str = "TGA", strand: str = "+", shift: int = 0) -> str:
    nt = back_translate(BODY + "*" + FLANK + "*", {P: stop, P + 1 + len(FLANK): "TAA"})
    assert len(nt) == N
    nt = nt if strand == "+" else rc(nt)
    return nt[shift:] + nt[:shift]


def _window(strand: str = "+", shift: int = 0) -> tuple[int, int]:
    s, e = 3 * WIN[0], 3 * WIN[1]
    if strand == "-":
        s, e = N - e, N - s
    a = (s - shift) % N
    return a, a + e - s


def _stop(strand: str = "+") -> int:
    return 3 * P if strand == "+" else N - 3 * P - 3


def _align(
    pep: str,
    n_diff: int,
    strand: str = "+",
    shift: int = 0,
    sites: frozenset[tuple[str, str, int]] | None = None,
) -> Any:
    topo = "circular" if shift else "linear"
    rep = Replicon("c1", N, topo)
    s, e = _window(strand, shift)
    if sites is None:
        sites = frozenset({("c1", strand, (_stop(strand) - shift) % N)})
    return align_hit(pep, "c1", s, e, strand, n_diff, None, _nt("TGA", strand, shift), rep, sites)


def _diffs(a: Any) -> list[tuple[str, tuple[int, ...], str, str]]:
    return [
        (o.kind, tuple(c % N for c in o.codons), o.ref, o.alt)
        for o in a.ops
        if o.kind not in ("match", "il")
    ]


# ------------------------------------------------------------ must-fix 4: codon past the origin


@pytest.mark.parametrize(("strand", "shift"), [("+", 40), ("-", 30)])
def test_stop_past_the_origin_is_keyed_by_the_wrapped_start(strand: str, shift: int) -> None:
    """``+``, shift 40: window ``[83, 125)`` unwrapped, the codon at 101 (wrapped 11). ``-``,
    shift 30: window ``[75, 117)``, the codon at 96 (wrapped 6). The site set holds only the
    wrapped start."""
    s, _e = _window(strand, shift)
    stop = (_stop(strand) - shift) % N
    assert s < N and s + (stop - s) % N >= N  # the codon lies past the origin, unwrapped
    a = _align("WQLNEFWVPAGHTL", 1, strand, shift)
    assert _diffs(a) == [("substitution", (stop,), "*", "W")]


def test_stop_past_the_origin_needs_the_site() -> None:
    """Control: without the site the ``*`` still raises."""
    with pytest.raises(ValueError, match=r"\*"):
        _align("WQLNEFWVPAGHTL", 1, "+", 40, sites=frozenset())


# ------------------------------------------------------------ must-fix 5: only move 0 at the cell


def test_cell_is_one_substitution_of_cost_one() -> None:
    a = _align("WQLNEFWVPAGHTL", 1)
    assert _diffs(a) == [("substitution", (51,), "*", "W")]
    (op,) = [o for o in a.ops if o.kind == "substitution"]
    assert op.pep[1] - op.pep[0] == 1 and len(op.codons) == 1


@pytest.mark.parametrize("strand", ["+", "-"])
def test_no_block_shape_over_the_cell(strand: str) -> None:
    """``WQLNEW`` + ``VPAGHTL`` (13 aa) against ``WQLNEF*VPAGHTL`` (14 cells): the only cost-1
    shape would read ``W`` for the two cells ``F *`` (a 1-to-2 block over the cell). Rule 10
    allows no such move, and no other alignment lies in the band."""
    with pytest.raises(ValueError, match="none within the band"):
        _align("WQLNEWVPAGHTL", 1, strand)


def test_swap_shape_into_the_cell_is_two_substitutions() -> None:
    """``...EWF VPAG...`` against ``...EF* VPAG...``: ``WF`` / ``F*`` looks like an adjacent pair,
    but the cell takes no swap; it is two substitutions (cost 2)."""
    a = _align("WQLNEWFVPAGHTL", 2)
    assert _diffs(a) == [("substitution", (48,), "F", "W"), ("substitution", (51,), "*", "F")]


def test_swap_right_after_the_cell_is_allowed() -> None:
    """``*VP`` read as ``WPV``: the cell is a substitution, the next pair a swap."""
    a = _align("WQLNEFWPVAGHTL", 2)
    assert _diffs(a) == [("substitution", (51,), "*", "W"), ("swap", (54, 57), "VP", "PV")]


# ---------------------------------------------------------------- statistical review item 4


def _run_sites(alt: str, stop: str) -> pd.DataFrame:
    from gippa.steps.variants import run_variants

    s, e = _window()
    pep = "WQLNEF" + alt + "VPAGHTL"
    hits = pd.DataFrame(
        [
            {
                "peptide": pep,
                "peptide_len": len(pep),
                "source": "denovo",
                "contig": "c1",
                "start": s,
                "end": e,
                "strand": "+",
                "n_diff": 1,
                "diff_types": "substitution",
                "placement_q": 0.0,
                "uniqueness": "unique",
            }
        ]
    )
    reps = {"c1": Replicon("c1", N, "linear")}
    gm = model(reps, [gene("G1", "c1", 0, 54, "+")])
    peptides = pd.DataFrame({"peptide": [pep], "search_db": "same"})
    cfg = GippaConfig.model_validate({"genome_relation": "same", "sample_name": "S1"})
    _pile, table, _records = run_variants(hits, None, peptides, gm, {"c1": _nt(stop)}, cfg)
    return table


@pytest.mark.parametrize(("alt", "stop"), [("W", "TGA"), ("Q", "TAA"), ("U", "TGA"), ("O", "TAG")])
def test_stop_sites_pool_in_the_substitution_class(alt: str, stop: str) -> None:
    table = _run_sites(alt, stop)
    assert len(table) == 1, table
    (r,) = table.to_dict("records")
    assert (r["ref_aa"], r["alt_aa"]) == ("*", alt)
    assert site_class(r["alt_type"], r["flags"]) == ("substitution", False)
