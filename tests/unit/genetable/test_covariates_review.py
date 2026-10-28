"""M4-10 review gaps (docs/specs/M4-10-protein-covariates.md, "Lead decisions after code review").

Item 1: TM boundary (window mean exactly 1.6 counts), internal X kept in place (0 in TM windows,
inside peptides), ambiguous codons -> X, lower-case ``nt`` upper-cased, ``missed < 0`` raises.
Item 2: an internal ``*`` (after stripping leading/trailing ones) raises ValueError.
"""

from __future__ import annotations

import dataclasses

import pytest


def _cov():
    from gippa.genetable import covariates as cov

    return cov


# ------------------------------------------------------------ item 1: TM boundary


def test_tm_window_mean_exactly_threshold_counts() -> None:
    # KD: L 3.8 + W -0.9 + Y -1.3 + 16 * A 1.8 = 30.4 = 1.6 * 19 exactly (mean >= 1.6 counts).
    assert _cov().tm_segments("LWY" + "A" * 16) == 1


def test_tm_window_just_below_threshold_does_not_count() -> None:
    # One A -> G (KD -0.4): sum 28.2 < 30.4.
    assert _cov().tm_segments("LWY" + "A" * 15 + "G") == 0


# ------------------------------------------------------------ item 1: internal X


def test_internal_x_scores_zero_in_tm_window() -> None:
    # 19-residue window with X in place: 18 * 3.8 / 19 = 3.6 >= 1.6 -> one segment. Removing X
    # first would leave 18 residues (no window at all).
    c = _cov().covariates("L" * 9 + "X" + "L" * 9, "trypsin", (7, 30))
    assert c.tm_segments == 1
    assert c.length_aa == 19


def test_internal_x_zero_score_can_break_boundary() -> None:
    # X replaces an A (1.8 -> 0): sum 28.6 < 30.4, so the boundary window no longer counts.
    c = _cov().covariates("LWYX" + "A" * 15, "trypsin", (7, 30))
    assert c.tm_segments == 0


def test_internal_x_stays_inside_peptides() -> None:
    # Peptides as written: AAAXAAK (7) and AAAAAAAK (8) -> 2 observable in [7, 30].
    # Dropping X first would give AAAAAK (6) -> only 1.
    c = _cov().covariates("AAAXAAKAAAAAAAK", "trypsin", (7, 30))
    assert c.n_observable_peptides == 2
    assert c.length_aa == 15


# ------------------------------------------------------------ item 1: translation


def test_ambiguous_codon_translates_to_x() -> None:
    assert _cov().protein_sequence(None, "ATGNNNGCT", "pseudogene") == "MXA"


def test_lower_case_nt_is_upper_cased() -> None:
    cov = _cov()
    assert cov.protein_sequence(None, "atggctgct", "pseudogene") == "MAA"
    # Lower-case stop codon still ends translation.
    assert cov.protein_sequence(None, "gtggcttaagct", "pseudogene") == "MA"
    assert cov.protein_sequence(None, "aTgGcTtAa", "pseudogene") == "MA"


# ------------------------------------------------------------ item 1: missed < 0


@pytest.mark.parametrize("missed", [-1, -5])
def test_negative_missed_raises(missed: int) -> None:
    with pytest.raises(ValueError):
        _cov().digest("MAGLLKAAGPEPTIDER", "trypsin", missed)


# ------------------------------------------------------------ item 2: internal stop


@pytest.mark.parametrize("aa", ["MK*AAAAAAAAK", "MK*AAAAAAAAK*", "*MK*AAAAAAAAK", "M**K"])
def test_internal_stop_raises(aa: str) -> None:
    with pytest.raises(ValueError):
        _cov().covariates(aa, "trypsin", (7, 30))


@pytest.mark.parametrize("aa", ["MKAAAAAAAAK*", "*MKAAAAAAAAK", "**MKAAAAAAAAK**"])
def test_leading_trailing_stops_still_stripped(aa: str) -> None:
    cov = _cov()
    got = dataclasses.astuple(cov.covariates(aa, "trypsin", (7, 30)))
    want = dataclasses.astuple(cov.covariates("MKAAAAAAAAK", "trypsin", (7, 30)))
    assert got == want
    assert got[0] == 11
