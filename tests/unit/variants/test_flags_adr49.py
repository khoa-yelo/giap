"""ADR-0049 flag changes (Khoa's decision, 2026-10-23; local/handoff/adr0049_decision.md).

Supersedes M7-05 lead decision 1 ("swaps carry no ``mass_equivalent``") and adds the A6a 2->1
block reading of substitutions. GIAP_SPEC.md §7.5.3: ``mass_equivalent`` when an equal-mass
reading explains the ALT, blocks included.

Interfaces pinned here (test-author choices, reported to the lead):

* ``giap.variants.flags.flags_for(site, starts, cfg, alt_sources)`` (signature unchanged) flags
  every ``alt_type = swap`` row ``mass_equivalent`` (decision 2, A6b), whatever the tolerance.
* ``giap.variants.flags.left_block_read(peptide, q, ref, k, tol_da) -> bool`` (new, pure): the
  ALT peptide ``peptide`` reads its residue ``q`` at the site, which is residue ``k`` of the
  reference residues ``ref`` (transcription order). True iff the peptide also aligns as the 2->1
  block ``ref[k-1:k+1] -> peptide[q]`` with its prefix shifted left by one, i.e. ``k >= 1``,
  ``mass_equivalent(ref[k-1:k+1], peptide[q], tol_da)`` and either ``q == 0`` or
  ``k - 1 - q >= 0`` and ``peptide[:q] == ref[k-1-q:k-1]`` (I and L folded, as in the mapper).
  Only the left form: never the right-hand ``ref[k:k+2]`` form, never a prefix that is not
  shifted. The caller (the variants step) supplies ``q``, ``ref`` and ``k``; the step-level
  behaviour is pinned in tests/unit/steps/test_variants_step_adr49.py.
* ``giap.variants.flags.identifiable(ref_aa, alt_aa, tol_da) -> bool`` (new, pure): the G1
  identifiable-truth rule (decision 4) for the M7-13a gate helpers, exactly
  ``not mass_equivalent(ref_aa, alt_aa, tol_da)``, so the gate and the call flag cannot drift.

Masses come from GIAP's own table (``giap.mapping.masses.RESIDUE_MASS``), as the decision asks.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any

import pytest

from giap.config.model import GiapConfig
from giap.mapping.masses import RESIDUE_MASS

ME = "mass_equivalent"
TOL = GiapConfig().variants.mass_equivalent_tol_da  # 0.05 Da, the configured default


def fl() -> ModuleType:
    """The module under test (imported per test, so each test fails on its own)."""
    return importlib.import_module("giap.variants.flags")


def cfg(tol: float | None = None, relation: str = "related") -> GiapConfig:
    over: dict[str, Any] = {"genome_relation": relation}
    if tol is not None:
        over["variants"] = {"mass_equivalent_tol_da": tol}
    return GiapConfig.model_validate(over)


def site(ref_codon: str, ref_aa: str, alt_aa: str, alt_type: str) -> dict[str, Any]:
    """One ``giap.variants.pileup.sites`` row (0-based ``codon_start``)."""
    return {
        "variant_id": "VAR_1", "contig": "c1", "codon_start": 99, "strand": "+",
        "gene": "G1", "aa_pos": 2, "ref_codon": ref_codon, "ref_aa": ref_aa, "alt_aa": alt_aa,
        "alt_type": alt_type, "call": "ALT", "confidence": "high", "n_ref": 0, "n_alt": 2,
        "source": "denovo",
    }  # fmt: skip


def test_the_package_mass_table_reads_ag_as_q() -> None:
    """The A6a examples rest on Q = A + G (0.00001 Da with the B.5 masses) and on G vs Q being
    no equal-mass reading on its own (so the substitution G>Q is not flagged by rule 1). A>Q is
    (A + Carbamidomethyl), so the examples always put the site on G."""
    assert abs(RESIDUE_MASS["A"] + RESIDUE_MASS["G"] - RESIDUE_MASS["Q"]) < 1e-4
    assert fl().mass_equivalent("AG", "Q", TOL)
    assert fl().mass_equivalent("GA", "Q", TOL)  # composition: the right-hand pair weighs as Q too
    assert not fl().mass_equivalent("G", "Q", TOL)


# --------------------------------------------------------------------------- swaps (A6b)

SWAPS = [
    # (ref_codon, ref_aa, alt_aa): the two codons of the pair, the first residue, the swapped pair
    ("AAAGCT", "K", "AK"),
    ("GGTACT", "G", "TG"),
    ("TCTGTT", "S", "VS"),
    ("TGGCAT", "W", "HW"),
    ("GAACGT", "E", "RE"),
]


@pytest.mark.parametrize(("ref_codon", "ref_aa", "alt"), SWAPS)
def test_every_swap_row_is_mass_equivalent(ref_codon: str, ref_aa: str, alt: str) -> None:
    """A6b: ``flags_for`` flags every ``swap`` row (75 of 75 in the statistician's data)."""
    got = fl().flags_for(site(ref_codon, ref_aa, alt, "swap"), frozenset(), cfg(), set())
    assert ME in got, got
    assert got[0] == ME, f"§9.4 order puts mass_equivalent first: {got}"


@pytest.mark.parametrize("tol", [0.0, 0.01, 0.05])
@pytest.mark.parametrize("relation", ["related", "same"])
def test_swap_flag_holds_at_any_tolerance_and_relation(tol: float, relation: str) -> None:
    got = fl().flags_for(site("TCTGTT", "S", "VS", "swap"), frozenset(), cfg(tol, relation), set())
    assert ME in got, got
    assert "snv" not in got and "multi_nt" not in got, "swaps still carry no nt class"


def test_swapped_pair_is_mass_equivalent_to_the_reference_pair() -> None:
    for ref_codon, _, alt in SWAPS:
        pair = alt[1] + alt[0]
        assert fl().mass_equivalent(pair, alt, 0.0), (ref_codon, pair, alt)


# --------------------------------------------------------------------------- A6a: left 2->1 block

# The decision's examples on reference residues (transcription order). ``k`` is the site.
#   N-terminal: ...WAGSVATEPIK..., ALT peptide QSVATEPIK reads G>Q at q = 0; AG -> Q.
#   internal:   ...KAAGSVATAEPIK..., ALT peptide AQSVATAEPIK reads G>Q at q = 1; A + AG -> A + Q.
REF_NT = "DWAGSVATEPIKR"
K_NT = 3  # REF_NT[3] == "G"
REF_IN = "DKAAGSVATAEPIKR"
K_IN = 4  # REF_IN[4] == "G"


def test_example_references_are_built_as_described() -> None:
    assert REF_NT[K_NT - 1 : K_NT + 1] == "AG" and REF_IN[K_IN - 2 : K_IN + 1] == "AAG"


def test_n_terminal_ag_to_q_reads_as_a_block() -> None:
    """Decision 2: the N-terminal AG -> Q read at ALT-peptide index 0 is a block."""
    assert fl().left_block_read("QSVATEPIK", 0, REF_NT, K_NT, TOL) is True


def test_internal_ag_to_q_reads_as_a_block_with_the_prefix_shifted() -> None:
    """``AQSVATAEPIK`` over ``...AAG...``: prefix ``A`` equals ``ref[k-2]``, so the peptide also
    aligns as ``A`` + (``AG`` -> ``Q``) + suffix."""
    assert fl().left_block_read("AQSVATAEPIK", 1, REF_IN, K_IN, TOL) is True


def test_longer_shifted_prefix_reads_as_a_block() -> None:
    # ref ...AAAG S..., peptide AAQS... (q = 2): it reads G>Q (prefix AA == ref[k-2:k]) and the
    # shifted prefix AA == ref[k-3:k-1] holds too
    ref = "DAAAGSWK"
    k = 4
    assert ref[k] == "G"
    assert fl().left_block_read("AAQSWK", 2, ref, k, TOL) is True


def test_prefix_that_is_not_shifted_does_not_read_as_a_block() -> None:
    """Looser form: ``ref[k-1:k+1]`` is AG -> Q and the peptide reads ``A`` before the site, but
    ``ref[k-2]`` is ``K``, not ``A``: the block alignment would need a second difference."""
    ref = "DKAGSVATAEPIK"
    k = 3
    assert ref[k - 1 : k + 1] == "AG"
    assert fl().left_block_read("AQSVATAEPIK", 1, ref, k, TOL) is False


def test_partly_matching_prefix_does_not_read_as_a_block() -> None:
    # q = 2: the peptide reads G>Q (prefix AA == ref[k-2:k]) but ref[k-3:k-1] is "SA", not "AA"
    ref = "DSAAGWK"
    k = 4
    assert fl().left_block_read("AAQWK", 2, ref, k, TOL) is False


def test_right_hand_block_does_not_read_as_a_block() -> None:
    """The right-hand form ``ref[k:k+2] -> alt`` (GA -> Q, suffix shifted) never flags, even when
    it is mass-equivalent and the suffix shift holds: the C-terminal ALT of ``DFHYWQ`` over
    ``...DFHYWGAM...`` and the internal ALT of ``RHYWQAA`` over ``...RHYWGAAAM...``."""
    ref = "PDFHYWGAMK"
    k = 6
    assert ref[k : k + 2] == "GA"
    assert fl().left_block_read("DFHYWQ", 5, ref, k, TOL) is False
    ref2 = "PRHYWGAAAMK"
    k2 = 5
    assert ref2[k2 : k2 + 2] == "GA" and ref2[k2 + 1 : k2 + 3] == ref2[k2 + 2 : k2 + 4]
    assert fl().left_block_read("RHYWQAA", 4, ref2, k2, TOL) is False


def test_left_block_needs_a_residue_before_the_site() -> None:
    assert fl().left_block_read("QSVATEPIK", 0, "GSVATEPIK", 0, TOL) is False


def test_left_block_needs_the_shifted_prefix_inside_the_reference() -> None:
    # q = 1 but only one residue before the site: ref[k-2] does not exist
    assert fl().left_block_read("AQSV", 1, "AGSV", 1, TOL) is False


def test_left_block_needs_the_pair_to_be_mass_equivalent() -> None:
    # SG -> Q: 144.05 vs 128.06, no modification delta closes it
    assert fl().left_block_read("QSVATEPIK", 0, "DWSGSVATEPIK", 3, TOL) is False


def test_left_block_folds_i_and_l_in_the_prefix() -> None:
    # LG + Oxidation weighs as W (0.021 Da). ref ...ILG..., peptide LW...: it reads G>W after L,
    # and the shifted prefix L equals ref[k-2] = I only with I and L folded.
    ref = "DILGSWK"
    k = 3
    assert fl().mass_equivalent("LG", "W", TOL)
    assert fl().left_block_read("LWSWK", 1, ref, k, TOL) is True


def test_left_block_follows_the_configured_tolerance() -> None:
    # AG -> K: K - (A + G) = 0.036 Da, inside 0.05 and outside 0.02
    ref = "DWAGSW"
    assert fl().left_block_read("KSW", 0, ref, 3, 0.05) is True
    assert fl().left_block_read("KSW", 0, ref, 3, 0.02) is False


def test_left_block_rejects_non_residue_reference_letters() -> None:
    assert fl().left_block_read("QSW", 0, "D*GSW", 2, TOL) is False
    assert fl().left_block_read("QSW", 0, "DXGSW", 2, TOL) is False


# --------------------------------------------------------------------------- identifiable truth


@pytest.mark.parametrize(
    ("ref", "alt", "want"),
    [
        ("G", "Q", True),  # 71.04 Da apart, no B.5 modification closes it
        ("A", "W", True),
        ("N", "D", False),  # deamidation
        ("Q", "E", False),  # deamidation
        ("K", "Q", False),  # 0.036 Da
        ("V", "L", False),  # methylation
        ("D", "E", False),  # methylation
        ("S", "T", False),  # methylation
        ("G", "A", False),  # methylation
        ("A", "S", False),  # oxidation
        ("L", "I", False),  # isobaric
    ],
)
def test_identifiable_examples(ref: str, alt: str, want: bool) -> None:
    """Decision 4: truth that is mass-equivalent to the reference leaves the G1 recall
    denominator (the statistician's examples V>L, N>K-style, D>E, S>T, G>A)."""
    assert fl().identifiable(ref, alt, TOL) is want


def test_identifiable_uses_the_tolerance_it_is_given() -> None:
    assert fl().identifiable("K", "Q", 0.05) is False
    assert fl().identifiable("K", "Q", 0.02) is True
