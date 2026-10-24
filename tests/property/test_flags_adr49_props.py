"""ADR-0049 flag properties (local/handoff/adr0049_decision.md, decisions 2 and 4; interfaces in
tests/unit/variants/test_flags_adr49.py).

* Swaps: any adjacent pair and its swap have the same composition, so ``mass_equivalent`` holds
  at every tolerance, and ``flags_for`` flags every ``swap`` row ``mass_equivalent``.
* ``left_block_read`` equals an independent oracle of the left-form rule, depends only on
  ``ref[k-1-q:k+1]`` and ``peptide[:q+1]`` (so no right-hand or suffix form can enter), and at
  ``q = 0`` reduces to ``mass_equivalent(ref[k-1:k+1], peptide[0])``.
* ``identifiable`` is exactly ``not mass_equivalent``.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config.model import GiapConfig

RES = "ACDEFGHIKLMNPQRSTVWY"
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
TOLS = st.sampled_from([0.0, 0.005, 0.02, 0.05, 0.1])


def fl() -> ModuleType:
    return importlib.import_module("giap.variants.flags")


def fold(s: str) -> str:
    return s.replace("I", "L")


@st.composite
def swap_pairs(draw: st.DrawFn) -> tuple[str, str]:
    a = draw(st.sampled_from(RES))
    b = draw(st.sampled_from(RES).filter(lambda x: fold(x) != fold(a)))
    return a, b


@settings(max_examples=300, deadline=None)
@given(pair=swap_pairs(), tol=TOLS)
def test_a_swapped_pair_is_mass_equivalent_to_its_reference(
    pair: tuple[str, str], tol: float
) -> None:
    a, b = pair
    assert fl().mass_equivalent(a + b, b + a, tol)


@settings(max_examples=300, deadline=None)
@given(
    pair=swap_pairs(),
    tol=TOLS,
    relation=st.sampled_from(["related", "same"]),
    strand=st.sampled_from(["+", "-"]),
    n_alt=st.integers(1, 5),
    sources=st.frozensets(
        st.tuples(st.sampled_from(["denovo", "database"]), st.sampled_from(["same", "other"])),
        max_size=3,
    ),
)
def test_flags_for_flags_every_swap_row(
    pair: tuple[str, str],
    tol: float,
    relation: str,
    strand: str,
    n_alt: int,
    sources: frozenset[tuple[str, str]],
) -> None:
    a, b = pair
    row: dict[str, Any] = {
        "variant_id": "VAR_1", "contig": "c1", "codon_start": 300, "strand": strand,
        "gene": "G1", "aa_pos": 7, "ref_codon": CODON[a] + CODON[b], "ref_aa": a,
        "alt_aa": b + a, "alt_type": "swap", "call": "ALT",
        "confidence": "high" if n_alt >= 2 else "candidate", "n_ref": 0, "n_alt": n_alt,
        "source": "denovo",
    }  # fmt: skip
    config = GiapConfig.model_validate(
        {"genome_relation": relation, "variants": {"mass_equivalent_tol_da": tol}}
    )
    got = fl().flags_for(row, frozenset(), config, sources)
    assert "mass_equivalent" in got, (row, got)
    assert row["confidence"] == ("high" if n_alt >= 2 else "candidate"), "flags never demote"


# --------------------------------------------------------------------------- left_block_read


def oracle(pep: str, q: int, ref: str, k: int, tol: float) -> bool:
    if k < 1 or not 0 <= q < len(pep) or k >= len(ref):
        return False
    if not fl().mass_equivalent(ref[k - 1 : k + 1], pep[q], tol):
        return False
    return q == 0 or (k - 1 - q >= 0 and fold(pep[:q]) == fold(ref[k - 1 - q : k - 1]))


@st.composite
def reads(draw: st.DrawFn) -> tuple[str, int, str, int]:
    """A substitution read of ``ref`` at ``k`` (peptide residue ``q``), biased towards the cases
    that matter: an equal-mass pair before the site and runs of the residue before it."""
    alpha = draw(st.sampled_from(["AGQ", "AGQKN", "ADEGLQSWK", RES]))
    ref = draw(st.text(alpha, min_size=3, max_size=14))
    k = draw(st.integers(0, len(ref) - 1))
    q = draw(st.integers(0, k))
    n_suf = draw(st.integers(0, len(ref) - 1 - k))
    alt = draw(st.sampled_from("QKNWR" + RES))
    pep = ref[k - q : k] + alt + ref[k + 1 : k + 1 + n_suf]
    return pep, q, ref, k


@settings(max_examples=500, deadline=None)
@given(r=reads(), tol=TOLS)
def test_left_block_read_matches_the_oracle(r: tuple[str, int, str, int], tol: float) -> None:
    pep, q, ref, k = r
    assert fl().left_block_read(pep, q, ref, k, tol) is oracle(pep, q, ref, k, tol), r


@settings(max_examples=300, deadline=None)
@given(
    r=reads(),
    tol=TOLS,
    ref_tail=st.text(RES, max_size=4),
    pep_tail=st.text(RES, max_size=4),
    ref_head=st.text(RES, max_size=4),
)
def test_left_block_read_ignores_everything_right_of_the_site(
    r: tuple[str, int, str, int], tol: float, ref_tail: str, pep_tail: str, ref_head: str
) -> None:
    """Only ``ref[k-1-q:k+1]`` and ``peptide[:q+1]`` matter: replacing what lies right of the
    site (where a right-hand ``ref[k:k+2]`` block would live) changes nothing, and neither does
    reference sequence left of the shifted prefix."""
    pep, q, ref, k = r
    if k - 1 - q < 0:
        return  # the pair or the shifted prefix runs off the left end; a new head would supply it
    base = fl().left_block_read(pep, q, ref, k, tol)
    lo = k - 1 - q
    ref2 = ref_head + ref[lo : k + 1] + ref_tail
    k2 = len(ref_head) + (k - lo)
    pep2 = pep[: q + 1] + pep_tail
    assert fl().left_block_read(pep2, q, ref2, k2, tol) is base, (r, ref2, k2, pep2)


@settings(max_examples=300, deadline=None)
@given(r=reads(), tol=TOLS)
def test_at_peptide_index_0_the_rule_is_the_pair_mass(
    r: tuple[str, int, str, int], tol: float
) -> None:
    pep, _, ref, k = r
    want = k >= 1 and fl().mass_equivalent(ref[k - 1 : k + 1], pep[0], tol)
    assert fl().left_block_read(pep, 0, ref, k, tol) is want


# --------------------------------------------------------------------------- identifiable


@settings(max_examples=500, deadline=None)
@given(
    ref=st.text(RES + "J", min_size=1, max_size=2),
    alt=st.text(RES + "J", min_size=1, max_size=2),
    tol=TOLS,
)
def test_identifiable_is_not_mass_equivalent(ref: str, alt: str, tol: float) -> None:
    assert fl().identifiable(ref, alt, tol) is (not fl().mass_equivalent(ref, alt, tol))
