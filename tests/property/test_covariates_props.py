"""Hypothesis properties of the M4-10 protein covariates (docs/specs/M4-10-protein-covariates.md).

Acceptance: ``"".join(digest(s, e)) == s``; ``n_observable_peptides`` <= number of fragments;
``tm_segments`` in [0, 10] and 0 for sequences of D, E, K, R only. Further invariants checked
against the independent oracle in tests/helpers/covariates_oracle.py.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import covariates_oracle as oracle

SETTINGS = settings(max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow])
ENZYMES = st.sampled_from(["trypsin_p", "trypsin"])
# Biased toward cleavage sites and P so fragments are short and K/R-P pairs are common.
proteins = st.text(alphabet=st.sampled_from(list(oracle.STANDARD_AA) + ["K", "R", "P"] * 3))
hydrophilic = st.text(alphabet=st.sampled_from("DEKR"), max_size=400)
tm_alphabet = st.sampled_from([*"LIVFAKRDE", "X"])
tm_seqs = st.text(alphabet=tm_alphabet, max_size=300)


def _cov():
    from giap.genetable import covariates as cov

    return cov


@SETTINGS
@given(proteins, ENZYMES)
def test_digest_concatenates_back(s: str, enzyme: str) -> None:
    frags = _cov().digest(s, enzyme)  # type: ignore[arg-type]
    assert "".join(frags) == s
    assert all(frags)  # no empty fragments


@SETTINGS
@given(proteins, ENZYMES, st.integers(min_value=0, max_value=3))
def test_digest_matches_oracle(s: str, enzyme: str, missed: int) -> None:
    got = _cov().digest(s, enzyme, missed)  # type: ignore[arg-type]
    exp = oracle.digest(s, enzyme, missed)
    if missed == 0:
        assert got == exp
    else:
        assert sorted(got) == sorted(exp)


@SETTINGS
@given(proteins)
def test_trypsin_fragments_are_unions_of_trypsin_p_fragments(s: str) -> None:
    cov = _cov()
    fine = cov.digest(s, "trypsin_p")
    coarse = cov.digest(s, "trypsin")
    cuts_fine = {sum(len(f) for f in fine[: i + 1]) for i in range(len(fine))}
    cuts_coarse = {sum(len(f) for f in coarse[: i + 1]) for i in range(len(coarse))}
    assert cuts_coarse <= cuts_fine
    assert len(coarse) <= len(fine)


@SETTINGS
@given(proteins, ENZYMES)
def test_observable_bounded_by_fragments_and_matches_oracle(s: str, enzyme: str) -> None:
    cov = _cov()
    n = cov.n_observable_peptides(s, enzyme)  # type: ignore[arg-type]
    assert 0 <= n <= len(cov.digest(s, enzyme))  # type: ignore[arg-type]
    assert n == oracle.n_observable(s, enzyme)


@SETTINGS
@given(proteins, ENZYMES, st.integers(1, 40), st.integers(0, 40))
def test_observable_monotone_in_length_window(s: str, enzyme: str, lo: int, width: int) -> None:
    cov = _cov()
    inner = cov.n_observable_peptides(s, enzyme, (lo, lo + width))  # type: ignore[arg-type]
    outer = cov.n_observable_peptides(s, enzyme, (max(1, lo - 1), lo + width + 1))  # type: ignore[arg-type]
    assert inner <= outer


@SETTINGS
@given(proteins, ENZYMES)
def test_observable_counts_distinct_peptides(s: str, enzyme: str) -> None:
    # A tandem duplicate of a sequence ending at a cleavage site yields the same distinct set.
    cov = _cov()
    unit = s + "K"
    if enzyme == "trypsin" and unit.startswith("P"):
        unit = "A" + unit
    assert cov.n_observable_peptides(unit * 2, enzyme) == cov.n_observable_peptides(unit, enzyme)  # type: ignore[arg-type]


@SETTINGS
@given(tm_seqs)
def test_tm_bounds_and_oracle(s: str) -> None:
    n = _cov().tm_segments(s)
    assert 0 <= n <= 10
    assert n == oracle.tm_segments(s)


@SETTINGS
@given(hydrophilic)
def test_tm_zero_for_hydrophilic(s: str) -> None:
    assert _cov().tm_segments(s) == 0


@SETTINGS
@given(
    st.text(alphabet=tm_alphabet, min_size=19, max_size=200),
    st.text(alphabet=tm_alphabet, min_size=19, max_size=200),
    st.integers(min_value=19, max_value=40),
)
def test_tm_additive_across_arginine_spacer(a: str, b: str, k: int) -> None:
    # R (-4.5) is the most hydrophilic residue: windows reaching from A (or B) into the spacer
    # never pass unless the adjacent all-A (all-B) window does; spacer-only windows never pass.
    cov = _cov()
    big = 10**6
    joined = a + "R" * k + b
    total = cov.tm_segments(a, cap=big) + cov.tm_segments(b, cap=big)
    assert cov.tm_segments(joined, cap=big) == total
    assert cov.tm_segments(joined) == min(10, total)


@SETTINGS
@given(st.integers(min_value=19, max_value=400))
def test_tm_single_hydrophobic_run(n: int) -> None:
    assert _cov().tm_segments("L" * n) == 1


@SETTINGS
@given(st.text(alphabet=st.sampled_from(list(oracle.STANDARD_AA) + list("XBZJUO")), max_size=120))
def test_covariates_bounds_and_biopython(s: str) -> None:
    c = _cov().covariates(s, "trypsin_p", (7, 30))
    g, p, m = oracle.biopython_props(s)
    if math.isnan(g):  # nothing left after cleaning; length_aa not pinned here (spec ambiguity)
        assert math.isnan(c.gravy) and math.isnan(c.pi) and math.isnan(c.mw_da)
        return
    assert c.length_aa == len(s)
    assert -4.5 <= c.gravy <= 4.5 and 0.0 <= c.pi <= 14.0 and c.mw_da > 0.0
    assert abs(c.gravy - g) <= 1e-9 and abs(c.pi - p) <= 1e-9 and abs(c.mw_da - m) <= 1e-9
    assert 0 <= c.tm_segments <= 10


@SETTINGS
@given(
    st.lists(
        st.sampled_from(sorted(k for k, v in oracle.TABLE11.items() if v != "*")),
        min_size=1,
        max_size=60,
    )
)
def test_pseudogene_translation_matches_oracle(codons: list[str]) -> None:
    nt = "".join(codons) + "TGA" + "GCT" * 3
    got = _cov().protein_sequence(None, nt, "pseudogene")
    assert got == oracle.translate_pseudogene(nt)
    assert len(got) == len(codons) and got[0] == "M"
