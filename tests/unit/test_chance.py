"""Chance-hit model E(L, m) and the ranking score E_i (issue M3-09; GIAP_SPEC.md §7.2.5, G3,
Appendix B.2).

E(L, m) = N * sum_{j<=m} C(L, j) (1 - s)^j s^(L - j), s = sum f_a^2, N = six-frame positions.
E_i = N * P(<= m mismatches) with per-position match probability f_{a_i} (Poisson-binomial).

Conventions pinned here:
- ``ResidueFreqs.f`` is indexed by residue code (``giap.mapping.sixframe.ALPHABET[:19]``, I -> L);
  it sums to 1 over residues and excludes ``*`` and ``X``.
- ``n_positions`` = sum of ``len(frame.aa)`` over the target frames (stops and X included, as the
  spec states), with circular-extension positions excluded.
- a query residue outside the 19 codes (``NOMATCH``, e.g. ``U``) has f = 0.
- ``m > L`` is a ``ValueError``.
"""

from __future__ import annotations

import dataclasses
import itertools
import math

import numpy as np
import pytest

from giap.mapping.chance import (
    ResidueFreqs,
    expected_hits,
    expected_hits_uniform,
    residue_freqs,
)
from giap.mapping.sixframe import ALPHABET, Frame, SixFrameIndex, build_sixframe, encode_target
from giap.model.coords import Replicon

S = 0.069
N7 = 1e7
RES = ALPHABET[:19]


def _closed_form(L: int, m: int, s: float = S, n: float = N7) -> float:
    return n * sum(math.comb(L, j) * (1 - s) ** j * s ** (L - j) for j in range(m + 1))


def _freqs(mapping: dict[str, float], n_positions: int = 10**7) -> ResidueFreqs:
    """19-vector with the given letters; the rest share the remaining mass evenly."""
    f = np.zeros(19)
    for a, v in mapping.items():
        f[RES.index(a)] = v
    rest = [i for i in range(19) if RES[i] not in mapping]
    f[rest] = (1.0 - sum(mapping.values())) / len(rest)
    return ResidueFreqs(f=f, n_positions=n_positions)


# ---------------------------------------------------------------- B.2 table (examples 1-2)

# Example 1: six significant digits as printed in the M3-09 spec.
SPEC_6SIG = {
    (7, 0): 0.0744635,
    (8, 1): 0.559742,
    (9, 2): 2.36692,
    (10, 1): 0.00332505,
    (12, 2): 0.00141835,
    (14, 1): 1.05295e-7,
}

# Appendix B.2 exactly as printed (value, significant digits shown).
B2 = {
    6: [(1.1, 2), (88, 2), (3.0e3, 2)],
    7: [(0.074, 2), (7.1, 2), (2.9e2, 2)],
    8: [(0.0051, 2), (0.56, 2), (27, 2)],
    9: [(3.5e-4, 2), (0.043, 2), (2.4, 2)],
    10: [(2.4e-5, 2), (0.0033, 2), (0.20, 2)],
    12: [(1.2e-7, 2), (1.9e-5, 2), (0.0014, 2)],
    14: [(5.5e-10, 2), (1.1e-7, 2), (9.3e-6, 2)],
}


def _round_sig(x: float, sig: int) -> float:
    return float(f"{x:.{sig - 1}e}")


@pytest.mark.parametrize(("L", "m"), sorted(SPEC_6SIG))
def test_example1_spec_values(L: int, m: int) -> None:
    got = expected_hits_uniform(L, m, S, N7)
    assert isinstance(got, float)
    assert got == pytest.approx(SPEC_6SIG[(L, m)], rel=1e-5)  # printed to 6 sig. digits


@pytest.mark.parametrize("L", [6, 7, 8, 9, 10, 12, 14])
@pytest.mark.parametrize("m", [0, 1, 2])
def test_example1_exact_formula(L: int, m: int) -> None:
    assert expected_hits_uniform(L, m, S, N7) == pytest.approx(_closed_form(L, m), rel=1e-9)


@pytest.mark.parametrize("L", sorted(B2))
@pytest.mark.parametrize("m", [0, 1, 2])
def test_examples1_2_b2_table_rounds_to_printed(L: int, m: int) -> None:
    value, sig = B2[L][m]
    assert _round_sig(expected_hits_uniform(L, m, S, N7), sig) == pytest.approx(value, rel=1e-12)


def test_example2_l6() -> None:
    got = [expected_hits_uniform(6, m, S, N7) for m in range(3)]
    assert [_round_sig(x, 2) for x in got] == [1.1, 88.0, 3.0e3]


def test_uniform_is_linear_in_n_positions() -> None:
    a = expected_hits_uniform(9, 1, S, 1e6)
    b = expected_hits_uniform(9, 1, S, 1e7)
    assert b == pytest.approx(10 * a, rel=1e-12)
    assert expected_hits_uniform(9, 1, S, 0) == 0.0


# ---------------------------------------------------------------- monotonicity


def test_uniform_monotone_in_m_and_l() -> None:
    for L in range(1, 30):
        vals = [expected_hits_uniform(L, m, S, N7) for m in range(L + 1)]
        assert all(x < y for x, y in itertools.pairwise(vals))
        assert vals[-1] == pytest.approx(N7, rel=1e-12)  # m = L: every window
    for m in range(3):
        vals = [expected_hits_uniform(L, m, S, N7) for L in range(max(m, 1) + 1, 30)]
        assert all(x > y for x, y in itertools.pairwise(vals))


def test_uniform_monotone_in_s() -> None:
    vals = [expected_hits_uniform(10, 1, s, N7) for s in (0.03, 0.05, 0.069, 0.1, 0.3)]
    assert all(x < y for x, y in itertools.pairwise(vals))


def test_ranking_score_monotone() -> None:
    fr = _freqs({"A": 0.1, "C": 0.05})
    # more mismatches allowed -> more chance hits
    assert (
        expected_hits("ACACAC", 0, fr)
        < expected_hits("ACACAC", 1, fr)
        < expected_hits("ACACAC", 2, fr)
    )
    # longer peptide -> fewer chance hits
    assert expected_hits("ACACACA", 1, fr) < expected_hits("ACACAC", 1, fr)
    # a residue with higher genome frequency -> more chance hits
    assert expected_hits("AAAAAA", 1, fr) > expected_hits("CCCCCC", 1, fr)


# ---------------------------------------------------------------- E_i (examples 5, 6, 8)


def test_example5_heterogeneous_query() -> None:
    fr = _freqs(
        {"A": 0.1, "C": 0.1, "D": 0.1, "E": 0.1, "F": 0.05, "G": 0.05, "H": 0.05, "K": 0.05}
    )
    q = "ACDEFGHK"
    assert expected_hits(q, 0, fr) == pytest.approx(0.00625, rel=1e-9)
    assert expected_hits(q, 1, fr) == pytest.approx(0.70625, rel=1e-9)
    assert expected_hits(q, 2, fr) == pytest.approx(34.38125, rel=1e-9)


def test_example5_order_of_residues_irrelevant() -> None:
    fr = _freqs(
        {"A": 0.1, "C": 0.1, "D": 0.1, "E": 0.1, "F": 0.05, "G": 0.05, "H": 0.05, "K": 0.05}
    )
    for m in range(3):
        assert expected_hits("KHGFEDCA", m, fr) == pytest.approx(
            expected_hits("ACDEFGHK", m, fr), rel=1e-12
        )


def test_example6_constant_f_equals_uniform() -> None:
    # a query of one letter with f = 0.069 everywhere
    fr = _freqs({"W": S})
    got = expected_hits("WWWWWWWW", 1, fr)
    assert got == pytest.approx(expected_hits_uniform(8, 1, S, 1e7), rel=1e-9)
    assert got == pytest.approx(0.559742, rel=1e-5)


def test_ranking_score_scales_with_n_positions() -> None:
    f = _freqs({"A": 0.1}).f
    a = expected_hits("AAAAAAAA", 1, ResidueFreqs(f=f, n_positions=10**6))
    b = expected_hits("AAAAAAAA", 1, ResidueFreqs(f=f, n_positions=10**7))
    assert b == pytest.approx(10 * a, rel=1e-12)


def test_example8_nomatch_residue_has_zero_frequency() -> None:
    fr = _freqs({"A": 0.1})
    assert expected_hits("AAAUAAA", 0, fr) == 0.0
    # with one mismatch allowed U costs exactly one: N * P(all 6 A match) * 1
    assert expected_hits("AAAUAAA", 1, fr) == pytest.approx(1e7 * 0.1**6, rel=1e-9)
    assert expected_hits("AAUUAAA", 1, fr) == 0.0


def test_isoleucine_is_leucine() -> None:
    fr = _freqs({"L": 0.1, "K": 0.05})
    for m in range(3):
        assert expected_hits("KILIK", m, fr) == expected_hits("KLLLK", m, fr)


# ---------------------------------------------------------------- edge cases


def test_m_greater_than_l_is_value_error() -> None:
    fr = _freqs({"A": 0.1})
    with pytest.raises(ValueError):
        expected_hits_uniform(5, 6, S, N7)
    with pytest.raises(ValueError):
        expected_hits("AAAAA", 6, fr)


def test_m_equal_l_counts_every_position() -> None:
    fr = _freqs({"A": 0.1})
    assert expected_hits("AAAAA", 5, fr) == pytest.approx(1e7, rel=1e-12)
    assert expected_hits_uniform(5, 5, S, N7) == pytest.approx(N7, rel=1e-12)


def test_tiny_values_positive_and_finite() -> None:
    e = expected_hits_uniform(40, 0, S, N7)
    assert e > 0.0 and math.isfinite(e)
    assert e == pytest.approx(N7 * S**40, rel=1e-9)
    fr = _freqs({"A": 0.5})
    ei = expected_hits("A" * 40, 0, fr)  # 1e7 * 0.5**40 ~ 9.1e-6
    assert ei == pytest.approx(1e7 * 0.5**40, rel=1e-9)
    tiny = expected_hits("W" * 40, 0, _freqs({"W": 0.01}))  # 1e7 * 1e-80
    assert tiny > 0.0
    assert tiny == pytest.approx(1e7 * 0.01**40, rel=1e-9)


def test_results_are_python_floats() -> None:
    fr = _freqs({"A": 0.1})
    assert isinstance(expected_hits("AAAAAAA", 1, fr), float)
    assert isinstance(expected_hits_uniform(7, 1, S, N7), float)


# ---------------------------------------------------------------- ResidueFreqs / residue_freqs


def test_residue_freqs_is_frozen_with_s_property() -> None:
    fr = _freqs({"A": 0.2, "L": 0.3})
    assert dataclasses.is_dataclass(fr)
    with pytest.raises(dataclasses.FrozenInstanceError):
        fr.n_positions = 5  # type: ignore[misc]
    assert fr.s == pytest.approx(float((fr.f**2).sum()), rel=1e-12)
    assert isinstance(fr.s, float)


def _linear_index(seqs: dict[str, str]) -> SixFrameIndex:
    """One frame per contig (offset 0, + strand), linear replicons of 3 * len(aa) nt."""
    reps = {c: Replicon(c, 3 * len(aa), "linear") for c, aa in seqs.items()}
    frames = tuple(Frame(c, "+", 0, encode_target(aa), 3 * len(aa)) for c, aa in seqs.items())
    return SixFrameIndex(replicons=reps, frames=frames)


def test_example7_residue_freqs_merges_i_into_l() -> None:
    fr = residue_freqs(_linear_index({"c1": "AAAL", "c2": "ILLK"}))
    assert fr.f.shape == (19,)
    assert fr.f[RES.index("A")] == pytest.approx(3 / 8)
    assert fr.f[RES.index("L")] == pytest.approx(4 / 8)
    assert fr.f[RES.index("K")] == pytest.approx(1 / 8)
    assert fr.f.sum() == pytest.approx(1.0, abs=1e-12)
    assert fr.n_positions == 8
    assert isinstance(fr.n_positions, int)
    assert fr.s == pytest.approx((3 / 8) ** 2 + (4 / 8) ** 2 + (1 / 8) ** 2)


def test_residue_freqs_excludes_stop_and_x_from_f_but_not_from_n() -> None:
    fr = residue_freqs(_linear_index({"c1": "AA*X", "c2": "KKAA"}))
    assert fr.f.sum() == pytest.approx(1.0, abs=1e-12)
    assert fr.f[RES.index("A")] == pytest.approx(4 / 6)
    assert fr.f[RES.index("K")] == pytest.approx(2 / 6)
    assert fr.n_positions == 8  # N = sum len(frame.aa)


def test_residue_freqs_matches_manual_count_on_sixframe() -> None:
    rng = np.random.default_rng(7)
    seq = "".join(rng.choice(list("ACGT"), size=601))
    reps = {"chr": Replicon("chr", 601, "linear")}
    idx = build_sixframe({"chr": seq}, reps)
    fr = residue_freqs(idx)
    codes = np.concatenate([fr_.aa for fr_ in idx.frames])
    counts = np.bincount(codes[codes < 19], minlength=19).astype(float)
    np.testing.assert_allclose(fr.f, counts / counts.sum(), rtol=1e-12)
    assert fr.n_positions == sum(len(x.aa) for x in idx.frames)


def test_circular_extension_not_counted() -> None:
    rng = np.random.default_rng(11)
    seq = "".join(rng.choice(list("ACGT"), size=300))
    lin = residue_freqs(build_sixframe({"c": seq}, {"c": Replicon("c", 300, "linear")}))
    circ_idx = build_sixframe({"c": seq}, {"c": Replicon("c", 300, "circular")})
    circ = residue_freqs(circ_idx)
    with_ext = sum(len(x.aa) for x in circ_idx.frames)
    assert with_ext > 2 * 300  # sanity: the 120-nt wrap adds ~ 6 * 40 positions
    # at most one codon start per nt per strand; at least the linear count
    assert lin.n_positions <= circ.n_positions <= 2 * 300
    assert circ.n_positions < with_ext


def test_deterministic() -> None:
    rng_seq = np.random.default_rng(3)
    seq = "".join(rng_seq.choice(list("ACGT"), size=900))
    reps = {"c": Replicon("c", 900, "circular")}
    a = residue_freqs(build_sixframe({"c": seq}, reps))
    b = residue_freqs(build_sixframe({"c": seq}, reps))
    np.testing.assert_array_equal(a.f, b.f)
    assert a.n_positions == b.n_positions
    assert expected_hits("PEPTIDEK", 1, a) == expected_hits("PEPTIDEK", 1, b)
