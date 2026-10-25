"""The pure estimator API of ``validation/pilot/entrapment.py`` (M10-07).

Spec: docs/specs/M10-07-free-entrapment.md, "Estimates", examples 3-7, "Edge cases" (r undefined,
empty input, above one), acceptance criterion 3, and the binding "Statistician decisions": S3
(``ucb95``, its table and the brute-force criterion S-2, the ordering criterion S-3), S6 (labels
and boundaries in ``fractions.Fraction``), and the pure API

    combined(n_t, n_e, n_k12_only, n_shared) -> Fraction | None
    lower(n_t, n_e) -> Fraction | None
    ucb95(n_t, n_e, n_k12_only, n_shared, level=0.05) -> float | None
    assess(combined, lower, alpha=Fraction(1, 100)) -> str

Interface fixed here (the spec leaves it open): ``assess(None, None)`` is ``"not_evaluable"``
(the empty-input label); an entrapment discovery with ``n_k12_only = 0`` raises
``AssertionError`` or ``ValueError`` ("an assertion guards it").

S-3 as written ("``ucb95`` >= ``fdp_combined`` whenever both are defined and N_E >= 1") is false
in two corners of the input space: ``fdp_combined`` > 1 (the bound inverts over M <= N, so it is
at most 1; the spec's own per-run WT_3 row has 3.0 and 1.0), and pi = k/(k+s) close to 1 (for
example N_E = 1, k = 39, s = 1). The tests below check the true statement
``ucb95 >= min(fdp_combined, 1)`` for pi <= 0.9 (r <= 9) and pin both corners; the wording is
reported to the statistician.
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Any

import numpy as np
import pytest
from helpers.entrapment_m10_07 import load_entrapment
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.stats import binom

from giap.rng import for_step

ALPHA = Fraction(1, 100)
SEED = 20261014 + 707  # seeds.calibration_base + 707 (the attack's offset), read-only here
LABELS = {"consistent_with_control", "evidence_of_failure", "inconclusive"}


def E() -> Any:
    return load_entrapment()


def brute_ucb95(n_t: int, n_e: int, k: int, s: int, level: float = 0.05) -> float:
    """S3: ``max{M in N_E..N : BinomCDF(N_E; M, pi) > level} / N``, by a full scan."""
    n = n_t + n_e
    m = np.arange(n_e, n + 1)
    ok = m[binom.cdf(n_e, m, k / (k + s)) > level]
    return float(ok.max()) / n


def _combined_ref(n_t: int, n_e: int, k: int, s: int) -> Fraction:
    return Fraction(n_e * (k + s), k * (n_t + n_e))


# ------------------------------------------------------------------ examples 3-7 (and S3 table)

#: (N_T, N_E, k12_only, shared, fdp_combined, fdp_lower, ucb95, label)
EXAMPLES = {
    "3": (4950, 20, 25000, 100000, Fraction(100, 4970), Fraction(20, 4970), 140 / 4970,
          "inconclusive"),
    "4": (9900, 10, 50000, 100000, Fraction(30, 9910), Fraction(10, 9910), 47 / 9910,
          "consistent_with_control"),
    "5": (900, 100, 1, 1, Fraction(1, 5), Fraction(1, 10), 0.225, "evidence_of_failure"),
    "6": (2990, 10, 1, 2, Fraction(1, 100), Fraction(10, 3000), 47 / 3000,
          "consistent_with_control"),
    "6-": (2989, 10, 1, 2, Fraction(30, 2999), Fraction(10, 2999), 47 / 2999, "inconclusive"),
    "7": (5000, 0, 2, 5, Fraction(0), Fraction(0), 8 / 5000, "consistent_with_control"),
}  # fmt: skip


@pytest.mark.parametrize("name", sorted(EXAMPLES))
def test_examples_3_to_7(name: str) -> None:
    n_t, n_e, k, s, comb, low, ucb, label = EXAMPLES[name]
    c, lo = E().combined(n_t, n_e, k, s), E().lower(n_t, n_e)
    assert isinstance(c, Fraction) and isinstance(lo, Fraction)
    assert c == comb and lo == low  # exact rationals
    assert E().ucb95(n_t, n_e, k, s) == pytest.approx(ucb, abs=1e-12)
    assert E().assess(c, lo) == label
    assert E().assess(c, lo, alpha=ALPHA) == label


def test_example_3_and_4_decimals() -> None:
    assert float(E().combined(4950, 20, 25000, 100000)) == pytest.approx(0.0201207, abs=5e-8)
    assert float(E().lower(4950, 20)) == pytest.approx(0.0040241, abs=5e-8)
    assert float(E().combined(9900, 10, 50000, 100000)) == pytest.approx(0.0030272, abs=5e-8)
    assert float(E().lower(9900, 10)) == pytest.approx(0.0010091, abs=5e-8)
    assert float(E().combined(2989, 10, 1, 2)) == pytest.approx(0.0100033, abs=5e-8)


def test_example_6_boundary_is_exact_and_inclusive() -> None:
    """``fdp_combined`` = 0.01 exactly; ``<= alpha`` is inclusive (S6)."""
    c = E().combined(2990, 10, 1, 2)
    assert c == ALPHA
    assert E().assess(c, E().lower(2990, 10)) == "consistent_with_control"


def test_example_2_ucb95_rows() -> None:
    """S3 table: example 2 peptide (pi 1/3), precursor (pi 1/4) and per-run WT_3 rows."""
    assert E().ucb95(3, 1, 2, 4) == pytest.approx(1.0, abs=1e-12)
    assert E().ucb95(4, 1, 6, 18) == pytest.approx(1.0, abs=1e-12)
    assert E().ucb95(0, 1, 2, 4) == pytest.approx(1.0, abs=1e-12)


def test_combined_above_one_is_not_clipped() -> None:
    """Example 2 per run WT_3: N_T 0, N_E 1, r 0.5 gives 3.0 (reported as computed)."""
    assert E().combined(0, 1, 2, 4) == Fraction(3)
    assert E().lower(0, 1) == Fraction(1)


# ------------------------------------------------------------------ assess (S6)


def test_assess_labels_and_boundaries_in_fractions() -> None:
    a = E().assess
    tiny = Fraction(1, 10**12)
    assert a(ALPHA, ALPHA / 3) == "consistent_with_control"
    assert a(ALPHA + tiny, ALPHA / 3) == "inconclusive"
    assert a(Fraction(1, 2), ALPHA) == "inconclusive"  # lower == alpha is not > alpha
    assert a(Fraction(1, 2), ALPHA + tiny) == "evidence_of_failure"
    assert a(Fraction(0), Fraction(0)) == "consistent_with_control"


def test_assess_alpha_is_a_parameter() -> None:
    a = E().assess
    assert (
        a(Fraction(3, 100), Fraction(1, 100), alpha=Fraction(5, 100)) == "consistent_with_control"
    )
    assert a(Fraction(3, 100), Fraction(1, 100)) == "inconclusive"


def test_assess_empty_is_not_evaluable() -> None:
    assert E().assess(None, None) == "not_evaluable"


@settings(max_examples=300, derandomize=True, deadline=None)
@given(st.integers(0, 10_000), st.integers(0, 500), st.integers(1, 10**6), st.integers(1, 10**6))
def test_assess_reading_rules(n_t: int, n_e: int, k: int, s: int) -> None:
    """The three readings are exclusive and follow lower <= combined (S6)."""
    if n_t + n_e == 0:
        return
    c, lo = E().combined(n_t, n_e, k, s), E().lower(n_t, n_e)
    assert lo <= c
    label = E().assess(c, lo)
    assert label in LABELS
    if lo > ALPHA:
        assert label == "evidence_of_failure"
    elif c <= ALPHA:
        assert label == "consistent_with_control"
    else:
        assert label == "inconclusive"


# ------------------------------------------------------------------ None cases (edge cases)


def test_no_shared_peptides_gives_none_everywhere() -> None:
    assert E().combined(10, 2, 5, 0) is None
    assert E().ucb95(10, 2, 5, 0) is None


def test_no_entrapment_space_gives_none_combined_but_a_lower_bound() -> None:
    """``|k12_only| = 0``: 1/r undefined; ``fdp_combined`` and ``ucb95`` are None even when
    N_E = 0, and ``fdp_lower`` is still 0."""
    assert E().combined(10, 0, 0, 5) is None
    assert E().ucb95(10, 0, 0, 5) is None
    assert E().lower(10, 0) == Fraction(0)


def test_entrapment_discovery_without_entrapment_space_is_guarded() -> None:
    with pytest.raises((AssertionError, ValueError)):
        E().combined(10, 1, 0, 5)


def test_empty_discoveries_give_none() -> None:
    assert E().combined(0, 0, 3, 5) is None
    assert E().lower(0, 0) is None
    assert E().ucb95(0, 0, 3, 5) is None


# ------------------------------------------------------------------ criterion 3 and S-2 (seeded)


def _seeded_inputs(name: str, n: int = 200) -> list[tuple[int, int, int, int]]:
    rng = for_step(name, SEED)
    out = []
    for _ in range(n):
        n_e = int(rng.integers(0, 301))
        n_t = int(rng.integers(0 if n_e else 1, 3001))
        k, s = (int(x) for x in rng.integers(1, 100_001, size=2))
        out.append((n_t, n_e, k, s))
    return out


def test_criterion_3_combined_equals_the_fraction_formula() -> None:
    for n_t, n_e, k, s in _seeded_inputs("test_m10_07_formula"):
        got = E().combined(n_t, n_e, k, s)
        assert got == _combined_ref(n_t, n_e, k, s), (n_t, n_e, k, s)


def test_s2_ucb95_equals_a_brute_force_scan() -> None:
    for n_t, n_e, k, s in _seeded_inputs("test_m10_07_ucb95"):
        got = E().ucb95(n_t, n_e, k, s)
        assert got is not None
        assert got == pytest.approx(brute_ucb95(n_t, n_e, k, s), abs=1e-12), (n_t, n_e, k, s)


def test_s2_brute_force_on_small_inputs_with_large_pi() -> None:
    """Small N and every pi of a grid (including pi near 1), where a bisection is easy to get
    wrong at the ends of {N_E, ..., N}."""
    for n_t in (0, 1, 3, 17):
        for n_e in (0, 1, 2, 9):
            if n_t + n_e == 0:
                continue
            for k, s in ((1, 39), (1, 3), (1, 1), (3, 1), (39, 1)):
                assert E().ucb95(n_t, n_e, k, s) == pytest.approx(
                    brute_ucb95(n_t, n_e, k, s), abs=1e-12
                ), (n_t, n_e, k, s)


def test_ucb95_level_is_a_parameter() -> None:
    assert E().ucb95(4950, 20, 25000, 100000, level=0.05) == pytest.approx(140 / 4970, abs=1e-12)
    want = brute_ucb95(4950, 20, 25000, 100000, level=0.2)
    assert E().ucb95(4950, 20, 25000, 100000, level=0.2) == pytest.approx(want, abs=1e-12)
    assert want < 140 / 4970


# ------------------------------------------------------------------ S-3 (properties)


@settings(max_examples=300, derandomize=True, deadline=None)
@given(st.integers(0, 3000), st.integers(1, 300), st.integers(1, 10**5), st.integers(1, 10**5))
def test_s3_ucb95_is_at_least_the_capped_point_estimate(n_t: int, n_e: int, k: int, s: int) -> None:
    c, u = E().combined(n_t, n_e, k, s), E().ucb95(n_t, n_e, k, s)
    assert u is not None and c is not None
    n = n_t + n_e
    assert round(u * n) >= min(math.floor(c * n), n)  # S-3 as reworded: every pi
    if k > 9 * s:  # pi > 0.9: the corner pinned below
        return
    assert u >= min(float(c), 1.0) - 1e-12
    assert 0.0 < u <= 1.0


@settings(max_examples=200, derandomize=True, deadline=None)
@given(st.integers(1, 5000), st.integers(1, 10**5), st.integers(1, 10**5))
def test_s3_with_no_entrapment_discovery_ucb95_is_non_negative(n_t: int, k: int, s: int) -> None:
    assert E().combined(n_t, 0, k, s) == 0
    u = E().ucb95(n_t, 0, k, s)
    assert u is not None and 0.0 <= u <= 1.0


def test_s3_corners_where_the_written_criterion_fails() -> None:
    """Pinned facts behind the module docstring's note on S-3."""
    # fdp_combined above one: the bound is capped at M = N
    assert float(E().combined(0, 1, 2, 4)) == 3.0
    assert E().ucb95(0, 1, 2, 4) == pytest.approx(1.0, abs=1e-12)
    # pi = 39/40: M_max = N_E = 1 < N_E / pi
    assert E().ucb95(1, 1, 39, 1) == pytest.approx(0.5, abs=1e-12)
    assert float(E().combined(1, 1, 39, 1)) > 0.5


@settings(max_examples=200, derandomize=True, deadline=None)
@given(st.integers(0, 2000), st.integers(0, 200), st.integers(1, 10**4), st.integers(1, 10**4))
def test_ucb95_is_a_rational_m_over_n(n_t: int, n_e: int, k: int, s: int) -> None:
    """ucb95 = M_max / N with N_E <= M_max <= N."""
    n = n_t + n_e
    if n == 0:
        return
    u = E().ucb95(n_t, n_e, k, s)
    assert u is not None
    m = round(u * n)
    assert abs(u - m / n) < 1e-12
    assert n_e <= m <= n


@settings(max_examples=200, derandomize=True, deadline=None)
@given(st.integers(1, 2000), st.integers(0, 200), st.integers(1, 10**4), st.integers(1, 10**4))
def test_ucb95_grows_with_n_e(n_t: int, n_e: int, k: int, s: int) -> None:
    """More entrapment discoveries among the same N never lower the bound."""
    a = E().ucb95(n_t, n_e, k, s)
    b = E().ucb95(n_t - 1, n_e + 1, k, s)
    assert a is not None and b is not None and b >= a - 1e-12
