"""Hypothesis properties of the Poisson-binomial DP and the chance-hit model (issue M3-09;
GIAP_SPEC.md §7.2.5, §8.3).

- pmf sums to 1 (<= 1e-12), is non-negative, equals the binomial pmf for constant p, and is
  invariant to permuting p; cdf(k) = sum pmf[:k+1]; sf(k) = sum pmf[k:];
- E_i equals the uniform E(L, m) when the query residue has a constant frequency;
  E_i is invariant to permuting the query, non-decreasing in m, non-increasing when a residue is
  appended, and linear in N.
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.mapping.chance import ResidueFreqs, expected_hits, expected_hits_uniform
from giap.mapping.sixframe import ALPHABET
from giap.stats.poisson_binomial import (
    poisson_binomial_cdf,
    poisson_binomial_pmf,
    poisson_binomial_sf,
)

FAST = settings(max_examples=200, deadline=None)
P = st.floats(0.0, 1.0, allow_nan=False)
PS = st.lists(P, min_size=0, max_size=60)
RES = ALPHABET[:19]


@FAST
@given(PS)
def test_pmf_is_a_distribution(p: list[float]) -> None:
    pmf = poisson_binomial_pmf(p)
    assert pmf.shape == (len(p) + 1,)
    assert np.all(pmf >= 0.0)
    assert abs(pmf.sum() - 1.0) <= 1e-12


@FAST
@given(st.integers(0, 80), P)
def test_constant_p_is_binomial(n: int, q: float) -> None:
    pmf = poisson_binomial_pmf([q] * n)
    ref = np.array([math.comb(n, j) * q**j * (1 - q) ** (n - j) for j in range(n + 1)])
    np.testing.assert_allclose(pmf, ref, rtol=1e-9, atol=1e-15)


@FAST
@given(PS, st.randoms(use_true_random=False))
def test_permutation_invariant(p: list[float], rnd: object) -> None:
    q = list(p)
    rnd.shuffle(q)  # type: ignore[attr-defined]
    np.testing.assert_allclose(
        poisson_binomial_pmf(q), poisson_binomial_pmf(p), rtol=1e-9, atol=1e-15
    )


@FAST
@given(PS, st.data())
def test_cdf_sf_agree_with_pmf(p: list[float], data: st.DataObject) -> None:
    k = data.draw(st.integers(0, len(p)))
    pmf = poisson_binomial_pmf(p)
    assert math.isclose(
        poisson_binomial_cdf(p, k), float(pmf[: k + 1].sum()), rel_tol=1e-9, abs_tol=1e-300
    )
    assert math.isclose(
        poisson_binomial_sf(p, k), float(pmf[k:].sum()), rel_tol=1e-9, abs_tol=1e-300
    )
    assert 0.0 <= poisson_binomial_cdf(p, k) <= 1.0 + 1e-12
    assert 0.0 <= poisson_binomial_sf(p, k) <= 1.0 + 1e-12


# --------------------------------------------------------------------------- E_i

F_VEC = (
    st.lists(st.floats(0.0, 1.0, allow_nan=False), min_size=19, max_size=19)
    .filter(lambda v: sum(v) > 1e-6)
    .map(lambda v: np.asarray(v) / sum(v))
)
QUERY = st.text(alphabet=RES, min_size=1, max_size=25)


@FAST
@given(F_VEC, QUERY, st.data())
def test_ei_permutation_invariant_and_monotone_in_m(
    f: np.ndarray, query: str, data: st.DataObject
) -> None:
    fr = ResidueFreqs(f=f, n_positions=10**7)
    m = data.draw(st.integers(0, len(query)))
    perm = "".join(data.draw(st.permutations(list(query))))
    e = expected_hits(query, m, fr)
    assert math.isclose(expected_hits(perm, m, fr), e, rel_tol=1e-9, abs_tol=1e-300)
    assert e >= 0.0 and math.isfinite(e)
    if m < len(query):
        assert expected_hits(query, m + 1, fr) >= e * (1 - 1e-12)


@FAST
@given(F_VEC, QUERY, st.sampled_from(RES), st.integers(0, 3))
def test_ei_non_increasing_when_residue_appended(
    f: np.ndarray, query: str, extra: str, m: int
) -> None:
    fr = ResidueFreqs(f=f, n_positions=10**7)
    m = min(m, len(query))
    assert expected_hits(query + extra, m, fr) <= expected_hits(query, m, fr) * (1 + 1e-12)


@FAST
@given(st.floats(0.01, 0.99), st.integers(1, 20), st.data())
def test_ei_constant_f_equals_uniform(fa: float, L: int, data: st.DataObject) -> None:
    m = data.draw(st.integers(0, L))
    f = np.full(19, (1 - fa) / 18)
    f[RES.index("A")] = fa
    fr = ResidueFreqs(f=f, n_positions=10**6)
    assert math.isclose(
        expected_hits("A" * L, m, fr),
        expected_hits_uniform(L, m, fa, 1e6),
        rel_tol=1e-9,
        abs_tol=1e-300,
    )


@FAST
@given(F_VEC, QUERY, st.integers(1, 10**9))
def test_ei_linear_in_n(f: np.ndarray, query: str, n: int) -> None:
    a = expected_hits(query, 1 if len(query) > 1 else 0, ResidueFreqs(f=f, n_positions=n))
    b = expected_hits(query, 1 if len(query) > 1 else 0, ResidueFreqs(f=f, n_positions=1))
    assert math.isclose(a, n * b, rel_tol=1e-9, abs_tol=1e-300)
