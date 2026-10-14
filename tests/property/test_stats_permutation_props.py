"""Hypothesis properties of permutation p-values and circular shifts (issue M4-03; GIAP_SPEC.md
§7.3.5(b), §8.2, §8.3).

- shifted sets keep their size and pairwise gaps (circular distance on circular replicons, plain
  difference on linear ones); linear shifts never wrap; linear K = n - span;
- ``perm_p`` in [1/(B+1), 1], never 0; two-sided ``perm_p`` invariant to adding a constant;
- ``choose_shifts`` returns a sorted, duplicate-free subset of ``valid``;
- ``besag_clifford`` p in [1/(max_B+1), 1] and B_used <= max_B.
"""

from __future__ import annotations

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.rng import for_step
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

FAST = settings(max_examples=200, deadline=None)
EPS = 1e-12
ALT = st.sampled_from(["two-sided", "greater", "less"])
# small integer grid makes ties frequent
VAL = st.one_of(st.integers(-5, 5).map(float), st.floats(-1e3, 1e3, allow_nan=False))
VAL_OR_NAN = st.one_of(VAL, st.just(float("nan")))


@st.composite
def replicon_and_set(draw: st.DrawFn) -> tuple[int, list[int]]:
    n = draw(st.integers(1, 60))
    idx = draw(st.lists(st.integers(0, n - 1), min_size=1, max_size=min(n, 8), unique=True))
    return n, idx


# ---------------------------------------------------------------------------- shifts


@FAST
@given(replicon_and_set())
def test_circular_shifts_are_all_nonzero_k(case: tuple[int, list[int]]) -> None:
    n, idx = case
    k = circular_shifts(n, idx, circular=True)
    assert k.dtype == np.int64
    assert list(k) == list(range(1, n))


@FAST
@given(replicon_and_set())
def test_linear_shifts_subset_sorted_count_is_n_minus_span(case: tuple[int, list[int]]) -> None:
    n, idx = case
    k = [int(x) for x in circular_shifts(n, idx, circular=False)]
    assert k == sorted(set(k))
    assert set(k) <= set(range(1, n))
    assert len(k) == n - (max(idx) - min(idx) + 1)


@FAST
@given(replicon_and_set())
def test_circular_shift_keeps_size_and_circular_gaps(case: tuple[int, list[int]]) -> None:
    n, idx = case
    a = np.asarray(idx)
    for k in circular_shifts(n, idx, circular=True):
        s = apply_shift(idx, int(k), n)
        assert s.dtype == np.int64
        assert len(set(s.tolist())) == len(idx)
        assert ((s >= 0) & (s < n)).all()
        np.testing.assert_array_equal((s[:, None] - s[None, :]) % n, (a[:, None] - a[None, :]) % n)


@FAST
@given(replicon_and_set())
def test_linear_shift_never_wraps_and_keeps_plain_gaps(case: tuple[int, list[int]]) -> None:
    n, idx = case
    a = np.asarray(idx)
    for k in circular_shifts(n, idx, circular=False):
        s = apply_shift(idx, int(k), n)
        d = s - a
        assert (d == k).all() or (d == k - n).all()  # one rigid translation, no split
        np.testing.assert_array_equal(s[:, None] - s[None, :], a[:, None] - a[None, :])
        assert len(set(s.tolist())) == len(idx)


@FAST
@given(
    st.lists(st.integers(1, 200), max_size=60, unique=True),
    st.integers(1, 300),
    st.integers(1, 300),
    st.integers(0, 2**31 - 1),
)
def test_choose_shifts_sorted_unique_subset(valid: list[int], n: int, B: int, seed: int) -> None:
    n_max = 100
    shifts, exhaustive = choose_shifts(valid, n_max, n, B, for_step("context", seed))
    s = [int(x) for x in shifts]
    assert shifts.dtype == np.int64
    assert s == sorted(set(s))
    assert set(s) <= set(valid)
    if n <= n_max:
        assert bool(exhaustive) is True
        assert s == sorted(valid)
    else:
        assert len(s) == min(B, len(valid))
        if valid:  # lead decision 3: exhaustive whenever every valid shift is used
            assert bool(exhaustive) is (len(valid) <= B)


# ---------------------------------------------------------------------------- perm_p


@FAST
@given(VAL, st.lists(VAL_OR_NAN, min_size=1, max_size=50), ALT)
def test_perm_p_in_floor_one_never_zero(obs: float, null: list[float], alt: str) -> None:
    B = int(np.isfinite(null).sum())
    p = perm_p(obs, null, alt)  # type: ignore[arg-type]
    if B == 0:
        assert np.isnan(p)
        return
    assert 0 < p <= 1
    assert p >= 1 / (B + 1) - EPS
    # p is on the (1 + b)/(1 + B) grid
    b = p * (B + 1) - 1
    assert abs(b - round(b)) < 1e-9


@FAST
@given(
    st.integers(-20, 20),
    st.lists(st.integers(-20, 20), min_size=1, max_size=40),
    st.integers(-1000, 1000),
)
def test_two_sided_invariant_to_constant_shift(obs: int, null: list[int], c: int) -> None:
    p0 = perm_p(float(obs), np.asarray(null, dtype=float))
    p1 = perm_p(float(obs + c), np.asarray(null, dtype=float) + c)
    assert abs(p0 - p1) < 1e-12


@FAST
@given(VAL, st.lists(VAL, min_size=1, max_size=40))
def test_greater_and_less_mirror_under_negation(obs: float, null: list[float]) -> None:
    g = perm_p(obs, null, "greater")
    ell = perm_p(-obs, [-x for x in null], "less")
    assert abs(g - ell) < 1e-12


@FAST
@given(VAL, st.lists(VAL, min_size=1, max_size=40))
def test_greater_monotone_in_obs(obs: float, null: list[float]) -> None:
    assert perm_p(obs + 1.0, null, "greater") <= perm_p(obs, null, "greater") + 1e-12


# ---------------------------------------------------------------------------- besag_clifford


@FAST
@given(
    st.floats(-3, 3, allow_nan=False),
    st.integers(1, 10),
    st.integers(0, 200),
    st.integers(0, 2**31 - 1),
)
def test_besag_clifford_bounds(obs: float, h: int, extra: int, seed: int) -> None:
    max_B = h + extra

    def fn(rng: np.random.Generator) -> float:
        return float(rng.normal())

    p, b_used = besag_clifford(fn, obs, h, max_B, for_step("bc", seed))
    assert 0 < p <= 1
    assert p >= 1 / (max_B + 1) - EPS
    assert 1 <= b_used <= max_B
    assert b_used >= h or b_used == max_B
