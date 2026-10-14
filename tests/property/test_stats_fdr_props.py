"""Hypothesis properties of the FDR engine (issue M4-01; GIAP_SPEC.md §8.3).

- ``bh``: q >= p, q <= 1, q monotone in p, NaN positions preserved, invariant to input order;
- ``tda_qvalues``: q non-increasing in score, adding a decoy never lowers any q, invariant to
  input order, q in [0, 1];
- ``monotone_q``: non-decreasing from strict to loose, <= the clipped input;
- ``picked_fdr`` with every decoy absent equals ``tda_qvalues(target, [])``.
"""

from __future__ import annotations

import numpy as np
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.stats.multitest import bh
from giap.stats.tda import monotone_q, picked_fdr, tda_qvalues

FAST = settings(max_examples=200, deadline=None)
EPS = 1e-12

P = st.floats(0.0, 1.0, allow_nan=False)
P_OR_NAN = st.one_of(P, st.just(float("nan")))
# small integer grid makes ties frequent
SCORE = st.one_of(st.integers(-5, 5).map(float), st.floats(-10, 10, allow_nan=False))
SCORE_OR_NAN = st.one_of(SCORE, st.just(float("nan")))


def _perm(n: int) -> st.SearchStrategy[list[int]]:
    return st.permutations(list(range(n)))


# ---------------------------------------------------------------------------- bh


@FAST
@given(st.lists(P_OR_NAN, max_size=60))
def test_bh_q_at_least_p_at_most_one_nan_preserved(p: list[float]) -> None:
    pa = np.asarray(p, dtype=float)
    q = bh(pa)
    assert q.shape == pa.shape
    np.testing.assert_array_equal(np.isnan(q), np.isnan(pa))
    ok = ~np.isnan(pa)
    assert np.all(q[ok] >= pa[ok] - EPS)
    assert np.all(q[ok] <= 1.0)


@FAST
@given(st.lists(P_OR_NAN, max_size=60))
def test_bh_monotone_in_p(p: list[float]) -> None:
    pa = np.asarray(p, dtype=float)
    q = bh(pa)
    ok = ~np.isnan(pa)
    ps, qs = pa[ok], q[ok]
    order = np.argsort(ps, kind="stable")
    assert np.all(np.diff(qs[order]) >= -EPS)
    # equal p -> equal q
    for v in np.unique(ps):
        assert np.ptp(qs[ps == v]) <= EPS


@FAST
@given(
    st.lists(P_OR_NAN, min_size=1, max_size=60).flatmap(
        lambda p: st.tuples(st.just(p), _perm(len(p)))
    )
)
def test_bh_invariant_to_input_order(args: tuple[list[float], list[int]]) -> None:
    p, perm = args
    pa = np.asarray(p, dtype=float)
    idx = np.asarray(perm)
    np.testing.assert_allclose(bh(pa[idx]), bh(pa)[idx], rtol=1e-12, atol=0, equal_nan=True)


# ---------------------------------------------------------------------------- monotone_q


@FAST
@given(st.lists(st.floats(0.0, 5.0), max_size=40))
def test_monotone_q_nondecreasing_and_bounded(fdr: list[float]) -> None:
    f = np.asarray(fdr, dtype=float)
    q = monotone_q(f)
    assert q.shape == f.shape
    assert np.all(np.diff(q) >= -EPS)
    assert np.all((q >= 0.0) & (q <= 1.0))
    assert np.all(q <= np.clip(f, 0.0, 1.0) + EPS)


# ---------------------------------------------------------------------------- tda_qvalues

TDA_ARGS = st.tuples(
    st.lists(SCORE_OR_NAN, max_size=40),
    st.lists(SCORE_OR_NAN, max_size=40),
    st.sampled_from([1.0, 0.5, 0.123]),
    st.booleans(),
)


@FAST
@given(TDA_ARGS)
def test_tda_q_nonincreasing_in_score(args: tuple[list[float], list[float], float, bool]) -> None:
    t, d, pi, plus_one = args
    ta = np.asarray(t, dtype=float)
    q = tda_qvalues(ta, d, pi=pi, plus_one=plus_one)
    assert q.shape == ta.shape
    np.testing.assert_array_equal(np.isnan(q), np.isnan(ta))
    ok = ~np.isnan(ta)
    s, qs = ta[ok], q[ok]
    assert np.all((qs >= 0.0) & (qs <= 1.0))
    order = np.argsort(-s, kind="stable")  # best first
    assert np.all(np.diff(qs[order]) >= -EPS)
    for v in np.unique(s):
        assert np.ptp(qs[s == v]) <= EPS


@FAST
@given(TDA_ARGS, SCORE)
def test_tda_adding_a_decoy_never_lowers_q(
    args: tuple[list[float], list[float], float, bool], extra: float
) -> None:
    t, d, pi, plus_one = args
    before = tda_qvalues(t, d, pi=pi, plus_one=plus_one)
    after = tda_qvalues(t, [*d, extra], pi=pi, plus_one=plus_one)
    np.testing.assert_array_equal(np.isnan(after), np.isnan(before))
    ok = ~np.isnan(before)
    assert np.all(after[ok] >= before[ok] - EPS)


@FAST
@given(
    st.lists(SCORE_OR_NAN, min_size=1, max_size=40).flatmap(
        lambda t: st.tuples(st.just(t), _perm(len(t)))
    ),
    st.lists(SCORE_OR_NAN, max_size=40).flatmap(lambda d: st.tuples(st.just(d), _perm(len(d)))),
)
def test_tda_invariant_to_input_order(
    targets: tuple[list[float], list[int]], decoys: tuple[list[float], list[int]]
) -> None:
    t, tp = targets
    d, dp = decoys
    ta, da = np.asarray(t, dtype=float), np.asarray(d, dtype=float)
    ti = np.asarray(tp, dtype=int)
    di = np.asarray(dp, dtype=int)
    q = tda_qvalues(ta, da)
    q2 = tda_qvalues(ta[ti], da[di] if da.size else da)
    np.testing.assert_allclose(q2, q[ti], rtol=1e-12, atol=0, equal_nan=True)


# ---------------------------------------------------------------------------- picked_fdr


@FAST
@given(st.lists(SCORE_OR_NAN, max_size=40))
def test_picked_without_decoys_equals_tda_without_decoys(t: list[float]) -> None:
    ta = np.asarray(t, dtype=float)
    res = picked_fdr(ta, np.full(ta.shape, np.nan))
    expected = tda_qvalues(ta, [])
    np.testing.assert_allclose(res.q, expected, rtol=1e-12, atol=0, equal_nan=True)
    assert not res.decoy_wins.any()
    assert res.n_decoy_wins == 0
    assert res.n_target_wins == int(np.sum(~np.isnan(ta)))


@FAST
@given(st.lists(st.tuples(SCORE_OR_NAN, SCORE_OR_NAN), max_size=40))
def test_picked_winner_counts_and_nan_pattern(pairs: list[tuple[float, float]]) -> None:
    t = np.asarray([a for a, _ in pairs], dtype=float)
    d = np.asarray([b for _, b in pairs], dtype=float)
    res = picked_fdr(t, d)
    has_t, has_d = ~np.isnan(t), ~np.isnan(d)
    expect_decoy = has_d & (~has_t | (d >= np.where(has_t, t, -np.inf)))
    np.testing.assert_array_equal(res.decoy_wins, expect_decoy)
    target_wins = has_t & ~expect_decoy
    np.testing.assert_array_equal(~np.isnan(res.q), target_wins)
    assert res.n_target_wins == int(target_wins.sum())
    assert res.n_decoy_wins == int(expect_decoy.sum())
    assert np.all((res.q[target_wins] >= 0.0) & (res.q[target_wins] <= 1.0))
