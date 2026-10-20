"""Moran's I, join counts, label and TU-block permutations, ``lag_test``, ``choose_block_len``
(issue M5-02; GIAP_SPEC.md §7.3.5(c), §7.3.5(d) step 5, §7.3.6, §8.3, Appendix B.4).

Spec examples 1-10 and 13 (docs/specs/M5-02-stats-moran-joincount.md), the edge cases, and the
lead decisions 4-8:

- decision 6: a constant null gives ``p_value`` NaN even when the observation differs (``perm_p``
  alone would give 1/(B + 1));
- decision 7 as amended by decision 17 (ADR-0020 Amendment 1): ``join_count(d, lag, circular,
  p=None)`` and ``join_count_matrix(D, lags, circular, P=None)`` return the Pearson residual sum
  Σ r_i r_j, r = (d - p)/sqrt(p(1 - p)) with p unclipped, over the valid lag-k pairs; genes with
  p in {0, 1} or NaN (and NaN entries of a float D) are missing and join no pair. Over the
  valid lag-k pairs when p is given; p moves with d under both nulls; all p equal and all d equal
  give a constant null, hence NaN.
- ``choose_block_len`` returns the spec's l-hat only (examples 11-13, decision 10);
- decision 10: ``block_len_prime(l_hat, T)`` = min(20, 2 l-hat + 1, T) (ADR-0021), shared by M5-06
  and M5-08;
- decision 11: ``lag_test`` with NaN obs: ``n_valid`` = count of finite null values, every other
  field except ``value`` NaN;
- decision 12: the < 5 s performance target is a benchmark (benchmarks/), not a test here.

Examples 11 and 12 and decision 8 are statistical claims, in
tests/calibration/test_engine_calibration.py.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest

from giap.stats.autocorr import (
    LagTest,
    block_len_prime,
    choose_block_len,
    join_count,
    join_count_matrix,
    label_permutations,
    lag_pairs,
    lag_test,
    max_lag,
    morans_i,
    morans_i_matrix,
    tu_block_permutations,
    tu_blocks,
)
from giap.stats.permutation import perm_p

NAN = float("nan")


def _rng(seed: int = 0) -> np.random.Generator:
    return np.random.default_rng(seed)


def _isnan(v: float) -> bool:
    return isinstance(v, float) and math.isnan(v)


# ----------------------------------------------------------------------------- lags


@pytest.mark.parametrize(
    ("n", "circular", "expected"),
    [
        (0, True, 0),
        (1, True, 0),
        (2, True, 0),
        (3, True, 1),
        (6, True, 2),
        (7, True, 3),
        (10, True, 4),
        (0, False, 0),
        (1, False, 0),
        (2, False, 1),
        (6, False, 5),
    ],
)
def test_max_lag(n: int, circular: bool, expected: int) -> None:
    assert max_lag(n, circular) == expected


def _pairs(n: int, lag: int, circular: bool) -> set[tuple[int, int]]:
    i, j = lag_pairs(n, lag, circular)
    assert i.dtype == np.int64 and j.dtype == np.int64
    assert i.shape == j.shape
    return set(zip(i.tolist(), j.tolist(), strict=True))


def test_lag_pairs_circular() -> None:
    assert _pairs(5, 2, True) == {(0, 2), (1, 3), (2, 4), (3, 0), (4, 1)}
    assert _pairs(5, 1, True) == {(i, (i + 1) % 5) for i in range(5)}


def test_lag_pairs_linear() -> None:
    assert _pairs(5, 2, False) == {(0, 2), (1, 3), (2, 4)}
    assert _pairs(5, 4, False) == {(0, 4)}


@pytest.mark.parametrize(
    ("n", "lag", "circular"),
    [(6, 3, True), (6, 0, True), (6, -1, True), (5, 5, False), (5, 0, False), (1, 1, False),
     (0, 1, True), (2, 1, True)],
)  # fmt: skip
def test_lag_pairs_empty_for_invalid_lag(n: int, lag: int, circular: bool) -> None:
    assert _pairs(n, lag, circular) == set()


# ----------------------------------------------------------------------------- Moran's I


def test_example1_alternating_series_circular() -> None:
    """B.4 / §7.3.6: x = +1, -1, ... (n = 10), circular: I_1 = -1, I_2 = +1."""
    x = np.array([1.0, -1.0] * 5)
    assert morans_i(x, 1, True) == pytest.approx(-1.0, abs=1e-12)
    assert morans_i(x, 2, True) == pytest.approx(1.0, abs=1e-12)


def test_example2_linear_and_circular() -> None:
    x = [1.0, 2.0, 3.0, 4.0]
    assert morans_i(x, 1, False) == pytest.approx(1.0 / 3.0, abs=1e-9)  # W = 3
    assert morans_i(x, 1, True) == pytest.approx(-0.2, abs=1e-9)  # W = 4


def test_example3_missing_values() -> None:
    x = [1.0, NAN, 3.0, 1.0, 3.0]
    assert morans_i(x, 1, False) == pytest.approx(-1.0, abs=1e-9)  # W = 2
    assert morans_i(x, 2, False) == pytest.approx(0.0, abs=1e-9)  # W = 2
    assert morans_i(x, 2, True) == pytest.approx(1.0 / 3.0, abs=1e-9)  # W = 3


def test_example4_no_variance_is_nan() -> None:
    assert _isnan(morans_i([2.0] * 5, 1, False))
    assert _isnan(morans_i([2.0] * 5, 1, True))


def test_no_variance_relative_threshold_is_nan() -> None:
    """Σ z² ≤ 1e-24·max(1, max x²): a large constant with rounding noise has no variance."""
    x = np.full(8, 1e6) + np.array([0, 1e-16, 0, 0, 0, 0, 0, 0]) * 1e6
    assert _isnan(morans_i(x, 1, True))


def test_example5_no_pair_with_values_is_nan() -> None:
    assert _isnan(morans_i([1.0, NAN, NAN, 2.0], 1, False))  # W = 0


def test_example6_invalid_circular_lag() -> None:
    x = np.arange(6, dtype=float)
    d = np.array([True, False, True, True, False, True])
    assert max_lag(6, True) == 2
    assert _isnan(morans_i(x, 3, True))
    assert _isnan(join_count(d, 3, True))
    i, j = lag_pairs(6, 3, True)
    assert i.size == 0 and j.size == 0


@pytest.mark.parametrize("x", [[], [1.0], [1.0, NAN, NAN], [NAN, NAN, NAN, NAN]])
def test_fewer_than_two_values_is_nan(x: list[float]) -> None:
    for circular in (True, False):
        assert _isnan(morans_i(x, 1, circular))


def test_morans_i_returns_float() -> None:
    v = morans_i([1.0, 2.0, 3.0, 4.0], 1, False)
    assert isinstance(v, float)


# ----------------------------------------------------------------------------- join counts


def test_example7_join_counts() -> None:
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    assert join_count(d, 1, True) == 2.0
    assert join_count(d, 1, False) == 2.0
    assert join_count(d, 2, True) == 2.0
    assert join_count(d, 2, False) == 1.0
    assert isinstance(join_count(d, 1, True), float)


def test_join_count_zero_is_a_value() -> None:
    d = np.zeros(6, dtype=bool)
    assert join_count(d, 1, True) == 0.0
    assert join_count(d, 2, False) == 0.0


def test_join_count_invalid_lag_is_nan() -> None:
    d = np.ones(5, dtype=bool)
    assert _isnan(join_count(d, 0, True))
    assert _isnan(join_count(d, 3, True))  # max_lag(5, circular) = 2
    assert _isnan(join_count(d, 5, False))


def test_residual_join_count_hand_values() -> None:
    """Decision 17 / ADR-0020 Amendment 1: Σ r_i r_j over the valid lag-k pairs,
    r = (d - p)/sqrt(p(1 - p)) = [1/3, 2, -1/3, 1, sqrt(3/7), -sqrt(3/7)] here."""
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    p = np.array([0.9, 0.2, 0.1, 0.5, 0.7, 0.3])
    q = math.sqrt(3.0 / 7.0)
    assert join_count(d, 1, False, p=p) == pytest.approx(-1 / 3 + q - 3 / 7, abs=1e-12)
    assert join_count(d, 1, False, p=p) == pytest.approx(-0.10725109119678478, abs=1e-12)
    assert join_count(d, 1, True, p=p) == pytest.approx(-0.3254689814327771, abs=1e-12)
    assert join_count(d, 2, False, p=p) == pytest.approx(1.0160173279449194, abs=1e-12)
    assert join_count(d, 2, True, p=p) == pytest.approx(-0.07507212323504242, abs=1e-12)
    half = np.full(6, 0.5)  # r = +-1
    assert join_count(d, 1, True, p=half) == pytest.approx(-2.0, abs=1e-12)
    assert join_count(d, 1, False, p=half) == pytest.approx(-1.0, abs=1e-12)
    assert _isnan(join_count(d, 3, True, p=p))


def test_residual_is_not_the_raw_residual() -> None:
    """The amended statistic differs from the superseded Σ (d - p)(d - p) (0.01 here)."""
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    p = np.array([0.9, 0.2, 0.1, 0.5, 0.7, 0.3])
    assert abs(join_count(d, 1, False, p=p) - 0.01) > 0.1


@pytest.mark.parametrize("bad", [0.0, 1.0, NAN])
def test_residual_genes_with_p_zero_one_or_nan_are_missing(bad: float) -> None:
    """Genes 1 and 3 get p in {0, 1} or NaN: no residual, they join no pair (p unclipped).
    r = [1/3, -, -1/3, -, sqrt(3/7), -sqrt(3/7)]."""
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    p = np.array([0.9, bad, 0.1, bad, 0.7, 0.3])
    if bad == 1.0:
        p[1] = 0.0  # mix p = 0 and p = 1
    q = math.sqrt(3.0 / 7.0)
    assert join_count(d, 1, False, p=p) == pytest.approx(-3 / 7, abs=1e-12)
    assert join_count(d, 1, True, p=p) == pytest.approx(-3 / 7 - q / 3, abs=1e-12)
    assert join_count(d, 2, False, p=p) == pytest.approx(-1 / 9 - q / 3, abs=1e-12)
    assert join_count(d, 2, True, p=p) == pytest.approx(-1 / 9, abs=1e-12)
    for lag, circular in [(1, True), (1, False), (2, True), (2, False)]:
        v = join_count(d, lag, circular, p=p)
        assert math.isfinite(v)
        P = np.vstack([p, p])
        np.testing.assert_allclose(
            join_count_matrix(np.vstack([d, d]), [lag], circular, P=P)[:, 0], [v, v], atol=1e-12
        )


def test_residual_p_unclipped_near_bounds() -> None:
    """p close to 0 or 1 is used as is (no floor): r = (d - p)/sqrt(p(1 - p))."""
    d = np.array([True, True, False, False])
    p = np.array([1e-6, 0.5, 1 - 1e-6, 0.5])
    r = (d - p) / np.sqrt(p * (1 - p))
    assert join_count(d, 1, False, p=p) == pytest.approx(float((r[:-1] * r[1:]).sum()), rel=1e-12)


def test_residual_all_genes_missing_is_zero() -> None:
    """Every p in {0, 1}: no pair joins; a join count is a sum, so 0 (it exists at valid lags)."""
    d = np.array([1, 0, 1, 1, 0], dtype=bool)
    p = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
    assert join_count(d, 1, True, p=p) == 0.0
    assert join_count(d, 2, False, p=p) == 0.0
    assert _isnan(join_count(d, 3, True, p=p))


def test_nan_in_float_d_with_p_is_missing() -> None:
    """Code re-review: a NaN in a float D together with P is a missing gene (joins no pair)."""
    d = np.array([1.0, NAN, 0.0, 1.0, 1.0, 0.0])
    p = np.array([0.9, 0.2, 0.1, 0.5, 0.7, 0.3])
    r = np.array([1 / 3, 0.0, -1 / 3, 1.0, math.sqrt(3 / 7), -math.sqrt(3 / 7)])
    for k in (1, 2):
        assert join_count(d, k, False, p=p) == pytest.approx(
            float((r[:-k] * r[k:]).sum()), abs=1e-12
        )
        assert join_count(d, k, True, p=p) == pytest.approx(
            float((r * np.roll(r, -k)).sum()), abs=1e-12
        )
    got = join_count_matrix(np.vstack([d, d[::-1]]), [1, 2], True, P=np.vstack([p, p[::-1]]))
    assert np.isfinite(got).all()
    np.testing.assert_allclose(got[0], [join_count(d, 1, True, p=p), join_count(d, 2, True, p=p)])


def test_join_count_matrix_residual_p_moves_with_d() -> None:
    """P = p[orders] (the caller permutes p with d); row b = scalar residual of (D[b], P[b])."""
    rng = _rng(3)
    n = 30
    p = rng.uniform(0.05, 0.95, size=n)
    d = rng.uniform(size=n) < p
    orders = label_permutations(n, 25, rng)
    D, P = d[orders], p[orders]
    lags = [1, 2, 5, 14, 15, 40]
    for circular in (True, False):
        got = join_count_matrix(D, lags, circular, P=P)
        assert got.shape == (25, len(lags))
        want = np.array([[join_count(D[b], k, circular, p=P[b]) for k in lags] for b in range(25)])
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-12, equal_nan=True)
        raw = join_count_matrix(D, lags, circular)
        want_raw = np.array([[join_count(D[b], k, circular) for k in lags] for b in range(25)])
        np.testing.assert_array_equal(raw, want_raw)


@pytest.mark.parametrize("all_detected", [True, False])
def test_all_equal_d_gives_constant_null_and_nan_p(all_detected: bool) -> None:
    """Edge case + decision 7: all (or no) genes detected -> BB constant under both nulls -> p NaN;
    with every p equal the residual statistic is constant too."""
    n = 40
    rng = _rng(4)
    d = np.full(n, all_detected)
    p = np.full(n, 0.3)
    labels = [f"TU{i // 3}" for i in range(n)]
    for circular in (True, False):
        nulls = {
            "label": label_permutations(n, 99, rng),
            "tu": tu_block_permutations(tu_blocks(labels, circular), 99, rng),
        }
        for orders in nulls.values():
            for P in (None, p[orders]):
                obs = join_count(d, 1, circular, p=None if P is None else p)
                null = join_count_matrix(d[orders], [1], circular, P=P)[:, 0]
                res = lag_test(obs, null)
                assert _isnan(res.p_value)
                assert res.value == obs
                assert res.n_valid == 99


def test_residual_with_varying_p_is_not_degenerate_when_d_constant() -> None:
    """All d equal but p varying: the Pearson residual varies, so the null is not constant."""
    n = 40
    rng = _rng(5)
    d = np.ones(n, dtype=bool)
    p = rng.uniform(0.1, 0.9, size=n)
    orders = label_permutations(n, 99, rng)
    obs = join_count(d, 1, True, p=p)
    null = join_count_matrix(d[orders], [1], True, P=p[orders])[:, 0]
    assert np.ptp(null) > 0
    res = lag_test(obs, null)
    assert 0 < res.p_value <= 1


# ----------------------------------------------------------------------------- matrices


def test_morans_i_matrix_equals_scalar_rows() -> None:
    rng = _rng(6)
    n = 25
    x = rng.normal(size=n)
    x[[3, 11, 12]] = NAN
    orders = label_permutations(n, 30, rng)
    X = x[orders]
    lags = [1, 2, 3, 12, 13, 24, 25]
    for circular in (True, False):
        got = morans_i_matrix(X, lags, circular)
        assert got.shape == (30, len(lags))
        assert got.dtype == np.float64
        want = np.array([[morans_i(X[b], k, circular) for k in lags] for b in range(30)])
        np.testing.assert_allclose(got, want, rtol=1e-10, atol=1e-12, equal_nan=True)


def test_matrices_large_b_independent_of_chunking() -> None:
    """Many rows (an implementation chunks): every row still equals the scalar function."""
    rng = _rng(7)
    n, B = 60, 2_500
    x = rng.normal(size=n)
    d = rng.uniform(size=n) < 0.4
    orders = label_permutations(n, B, rng)
    lags = list(range(1, 21))
    M = morans_i_matrix(x[orders], lags, True)
    J = join_count_matrix(d[orders], lags, False)
    assert M.shape == J.shape == (B, 20)
    for b in (0, 1, 999, 1000, 1001, 1777, B - 1):
        np.testing.assert_allclose(
            M[b], [morans_i(x[orders[b]], k, True) for k in lags], rtol=1e-10, atol=1e-12
        )
        np.testing.assert_array_equal(J[b], [join_count(d[orders[b]], k, False) for k in lags])


# ----------------------------------------------------------------------------- permutations


def test_label_permutations_rows_are_permutations() -> None:
    out = label_permutations(12, 50, _rng(8))
    assert out.shape == (50, 12) and out.dtype == np.int64
    for row in out:
        assert sorted(row.tolist()) == list(range(12))


def test_label_permutations_deterministic_in_rng() -> None:
    a = label_permutations(20, 10, _rng(9))
    b = label_permutations(20, 10, _rng(9))
    np.testing.assert_array_equal(a, b)


def test_label_permutations_uniform() -> None:
    from scipy.stats import chisquare

    out = label_permutations(3, 6_000, _rng(10))
    keys, counts = np.unique(out, axis=0, return_counts=True)
    assert len(keys) == 6
    assert chisquare(counts).pvalue > 1e-4


def test_label_permutations_movable_keeps_fixed_positions() -> None:
    movable = np.array([True, False, True, False, True, True, False])
    out = label_permutations(7, 400, _rng(11), movable=movable)
    assert out.shape == (400, 7)
    fixed = np.flatnonzero(~movable)
    mov = np.flatnonzero(movable)
    np.testing.assert_array_equal(out[:, fixed], np.broadcast_to(fixed, (400, fixed.size)))
    for row in out:
        assert sorted(row[mov].tolist()) == mov.tolist()
    assert len(np.unique(out[:, mov], axis=0)) == 24  # all 4! arrangements appear


@pytest.mark.parametrize("n", [0, 1])
def test_permutations_of_length_zero_or_one_are_identity(n: int) -> None:
    lab = label_permutations(n, 5, _rng(12))
    assert lab.shape == (5, n)
    np.testing.assert_array_equal(lab, np.zeros((5, n), dtype=np.int64))
    blocks = tu_blocks(["A"] * n, True)
    tub = tu_block_permutations(blocks, 5, _rng(12))
    assert tub.shape == (5, n)
    np.testing.assert_array_equal(tub, np.zeros((5, n), dtype=np.int64))
    assert _isnan(morans_i(np.ones(n), 1, True))
    assert _isnan(join_count(np.ones(n, dtype=bool), 1, True))


def test_example8_tu_blocks() -> None:
    labels = ["A", "A", "B", "C", "C", "A"]
    np.testing.assert_array_equal(tu_blocks(labels, False), [0, 0, 1, 2, 2, 3])
    np.testing.assert_array_equal(tu_blocks(labels, True), [0, 0, 1, 2, 2, 0])
    assert tu_blocks(labels, True).dtype == np.int64


def test_tu_blocks_interrupted_tu_forms_several_blocks() -> None:
    np.testing.assert_array_equal(tu_blocks(["A", "B", "A", "C"], False), [0, 1, 2, 3])
    np.testing.assert_array_equal(tu_blocks(["A", "B", "A", "C"], True), [0, 1, 2, 3])
    np.testing.assert_array_equal(tu_blocks(["A", "B", "A"], True), [0, 1, 0])


def test_tu_blocks_single_tu() -> None:
    for circular in (True, False):
        np.testing.assert_array_equal(tu_blocks(["A"] * 5, circular), [0] * 5)
    assert tu_blocks([], True).shape == (0,)


def test_tu_block_permutations_exact_rows_with_wrapped_block() -> None:
    """Blocks in random order, internal order kept; the origin-crossing block is emitted as its
    end-part followed by its start-part: block 0 = positions [5, 0, 1]."""
    blocks = tu_blocks(["A", "A", "B", "C", "C", "A"], True)
    out = tu_block_permutations(blocks, 600, _rng(13))
    assert out.shape == (600, 6) and out.dtype == np.int64
    parts = {0: (5, 0, 1), 1: (2,), 2: (3, 4)}
    allowed = {
        tuple(i for b in order for i in parts[b])
        for order in [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    }
    rows = {tuple(r) for r in out.tolist()}
    assert rows == allowed


def test_tu_block_permutations_linear_exact_rows() -> None:
    blocks = tu_blocks(["A", "A", "B", "C", "C", "A"], False)
    out = tu_block_permutations(blocks, 2_400, _rng(14))
    parts = [(0, 1), (2,), (3, 4), (5,)]
    rows = {tuple(r) for r in out.tolist()}
    assert len(rows) == 24
    for r in rows:
        seen = []
        k = 0
        while k < len(r):
            part = next(p for p in parts if p[0] == r[k])
            assert r[k : k + len(part)] == part
            seen.append(part)
            k += len(part)
        assert sorted(seen) == sorted(parts)


def test_tu_block_permutations_uniform_over_block_orders() -> None:
    from scipy.stats import chisquare

    blocks = tu_blocks(["A", "B", "B", "C"], False)
    out = tu_block_permutations(blocks, 6_000, _rng(15))
    keys, counts = np.unique(out, axis=0, return_counts=True)
    assert len(keys) == 6
    assert chisquare(counts).pvalue > 1e-4


def test_tu_block_permutations_two_wrapped_parts_example() -> None:
    blocks = tu_blocks(["A", "B", "A"], True)  # [0, 1, 0]: block 0 = [2, 0]
    out = tu_block_permutations(blocks, 200, _rng(16))
    assert {tuple(r) for r in out.tolist()} == {(2, 0, 1), (1, 2, 0)}


def test_one_tu_replicon_tu_null_is_identity_and_degenerate() -> None:
    """Edge case: one TU -> one block -> identity rows -> lag_test degenerate -> p NaN."""
    rng = _rng(17)
    x = rng.normal(size=9)
    blocks = tu_blocks(["T1"] * 9, True)
    out = tu_block_permutations(blocks, 50, rng)
    np.testing.assert_array_equal(out, np.broadcast_to(np.arange(9), (50, 9)))
    null = morans_i_matrix(x[out], [1], True)[:, 0]
    res = lag_test(morans_i(x, 1, True), null)
    assert _isnan(res.p_value)


def test_tu_block_permutations_deterministic_in_rng() -> None:
    blocks = tu_blocks([f"T{i // 2}" for i in range(20)], True)
    np.testing.assert_array_equal(
        tu_block_permutations(blocks, 7, _rng(18)), tu_block_permutations(blocks, 7, _rng(18))
    )


# ----------------------------------------------------------------------------- lag_test


def test_example9_lag_test() -> None:
    res = lag_test(3.0, [0.0, 0.0, 0.0, 1.0])
    assert isinstance(res, LagTest)
    assert res.value == 3.0
    assert res.p_value == pytest.approx(0.2, abs=1e-12)
    assert res.null_mean == pytest.approx(0.25, abs=1e-12)
    assert res.null_lo == pytest.approx(0.0, abs=1e-12)
    assert res.null_hi == pytest.approx(0.925, abs=1e-12)
    assert res.n_valid == 4


def test_example10_degenerate_null_equal_to_obs() -> None:
    res = lag_test(0.0, [0.0, 0.0, 0.0])
    assert _isnan(res.p_value)
    assert res.null_mean == 0.0
    assert res.n_valid == 3


def test_decision6_constant_null_differing_obs_is_nan() -> None:
    """Decision 6: the NA rule tests "null constant"; perm_p alone gives 1/(B + 1) here."""
    assert perm_p(1.0, [0.0, 0.0, 0.0]) == pytest.approx(0.25)
    res = lag_test(1.0, [0.0, 0.0, 0.0])
    assert _isnan(res.p_value)
    assert res.value == 1.0
    assert res.null_mean == 0.0
    assert res.null_lo == 0.0 and res.null_hi == 0.0
    assert res.n_valid == 3


def test_lag_test_drops_nan_null_entries() -> None:
    res = lag_test(3.0, [0.0, NAN, 0.0, 0.0, NAN, 1.0])
    assert res.n_valid == 4
    assert res.p_value == pytest.approx(0.2, abs=1e-12)
    assert res.null_mean == pytest.approx(0.25, abs=1e-12)
    assert res.null_hi == pytest.approx(0.925, abs=1e-12)


def test_lag_test_obs_nan() -> None:
    """Decision 11: n_valid counts the finite null values; every other field but value is NaN."""
    for null, n_valid in (([0.0, 1.0, 2.0], 3), ([0.0, NAN, 2.0, NAN], 2), ([NAN], 0), ([], 0)):
        res = lag_test(NAN, null)
        assert _isnan(res.value)
        assert res.n_valid == n_valid
        assert isinstance(res.n_valid, int)
        for v in (res.p_value, res.null_mean, res.null_lo, res.null_hi):
            assert _isnan(v)


def test_lag_test_no_valid_null() -> None:
    for null in ([], [NAN, NAN]):
        res = lag_test(2.0, null)
        assert res.value == 2.0
        assert res.n_valid == 0
        for v in (res.p_value, res.null_mean, res.null_lo, res.null_hi):
            assert _isnan(v)


def test_lag_test_matches_perm_p_and_quantiles() -> None:
    rng = _rng(19)
    null = rng.normal(size=199)
    for obs in (-3.0, -0.1, 0.0, 0.4, 2.5):
        res = lag_test(obs, null)
        assert res.p_value == perm_p(obs, null, "two-sided")
        assert res.null_mean == pytest.approx(float(null.mean()), abs=1e-12)
        assert res.null_lo == pytest.approx(float(np.quantile(null, 0.025)), abs=1e-12)
        assert res.null_hi == pytest.approx(float(np.quantile(null, 0.975)), abs=1e-12)
        assert res.n_valid == 199


def test_lag_test_is_frozen() -> None:
    res = lag_test(3.0, [0.0, 0.0, 0.0, 1.0])
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.p_value = 0.5  # type: ignore[misc]
    assert {f.name for f in dataclasses.fields(LagTest)} == {
        "value", "null_mean", "null_lo", "null_hi", "p_value", "n_valid"
    }  # fmt: skip


# ----------------------------------------------------------------------------- choose_block_len


def test_example13_constant_z_gives_one() -> None:
    for circular in (True, False):
        assert choose_block_len(np.full(50, 1.5), circular, 99, _rng(20)) == 1


@pytest.mark.parametrize("n", [0, 1, 2])
def test_fewer_than_three_values_gives_one(n: int) -> None:
    for circular in (True, False):
        assert choose_block_len(np.arange(n, dtype=float), circular, 99, _rng(21)) == 1


def test_all_lags_significant_bounded_by_max_len() -> None:
    """Alternating series, n = 200 circular: |I_k| = 1 at every lag -> l = min(20, 100, 200)."""
    z = np.array([1.0, -1.0] * 100)
    assert choose_block_len(z, True, 199, _rng(22)) == 20
    assert choose_block_len(z, True, 199, _rng(22), max_len=3) == 3
    z60 = np.array([1.0, -1.0] * 30)
    assert choose_block_len(z60, False, 199, _rng(22)) == 20


def test_all_lags_significant_bounded_by_max_lag_plus_one() -> None:
    """n = 30 circular: max_lag = 14, every lag has |I_k| = 1 -> l = min(20, 15, 30) = 15."""
    z = np.array([1.0, -1.0] * 15)
    assert choose_block_len(z, True, 199, _rng(23)) == 15


def test_first_non_significant_lag_is_returned() -> None:
    """MA(1) z_i = e_i + 0.2 e_(i-1), n = 4,000: I_1 ≈ 0.19 (≈ 12 SE, significant), true lag-2
    autocorrelation 0 -> l = 2 unless lag 2 is a chance rejection (Bartlett: rate ≈ 0.06)."""
    hits = 0
    for s in range(20):
        rng = _rng(100 + s)
        e = rng.normal(size=4_001)
        z = e[1:] + 0.2 * e[:-1]
        hits += choose_block_len(z, True, 199, rng) == 2
    assert hits >= 15


def test_alpha_is_the_threshold() -> None:
    """With B = 19 the smallest p is 0.05: p > alpha = 0.01 at lag 1 for any series -> l = 1;
    at alpha = 0.05 the alternating series is significant at every lag."""
    z = np.array([1.0, -1.0] * 100)
    assert choose_block_len(z, True, 19, _rng(24), alpha=0.01) == 1
    assert choose_block_len(z, True, 19, _rng(24), alpha=0.05) == 20


def test_choose_block_len_deterministic_and_returns_int() -> None:
    z = np.repeat(_rng(25).normal(size=300), 3)
    a = choose_block_len(z, True, 99, _rng(26))
    b = choose_block_len(z, True, 99, _rng(26))
    assert a == b
    assert isinstance(a, int)


def test_choose_block_len_draws_fixed_by_b_and_n_only() -> None:
    """One permutation matrix for all lags: the generator ends in the same state whatever z is."""
    n, B = 120, 49
    states = []
    for z in (
        np.zeros(n),
        _rng(27).normal(size=n),
        np.array([1.0, -1.0] * (n // 2)),
        np.repeat(_rng(28).normal(size=n // 4), 4),
    ):
        rng = _rng(29)
        choose_block_len(z, True, B, rng)
        states.append(float(rng.random()))
    assert len(set(states)) == 1


# ----------------------------------------------------------------------------- block_len_prime


@pytest.mark.parametrize(
    ("l_hat", "T", "expected"),
    [
        (1, 1_000, 3),
        (1, 3, 3),
        (1, 2, 2),
        (8, 1_000, 17),
        (8, 17, 17),
        (8, 12, 12),
        (10, 1_000, 20),  # 2*10 + 1 = 21 capped at 20
        (10, 25, 20),
        (10, 15, 15),
        (20, 1_000, 20),
        (20, 20, 20),
        (20, 7, 7),
        (1, 1, 1),
        (8, 1, 1),
        (20, 1, 1),
    ],
)
def test_block_len_prime(l_hat: int, T: int, expected: int) -> None:
    """Decision 10 / ADR-0021: l' = min(20, 2 l-hat + 1, T)."""
    got = block_len_prime(l_hat, T)
    assert got == expected
    assert isinstance(got, int)


# ----------------------------------------------------------------------------- review of #116


def test_two_tu_blocks_circular_rotation_null_gives_nan_p() -> None:
    """Decision 14: on a circular replicon with 2 TU blocks the only block orders are rotations, so
    Moran's I and the (residual) join count are constant under the TU null up to rounding
    (spread ~1 ulp); "constant" is spread <= 1e-12 * max(1, max |null|) -> p NaN."""
    for s in range(60):
        rng = _rng(1_000 + s)
        n = int(rng.integers(12, 200))
        m = int(rng.integers(1, n))
        labels = ["A"] * m + ["B"] * (n - m)
        orders = tu_block_permutations(tu_blocks(labels, True), 40, rng)
        x = rng.normal(size=n) * 3.7 + 11.0
        p = rng.uniform(0.05, 0.95, size=n)
        d = rng.uniform(size=n) < p
        for k in (1, 2, 5):
            if k > max_lag(n, True):
                continue
            cases = {
                "moran": (morans_i(x, k, True), morans_i_matrix(x[orders], [k], True)[:, 0]),
                "bb": (join_count(d, k, True), join_count_matrix(d[orders], [k], True)[:, 0]),
                "bb_resid": (
                    join_count(d, k, True, p=p),
                    join_count_matrix(d[orders], [k], True, P=p[orders])[:, 0],
                ),
            }
            for name, (obs, null) in cases.items():
                assert np.ptp(null) <= 1e-12 * max(1.0, float(np.max(np.abs(null)))), name
                res = lag_test(obs, null)
                assert _isnan(res.p_value), (s, n, m, k, name, float(np.ptp(null)))


def test_lag_test_near_constant_null_is_degenerate() -> None:
    """Decision 14 at the tie scale: a 1-ulp spread is constant, a 1e-9 relative spread is not."""
    base = 0.123456789
    tiny = np.array([base, np.nextafter(base, 1.0), base, base])
    assert _isnan(lag_test(0.5, tiny).p_value)
    spread = np.array([base, base * (1 + 1e-9), base, base])
    assert not _isnan(lag_test(0.5, spread).p_value)


def test_lag_test_near_constant_null_large_magnitude() -> None:
    """Decision 14 tolerance scales with max |null|: at 1e6 a spread of 1e-7 (<= 1e-6) is constant,
    a spread of 1e-5 is not."""
    c = 1e6
    assert _isnan(lag_test(c + 1.0, [c, c + 1e-7, c, c]).p_value)
    assert _isnan(lag_test(-c - 1.0, [-c, -c - 1e-7, -c]).p_value)
    res = lag_test(c + 1.0, [c, c + 1e-5, c, c])
    assert not _isnan(res.p_value)
    assert res.p_value == perm_p(c + 1.0, [c, c + 1e-5, c, c], "two-sided")


@pytest.mark.parametrize(
    "null", [[0.0, 0.5, math.inf, 2.0], [-math.inf, 1.0, 3.0], [0.0, 0.0, math.inf]]
)
def test_lag_test_p_is_perm_p_with_inf_and_finite(null: list[float]) -> None:
    """Code re-review: a null mixing inf with finite values is not constant; p = perm_p."""
    for obs in (1.0, 10.0):
        res = lag_test(obs, null)
        assert res.n_valid == len(null)
        assert res.p_value == perm_p(obs, null, "two-sided")


def test_lag_test_inf_null_entries_count_in_n_valid() -> None:
    """Decision 13: n_valid counts the non-NaN null values; +-inf are valid (as in perm_p)."""
    inf = math.inf
    null = [0.0, inf, 1.0, NAN, -inf, 2.0]
    res = lag_test(1.0, null)
    assert res.n_valid == 5
    assert res.p_value == perm_p(1.0, null, "two-sided")
    res_nan_obs = lag_test(NAN, null)
    assert res_nan_obs.n_valid == 5
    assert _isnan(res_nan_obs.p_value)


def test_no_variance_threshold_boundary() -> None:
    """Σ z² <= 1e-24 · max(1, max x²): at x ≈ 1e6 the bound is ≈ 1e-12. One value off by δ gives
    Σ z² = δ² · 5/6: δ = 2e-6 (3.3e-12) has variance, δ = 5e-7 (2.1e-13) has none."""
    c = 1e6
    above = np.array([c, c, c, c, c, c + 2e-6])
    below = np.array([c, c, c, c, c, c + 5e-7])
    assert not _isnan(morans_i(above, 1, False))
    assert _isnan(morans_i(below, 1, False))
    small = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 2e-12])  # max(1, max x²) = 1: 3.3e-24 > 1e-24
    assert not _isnan(morans_i(small, 1, False))
    assert _isnan(morans_i(small * 0.5, 1, False))  # 8.3e-25 <= 1e-24


def test_label_permutations_two_movable_positions_swap() -> None:
    out = label_permutations(5, 200, _rng(30), movable=[False, True, False, False, True])
    assert {tuple(r) for r in out.tolist()} == {(0, 1, 2, 3, 4), (0, 4, 2, 3, 1)}
    out2 = label_permutations(2, 200, _rng(31))
    assert {tuple(r) for r in out2.tolist()} == {(0, 1), (1, 0)}


def _manual_block_len(z: np.ndarray, orders: np.ndarray, circular: bool) -> int:
    n = z.size
    top = min(20, max_lag(n, circular))
    lags = list(range(1, top + 1))
    obs = morans_i_matrix(z[None, :], lags, circular)[0]
    null = morans_i_matrix(z[orders], lags, circular)
    for c, k in enumerate(lags):
        p = lag_test(float(obs[c]), null[:, c]).p_value
        if not p <= 0.05:
            return k
    return min(20, max_lag(n, circular) + 1, n)


def test_choose_block_len_moves_nans_with_values() -> None:
    """Decision 16: choose_block_len permutes all n positions (NaNs move with their values, no
    ``movable``), one matrix ``label_permutations(n, B, rng)`` for every lag."""
    from scipy.signal import lfilter

    differs = 0
    for s in range(10):
        r = _rng(200 + s)
        n = 300
        z = np.asarray(lfilter([1.0], [1.0, -0.8], r.normal(size=n)))
        z[np.arange(n) % 3 == 1] = NAN
        want = _manual_block_len(z, label_permutations(n, 199, _rng(s)), True)
        assert choose_block_len(z, True, 199, _rng(s)) == want, s
        fixed = label_permutations(n, 199, _rng(s), movable=~np.isnan(z))
        differs += _manual_block_len(z, fixed, True) != want
    assert differs >= 1  # the test tells the two nulls apart


# ----------------------------------------------------------------------------- ValueErrors


@pytest.mark.parametrize("X", [np.zeros(5), np.zeros((2, 3, 4)), 1.0])
def test_morans_i_matrix_requires_2d(X: object) -> None:
    with pytest.raises(ValueError):
        morans_i_matrix(X, [1], True)  # type: ignore[arg-type]


@pytest.mark.parametrize("D", [np.zeros(5, dtype=bool), np.zeros((2, 3, 4), dtype=bool)])
def test_join_count_matrix_requires_2d(D: np.ndarray) -> None:
    with pytest.raises(ValueError):
        join_count_matrix(D, [1], True)


@pytest.mark.parametrize(
    "P", [np.zeros((1, 6)), np.zeros(6), np.zeros((4, 5)), np.zeros((4, 7)), np.zeros((4, 6, 1))]
)
def test_join_count_matrix_p_shape_must_match_d(P: np.ndarray) -> None:
    D = np.zeros((4, 6), dtype=bool)
    with pytest.raises(ValueError):
        join_count_matrix(D, [1], True, P=P)


@pytest.mark.parametrize(("n", "B"), [(-1, 5), (5, -1), (-2, -2)])
def test_label_permutations_rejects_negative_n_or_b(n: int, B: int) -> None:
    with pytest.raises(ValueError):
        label_permutations(n, B, _rng(40))


@pytest.mark.parametrize("length", [0, 4, 6, 10])
def test_label_permutations_rejects_wrong_movable_length(length: int) -> None:
    with pytest.raises(ValueError):
        label_permutations(5, 3, _rng(41), movable=np.ones(length, dtype=bool))


def test_tu_block_permutations_rejects_negative_b() -> None:
    with pytest.raises(ValueError):
        tu_block_permutations(np.array([0, 0, 1]), -1, _rng(42))


@pytest.mark.parametrize(
    "blocks",
    [
        [0, 1, 2, 1],  # one gap in a block other than the one holding position 0
        [0, 0, 1, 0, 2],  # block of position 0 with a gap that does not end at n - 1
        [0, 1, 0, 1],  # two blocks with gaps
        [0, 1, 0, 2, 0],  # block 0 with two gaps
        [1, 0, 1, 0],
    ],
)
def test_tu_block_permutations_rejects_non_contiguous_blocks(blocks: list[int]) -> None:
    with pytest.raises(ValueError):
        tu_block_permutations(np.array(blocks), 3, _rng(43))


@pytest.mark.parametrize("max_len", [0, -1])
def test_choose_block_len_rejects_max_len_below_one(max_len: int) -> None:
    with pytest.raises(ValueError):
        choose_block_len(_rng(44).normal(size=50), True, 19, _rng(45), max_len=max_len)


@pytest.mark.parametrize("B", [0, -1])
def test_choose_block_len_rejects_b_below_one(B: int) -> None:
    """Decision 16: B = 0 raises (it would silently give l-hat = 1)."""
    with pytest.raises(ValueError):
        choose_block_len(_rng(46).normal(size=50), True, B, _rng(47))


@pytest.mark.parametrize(("l_hat", "T"), [(0, 10), (-1, 10), (5, 0), (5, -3), (0, 0)])
def test_block_len_prime_rejects_values_below_one(l_hat: int, T: int) -> None:
    with pytest.raises(ValueError):
        block_len_prime(l_hat, T)


# ----------------------------------------------------------------------------- decision 21


@pytest.mark.parametrize("bad", [-0.1, 1.5, -math.inf, math.inf])
def test_join_count_rejects_p_outside_unit_interval(bad: float) -> None:
    """Decision 21: a non-NaN p outside [0, 1] is an upstream bug and raises, scalar and matrix."""
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    p = np.array([0.9, 0.2, bad, 0.5, 0.7, 0.3])
    for circular in (True, False):
        with pytest.raises(ValueError):
            join_count(d, 1, circular, p=p)
        with pytest.raises(ValueError):
            join_count_matrix(np.vstack([d, d]), [1, 2], circular, P=np.vstack([p, p[::-1]]))
        good = np.array([0.9, 0.2, 0.1, 0.5, 0.7, 0.3])
        with pytest.raises(ValueError):  # one bad row is enough
            join_count_matrix(np.vstack([d, d]), [1], circular, P=np.vstack([good, p]))
    with pytest.raises(ValueError):  # also at an invalid lag (the check is on the input)
        join_count(d, 3, True, p=p)


def test_join_count_p_zero_one_nan_do_not_raise() -> None:
    """Decision 21: p in {0, 1} and NaN stay missing (no error), scalar and matrix."""
    d = np.array([1, 1, 0, 1, 1, 0], dtype=bool)
    p = np.array([0.0, 1.0, math.nan, 0.5, 0.7, 0.3])
    r = np.array([0.0, 0.0, 0.0, 1.0, math.sqrt(3 / 7), -math.sqrt(3 / 7)])
    for circular in (True, False):
        v = join_count(d, 1, circular, p=p)
        want = float((r * np.roll(r, -1)).sum() if circular else (r[:-1] * r[1:]).sum())
        assert v == pytest.approx(want, abs=1e-12)
        got = join_count_matrix(np.vstack([d, d]), [1], circular, P=np.vstack([p, p]))
        np.testing.assert_allclose(got[:, 0], [v, v], atol=1e-12)
