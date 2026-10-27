"""LOESS detrending with circular padding (issue M5-01; GIPPA_SPEC.md §7.3.5(c), §7.3.5(d) step 3).

Lead decisions 4-7 (2026-10-20): k = 3 on equal spacing gives fit == y (documented degeneracy);
y / pos length mismatch -> ValueError; ``loess_matrix(pos, frac, *, circular, period=None,
min_points=3)`` is the n x n smoother with ``loess_detrend(y, ...).fit == loess_matrix(...) @ y``;
no calibration test (deterministic helpers, no statistical claim).

``loess_span_points(n, frac, min_points=3)`` = min(n, max(min_points, ceil(frac * n - 1e-9))).
``loess_detrend(y, pos, frac, *, circular, period=None, min_points=3) -> LoessFit(fit, resid, k)``:
plain local-linear LOESS (tricube, k nearest by ``pos``, no robustness iterations: lead decision
2, ``it = 0``), numerically equal (<= 1e-10) to statsmodels ``lowess(y_pad, pos_pad,
frac=k/n_pad, it=0, delta=0, return_sorted=False)`` on the original points. Circular padding
prepends the last min(k, n) points at ``pos - period`` and appends the first min(k, n) at
``pos + period``. Lead decision 1 / ADR-0020: ``pos`` may be the ordinals 0 .. n-1, so the smoother
is one fixed linear operator (a matrix depending only on ``pos``, ``frac``, ``circular``,
``period``) applied to any permutation of the values.

When the estimator does not exist:
- n = 0 -> empty fit and resid, k = 0.
- 1 <= n < min_points -> fit = mean(y) everywhere (residuals centred), k = n.
- all y equal -> the local fits exist; residuals 0.
- ``pos`` not strictly increasing, non-finite y or pos, frac outside (0, 1], circular without a
  period, or circular with pos outside [0, period) -> ValueError.

The oracle is statsmodels (BSD-3) on the padded series, written here from the spec text.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from statsmodels.nonparametric.smoothers_lowess import lowess

from gippa.rng import for_step
from gippa.stats.smooth import LoessFit, loess_detrend, loess_matrix, loess_span_points

FAST = settings(max_examples=60, deadline=None)


def span(n: int, frac: float, min_points: int = 3) -> int:
    return min(n, max(min_points, math.ceil(frac * n - 1e-9)))


def oracle_fit(
    y: np.ndarray, pos: np.ndarray, frac: float, circular: bool, period: float | None = None
) -> np.ndarray:
    n = y.size
    if n < 3:
        return np.full(n, y.mean())
    k = span(n, frac)
    if circular:
        assert period is not None
        m = min(k, n)
        pp = np.concatenate([pos[-m:] - period, pos, pos[:m] + period])
        yp = np.concatenate([y[-m:], y, y[:m]])
        out = lowess(yp, pp, frac=k / pp.size, it=0, delta=0, return_sorted=False)
        return np.asarray(out[m : m + n], float)
    return np.asarray(lowess(y, pos, frac=k / n, it=0, delta=0, return_sorted=False), float)


def gapped_positions(n: int, seed: int, max_gap: int = 3) -> tuple[np.ndarray, float]:
    rng = for_step("test.smooth.positions", seed)
    gaps = rng.integers(1, max_gap + 1, size=n)
    pos = np.concatenate([[0], np.cumsum(gaps[:-1])]).astype(float)
    return pos, float(pos[-1] + gaps[-1])


# ---------------------------------------------------------------------------- example 7: span


@pytest.mark.parametrize(
    ("n", "frac", "k"),
    [(30, 0.10, 3), (25, 0.10, 3), (10, 0.10, 3), (1000, 0.10, 100)],
)
def test_example_7_span_points(n: int, frac: float, k: int) -> None:
    assert loess_span_points(n, frac) == k


def test_span_points_rounding_and_floor() -> None:
    assert loess_span_points(30, 0.1) == 3  # 0.1 * 30 is 3.0000000000000004: not 4
    assert loess_span_points(101, 0.1) == 11  # ceil(10.1)
    assert loess_span_points(45, 0.1) == 5  # ceil(4.5)
    assert loess_span_points(4500, 0.1) == 450
    assert loess_span_points(100, 0.01, min_points=5) == 5
    assert loess_span_points(2, 0.1) == 2  # capped at n
    assert loess_span_points(1, 0.5) == 1
    assert loess_span_points(0, 0.1) == 0
    assert loess_span_points(10, 1.0) == 10


@FAST
@given(st.integers(0, 5000), st.floats(0.001, 1.0), st.integers(1, 10))
def test_span_points_formula(n: int, frac: float, min_points: int) -> None:
    k = loess_span_points(n, frac, min_points)
    assert isinstance(k, int)
    assert k == span(n, frac, min_points)
    assert 0 <= k <= n


# ---------------------------------------------------------------------------- examples 8-11


def test_example_8_affine_exact_linear() -> None:
    pos = np.arange(100.0)
    res = loess_detrend(2 * pos + 1, pos, 0.10, circular=False)
    assert isinstance(res, LoessFit)
    assert res.k == 10
    assert np.max(np.abs(res.resid)) < 1e-9


def test_example_9_cosine_circular_vs_linear() -> None:
    pos = np.arange(100.0)
    y = np.cos(2 * np.pi * pos / 100)
    circ = loess_detrend(y, pos, 0.10, circular=True, period=100.0)
    lin = loess_detrend(y, pos, 0.10, circular=False)
    assert circ.k == lin.k == 10
    assert abs(np.max(np.abs(circ.resid)) - 0.00711) <= 1e-5
    assert abs(np.max(np.abs(lin.resid)) - 0.00977) <= 1e-5
    # linear: largest at the ends; padding removes the end effect
    assert int(np.argmax(np.abs(lin.resid))) in {0, 99}
    assert abs(circ.resid[0]) < abs(lin.resid[0])
    assert abs(circ.resid[99]) < abs(lin.resid[99])
    # circular residuals are periodic: as large at the wrap as at the opposite extremum (pos 50)
    assert abs(abs(circ.resid[0]) - abs(circ.resid[50])) <= 1e-9


@pytest.mark.parametrize("circular", [False, True])
def test_example_10_two_points_mean(circular: bool) -> None:
    res = loess_detrend(
        [1.0, 3.0], [0, 1], 0.10, circular=circular, period=2.0 if circular else None
    )
    np.testing.assert_allclose(res.fit, [2.0, 2.0], atol=1e-15)
    np.testing.assert_allclose(res.resid, [-1.0, 1.0], atol=1e-15)
    assert res.k == 2


@pytest.mark.parametrize("circular", [False, True])
def test_example_11_constant(circular: bool) -> None:
    pos = np.arange(50.0)
    res = loess_detrend(np.full(50, 5.0), pos, 0.10, circular=circular, period=50.0)
    assert np.max(np.abs(res.resid)) <= 1e-12
    np.testing.assert_allclose(res.fit, 5.0, atol=1e-12)


# ---------------------------------------------------------------------------- degenerate inputs


def test_empty() -> None:
    for circular in (False, True):
        res = loess_detrend([], [], 0.10, circular=circular, period=10.0)
        assert res.fit.shape == (0,)
        assert res.resid.shape == (0,)
        assert res.fit.dtype == np.float64
        assert res.resid.dtype == np.float64
        assert res.k == 0


def test_single_point() -> None:
    res = loess_detrend([4.0], [0], 0.5, circular=False)
    assert res.fit.tolist() == [4.0]
    assert res.resid.tolist() == [0.0]
    assert res.k == 1


def test_below_custom_min_points_is_mean() -> None:
    y = np.array([1.0, 2.0, 6.0, 7.0])
    res = loess_detrend(y, np.arange(4), 0.5, circular=False, min_points=5)
    np.testing.assert_allclose(res.fit, 4.0, atol=1e-15)
    np.testing.assert_allclose(res.resid, y - 4.0, atol=1e-15)
    assert res.k == 4


def test_three_point_floor_matches_statsmodels() -> None:
    # n = 30, frac 0.10 -> k = 3 (the floor), not 4
    y = for_step("test.smooth.floor", 1).normal(size=30)
    pos = np.arange(30.0)
    for circular in (False, True):
        res = loess_detrend(y, pos, 0.10, circular=circular, period=30.0)
        assert res.k == 3
        np.testing.assert_allclose(
            res.fit, oracle_fit(y, pos, 0.10, circular, 30.0), rtol=0, atol=1e-10
        )


@pytest.mark.parametrize("n", [3, 10, 30])
@pytest.mark.parametrize("spacing", [1.0, 2.5])
def test_three_point_floor_on_equal_spacing_is_identity(n: int, spacing: float) -> None:
    # Lead decision 4: k = 3 on equally spaced positions gives both neighbours tricube weight 0
    # (they sit at the radius), so fit == y and every residual is 0. Documented degeneracy, not an
    # error: no ValueError, no NaN; consumers (M5-05, M5-06) flag k < 4 as `few_values`.
    y = for_step("test.smooth.floor", n).normal(size=n)
    pos = spacing * np.arange(float(n))
    period = spacing * n
    for circular in (False, True):
        res = loess_detrend(y, pos, 0.10, circular=circular, period=period)
        assert res.k == 3
        np.testing.assert_allclose(res.fit, y, rtol=0, atol=1e-12)
        assert np.max(np.abs(res.resid)) <= 1e-12
        H = loess_matrix(pos, 0.10, circular=circular, period=period)
        np.testing.assert_allclose(H, np.eye(n), rtol=0, atol=1e-12)


def test_three_point_floor_with_gaps_is_identity_too() -> None:
    # k = 3 is the identity for any strictly increasing positions, not only equal spacing: the
    # farthest of the 3 neighbours sits at the radius (weight 0) and the local line through the
    # remaining two points (one is the point itself) returns y_i. Checked against statsmodels.
    pos = np.array([0.0, 1.0, 3.0, 4.0, 7.0, 8.0, 9.0, 12.0, 13.0, 15.0])
    y = for_step("test.smooth.floor.gaps", 0).normal(size=10)
    for circular in (False, True):
        res = loess_detrend(y, pos, 0.10, circular=circular, period=17.0)
        assert res.k == 3
        np.testing.assert_allclose(res.fit, y, rtol=0, atol=1e-12)
        np.testing.assert_allclose(
            res.fit, oracle_fit(y, pos, 0.10, circular, 17.0), rtol=0, atol=1e-10
        )


def test_three_points_exact() -> None:
    y = np.array([1.0, 5.0, 2.0])
    res = loess_detrend(y, [0, 1, 2], 0.10, circular=False)
    assert res.k == 3
    np.testing.assert_allclose(res.fit, oracle_fit(y, np.arange(3.0), 0.1, False), atol=1e-10)


@pytest.mark.parametrize("n", [3, 4, 5, 7])
def test_circular_k_equals_n_tiny_replicon(n: int) -> None:
    # frac = 1 -> k = n: padding duplicates every point once on each side
    pos, period = gapped_positions(n, seed=n)
    y = for_step("test.smooth.tiny", n).normal(size=n)
    res = loess_detrend(y, pos, 1.0, circular=True, period=period)
    assert res.k == n
    np.testing.assert_allclose(res.fit, oracle_fit(y, pos, 1.0, True, period), rtol=0, atol=1e-10)


def test_resid_is_y_minus_fit_and_dtypes() -> None:
    pos, period = gapped_positions(40, seed=5)
    y = for_step("test.smooth.resid", 5).normal(size=40)
    for circular in (False, True):
        res = loess_detrend(y, pos, 0.2, circular=circular, period=period)
        assert res.fit.dtype == np.float64
        assert res.resid.dtype == np.float64
        assert res.fit.shape == res.resid.shape == (40,)
        np.testing.assert_allclose(res.resid, y - res.fit, rtol=0, atol=1e-12)
        assert res.k == loess_span_points(40, 0.2)


def test_loessfit_is_frozen() -> None:
    res = loess_detrend([1.0, 2.0, 3.0], [0, 1, 2], 0.5, circular=False)
    with pytest.raises(dataclasses.FrozenInstanceError):
        res.k = 7  # type: ignore[misc]


# ---------------------------------------------------------------------------- validation


@pytest.mark.parametrize(
    "pos",
    [[0, 1, 1, 2, 3], [0, 2, 1, 3, 4], [4, 3, 2, 1, 0]],
    ids=["tie", "unsorted", "decreasing"],
)
def test_pos_not_strictly_increasing_raises(pos: list[int]) -> None:
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0, 3.0, 4.0, 5.0], pos, 0.5, circular=False)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_raises(bad: float) -> None:
    with pytest.raises(ValueError):
        loess_detrend([1.0, bad, 3.0, 4.0], [0, 1, 2, 3], 0.5, circular=False)
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0, 3.0, 4.0], [0, 1, 2, bad], 0.5, circular=False)


@pytest.mark.parametrize("frac", [0.0, -0.1, 1.0000001, 2.0, float("nan")])
def test_frac_out_of_range_raises(frac: float) -> None:
    with pytest.raises(ValueError):
        loess_detrend(np.arange(10.0), np.arange(10), frac, circular=False)


def test_frac_one_allowed() -> None:
    res = loess_detrend(np.arange(10.0) ** 2, np.arange(10), 1.0, circular=False)
    assert res.k == 10


def test_circular_requires_period() -> None:
    with pytest.raises(ValueError):
        loess_detrend(np.arange(10.0), np.arange(10), 0.3, circular=True)


@pytest.mark.parametrize(
    ("pos", "period"),
    [([0, 1, 2, 10], 10.0), ([-1, 1, 2, 3], 10.0), ([0, 1, 2, 3], 3.0)],
    ids=["equals-period", "negative", "beyond-period"],
)
def test_circular_pos_outside_period_raises(pos: list[int], period: float) -> None:
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0, 3.0, 4.0], pos, 0.5, circular=True, period=period)


def test_length_mismatch_raises() -> None:
    # lead decision 5
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0, 3.0, 4.0], [0, 1, 2], 0.5, circular=False)
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0, 3.0], [0, 1, 2, 3], 0.5, circular=True, period=4.0)
    with pytest.raises(ValueError):
        loess_detrend([1.0, 2.0], [0, 1, 2], 0.5, circular=False)  # below min_points too
    with pytest.raises(ValueError):
        loess_detrend([], [0], 0.5, circular=False)


# ---------------------------------------------------------------------------- statsmodels


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("frac", [0.05, 0.10, 0.30, 0.50])
@pytest.mark.parametrize("n", [3, 4, 5, 10, 57, 300])
def test_equals_statsmodels_grid(n: int, frac: float, circular: bool) -> None:
    pos, period = gapped_positions(n, seed=1000 * n + int(100 * frac))
    y = for_step("test.smooth.grid", n).normal(size=n) + np.sin(pos / 7.0)
    res = loess_detrend(y, pos, frac, circular=circular, period=period)
    assert res.k == span(n, frac)
    np.testing.assert_allclose(
        res.fit, oracle_fit(y, pos, frac, circular, period), rtol=0, atol=1e-10
    )


@FAST
@given(
    n=st.integers(3, 300),
    frac=st.floats(0.05, 0.5),
    circular=st.booleans(),
    seed=st.integers(0, 2**31 - 1),
    max_gap=st.integers(1, 4),
)
def test_equals_statsmodels_random(
    n: int, frac: float, circular: bool, seed: int, max_gap: int
) -> None:
    pos, period = gapped_positions(n, seed, max_gap)
    y = for_step("test.smooth.random", seed).normal(0.0, 2.0, size=n)
    res = loess_detrend(y, pos, frac, circular=circular, period=period)
    np.testing.assert_allclose(
        res.fit, oracle_fit(y, pos, frac, circular, period), rtol=0, atol=1e-10
    )


def test_equals_statsmodels_real_positions() -> None:
    rng = for_step("test.smooth.real", 0)
    pos = np.cumsum(rng.uniform(0.1, 5.0, size=120))
    period = float(pos[-1] + 1.0)
    pos = pos - pos[0]
    y = rng.normal(size=120)
    for circular in (False, True):
        res = loess_detrend(y, pos, 0.1, circular=circular, period=period)
        np.testing.assert_allclose(
            res.fit, oracle_fit(y, pos, 0.1, circular, period), rtol=0, atol=1e-10
        )


# ---------------------------------------------------------------------------- ordinal interface


@pytest.mark.parametrize("circular", [False, True])
def test_ordinal_positions_interface(circular: bool) -> None:
    # ADR-0020: positions 0 .. n-1 (ordinal among genes with values), period n; int input works
    n = 200
    y = for_step("test.smooth.ordinal", 0).normal(size=n)
    res_int = loess_detrend(y, list(range(n)), 0.10, circular=circular, period=n)
    res_flt = loess_detrend(y, np.arange(n, dtype=float), 0.10, circular=circular, period=float(n))
    np.testing.assert_array_equal(res_int.fit, res_flt.fit)
    np.testing.assert_allclose(
        res_int.fit, oracle_fit(y, np.arange(n, dtype=float), 0.10, circular, float(n)), atol=1e-10
    )


def test_gaps_in_positions_are_respected() -> None:
    n = 60
    y = for_step("test.smooth.gaps", 0).normal(size=n)
    gapped = np.concatenate([np.arange(30.0), np.arange(30.0) + 100.0])
    ordinal = np.arange(float(n))
    a = loess_detrend(y, gapped, 0.2, circular=False)
    b = loess_detrend(y, ordinal, 0.2, circular=False)
    np.testing.assert_allclose(a.fit, oracle_fit(y, gapped, 0.2, False), atol=1e-10)
    assert np.max(np.abs(a.fit - b.fit)) > 1e-3


def test_translation_and_scale_invariant_linear() -> None:
    n = 80
    y = for_step("test.smooth.affinepos", 0).normal(size=n)
    pos, _ = gapped_positions(n, seed=11)
    base = loess_detrend(y, pos, 0.15, circular=False).fit
    np.testing.assert_allclose(
        loess_detrend(y, pos + 1000, 0.15, circular=False).fit, base, atol=1e-9
    )
    np.testing.assert_allclose(
        loess_detrend(y, 10 * pos, 0.15, circular=False).fit, base, atol=1e-9
    )


# ---------------------------------------------------------------------------- linear operator


def smoother_matrix(
    pos: np.ndarray, frac: float, circular: bool, period: float | None
) -> np.ndarray:
    n = pos.size
    eye = np.eye(n)
    return np.column_stack(
        [loess_detrend(eye[j], pos, frac, circular=circular, period=period).fit for j in range(n)]
    )


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("n", [2, 3, 12, 45])
def test_smoother_is_one_linear_operator(n: int, circular: bool) -> None:
    pos, period = gapped_positions(n, seed=n)
    if n == 45:
        pos, period = np.arange(45.0), 45.0  # ordinals, as ADR-0020 uses them
    H = smoother_matrix(pos, 0.10, circular, period)
    # rows sum to 1: constants are reproduced
    np.testing.assert_allclose(H.sum(axis=1), 1.0, atol=1e-10)
    rng = for_step("test.smooth.linear", n)
    for _ in range(5):
        y = rng.normal(size=n)
        np.testing.assert_allclose(
            loess_detrend(y, pos, 0.10, circular=circular, period=period).fit, H @ y, atol=1e-10
        )
        # the same matrix serves every permutation of the values (it depends only on pos)
        yp = y[rng.permutation(n)]
        np.testing.assert_allclose(
            loess_detrend(yp, pos, 0.10, circular=circular, period=period).fit, H @ yp, atol=1e-10
        )
    # additivity and homogeneity directly
    a, b = rng.normal(size=n), rng.normal(size=n)
    fa = loess_detrend(a, pos, 0.10, circular=circular, period=period).fit
    fb = loess_detrend(b, pos, 0.10, circular=circular, period=period).fit
    fab = loess_detrend(2.5 * a - 0.7 * b, pos, 0.10, circular=circular, period=period).fit
    np.testing.assert_allclose(fab, 2.5 * fa - 0.7 * fb, atol=1e-10)


def test_below_min_points_operator_is_mean() -> None:
    H = smoother_matrix(np.array([0.0, 1.0]), 0.1, False, None)
    np.testing.assert_allclose(H, np.full((2, 2), 0.5), atol=1e-15)


def test_linear_operator_reproduces_positions() -> None:
    pos, _ = gapped_positions(50, seed=3)
    H = smoother_matrix(pos, 0.10, False, None)
    np.testing.assert_allclose(H @ pos, pos, atol=1e-9)


# ---------------------------------------------------------------------------- purity


def test_deterministic_and_inputs_not_mutated() -> None:
    pos, period = gapped_positions(70, seed=9)
    y = for_step("test.smooth.pure", 9).normal(size=70)
    y0, p0 = y.copy(), pos.copy()
    r1 = loess_detrend(y, pos, 0.1, circular=True, period=period)
    r2 = loess_detrend(y, pos, 0.1, circular=True, period=period)
    np.testing.assert_array_equal(r1.fit, r2.fit)
    np.testing.assert_array_equal(y, y0)
    np.testing.assert_array_equal(pos, p0)


# ------------------------------------------------------------------- loess_matrix (decision 6)


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("frac", [0.05, 0.10, 0.30, 1.0])
@pytest.mark.parametrize("n", [1, 2, 3, 4, 11, 60, 200])
def test_loess_matrix_matches_detrend(n: int, frac: float, circular: bool) -> None:
    pos, period = gapped_positions(n, seed=7 * n + int(100 * frac))
    H = loess_matrix(pos, frac, circular=circular, period=period)
    assert isinstance(H, np.ndarray)
    assert H.dtype == np.float64
    assert H.shape == (n, n)
    rng = for_step("test.smooth.matrix", n)
    for _ in range(3):
        y = rng.normal(0.0, 3.0, size=n)
        fit = loess_detrend(y, pos, frac, circular=circular, period=period).fit
        np.testing.assert_allclose(H @ y, fit, rtol=0, atol=1e-12)
    np.testing.assert_allclose(H.sum(axis=1), 1.0, rtol=0, atol=1e-12)


@pytest.mark.parametrize("circular", [False, True])
def test_loess_matrix_on_batch_matches_per_column(circular: bool) -> None:
    # ADR-0020: one matrix applied to B permuted replicates (columns of an (n, B) batch)
    n, B = 120, 25
    pos, period = np.arange(float(n)), float(n)
    H = loess_matrix(pos, 0.10, circular=circular, period=period)
    rng = for_step("test.smooth.batch", 0)
    y = rng.normal(size=n)
    Y = np.column_stack([y[rng.permutation(n)] for _ in range(B)])
    fits = H @ Y
    assert fits.shape == (n, B)
    for b in range(B):
        ref = loess_detrend(Y[:, b], pos, 0.10, circular=circular, period=period).fit
        np.testing.assert_allclose(fits[:, b], ref, rtol=0, atol=1e-12)


def test_loess_matrix_ordinals_and_gaps() -> None:
    pos = np.concatenate([np.arange(30.0), np.arange(30.0) + 100.0])
    H_gap = loess_matrix(pos, 0.2, circular=False)
    H_ord = loess_matrix(np.arange(60.0), 0.2, circular=False)
    assert np.max(np.abs(H_gap - H_ord)) > 1e-3
    np.testing.assert_allclose(H_gap @ pos, pos, atol=1e-9)  # local linear reproduces pos
    # int positions accepted, identical to float
    np.testing.assert_array_equal(
        loess_matrix(list(range(60)), 0.2, circular=True, period=60),
        loess_matrix(np.arange(60.0), 0.2, circular=True, period=60.0),
    )


def test_loess_matrix_degenerate() -> None:
    assert loess_matrix([], 0.1, circular=False).shape == (0, 0)
    assert loess_matrix([], 0.1, circular=True, period=5.0).shape == (0, 0)
    np.testing.assert_allclose(loess_matrix([0, 1], 0.1, circular=False), np.full((2, 2), 0.5))
    np.testing.assert_allclose(
        loess_matrix([0, 1, 2, 3], 0.5, circular=False, min_points=5), np.full((4, 4), 0.25)
    )
    np.testing.assert_array_equal(loess_matrix([3], 0.5, circular=True, period=4.0), [[1.0]])


def test_loess_matrix_matches_statsmodels_unit_vectors() -> None:
    pos, period = gapped_positions(40, seed=21)
    for circular in (False, True):
        H = loess_matrix(pos, 0.15, circular=circular, period=period)
        eye = np.eye(40)
        ref = np.column_stack([oracle_fit(eye[j], pos, 0.15, circular, period) for j in range(40)])
        np.testing.assert_allclose(H, ref, rtol=0, atol=1e-10)


@pytest.mark.parametrize(
    ("pos", "kw"),
    [
        ([0, 1, 1, 2], {"circular": False}),
        ([0, 1, float("nan"), 3], {"circular": False}),
        ([0, 1, 2, 3], {"circular": True}),
        ([0, 1, 2, 4], {"circular": True, "period": 4.0}),
    ],
    ids=["tie", "nan", "no-period", "outside-period"],
)
def test_loess_matrix_validation(pos: list[float], kw: dict[str, object]) -> None:
    # assumed: the pos / frac / period rules of loess_detrend apply unchanged (decision 6 gives
    # the same signature without y)
    with pytest.raises(ValueError):
        loess_matrix(pos, 0.5, **kw)  # type: ignore[arg-type]


@pytest.mark.parametrize("frac", [0.0, 1.5, float("nan")])
def test_loess_matrix_frac_out_of_range(frac: float) -> None:
    with pytest.raises(ValueError):
        loess_matrix(np.arange(10.0), frac, circular=False)


def test_loess_matrix_deterministic_and_pos_not_mutated() -> None:
    pos, period = gapped_positions(50, seed=4)
    p0 = pos.copy()
    np.testing.assert_array_equal(
        loess_matrix(pos, 0.1, circular=True, period=period),
        loess_matrix(pos, 0.1, circular=True, period=period),
    )
    np.testing.assert_array_equal(pos, p0)
