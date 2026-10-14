"""Benjamini-Hochberg q-values (issue M4-01; GIAP_SPEC.md §8.2, §8.3, Appendix B.1)."""

from __future__ import annotations

import numpy as np
import pytest
from statsmodels.stats.multitest import multipletests

from giap.stats.multitest import bh

NAN = float("nan")

B1_P = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
B1_Q = [0.01, 0.04, 0.084, 0.084, 0.084, 0.1, 0.105714, 0.216, 0.216, 0.216]


def test_example1_appendix_b1_exact_numbers() -> None:
    q = bh(B1_P)
    assert q.dtype == np.float64
    assert q.shape == (10,)
    np.testing.assert_allclose(q, B1_Q, rtol=0, atol=1e-6)


def test_appendix_b1_in_shuffled_order_keeps_input_order() -> None:
    perm = [7, 2, 9, 0, 4, 1, 8, 3, 6, 5]
    q = bh(np.asarray(B1_P)[perm])
    np.testing.assert_allclose(q, np.asarray(B1_Q)[perm], rtol=0, atol=1e-6)


def test_example2_nan_stays_nan_and_does_not_count_toward_m() -> None:
    q = bh([0.01, NAN, 0.04, 0.03])
    np.testing.assert_allclose(q, [0.03, NAN, 0.04, 0.04], rtol=0, atol=1e-12, equal_nan=True)


def test_example3_all_tied() -> None:
    np.testing.assert_allclose(bh([0.5, 0.5, 0.5]), [0.5, 0.5, 0.5], rtol=0, atol=1e-12)


def test_ties_get_equal_q() -> None:
    q = bh([0.02, 0.01, 0.02, 0.9])
    # sorted: 0.01, 0.02, 0.02, 0.9 -> raw 0.04, 0.04, 0.0267, 0.9 -> monotone 0.0267 x3, 0.9
    np.testing.assert_allclose(q, [0.08 / 3, 0.08 / 3, 0.08 / 3, 0.9], rtol=1e-12)
    assert q[0] == q[2]


def test_example3_empty_input() -> None:
    q = bh([])
    assert isinstance(q, np.ndarray)
    assert q.shape == (0,)
    assert q.dtype == np.float64


def test_all_nan_gives_all_nan() -> None:
    q = bh([NAN, NAN])
    assert q.shape == (2,)
    assert np.isnan(q).all()


def test_single_p_is_returned_unchanged() -> None:
    np.testing.assert_allclose(bh([0.3]), [0.3], rtol=0, atol=0)


def test_bounds_zero_and_one_are_valid() -> None:
    np.testing.assert_allclose(bh([0.0, 1.0]), [0.0, 1.0], rtol=0, atol=0)


def test_q_never_exceeds_one() -> None:
    q = bh([0.9, 0.95, 1.0, 0.99])
    assert np.all(q <= 1.0)


@pytest.mark.parametrize("bad", [[0.1, 1.0000001], [-1e-9, 0.5], [2.0], [-0.5]])
def test_p_outside_unit_interval_raises(bad: list[float]) -> None:
    with pytest.raises(ValueError):
        bh(bad)


def test_input_is_not_mutated() -> None:
    p = np.array([0.04, 0.01, NAN, 0.03])
    before = p.copy()
    bh(p)
    np.testing.assert_array_equal(p, before)


def test_accepts_list_tuple_and_array() -> None:
    a = bh([0.01, 0.04, 0.03])
    b = bh((0.01, 0.04, 0.03))
    c = bh(np.array([0.01, 0.04, 0.03]))
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(a, c)


def test_agrees_with_statsmodels_fdr_bh_on_100_random_vectors() -> None:
    rng = np.random.default_rng(20261014)
    for _ in range(100):
        n = int(rng.integers(1, 200))
        # mixture of nulls and small p values, with some exact ties
        p = np.where(rng.random(n) < 0.3, rng.beta(0.2, 5.0, n), rng.random(n))
        if n > 3:
            p[rng.integers(0, n, 3)] = p[0]
        expected = multipletests(p, method="fdr_bh")[1]
        np.testing.assert_allclose(bh(p), expected, rtol=1e-12, atol=0)


def test_agrees_with_statsmodels_after_dropping_nan() -> None:
    rng = np.random.default_rng(20261015)
    p = rng.random(50)
    p[rng.choice(50, 10, replace=False)] = np.nan
    ok = ~np.isnan(p)
    q = bh(p)
    assert np.isnan(q[~ok]).all()
    np.testing.assert_allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1], rtol=1e-12)
