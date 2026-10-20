"""Left-censored Blom normal scores (issue M5-01; GIAP_SPEC.md §7.3.5(d) step 1, §7.3.6, §8.3).

``blom_scores(x, censored_mask)``: censored entries (mask True) all get rank (n0 + 1) / 2, the
average of ranks 1 .. n0; uncensored entries get n0 + their average rank among the uncensored
values; score = ndtri((r - 3/8) / (n + 1/4)), output in input order.

When the estimator does not exist:
- n = 0 -> empty float64 array (no ranks to score).
- every entry censored -> every rank is (n + 1) / 2, so every score is exactly 0.0.
- an *uncensored* NaN or +-inf -> ValueError (the caller decides who is scored); the value of a
  censored entry is ignored (it may be NaN or inf).
- len(x) != len(censored_mask) -> ValueError.

The oracle below is written from the spec text with scipy.stats.rankdata (independent of
giap.stats.scores).
"""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.special import ndtri
from scipy.stats import rankdata

from giap.stats.scores import blom_scores

nan = float("nan")
F, T = False, True

FAST = settings(max_examples=150, deadline=None)


def oracle(x: list[float], mask: list[bool]) -> np.ndarray:
    xa = np.asarray(x, float)
    m = np.asarray(mask, bool)
    n = xa.size
    n0 = int(m.sum())
    r = np.empty(n)
    r[m] = (n0 + 1) / 2.0
    if (~m).any():
        r[~m] = n0 + rankdata(xa[~m], method="average")
    return np.asarray(ndtri((r - 0.375) / (n + 0.25)), float)


# ---------------------------------------------------------------------------- spec examples


def test_example_1_censored_tied_lowest() -> None:
    got = blom_scores([5, nan, 3, nan], [F, T, F, T])
    # ranks [4, 1.5, 3, 1.5], n = 4
    expected = ndtri((np.array([4, 1.5, 3, 1.5]) - 0.375) / 4.25)
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(got, [1.0491, -0.6289, 0.2993, -0.6289], rtol=0, atol=5e-5)


def test_example_2_ties_average_ranks() -> None:
    got = blom_scores([2, 2, 1], [F, F, F])
    expected = ndtri((np.array([2.5, 2.5, 1.0]) - 0.375) / 3.25)
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)
    np.testing.assert_allclose(got, [0.3957, 0.3957, -0.8694], rtol=0, atol=5e-5)


def test_example_3_all_censored_exactly_zero() -> None:
    got = blom_scores([nan, nan, nan, nan], [T, T, T, T])
    assert got.tolist() == [0.0, 0.0, 0.0, 0.0]


def test_example_4_single_uncensored_is_zero() -> None:
    got = blom_scores([7.5], [F])
    assert got.tolist() == [0.0]


def test_example_5_empty() -> None:
    got = blom_scores([], [])
    assert isinstance(got, np.ndarray)
    assert got.shape == (0,)
    assert got.dtype == np.float64


@pytest.mark.parametrize("bad", [nan, float("inf"), float("-inf")])
def test_example_6_uncensored_non_finite_raises(bad: float) -> None:
    with pytest.raises(ValueError):
        blom_scores([1.0, bad, 2.0], [F, F, F])


# ---------------------------------------------------------------------------- edge cases


def test_length_mismatch_raises() -> None:
    with pytest.raises(ValueError):
        blom_scores([1.0, 2.0, 3.0], [F, F])
    with pytest.raises(ValueError):
        blom_scores([], [T])


def test_single_censored_is_zero() -> None:
    assert blom_scores([nan], [T]).tolist() == [0.0]


@pytest.mark.parametrize("val", [nan, float("inf"), float("-inf"), 1e6, -1e6])
def test_censored_value_is_ignored(val: float) -> None:
    ref = blom_scores([5.0, nan, 3.0], [F, T, F])
    got = blom_scores([5.0, val, 3.0], [F, T, F])
    np.testing.assert_array_equal(got, ref)


def test_low_detected_still_ranks_above_censored() -> None:
    # censoring is "at the lowest rank", not "at its own value"
    got = blom_scores([-100.0, 50.0, 1.0], [F, T, F])
    # ranks [2, 1, 3]
    np.testing.assert_allclose(got, ndtri((np.array([2.0, 1.0, 3.0]) - 0.375) / 3.25), atol=1e-12)
    assert got[1] < got[0] < got[2]


def test_ties_in_detected_and_censored_block_averaged_separately() -> None:
    x = [4.0, nan, 4.0, 1.0, nan, nan, 9.0]
    m = [F, T, F, F, T, T, F]
    # n0 = 3 -> censored rank 2; detected [4, 4, 1, 9] -> ranks 2.5, 2.5, 1, 4 -> + 3
    ranks = np.array([5.5, 2.0, 5.5, 4.0, 2.0, 2.0, 7.0])
    np.testing.assert_allclose(blom_scores(x, m), ndtri((ranks - 0.375) / 7.25), atol=1e-12)


def test_all_detected_tied_is_zero() -> None:
    np.testing.assert_allclose(blom_scores([3.0] * 5, [F] * 5), np.zeros(5), atol=1e-15)


def test_output_in_input_order_and_dtype() -> None:
    x = [3.0, 1.0, 2.0]
    got = blom_scores(np.array(x), np.array([F, F, F]))
    assert got.dtype == np.float64
    assert got.shape == (3,)
    assert got[1] < got[2] < got[0]


def test_numpy_bool_mask_and_inputs_not_mutated() -> None:
    x = np.array([2.0, nan, 1.0])
    m = np.array([False, True, False])
    x0, m0 = x.copy(), m.copy()
    blom_scores(x, m)
    np.testing.assert_array_equal(x, x0)
    np.testing.assert_array_equal(m, m0)


def test_deterministic() -> None:
    x = [0.3, nan, -1.2, 4.0, 4.0]
    m = [F, T, F, F, F]
    np.testing.assert_array_equal(blom_scores(x, m), blom_scores(x, m))


# ---------------------------------------------------------------------------- properties

VALUES = st.lists(st.integers(-20, 20).map(float), min_size=0, max_size=40)


@FAST
@given(st.lists(st.floats(-1e6, 1e6, allow_nan=False), min_size=1, max_size=60, unique=True))
def test_symmetry_no_ties_no_censoring(x: list[float]) -> None:
    s = blom_scores(x, [False] * len(x))
    ordered = np.sort(s)
    np.testing.assert_allclose(ordered + ordered[::-1], 0.0, rtol=0, atol=1e-12)
    assert abs(float(s.sum())) <= 1e-10


@FAST
@given(st.data())
def test_matches_spec_oracle(data: st.DataObject) -> None:
    x = data.draw(VALUES)
    mask = data.draw(st.lists(st.booleans(), min_size=len(x), max_size=len(x)))
    xm = [nan if c else v for v, c in zip(x, mask, strict=True)]
    np.testing.assert_allclose(blom_scores(xm, mask), oracle(x, mask), rtol=0, atol=1e-12)


@FAST
@given(st.data())
def test_censored_below_every_detected(data: st.DataObject) -> None:
    x = data.draw(st.lists(st.floats(-1e3, 1e3), min_size=2, max_size=40))
    mask = data.draw(st.lists(st.booleans(), min_size=len(x), max_size=len(x)))
    m = np.asarray(mask, bool)
    s = blom_scores(x, mask)
    if m.any():
        assert np.all(s[m] == s[m][0])
    if m.any() and (~m).any():
        assert s[m].max() < s[~m].min()


@FAST
@given(st.data())
def test_monotone_transform_and_permutation_equivariance(data: st.DataObject) -> None:
    x = data.draw(VALUES)
    mask = data.draw(st.lists(st.booleans(), min_size=len(x), max_size=len(x)))
    s = blom_scores(x, mask)
    # strictly increasing transform of the values leaves scores unchanged
    np.testing.assert_array_equal(blom_scores([3.0 * v + 7.0 for v in x], mask), s)
    perm = data.draw(st.permutations(range(len(x))))
    xp = [x[i] for i in perm]
    mp = [mask[i] for i in perm]
    np.testing.assert_allclose(blom_scores(xp, mp), s[list(perm)], rtol=0, atol=1e-15)
