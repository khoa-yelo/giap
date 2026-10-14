"""Hodges-Lehmann two-sample shift (issue M4-04; GIAP_SPEC.md §7.3.5(b), §7.3.6, §8.3).

``hodges_lehmann(x, y)`` = median of all pairwise differences x_i - y_j; the median of an even
count is the mean of the two middle values. Exact on both paths: direct enumeration for
m*n <= 5e4 (_DIRECT_MAX, lead decision after review), bisection on sorted y above (it must equal
the direct median bit for bit). Empty input gives NaN; NaN or +-inf in the input raises
ValueError (callers drop missing abundances first); see test_robust_review.py.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from giap.rng import for_step
from giap.stats.robust import hodges_lehmann


def brute(x: object, y: object) -> float:
    return float(np.median(np.subtract.outer(np.asarray(x, float), np.asarray(y, float))))


# ---------------------------------------------------------------------------- spec examples


def test_example_1_hand_computed() -> None:
    # differences [2, 1, -1, 4, 3, 1] -> sorted [-1, 1, 1, 2, 3, 4] -> (1 + 2) / 2
    assert hodges_lehmann([10, 12], [8, 9, 11]) == 1.5


def test_example_2_single_pair_and_empty() -> None:
    assert hodges_lehmann([5], [5]) == 0.0
    assert math.isnan(hodges_lehmann([], [1]))
    assert math.isnan(hodges_lehmann([1], []))
    assert math.isnan(hodges_lehmann([], []))


def test_example_3_large_path_equals_direct_median_exactly() -> None:
    rng = for_step("test.hodges_lehmann", 3)
    x = rng.normal(1.0, 1.0, 2000)
    y = rng.normal(0.0, 1.0, 3000)  # 6e6 pairs > 5e4: bisection path
    expected = float(np.median(np.subtract.outer(x, y)))
    got = hodges_lehmann(x, y)
    assert isinstance(got, float)
    assert got == expected
    assert 0.8 < got < 1.2


# ---------------------------------------------------------------------------- more exact cases


def test_odd_count_is_middle_difference() -> None:
    # differences of [1, 4, 6] - [0]: [1, 4, 6] -> 4
    assert hodges_lehmann([1, 4, 6], [0]) == 4.0
    # x = [0], y = [1, 2, 7]: differences [-1, -2, -7] -> -2
    assert hodges_lehmann([0], [1, 2, 7]) == -2.0


def test_hand_computed_three_by_two() -> None:
    # x = [3, 7, 8], y = [1, 5]: differences [2, -2, 6, 2, 7, 3] -> sorted [-2, 2, 2, 3, 6, 7]
    assert hodges_lehmann([3, 7, 8], [1, 5]) == 2.5


def test_all_ties_gives_zero() -> None:
    assert hodges_lehmann([3.0] * 4, [3.0] * 5) == 0.0


def test_constant_samples_give_their_difference() -> None:
    assert hodges_lehmann([2.5] * 3, [-1.0] * 2) == 3.5


def test_large_path_with_heavy_ties_is_exact() -> None:
    # 1500 * 1000 = 1.5e6 pairs on a small integer grid: many tied differences, and the two middle
    # order statistics can differ (even count).
    rng = for_step("test.hodges_lehmann.ties", 7)
    x = rng.integers(0, 6, 1500).astype(float)
    y = rng.integers(0, 4, 1000).astype(float)
    assert hodges_lehmann(x, y) == brute(x, y)


@pytest.mark.parametrize(("m", "n"), [(1000, 1000), (1001, 999), (1, 1_000_001)])
def test_threshold_boundary_matches_direct(m: int, n: int) -> None:
    rng = for_step("test.hodges_lehmann.boundary", m + n)
    x = rng.normal(0.0, 2.0, m)
    y = rng.normal(0.3, 1.0, n)
    assert hodges_lehmann(x, y) == brute(x, y)


def test_accepts_lists_tuples_and_int_arrays() -> None:
    assert hodges_lehmann((10, 12), np.array([8, 9, 11])) == 1.5
    assert hodges_lehmann(np.array([10, 12], dtype=np.int64), [8, 9, 11]) == 1.5


def test_does_not_mutate_inputs() -> None:
    x = np.array([5.0, 1.0, 3.0])
    y = np.array([2.0, 0.0])
    hodges_lehmann(x, y)
    np.testing.assert_array_equal(x, [5.0, 1.0, 3.0])
    np.testing.assert_array_equal(y, [2.0, 0.0])


# ---------------------------------------------------------------------------- errors


@pytest.mark.parametrize(
    ("x", "y"),
    [([1.0, float("nan")], [0.0]), ([1.0], [float("nan")]), ([float("nan")], [])],
)
def test_nan_input_raises(x: list[float], y: list[float]) -> None:
    with pytest.raises(ValueError):
        hodges_lehmann(x, y)
