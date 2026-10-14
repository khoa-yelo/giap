"""Lead decisions after the statistician's review of M4-03 (PR #23;
docs/specs/M4-03-stats-permutation-shifts.md, "Lead decisions").

1. Two-sided ``perm_p`` centres on the pooled mean (sum(null) + obs)/(B + 1), not the null mean.
2. +-inf null values are kept (they count in B) and count as extreme; NaN stays absent.
3. ``choose_shifts`` reports exhaustive=True whenever every valid shift is used.
"""

from __future__ import annotations

import numpy as np
import pytest

from giap.rng import for_step
from giap.stats.permutation import choose_shifts, perm_p

NAN = float("nan")
INF = float("inf")


def close(actual: object, expected: float, atol: float = 1e-12) -> None:
    np.testing.assert_allclose(float(actual), expected, rtol=0, atol=atol)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------- 1. pooled centring


def test_pooled_centre_differs_from_null_mean_small_example() -> None:
    # centre = (0 + 2 + 3)/3 = 5/3; |obs dev| = 4/3; null |dev| = 5/3, 1/3 -> b = 1 -> 2/3.
    # (Null-mean centring would give centre 1, |obs dev| 2, null |dev| 1, 1 -> 1/3.)
    close(perm_p(3.0, [0.0, 2.0]), 2 / 3)


def test_pooled_centre_negative_obs() -> None:
    # centre = (0 + 0 + 4 - 3)/4 = 0.25; |obs dev| = 3.25; null |dev| = .25, .25, 3.75 -> 2/4.
    # (Null-mean centring: centre 4/3, |obs dev| 13/3 beats every null -> 1/4.)
    close(perm_p(-3.0, [0.0, 0.0, 4.0]), 2 / 4)


def test_pooled_centre_example1_from_spec() -> None:
    close(perm_p(5.5, [1, 2, 3, 4, 5]), 2 / 6)


@pytest.mark.parametrize("shift", [-1e3, -2.5, 0.0, 7.0, 1e4])
def test_pooled_centre_translation_invariant(shift: float) -> None:
    close(perm_p(3.0 + shift, np.array([0.0, 2.0]) + shift), 2 / 3)
    close(perm_p(-3.0 + shift, np.array([0.0, 0.0, 4.0]) + shift), 2 / 4)


def test_two_sided_exact_validity_under_exchangeability_B19() -> None:
    """With a common centre, rotating which of 20 exchangeable values is 'observed' gives at most
    k of the 20 positions with p <= k/20 (exact validity). Null-mean centring breaks this."""
    rng = for_step("test.perm_review", 19)
    for _ in range(200):
        vals = rng.standard_t(3, size=20)
        ps = np.array([perm_p(vals[i], np.delete(vals, i)) for i in range(20)])
        for k in range(1, 21):
            assert np.count_nonzero(ps <= k / 20 + 1e-12) <= k


# ---------------------------------------------------------------------------- 2. +-inf kept


def test_inf_counts_as_extreme_greater() -> None:
    close(perm_p(5, [1, 2, INF], "greater"), 2 / 4)
    close(perm_p(5, [1, 2, -INF], "greater"), 1 / 4)  # kept in B, not an exceedance


def test_inf_counts_as_extreme_less() -> None:
    close(perm_p(5, [7, 8, -INF], "less"), 2 / 4)
    close(perm_p(5, [7, 8, INF], "less"), 1 / 4)  # kept in B, not an exceedance


@pytest.mark.parametrize("v", [INF, -INF])
def test_inf_counts_as_extreme_two_sided(v: float) -> None:
    close(perm_p(5, [1, 2, v]), 2 / 4)
    close(perm_p(-5, [-1, -2, v]), 2 / 4)


def test_both_infinities_two_sided() -> None:
    close(perm_p(5, [1, 2, INF, -INF]), 3 / 5)


def test_inf_counts_in_B_nan_does_not() -> None:
    close(perm_p(5, [1, 2, INF, NAN], "greater"), 2 / 4)
    close(perm_p(5, [NAN, INF], "greater"), 1.0)
    close(perm_p(5, [NAN, -INF], "greater"), 1 / 2)


# ---------------------------------------------------------------------------- NaN still absent


def test_nan_absent_from_counts_and_centre() -> None:
    close(perm_p(5, [1, 2, NAN], "greater"), 1 / 3)
    close(perm_p(-5, [1, 2, NAN], "less"), 1 / 3)
    close(perm_p(3.0, [0.0, NAN, 2.0]), 2 / 3)  # same as without the NaN
    close(perm_p(-3.0, [NAN, 0.0, 0.0, NAN, 4.0]), 2 / 4)


# ---------------------------------------------------------------------------- 3. exhaustive flag


@pytest.mark.parametrize("B", [7, 8, 1_000])
def test_random_branch_using_all_valid_shifts_is_exhaustive(B: int) -> None:
    valid = [1, 2, 3, 4, 5, 8, 9]
    shifts, exhaustive = choose_shifts(valid, 5, 10, B, for_step("context", 7))
    assert [int(x) for x in shifts] == valid
    assert bool(exhaustive) is True


def test_random_branch_with_fewer_than_all_shifts_is_not_exhaustive() -> None:
    valid = [1, 2, 3, 4, 5, 8, 9]
    shifts, exhaustive = choose_shifts(valid, 5, 10, 6, for_step("context", 7))
    assert len(shifts) == 6
    assert bool(exhaustive) is False


def test_threshold_branch_still_exhaustive() -> None:
    shifts, exhaustive = choose_shifts(np.arange(1, 10), 10, 10, 3, for_step("context", 7))
    assert len(shifts) == 9
    assert bool(exhaustive) is True
