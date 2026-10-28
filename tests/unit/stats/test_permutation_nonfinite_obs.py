"""``perm_p`` with a non-finite observed statistic (M5-03 lead decision 17, PR #118 review B1/N1).

Before the fix ``gippa.stats.permutation._tol`` raised ``TypeError`` for a non-finite scale, so
``perm_p(+-inf, ...)`` raised for the one-sided alternatives. Pinned semantics (the existing rule
"values within the tie tolerance count as exceedances" with infinities compared exactly):

- NaN obs -> NaN (every alternative).
- "greater": +inf obs -> (1 + #{null == +inf}) / (1 + B), i.e. the +1 floor 1/(B + 1) when the
  null has no +inf and ties counted when it has; -inf obs -> every null value is >= -inf, p = 1.
- "less": mirror image; -inf obs -> (1 + #{null == -inf}) / (1 + B); +inf obs -> 1.
- "two-sided": the centre is the mean of the finite values (an infinite obs is left out, as the
  docstring says), |obs - centre| = inf, so b = #{null is +-inf}:
  p = (1 + #{|null| = inf}) / (1 + B) for either sign of obs.
- NaN null entries stay absent (B counts the rest).
"""

from __future__ import annotations

import numpy as np
import pytest

from gippa.stats.permutation import perm_p

NAN = float("nan")
INF = float("inf")


def close(actual: object, expected: float) -> None:
    np.testing.assert_allclose(float(actual), expected, rtol=0, atol=1e-15)  # type: ignore[arg-type]


@pytest.mark.parametrize("alternative", ["greater", "less", "two-sided"])
def test_nan_obs_is_nan(alternative: str) -> None:
    assert np.isnan(perm_p(NAN, [1.0, 2.0, INF], alternative))  # type: ignore[arg-type]


def test_plus_inf_greater_floor_without_inf_in_null() -> None:
    close(perm_p(INF, [1.0, 2.0, 3.0], "greater"), 1 / 4)
    close(perm_p(INF, [1.0, -INF, 3.0], "greater"), 1 / 4)


def test_plus_inf_greater_counts_inf_ties() -> None:
    close(perm_p(INF, [1.0, INF, 3.0], "greater"), 2 / 4)
    close(perm_p(INF, [INF, INF, NAN], "greater"), 1.0)


def test_minus_inf_greater_is_one() -> None:
    close(perm_p(-INF, [1.0, 2.0, 3.0], "greater"), 1.0)
    close(perm_p(-INF, [-INF, 2.0], "greater"), 1.0)


def test_minus_inf_less_mirrors_greater() -> None:
    close(perm_p(-INF, [1.0, 2.0, 3.0], "less"), 1 / 4)
    close(perm_p(-INF, [1.0, -INF, INF], "less"), 2 / 4)


def test_plus_inf_less_is_one() -> None:
    close(perm_p(INF, [1.0, 2.0, 3.0], "less"), 1.0)


@pytest.mark.parametrize("obs", [INF, -INF])
def test_inf_two_sided_counts_infinite_null(obs: float) -> None:
    close(perm_p(obs, [1.0, 2.0, 3.0]), 1 / 4)
    close(perm_p(obs, [1.0, INF, 3.0]), 2 / 4)
    close(perm_p(obs, [-INF, INF, 3.0, NAN]), 3 / 4)


@pytest.mark.parametrize("alternative", ["greater", "less", "two-sided"])
def test_inf_obs_empty_or_all_nan_null_is_nan(alternative: str) -> None:
    assert np.isnan(perm_p(INF, [], alternative))  # type: ignore[arg-type]
    assert np.isnan(perm_p(-INF, [NAN], alternative))  # type: ignore[arg-type]


def test_finite_results_unchanged() -> None:
    # the tie tolerance on finite obs is untouched by the fix
    close(perm_p(2.0, [2.0 - 1e-14, 1.0], "greater"), 2 / 3)
    close(perm_p(-2.0, [-2.0 + 1e-14, 1.0], "less"), 2 / 3)
    close(perm_p(5.5, [1, 2, 3, 4, 5]), 2 / 6)
