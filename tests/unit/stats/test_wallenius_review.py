"""M8-01 review of #153 (B1): extreme odds fail closed instead of returning NaN or aborting."""

from __future__ import annotations

import math

import pytest

from giap.stats.wallenius import ODDS_RANGE, wallenius_tails


@pytest.mark.parametrize("odds", [1e-7, 1e-12, 1e8, 1e300])
def test_odds_outside_the_reliable_range_raise(odds: float) -> None:
    with pytest.raises(ValueError, match="outside"):
        wallenius_tails(0, 1000, 20, 100, odds)


@pytest.mark.parametrize("odds", list(ODDS_RANGE))
def test_the_range_bounds_give_finite_tails(odds: float) -> None:
    lo = max(0, 100 - (1000 - 20))
    for k in (lo, 5, 20):
        over, under = wallenius_tails(k, 1000, 20, 100, odds)
        assert math.isfinite(over) and math.isfinite(under)
        assert 0.0 <= over <= 1.0 and 0.0 <= under <= 1.0


def test_exact_zero_odds_and_single_point_support_skip_the_range_check() -> None:
    assert wallenius_tails(0, 1000, 20, 100, 0.0) == (1.0, 1.0)
    assert wallenius_tails(5, 20, 20, 5, 1e-12) == (1.0, 1.0)  # support [5, 5]
