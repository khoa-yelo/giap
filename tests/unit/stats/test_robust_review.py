"""Lead decisions after statistician review of PR #31 (docs/specs/M4-04-stats-bootstrap-hl.md).

1. ``hodges_lehmann`` raises ValueError on any non-finite input (NaN, +inf, -inf) on both the direct
   path and the bisection path; when pairwise differences overflow (values near the float64 limit)
   it terminates and either returns the exact pairwise median or raises ValueError, never loops.
2. The direct-enumeration threshold ``_DIRECT_MAX`` is m*n <= 5e4.
"""

from __future__ import annotations

import math
import threading
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

from gippa.rng import for_step
from gippa.stats import robust
from gippa.stats.robust import hodges_lehmann

TIME_LIMIT_S = 10.0
NAN, INF = float("nan"), float("inf")


def run_with_limit(fn: Callable[[], Any], limit: float = TIME_LIMIT_S) -> tuple[str, Any]:
    """Run fn in a daemon thread; ('ok', value) or ('raise', exc). Fails the test on timeout."""
    box: dict[str, Any] = {}

    def target() -> None:
        try:
            box["ok"] = fn()
        except BaseException as exc:
            box["raise"] = exc

    t = threading.Thread(target=target, daemon=True)
    t.start()
    t.join(limit)
    if t.is_alive():
        pytest.fail(f"hodges_lehmann did not terminate within {limit} s (infinite loop?)")
    return next(iter(box.items()))


def samples(m: int, n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = for_step("test.hodges_lehmann_review", seed)
    return rng.normal(1.0, 1.0, m), rng.normal(0.0, 1.0, n)


# ---------------------------------------------------------------------------- 2. threshold


def test_direct_max_is_5e4() -> None:
    assert robust._DIRECT_MAX == 50_000


# ---------------------------------------------------------------------------- 1. non-finite


SIZES = [pytest.param(4, 5, id="direct"), pytest.param(200, 300, id="bisection")]
BAD = [pytest.param(NAN, id="nan"), pytest.param(INF, id="+inf"), pytest.param(-INF, id="-inf")]


@pytest.mark.parametrize(("m", "n"), SIZES)
@pytest.mark.parametrize("bad", BAD)
@pytest.mark.parametrize("where", ["x", "y", "both"])
def test_non_finite_raises(m: int, n: int, bad: float, where: str) -> None:
    assert (m * n > 50_000) == (m == 200)  # the bisection case really is above the threshold
    x, y = samples(m, n, 1)
    if where in ("x", "both"):
        x[m // 2] = bad
    if where in ("y", "both"):
        y[0] = bad
    kind, val = run_with_limit(lambda: hodges_lehmann(x, y))
    assert kind == "raise", f"returned {val!r} instead of raising"
    assert isinstance(val, ValueError)


@pytest.mark.parametrize(("m", "n"), SIZES)
def test_opposite_infinities_raise(m: int, n: int) -> None:
    x, y = samples(m, n, 2)
    x[0], y[-1] = INF, -INF
    kind, val = run_with_limit(lambda: hodges_lehmann(x, y))
    assert kind == "raise" and isinstance(val, ValueError)


# ---------------------------------------------------------------------------- 1. overflow


def _overflowing(m: int, n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Finite inputs near +-1e308 whose differences overflow both ways; most differences small."""
    x, y = samples(m, n, seed)
    x[: m // 10], x[m // 10 : m // 5] = 1.5e308, -1.5e308
    y[: n // 10], y[n // 10 : n // 5] = 1.7e308, -1.7e308
    return x, y


@pytest.mark.parametrize(("m", "n", "seed"), [(200, 300, 3), (250, 401, 4), (4, 5, 5)])
def test_overflowing_differences_terminate(m: int, n: int, seed: int) -> None:
    x, y = _overflowing(m, n, seed)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    with np.errstate(over="ignore", invalid="ignore"):
        expected = float(np.median(np.subtract.outer(x, y)))
    kind, val = run_with_limit(lambda: hodges_lehmann(x, y))
    if kind == "raise":
        assert isinstance(val, ValueError)
    else:
        assert isinstance(val, float)
        assert val == expected or (math.isnan(val) and math.isnan(expected))


@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_all_differences_overflow_one_way_terminate(sign: float) -> None:
    # every x_i - y_j overflows to sign*inf; bounds of the bisection are infinite
    x = np.full(200, sign * 1.6e308) + np.arange(200) * 1e292
    y = np.full(300, -sign * 1.6e308) - np.arange(300) * 1e292
    with np.errstate(over="ignore"):
        expected = float(np.median(np.subtract.outer(x, y)))
    kind, val = run_with_limit(lambda: hodges_lehmann(x, y))
    if kind == "raise":
        assert isinstance(val, ValueError)
    else:
        assert val == expected
