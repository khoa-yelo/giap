"""Window scan, maximum-statistic null, FWER p-values and window merging (issue M5-03).

GIAP_SPEC.md §7.3.5(d) steps 4-6, §7.3.6 ("scan window enumeration with wrap-around"), §8.1, §8.3;
lead decisions 9-10 (2026-10-20): circular windows are capped at n - 1 TUs; windows merge only on
>= 1 shared TU (touching windows stay separate).

Interface assumed (docs/specs/M5-03-stats-scan.md):
- ``ScanResult(start, length, s, max_abs)``: int64 / int64 / float64 arrays and a float; windows
  ordered by length, then start; circular starts in [0, n).
- ``n_windows(n, wmin, wmax, circular)``; ``scan_max(z, wmin, wmax, circular)``;
  ``scan_null_max(Z, wmin, wmax, circular)`` ((B, n) -> (B,)); ``fwer_p(s_abs, null_max)``;
  ``Cluster(direction, best, members)``; ``merge_windows(res, significant, n, circular)``.

When the statistic does not exist:
- n < wmin, or a circular replicon with n - 1 < wmin -> no window: empty arrays, ``max_abs`` NaN.
- NaN or +-inf in z -> ValueError; wmin < 1 or wmin > wmax -> ValueError.
- ``fwer_p`` with no finite null maximum -> NaN.

The oracle below enumerates windows by direct summation, written from the spec text.
"""

from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import giap
from giap.stats.permutation import perm_p
from giap.stats.scan import (
    Cluster,
    ScanResult,
    fwer_p,
    merge_windows,
    n_windows,
    scan_max,
    scan_null_max,
)

SRC = Path(giap.__file__).resolve().parents[1]  # the giap under test (src/ in a checkout)
R2, R3 = math.sqrt(2.0), math.sqrt(3.0)


def oracle_windows(n: int, wmin: int, wmax: int, circular: bool) -> list[tuple[int, int]]:
    cap = n - 1 if circular else n
    out = []
    for w in range(wmin, min(wmax, cap) + 1):
        starts = range(n) if circular else range(n - w + 1)
        out.extend((s, w) for s in starts)
    return out


def oracle_s(z: np.ndarray, start: int, w: int) -> float:
    n = z.size
    return float(sum(z[(start + i) % n] for i in range(w)) / math.sqrt(w))


def make_res(windows: list[tuple[int, int, float]]) -> ScanResult:
    s = np.array([w[2] for w in windows], dtype=np.float64)
    return ScanResult(
        start=np.array([w[0] for w in windows], dtype=np.int64),
        length=np.array([w[1] for w in windows], dtype=np.int64),
        s=s,
        max_abs=float(np.abs(s).max()) if s.size else float("nan"),
    )


def window_of(res: ScanResult, start: int, length: int) -> int:
    hits = np.flatnonzero((res.start == start) & (res.length == length))
    assert hits.size == 1, (start, length)
    return int(hits[0])


# ---- examples 1-7: scan_max -------------------------------------------------------------------


def test_example1_linear_window_count_and_max() -> None:
    z = [0.0, 0.0, 3.0, 3.0, 0.0, 0.0]
    res = scan_max(z, 2, 3, circular=False)
    assert res.start.size == 9 == n_windows(6, 2, 3, circular=False)
    assert list(zip(res.start.tolist(), res.length.tolist(), strict=True)) == oracle_windows(
        6, 2, 3, False
    )
    assert res.max_abs == pytest.approx(6 / R2, abs=1e-6)
    assert res.max_abs == pytest.approx(4.242641, abs=1e-6)
    i = int(np.argmax(np.abs(res.s)))
    assert (int(res.start[i]), int(res.length[i])) == (2, 2)
    want = [oracle_s(np.asarray(z), s, w) for s, w in oracle_windows(6, 2, 3, False)]
    np.testing.assert_allclose(res.s, want, rtol=0, atol=1e-12)


def test_example2_circular_window_count() -> None:
    z = [0.0, 0.0, 3.0, 3.0, 0.0, 0.0]
    res = scan_max(z, 2, 3, circular=True)
    assert res.start.size == 12 == n_windows(6, 2, 3, circular=True)
    assert list(zip(res.start.tolist(), res.length.tolist(), strict=True)) == oracle_windows(
        6, 2, 3, True
    )
    assert res.max_abs == pytest.approx(6 / R2, abs=1e-12)


def test_example3_circular_wraps_to_position_0() -> None:
    z = np.array([3.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    res = scan_max(z, 2, 3, circular=True)
    i = window_of(res, 5, 2)
    assert res.s[i] == pytest.approx(4.242641, abs=1e-6)
    assert res.max_abs == pytest.approx(6 / R2, abs=1e-12)
    assert np.flatnonzero(np.isclose(np.abs(res.s), res.max_abs)).tolist() == [i]
    # the 3-windows (4, 3) and (5, 3) also wrap
    assert res.s[window_of(res, 4, 3)] == pytest.approx(6 / R3, abs=1e-12)
    assert res.s[window_of(res, 5, 3)] == pytest.approx(6 / R3, abs=1e-12)
    cl = merge_windows(res, np.abs(res.s) >= res.max_abs - 1e-9, 6, circular=True)
    assert [(c.direction, c.best) for c in cl] == [("high", i)]


def test_example4_linear_tie_best_is_start_0() -> None:
    z = np.array([3.0, 0.0, 0.0, 0.0, 0.0, 3.0])
    res = scan_max(z, 2, 3, circular=False)
    assert res.max_abs == pytest.approx(2.121320, abs=1e-6)
    a, b = window_of(res, 0, 2), window_of(res, 4, 2)
    assert res.s[a] == pytest.approx(3 / R2, abs=1e-12)
    assert res.s[b] == pytest.approx(3 / R2, abs=1e-12)
    # [0, 1] and [4, 5] share no TU: two islands; ordered by start, the best one first is start 0
    cl = merge_windows(res, np.abs(res.s) >= res.max_abs - 1e-9, 6, circular=False)
    assert [(c.direction, c.best, c.members) for c in cl] == [
        ("high", a, (a,)),
        ("high", b, (b,)),
    ]
    assert int(res.start[cl[0].best]) == 0


def test_example5_low_direction() -> None:
    z = np.array([0.0, -2.0, -2.0, 0.0])
    res = scan_max(z, 2, 2, circular=False)
    assert res.start.tolist() == [0, 1, 2]
    assert res.length.tolist() == [2, 2, 2]
    i = int(np.argmin(res.s))
    assert int(res.start[i]) == 1
    assert res.s[i] == pytest.approx(-2.828427, abs=1e-6)
    assert res.max_abs == pytest.approx(2.828427, abs=1e-6)
    cl = merge_windows(res, [False, True, False], 4, circular=False)
    assert cl == [Cluster(direction="low", best=1, members=(1,))]


@pytest.mark.parametrize("circular", [False, True])
def test_example6_no_window(circular: bool) -> None:
    res = scan_max([1.5], 2, 20, circular=circular)
    assert n_windows(1, 2, 20, circular) == 0
    assert res.start.size == res.length.size == res.s.size == 0
    assert res.start.dtype == np.int64 and res.length.dtype == np.int64
    assert res.s.dtype == np.float64
    assert math.isnan(res.max_abs)


def test_example7_circular_cap_n_minus_1() -> None:
    res = scan_max([1.0, 2.0, 3.0], 2, 20, circular=True)
    assert n_windows(3, 2, 20, circular=True) == 3
    assert res.length.tolist() == [2, 2, 2]
    assert res.start.tolist() == [0, 1, 2]
    np.testing.assert_allclose(res.s, [3 / R2, 5 / R2, 4 / R2], rtol=0, atol=1e-12)
    # linear: the whole replicon is a window
    lin = scan_max([1.0, 2.0, 3.0], 2, 20, circular=False)
    got = list(zip(lin.start.tolist(), lin.length.tolist(), strict=True))
    assert got == [(0, 2), (1, 2), (0, 3)]
    np.testing.assert_allclose(lin.s, [3 / R2, 5 / R2, 6 / R3], rtol=0, atol=1e-12)
    assert lin.max_abs == pytest.approx(5 / R2, abs=1e-12)


def test_circular_n_equal_wmin_has_no_window() -> None:
    # cap n - 1 < wmin: no window covers the whole circle (lead decision 9)
    res = scan_max([1.0, 2.0], 2, 20, circular=True)
    assert res.start.size == 0 and math.isnan(res.max_abs)
    assert n_windows(2, 2, 20, circular=True) == 0
    assert n_windows(2, 2, 20, circular=False) == 1


# ---- example 8: fwer_p --------------------------------------------------------------------------


def test_example8_fwer_p_ties_count() -> None:
    p = fwer_p([6.5, 9.0, 10.0], np.arange(1, 10, dtype=float))
    assert p.dtype == np.float64
    np.testing.assert_allclose(p, [0.4, 0.2, 0.1], rtol=0, atol=1e-15)


def test_fwer_p_equals_perm_p_greater() -> None:
    rng = np.random.default_rng(3)
    m = np.abs(rng.normal(size=99)) * 3
    s = np.concatenate([np.abs(rng.normal(size=20)) * 3, m[:5]])
    got = fwer_p(s, m)
    want = [perm_p(float(v), m, "greater") for v in s]
    np.testing.assert_array_equal(got, want)


def test_fwer_p_tie_tolerance_counts_near_ties() -> None:
    p = fwer_p([2.0], [2.0 - 1e-14, 1.0])
    assert p.tolist() == [pytest.approx(2 / 3)]


def test_fwer_p_drops_nan_null() -> None:
    p = fwer_p([5.0], [np.nan, 6.0, 1.0, np.nan])
    assert p.tolist() == [pytest.approx(2 / 3)]


@pytest.mark.parametrize("null", [[], [np.nan, np.nan]])
def test_fwer_p_no_valid_null_is_nan(null: list[float]) -> None:
    p = fwer_p([1.0, 2.0], null)
    assert p.shape == (2,)
    assert np.isnan(p).all()


def test_fwer_p_empty_windows() -> None:
    p = fwer_p([], [1.0, 2.0])
    assert p.shape == (0,) and p.dtype == np.float64


# ---- examples 9-11: merge_windows --------------------------------------------------------------


def test_example9_merge_by_direction() -> None:
    res = make_res([(2, 2, 4.2), (3, 2, 3.0), (5, 2, 2.9), (3, 2, -3.1)])
    cl = merge_windows(res, [True, True, True, True], 10, circular=False)
    # ordered by the best window's start: high (2, 2), low (3, 2), high (5, 2)
    assert cl == [
        Cluster(direction="high", best=0, members=(0, 1)),
        Cluster(direction="low", best=3, members=(3,)),
        Cluster(direction="high", best=2, members=(2,)),
    ]


def test_example10_chain_is_transitive() -> None:
    res = make_res([(0, 3, 3.0), (2, 3, 3.5), (4, 3, 3.2)])
    cl = merge_windows(res, np.ones(3, bool), 10, circular=False)
    # (0, 3) and (4, 3) share no TU but are linked through (2, 3)
    assert cl == [Cluster(direction="high", best=1, members=(0, 1, 2))]


def test_example11_circular_shared_position_0() -> None:
    res = make_res([(0, 2, 3.0), (9, 2, 3.3)])
    cl = merge_windows(res, [True, True], 10, circular=True)
    assert cl == [Cluster(direction="high", best=1, members=(0, 1))]


def test_touching_windows_stay_separate() -> None:
    # lead decision 10: [0, 2) and [2, 4) touch but share no TU
    res = make_res([(0, 2, 3.0), (2, 2, 3.3)])
    cl = merge_windows(res, [True, True], 10, circular=False)
    assert cl == [
        Cluster(direction="high", best=0, members=(0,)),
        Cluster(direction="high", best=1, members=(1,)),
    ]


def test_circular_touching_across_origin_stays_separate() -> None:
    res = make_res([(8, 2, 3.0), (0, 2, 3.3)])  # {8, 9} and {0, 1}
    cl = merge_windows(res, [True, True], 10, circular=True)
    assert [c.members for c in cl] == [(1,), (0,)]


def test_equal_best_start_high_before_low() -> None:
    # decision 12: the low window is listed first but the high cluster comes first
    res = make_res([(3, 2, -3.1), (3, 3, 3.0), (7, 2, 2.0)])
    cl = merge_windows(res, [True, True, True], 10, circular=False)
    assert cl == [
        Cluster(direction="high", best=1, members=(1,)),
        Cluster(direction="low", best=0, members=(0,)),
        Cluster(direction="high", best=2, members=(2,)),
    ]


def test_zero_s_significant_window_is_dropped() -> None:
    # decision 13: S = 0 has no direction; it joins no cluster, bridges none, raises no error
    res = make_res([(0, 3, 3.0), (2, 3, 0.0), (4, 3, 3.2), (6, 2, -0.0)])
    cl = merge_windows(res, [True, True, True, True], 10, circular=False)
    assert cl == [
        Cluster(direction="high", best=0, members=(0,)),
        Cluster(direction="high", best=2, members=(2,)),
    ]
    only_zero = make_res([(0, 2, 0.0), (1, 2, 0.0)])
    assert merge_windows(only_zero, [True, True], 10, circular=True) == []


def test_opposite_directions_never_merge() -> None:
    res = make_res([(0, 3, 3.0), (1, 3, -3.5)])
    cl = merge_windows(res, [True, True], 10, circular=False)
    assert cl == [
        Cluster(direction="high", best=0, members=(0,)),
        Cluster(direction="low", best=1, members=(1,)),
    ]


def test_not_significant_windows_ignored_and_do_not_bridge() -> None:
    res = make_res([(0, 3, 3.0), (2, 3, 1.0), (4, 3, 3.2)])
    cl = merge_windows(res, [True, False, True], 10, circular=False)
    assert [c.members for c in cl] == [(0,), (2,)]


def test_no_significant_window_gives_no_cluster() -> None:
    res = make_res([(0, 2, 1.0), (1, 2, -1.0)])
    assert merge_windows(res, [False, False], 10, circular=False) == []
    empty = scan_max([1.0], 2, 3, circular=False)
    assert merge_windows(empty, np.zeros(0, bool), 1, circular=False) == []


def test_best_tie_fewer_positions_then_smaller_start() -> None:
    res = make_res([(1, 3, 3.0), (0, 2, 3.0), (1, 2, -9.0), (1, 2, 3.0)])
    cl = merge_windows(res, [True, True, False, True], 10, circular=False)
    # equal |S|: fewer positions wins over (1, 3), then the smaller start (0, 2) over (1, 2)
    assert cl == [Cluster(direction="high", best=1, members=(0, 1, 3))]


def test_best_tie_circular_start_in_range() -> None:
    res = make_res([(9, 2, 3.0), (0, 2, 3.0)])
    cl = merge_windows(res, [True, True], 10, circular=True)
    assert cl == [Cluster(direction="high", best=1, members=(0, 1))]


def test_best_is_largest_abs_s_in_low_cluster() -> None:
    res = make_res([(0, 3, -3.0), (1, 3, -4.0), (2, 2, -3.5)])
    cl = merge_windows(res, [True, True, True], 10, circular=False)
    assert cl == [Cluster(direction="low", best=1, members=(0, 1, 2))]


def test_significant_length_mismatch_raises() -> None:
    res = make_res([(0, 2, 3.0), (1, 2, 3.0)])
    with pytest.raises(ValueError):
        merge_windows(res, [True], 10, circular=False)


# ---- edge cases ---------------------------------------------------------------------------------


@pytest.mark.parametrize("circular", [False, True])
def test_all_zero(circular: bool) -> None:
    res = scan_max(np.zeros(8), 2, 4, circular)
    assert res.start.size == n_windows(8, 2, 4, circular) > 0
    assert (res.s == 0).all()
    assert res.max_abs == 0.0
    null = scan_null_max(np.random.default_rng(0).normal(size=(5, 8)), 2, 4, circular)
    np.testing.assert_array_equal(fwer_p(np.abs(res.s), null), np.ones(res.s.size))
    np.testing.assert_array_equal(fwer_p(np.abs(res.s), np.zeros(5)), np.ones(res.s.size))


@pytest.mark.parametrize(("wmin", "wmax"), [(3, 2), (0, 5), (-1, 2), (0, 0)])
def test_bad_window_range_raises(wmin: int, wmax: int) -> None:
    with pytest.raises(ValueError):
        scan_max(np.ones(10), wmin, wmax, circular=False)
    with pytest.raises(ValueError):
        scan_max(np.ones(10), wmin, wmax, circular=True)
    with pytest.raises(ValueError):
        scan_null_max(np.ones((2, 10)), wmin, wmax, circular=False)
    with pytest.raises(ValueError):
        n_windows(10, wmin, wmax, circular=False)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_nonfinite_z_raises(bad: float) -> None:
    z = np.array([0.0, 1.0, bad, 2.0])
    with pytest.raises(ValueError):
        scan_max(z, 2, 3, circular=False)
    with pytest.raises(ValueError):
        scan_null_max(np.vstack([np.zeros(4), z]), 2, 3, circular=True)


def test_wmin_1_includes_single_positions() -> None:
    z = np.array([1.0, -4.0, 2.0])
    res = scan_max(z, 1, 1, circular=False)
    np.testing.assert_array_equal(res.s, z)
    assert res.max_abs == 4.0


def test_accepts_list_and_int_input() -> None:
    res = scan_max([1, 2, 3, 4], 2, 2, circular=False)
    assert res.s.dtype == np.float64
    np.testing.assert_allclose(res.s, [3 / R2, 5 / R2, 7 / R2], rtol=0, atol=1e-12)


def test_result_is_frozen() -> None:
    res = scan_max([1.0, 2.0, 3.0], 2, 2, circular=False)
    with pytest.raises((AttributeError, TypeError)):
        res.max_abs = 0.0  # type: ignore[misc]
    c = Cluster(direction="high", best=0, members=(0,))
    with pytest.raises((AttributeError, TypeError)):
        c.best = 1  # type: ignore[misc]


# ---- n_windows equals the window count ----------------------------------------------------------


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize(("wmin", "wmax"), [(1, 1), (1, 5), (2, 3), (2, 20), (5, 5), (7, 60)])
def test_n_windows_matches_enumeration(circular: bool, wmin: int, wmax: int) -> None:
    rng = np.random.default_rng(wmin * 100 + wmax)
    for n in range(1, 51):
        want = oracle_windows(n, wmin, wmax, circular)
        assert n_windows(n, wmin, wmax, circular) == len(want), n
        res = scan_max(rng.normal(size=n), wmin, wmax, circular)
        assert res.start.size == res.length.size == res.s.size == len(want), n
        got = list(zip(res.start.tolist(), res.length.tolist(), strict=True))
        assert got == want, n


def test_n_windows_zero_replicon() -> None:
    assert n_windows(0, 2, 20, circular=False) == 0
    assert n_windows(0, 2, 20, circular=True) == 0


# ---- scan_null_max ------------------------------------------------------------------------------


@pytest.mark.parametrize("circular", [False, True])
def test_scan_null_max_rowwise(circular: bool) -> None:
    Z = np.random.default_rng(11).normal(size=(40, 37))
    got = scan_null_max(Z, 2, 20, circular)
    assert got.shape == (40,) and got.dtype == np.float64
    want = [scan_max(row, 2, 20, circular).max_abs for row in Z]
    np.testing.assert_allclose(got, want, rtol=0, atol=1e-9)


def test_scan_null_max_no_window_is_nan() -> None:
    got = scan_null_max(np.ones((3, 1)), 2, 20, circular=False)
    assert got.shape == (3,) and np.isnan(got).all()


def test_scan_null_max_no_rows() -> None:
    got = scan_null_max(np.zeros((0, 10)), 2, 5, circular=True)
    assert got.shape == (0,) and got.dtype == np.float64


def test_planted_island_detected() -> None:
    rng = np.random.default_rng(2026)
    n = 400
    z = rng.normal(size=n)
    z[100:110] += 2.5
    Z = np.array([rng.permutation(z) for _ in range(199)])
    res = scan_max(z, 2, 20, circular=True)
    p = fwer_p(np.abs(res.s), scan_null_max(Z, 2, 20, circular=True))
    cl = merge_windows(res, p <= 0.05, n, circular=True)
    assert len(cl) >= 1
    top = max(cl, key=lambda c: abs(res.s[c.best]))
    assert top.direction == "high"
    b = top.best
    covered = {(int(res.start[b]) + i) % n for i in range(int(res.length[b]))}
    assert covered & set(range(100, 110))


# ---- numba on and off give identical results ---------------------------------------------------

_SCRIPT = """
import json, sys
import numpy as np
from giap.stats.scan import scan_max, scan_null_max
rng = np.random.default_rng(5)
out = []
for circular in (False, True):
    for n in (1, 2, 3, 7, 30, 101):
        z = rng.normal(size=n)
        r = scan_max(z, 2, 20, circular)
        Z = rng.normal(size=(6, n))
        out.append([r.start.tolist(), r.length.tolist(), r.s.tolist(), repr(r.max_abs),
                    [repr(v) for v in scan_null_max(Z, 2, 20, circular)]])
print(json.dumps(out))
"""


def _run(disable_jit: bool) -> str:
    env = dict(os.environ)
    env.pop("NUMBA_DISABLE_JIT", None)
    if disable_jit:
        env["NUMBA_DISABLE_JIT"] = "1"
    env["PYTHONPATH"] = os.pathsep.join(p for p in (str(SRC), env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT],
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-3000:]
    return proc.stdout


def test_numba_disabled_identical() -> None:
    assert _run(disable_jit=True) == _run(disable_jit=False)


# ---- decision 17: non-finite |S| in fwer_p -------------------------------------------------------


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_fwer_p_nonfinite_s_is_nan(bad: float) -> None:
    p = fwer_p([1.5, bad, 0.5], [1.0, 2.0, 3.0])
    assert p.shape == (3,)
    assert np.isnan(p[1])
    np.testing.assert_allclose(p[[0, 2]], [3 / 4, 1.0], rtol=0, atol=1e-15)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_fwer_p_nonfinite_s_is_nan_with_infinite_null(bad: float) -> None:
    p = fwer_p([bad], [1.0, np.inf, -np.inf, np.nan])
    assert p.shape == (1,) and np.isnan(p[0])


def test_fwer_p_nan_max_abs_of_windowless_replicon() -> None:
    res = scan_max([1.0], 2, 3, circular=False)
    assert math.isnan(res.max_abs)
    p = fwer_p([res.max_abs], [1.0, 2.0])
    assert p.shape == (1,) and np.isnan(p[0])
    p = fwer_p(np.abs(res.s), [1.0, 2.0])
    assert p.shape == (0,)


# ---- decision 18: input checks -----------------------------------------------------------------


@pytest.mark.parametrize("bad", [np.float64(1.0), np.ones((2, 2))], ids=["0d", "2d"])
def test_fwer_p_rejects_non_1d(bad: np.ndarray) -> None:
    with pytest.raises(ValueError):
        fwer_p(bad, [1.0, 2.0])
    with pytest.raises(ValueError):
        fwer_p([1.0, 2.0], bad)


def test_merge_windows_rejects_non_1d_significant() -> None:
    res = make_res([(0, 2, 3.0), (1, 2, 3.0)])
    with pytest.raises(ValueError):
        merge_windows(res, np.ones((1, 2), bool), 10, circular=False)
    with pytest.raises(ValueError):
        merge_windows(res, np.ones((2, 1), bool), 10, circular=False)
    with pytest.raises(ValueError):
        merge_windows(make_res([(0, 2, 3.0)]), np.bool_(True), 10, circular=False)


@pytest.mark.parametrize("circular", [False, True])
def test_merge_windows_rejects_negative_n(circular: bool) -> None:
    empty = scan_max(np.zeros(0), 2, 3, circular)
    with pytest.raises(ValueError):
        merge_windows(empty, np.zeros(0, bool), -1, circular)


@pytest.mark.parametrize("circular", [False, True])
def test_merge_windows_n0(circular: bool) -> None:
    res = scan_max(np.zeros(0), 2, 3, circular)
    assert res.s.size == 0 and math.isnan(res.max_abs)
    assert merge_windows(res, np.zeros(0, bool), 0, circular) == []
    assert merge_windows(res, [], 0, circular) == []


@pytest.mark.parametrize("bad", [np.float64(1.0), np.ones((2, 5))], ids=["0d", "2d"])
@pytest.mark.parametrize("circular", [False, True])
def test_scan_max_rejects_non_1d(bad: np.ndarray, circular: bool) -> None:
    with pytest.raises(ValueError):
        scan_max(bad, 1, 3, circular)


@pytest.mark.parametrize("bad", [np.ones(5), np.ones((2, 5, 1))], ids=["1d", "3d"])
@pytest.mark.parametrize("circular", [False, True])
def test_scan_null_max_rejects_non_2d(bad: np.ndarray, circular: bool) -> None:
    with pytest.raises(ValueError):
        scan_null_max(bad, 1, 3, circular)


# ---- decision 22: a negative finite |S| is a caller error --------------------------------------


@pytest.mark.parametrize("bad", [-5.0, -1e-300, -np.finfo(float).max])
def test_fwer_p_negative_s_abs_raises(bad: float) -> None:
    with pytest.raises(ValueError):
        fwer_p([1.0, bad], [1.0, 2.0])
    with pytest.raises(ValueError):
        fwer_p([bad], [np.nan])  # checked even when no null maximum is valid


def test_fwer_p_signed_s_raises() -> None:
    res = scan_max([0.0, -2.0, -2.0, 0.0], 2, 2, circular=False)
    with pytest.raises(ValueError):
        fwer_p(res.s, [1.0, 2.0])


def test_fwer_p_negative_zero_and_nonfinite_do_not_raise() -> None:
    p = fwer_p([-0.0, 0.0, np.nan, -np.inf, np.inf], [1.0, 2.0])
    np.testing.assert_array_equal(p[:2], [1.0, 1.0])
    assert np.isnan(p[2:]).all()  # decision 19: non-finite |S| -> NaN, not an error
