"""FWER calibration of the M5-03 window scan (GIPPA_SPEC.md §7.3.5(d) steps 4-6, §8.1, §8.4).

Claim (§8.1 "Islands ... maximum statistic | FWER | 0.05"): with an exchangeable null, the
maximum-statistic p-value controls the family-wise error over all windows of both signs, i.e.
P(any window has p_w <= alpha) <= alpha. Each replicate draws z iid N(0, 1) on n TUs, takes B = 199
uniform permutations of z as the null rows (exactly exchangeable with the observed series),
and returns min_w p_w = ``fwer_p(max |S_w|, scan_null_max(Z))``; the level criterion of
``tolerances.yaml`` is checked at every alpha in ``levels`` with R = ``replicates_fast.islands``
(fast arm) or ``replicates_full.islands`` = 500 (``calibration_full`` arm, lead decision 21).
The share of replicates with >= 1 island after ``merge_windows`` equals the same rejection rate.

Planted island (decision 21; statistician's review of PR #118, item 3): one true island of L TUs
shifted by ``eff`` SD at a random start, null rows = permutations of the observed series. A
cluster is false when none of its member windows overlaps the true island. Every replicate checks
the implication "false cluster => some significant window is outside the island" (so the
false-island rate is <= the window FWER), and the false-island indicator goes through the ``fwer``
criterion at every alpha in ``levels``.

Teeth: a per-window permutation p-value (each window against the null of its own |S_w|, no
maximum) is anti-conservative over the family and must be a confirmed fail.
"""

from __future__ import annotations

from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from gippa.stats.scan import ScanResult, fwer_p, merge_windows, scan_max, scan_null_max

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "islands"
B = 199
WMIN, WMAX = 2, 20


def _min_p(rng: np.random.Generator, n: int, circular: bool) -> float:
    z = rng.normal(size=n)
    Z = np.array([rng.permutation(z) for _ in range(B)])
    res = scan_max(z, WMIN, WMAX, circular)
    null = scan_null_max(Z, WMIN, WMAX, circular)
    assert null.shape == (B,) and np.isfinite(null).all()
    p = fwer_p(np.abs(res.s), null)
    assert p.shape == res.s.shape and ((p > 0) & (p <= 1)).all()
    pmin = float(fwer_p([res.max_abs], null)[0])
    assert pmin == float(p.min())
    # islands exist iff the best window is significant
    for alpha in C.load_tolerances().levels:
        assert bool(merge_windows(res, p <= alpha, n, circular)) == (pmin <= alpha)
    return pmin


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize(("n", "circular"), [(300, True), (300, False), (40, True)])
def test_scan_fwer_level(n: int, circular: bool, full: bool) -> None:
    r = C.replicates(SUITE, full=full)
    levels = C.load_tolerances().levels
    topo = "circular" if circular else "linear"
    calibrate_levels(
        SUITE,
        "level",
        levels,
        lambda rng: _min_p(rng, n, circular),
        r,
        stream=f"{SUITE}.stats.scan.{topo}.{n}",
    ).assert_ok()


def _covered(res: ScanResult, i: int, n: int) -> set[int]:
    return {(int(res.start[i]) + k) % n for k in range(int(res.length[i]))}


def _false_island(
    rng: np.random.Generator, n: int, circular: bool, length: int, eff: float
) -> np.ndarray:
    levels = C.load_tolerances().levels
    z = rng.normal(size=n)
    first = int(rng.integers(0, n if circular else n - length + 1))
    island_idx = (first + np.arange(length)) % n
    z[island_idx] += eff
    island = set(island_idx.tolist())
    Z = np.array([rng.permutation(z) for _ in range(B)])
    res = scan_max(z, WMIN, WMAX, circular)
    p = fwer_p(np.abs(res.s), scan_null_max(Z, WMIN, WMAX, circular))
    out = []
    for alpha in levels:
        sig = p <= alpha
        window_false = any(not (_covered(res, i, n) & island) for i in np.flatnonzero(sig))
        clusters = merge_windows(res, sig, n, circular)
        cluster_false = any(
            all(not (_covered(res, i, n) & island) for i in c.members) for c in clusters
        )
        assert window_false or not cluster_false  # merging never adds a false island
        out.append(float(cluster_false))
    return np.array(out)


PLANTED = [(300, True, 5, 1.5), (300, True, 20, 1.0), (300, True, 40, 1.0), (300, False, 20, 1.0)]


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize(("n", "circular", "length", "eff"), PLANTED)
def test_planted_island_false_island_rate(
    n: int, circular: bool, length: int, eff: float, full: bool
) -> None:
    topo = "circular" if circular else "linear"
    calibrate_levels(
        SUITE,
        "fwer",
        C.load_tolerances().levels,
        lambda rng: _false_island(rng, n, circular, length, eff),
        C.replicates(SUITE, full=full),
        stream=f"{SUITE}.stats.scan.planted.{topo}.{n}.{length}x{eff:g}",
    ).assert_ok()


@pytest.mark.calibration
def test_teeth_per_window_p_is_not_fwer(teeth_log: Path) -> None:
    """Without the maximum, the smallest per-window p-value is far below alpha too often."""
    n = 300

    def rep(rng: np.random.Generator) -> float:
        z = rng.normal(size=n)
        Z = np.array([rng.permutation(z) for _ in range(B)])
        res = scan_max(z, WMIN, WMAX, circular=True)
        # per-window: the null of each window is its own |S| across rows (no maximum)
        null_s = np.array([np.abs(scan_max(row, WMIN, WMAX, circular=True).s) for row in Z])
        p = (1 + (null_s >= np.abs(res.s)[None, :]).sum(axis=0)) / (1 + B)
        return float(p.min())

    rec = C.calibrate(
        SUITE,
        "level",
        0.05,
        rep,
        C.replicates(SUITE, full=False),
        stream=f"{SUITE}.stats.scan.teeth_per_window",
        log_dir=teeth_log,
    )
    assert rec["status"] == "fail", rec
