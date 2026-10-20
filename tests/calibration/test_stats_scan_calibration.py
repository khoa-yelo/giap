"""Fast FWER calibration of the M5-03 window scan (GIAP_SPEC.md §7.3.5(d) steps 4-6, §8.1, §8.4).

Claim (§8.1 "Islands ... maximum statistic | FWER | 0.05"): with an exchangeable null, the
maximum-statistic p-value controls the family-wise error over all windows of both signs, i.e.
P(any window has p_w <= alpha) <= alpha. Each replicate draws z iid N(0, 1) on n TUs, takes B = 199
uniform permutations of z as the null rows (exactly exchangeable with the observed series),
and returns min_w p_w = ``fwer_p(max |S_w|, scan_null_max(Z))``; the level criterion of
``tolerances.yaml`` is checked at every alpha in ``levels`` with R = ``replicates_fast.islands``.
The share of replicates with >= 1 island after ``merge_windows`` equals the same rejection rate.

Teeth: a per-window permutation p-value (each window against the null of its own |S_w|, no
maximum) is anti-conservative over the family and must be a confirmed fail.
"""

from __future__ import annotations

from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.stats.scan import fwer_p, merge_windows, scan_max, scan_null_max

pytestmark = pytest.mark.calibration

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


@pytest.mark.parametrize(("n", "circular"), [(300, True), (300, False), (40, True)])
def test_scan_fwer_level(n: int, circular: bool) -> None:
    r = C.replicates(SUITE, full=False)
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
