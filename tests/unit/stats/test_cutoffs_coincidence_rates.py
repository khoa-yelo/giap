"""Per-placement coincidence rates (PR #264 review, blocking item 1; lead decision of 2026-10-23).

The K-12 driver used 20 calls of ``coincidence_rate(R=R/20)`` and a batch-means Monte Carlo SE.
The lead replaced that design: part A returns the per-placement rates, so one call at ``R`` gives
c = their mean and the MC SE = their sample sd (ddof = 1) / sqrt(R). Placements are iid, so this is
the plain SE of a mean.

Interface pinned here (``giap.stats.cutoffs``):

- ``coincidence_rate(calls, confidence, windows, truth, replicon, threshold, *, rng, R=200,
  tol_bp=50, return_rates=False)``. With ``return_rates=False`` (the default) it is unchanged and
  returns a float. With ``return_rates=True`` it returns a float64 ``numpy.ndarray`` of shape
  ``(R,)``. Element r is the matched share of the kept calls (confidence >= threshold) in
  placement r, so each element is a multiple of 1 / n_kept in [0, 1].
- The draws are the same in both modes. For the same rng state,
  ``coincidence_rate(..., return_rates=True).mean()`` equals ``coincidence_rate(...)``, and both
  modes leave the generator in the same state.
- No call at or above the threshold gives ``np.full(R, nan)``, matching the float mode's NaN, and
  draws nothing.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats.cutoffs import (
    DEFAULT_SEED,
    STEP,
    Window,
    coincidence_rate,
    downstream_windows,
)

C = "chr1"
LIN = Replicon(C, 10_000, "linear")


def gene(tag: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(C, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def iv(start: int, end: int, strand: str) -> Interval:
    return Interval(C, start, end, strand)  # type: ignore[arg-type]


def rng() -> np.random.Generator:
    return for_step(STEP, DEFAULT_SEED)


# + windows: gA [950, 1500), gB [2350, 2950). - windows: gC [2400, 3000), h1 [4000, 5050),
# h2 [6000, 6550) (the same toy as test_cutoffs_terminators.py).
NULL_GENES = [
    gene("gA", 100, 1000, "+"),
    gene("gB", 1500, 2400, "+"),
    gene("gC", 2950, 4000, "-"),
    gene("h1", 5000, 6000, "-"),
    gene("h2", 6500, 7000, "-"),
]


def windows() -> list[Window]:
    return downstream_windows(NULL_GENES, LIN)


CALLS = [iv(1000, 1030, "+"), iv(5000, 5020, "-"), iv(1130, 1160, "+")]
CONF = [0.9, 0.8, 0.7]
TRUTH = [iv(1000, 1030, "+"), iv(5000, 5020, "-"), iv(2530, 2560, "+")]


def test_rates_shape_dtype_and_grid() -> None:
    rates = coincidence_rate(
        CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=rng(), R=40, return_rates=True
    )
    assert isinstance(rates, np.ndarray)
    assert rates.dtype == np.float64
    assert rates.shape == (40,)
    assert np.all((rates >= 0.0) & (rates <= 1.0))
    k = rates * 3  # three kept calls: every rate is k / 3
    assert np.allclose(k, np.round(k), atol=1e-12)


def test_rates_mean_is_the_float_mode_and_the_rng_ends_in_the_same_state() -> None:
    g1, g2 = rng(), rng()
    rates = coincidence_rate(
        CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=g1, R=40, return_rates=True
    )
    c = coincidence_rate(CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=g2, R=40)
    assert isinstance(c, float)
    assert float(np.mean(rates)) == pytest.approx(c, abs=1e-15)
    assert g1.random() == g2.random()  # same number of draws in both modes


def test_rates_default_is_the_float() -> None:
    """``return_rates`` defaults to False: existing callers keep the float."""
    c = coincidence_rate(
        CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=rng(), R=10, return_rates=False
    )
    assert isinstance(c, float)
    assert c == coincidence_rate(CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=rng(), R=10)


def test_rates_single_call_two_windows_are_zero_or_one() -> None:
    """One + call on its truth in gA: placement r lands in gA (rate 1) or gB (rate 0)."""
    w = windows()
    rates = coincidence_rate(
        [iv(1000, 1030, "+")], [0.9], w, [iv(1000, 1030, "+")], LIN, 0.5, rng=rng(), R=200,
        return_rates=True,
    )  # fmt: skip
    assert set(np.unique(rates)) == {0.0, 1.0}
    assert rates.mean() == pytest.approx(0.5, abs=0.12)


def test_rates_constant_when_only_one_window_fits() -> None:
    """Offset 560 + length 30 > 550: gA cannot hold the call; every placement is gB's match."""
    rates = coincidence_rate(
        [iv(2910, 2940, "+")], [0.9], windows(), [iv(2920, 2925, "+")], LIN, 0.5, rng=rng(), R=20,
        return_rates=True,
    )  # fmt: skip
    assert np.array_equal(rates, np.ones(20))


def test_rates_r_one() -> None:
    rates = coincidence_rate(
        CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=rng(), R=1, return_rates=True
    )
    assert rates.shape == (1,)


def test_rates_no_call_at_threshold_is_all_nan_and_draws_nothing() -> None:
    g, ref = rng(), rng()
    rates = coincidence_rate(
        CALLS, CONF, windows(), TRUTH, LIN, 0.95, rng=g, R=7, return_rates=True
    )
    assert rates.shape == (7,)
    assert np.all(np.isnan(rates))
    assert g.random() == ref.random()


def test_mc_se_from_one_call_is_sd_over_root_r() -> None:
    """The driver's c and MC SE come from one call; for the 3-call toy the SE is positive and
    shrinks like 1/sqrt(R) (R = 400 against R = 25: the ratio is about 4)."""
    se = {}
    for r in (25, 400):
        rates = coincidence_rate(
            CALLS, CONF, windows(), TRUTH, LIN, 0.0, rng=rng(), R=r, return_rates=True
        )
        se[r] = float(np.std(rates, ddof=1) / math.sqrt(r))
    assert se[25] > 0.0
    assert se[25] / se[400] == pytest.approx(4.0, rel=0.35)


@settings(max_examples=40, deadline=None)
@given(
    picks=st.lists(
        st.tuples(st.sampled_from([0, 1, 2, 3]), st.integers(0, 200), st.floats(0.0, 1.0)),
        min_size=1,
        max_size=6,
    ),
    t=st.floats(0.0, 1.0),
    r=st.integers(1, 30),
    seed=st.integers(0, 2**31 - 1),
)
def test_rates_property_mean_grid_and_draws(
    picks: list[tuple[int, int, float]], t: float, r: int, seed: int
) -> None:
    """For any eligible calls: the rates' mean is the float mode, each rate is k / n_kept and both
    modes consume the same draws."""
    w = windows()
    slots = [w[0], w[1], w[3], w[4]]  # gA, gB (+), h1, h2 (-)
    calls = []
    for slot, off, _ in picks:
        win = slots[slot]
        if win.strand == "+":
            calls.append(iv(win.lo + off, win.lo + off + 30, "+"))
        else:
            calls.append(iv(win.hi - off - 30, win.hi - off, "-"))
    conf = [p for _, _, p in picks]
    truth = [iv(1000, 1030, "+"), iv(2400, 2440, "+"), iv(4900, 5000, "-"), iv(6400, 6450, "-")]
    g1, g2 = for_step(STEP, seed), for_step(STEP, seed)
    rates = coincidence_rate(calls, conf, w, truth, LIN, t, rng=g1, R=r, return_rates=True)
    c = coincidence_rate(calls, conf, w, truth, LIN, t, rng=g2, R=r)
    assert rates.shape == (r,)
    n_kept = sum(p >= t for p in conf)
    if n_kept == 0:
        assert math.isnan(c) and np.all(np.isnan(rates))
    else:
        assert float(rates.mean()) == pytest.approx(c, abs=1e-12)
        k = rates * n_kept
        assert np.allclose(k, np.round(k), atol=1e-9)
    assert g1.random() == g2.random()
