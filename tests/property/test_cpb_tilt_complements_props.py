"""Hypothesis properties of the complement-keeping tilt (M8-04a review B1; docs/specs/M8-04 lead
decision 11).

- contract: for random q with entries within 1e-15..1e-9 of 0 and of 1 (plus mid values and exact
  0s and 1s), ``conditional_poisson_binomial`` on ``tilt_pair(q, N)`` (split into C and its
  complement, complements passed) equals the UNTILTED exact rational law of the raw q to 1e-11
  relative on p_over, p_under and mean, with the same support, and never raises; the same holds
  for a total one away from N;
- ``tilt_pair``: comp = expit(-(logit q + theta)) to 1e-12 relative against a 60-digit reference
  built from the exact rational q (never 1 - p); exact 0s and 1s kept with comp 1 or 0; every
  open entry has p, comp >= the smallest normal; p + comp = 1 within rounding;
- ``tilt`` is ``tilt_pair(...)[0]`` exactly when it does not raise, and raises exactly when an
  open entry of ``tilt_pair`` is 1.0 in floating point or clamped at the smallest normal from an
  underflow to 0.
"""

from __future__ import annotations

import numpy as np
from helpers.cpb_oracle import assert_tail_close, cpb_exact, rel_err
from helpers.cpb_pair_oracle import dec_rel_err, tilted_reference
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from gippa.stats.poisson_binomial import conditional_poisson_binomial, tilt, tilt_pair

CONTRACT = 1e-11
RTOL = 1e-12
TINY = float(np.finfo(np.float64).tiny)
EPS = float(np.finfo(np.float64).eps)
FAST = settings(max_examples=150, deadline=None)
SLOW = settings(max_examples=60, deadline=None)

MID = st.floats(0.01, 0.99, allow_nan=False)
NEAR1 = st.floats(9.0, 15.0).map(lambda u: 1.0 - 10.0**-u)
NEAR0 = st.floats(9.0, 15.0).map(lambda u: 10.0**-u)
Q = st.one_of(MID, MID, NEAR1, NEAR0, st.sampled_from([0.0, 1.0]))
# tilt_pair alone (no oracle DP): also the extremes that underflow or saturate
EXTREME = st.sampled_from([1e-300, 5e-324, 1.0 - 2.0**-53, 1.0 - 1e-15])
QX = st.one_of(MID, MID, NEAR1, NEAR0, EXTREME, st.sampled_from([0.0, 1.0]))


@st.composite
def universes(draw: st.DrawFn, max_t: int = 24) -> tuple[np.ndarray, int, int, int]:
    """(q, n_in, k, N) with S = N of positive probability under the raw q and k in its support."""
    q = np.asarray(draw(st.lists(Q, min_size=1, max_size=max_t)), dtype=float)
    n_in = draw(st.integers(0, q.size))
    ones, pos = int((q == 1).sum()), int((q > 0).sum())
    N = draw(st.integers(ones, pos))
    a_in, a_out = q[:n_in], q[n_in:]
    lo = max(int((a_in == 1).sum()), N - int((a_out > 0).sum()), 0)
    hi = min(int((a_in > 0).sum()), N - int((a_out == 1).sum()), N)
    assume(lo <= hi)
    k = draw(st.integers(lo, hi))
    return q, n_in, k, N


def _check(q: np.ndarray, n_in: int, k: int, N: int, total: int) -> None:
    p, c = tilt_pair(q, total)
    t = conditional_poisson_binomial(  # never raises: the raw law has a non-empty support
        p[:n_in], p[n_in:], k, N, comp_in=c[:n_in], comp_out=c[n_in:]
    )
    ex = cpb_exact(list(q[:n_in]), list(q[n_in:]), k, N)  # the exact law of the RAW q
    assert ex is not None
    assert (t.lo, t.hi) == (ex.lo, ex.hi)
    assert_tail_close(t.p_over, ex.p_over, rtol=CONTRACT)
    assert_tail_close(t.p_under, ex.p_under, rtol=CONTRACT)
    assert rel_err(t.mean, ex.mean) <= CONTRACT, (t.mean, float(ex.mean))


@SLOW
@given(universes())
def test_tilt_pair_to_N_gives_the_untilted_exact_law(u: tuple) -> None:
    q, n_in, k, N = u
    _check(q, n_in, k, N, N)


@SLOW
@given(universes(), st.sampled_from([-1, 1]))
def test_tilt_pair_to_a_nearby_total_gives_the_untilted_exact_law(u: tuple, d: int) -> None:
    q, n_in, k, N = u
    ones, pos = int((q == 1).sum()), int((q > 0).sum())
    total = N + d
    assume(ones < total < pos)  # interior: a finite theta
    _check(q, n_in, k, N, total)


@st.composite
def tilts(draw: st.DrawFn) -> tuple[np.ndarray, int]:
    q = np.asarray(draw(st.lists(QX, min_size=2, max_size=200)), dtype=float)
    ones, pos = int((q == 1).sum()), int((q > 0).sum())
    assume(pos - ones >= 2)
    m = (q > 0) & (q < 1)
    assume(np.any((q[m] >= 0.01) & (q[m] <= 0.99)))  # a well-conditioned entry anchors theta
    return q, draw(st.integers(ones + 1, pos - 1))


@FAST
@given(tilts())
def test_tilt_pair_matches_the_exact_tilt(qt: tuple) -> None:
    q, total = qt
    p, c = tilt_pair(q, total)
    assert p.shape == q.shape and c.shape == q.shape
    assert p.dtype == np.float64 and c.dtype == np.float64
    assert abs(p.sum() - total) <= 1e-9 * total
    np.testing.assert_array_equal(p[q == 0], 0.0)
    np.testing.assert_array_equal(c[q == 0], 1.0)
    np.testing.assert_array_equal(p[q == 1], 1.0)
    np.testing.assert_array_equal(c[q == 1], 0.0)
    m = (q > 0) & (q < 1)
    assert np.all(p[m] >= TINY) and np.all(c[m] >= TINY)
    assert np.all(p <= 1.0) and np.all(c <= 1.0)
    assert np.all(np.abs(p + c - 1.0) <= 4 * EPS)
    ref_p, ref_c = tilted_reference(q, p, c)
    for i in np.flatnonzero(m):
        rp, rc = ref_p[i], ref_c[i]
        assert rp is not None and rc is not None
        for got, ref in ((float(p[i]), rp), (float(c[i]), rc)):
            if ref >= 1e-300:
                assert dec_rel_err(got, ref) <= RTOL, (float(q[i]), got, float(ref))
            else:
                assert TINY <= got <= 1e-300


@FAST
@given(tilts())
def test_tilt_is_tilt_pair_or_raises_on_saturation(qt: tuple) -> None:
    q, total = qt
    p, c = tilt_pair(q, total)
    m = (q > 0) & (q < 1)
    # an open entry that a single array cannot keep open: 1.0, or an underflow to 0 (clamped)
    saturated = bool(np.any(p[m] == 1.0))
    clamped = bool(np.any(p[m] == TINY))
    try:
        got = tilt(q, total)
    except ValueError:
        assert saturated or clamped, (p[m].min(), c[m].min())
        return
    assert not saturated
    np.testing.assert_array_equal(got, p)
