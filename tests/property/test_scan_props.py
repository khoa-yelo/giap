"""Hypothesis properties of the window scan (issue M5-03; GIAP_SPEC.md §7.3.5(d) steps 4-6, §8.3).

- S on random z equals the direct sum over the window (wrapping on circular input) divided by
  sqrt(n_w), within 1e-9 * max(1, sum |z|); windows are exactly the enumerated set, in order.
- Circular ``max_abs`` is invariant under rotation of z; linear ``max_abs`` under reversal; both
  under negation (the scan is two-sided).
- ``scan_null_max`` equals row-wise ``scan_max(...).max_abs``.
- ``fwer_p`` is monotone non-increasing in |S|, lies in [1/(1+B), 1] and equals the +1 rule.
- ``merge_windows``: clusters partition the significant windows; members share the direction of
  the cluster; same-direction windows sharing a TU are in one cluster; each cluster is connected;
  the best window maximises |S| among its members; clusters are ordered by the best start.
"""

from __future__ import annotations

import math

import numpy as np
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.stats.scan import fwer_p, merge_windows, n_windows, scan_max, scan_null_max

SETTINGS = settings(max_examples=100, deadline=None, suppress_health_check=[HealthCheck.too_slow])

VAL = st.floats(-50.0, 50.0, allow_nan=False, allow_infinity=False)


@st.composite
def scan_case(draw: st.DrawFn, max_n: int = 60) -> tuple[np.ndarray, int, int, bool]:
    n = draw(st.integers(1, max_n))
    z = np.array(draw(st.lists(VAL, min_size=n, max_size=n)), dtype=float)
    wmin = draw(st.integers(1, 8))
    wmax = draw(st.integers(wmin, 25))
    circular = draw(st.booleans())
    return z, wmin, wmax, circular


def positions(start: int, length: int, n: int, circular: bool) -> set[int]:
    if circular:
        return {(start + i) % n for i in range(length)}
    return set(range(start, start + length))


@SETTINGS
@given(scan_case())
def test_s_equals_direct_sum(case: tuple[np.ndarray, int, int, bool]) -> None:
    z, wmin, wmax, circular = case
    n = z.size
    res = scan_max(z, wmin, wmax, circular)
    cap = n - 1 if circular else n
    want = [
        (s, w)
        for w in range(wmin, min(wmax, cap) + 1)
        for s in (range(n) if circular else range(n - w + 1))
    ]
    assert list(zip(res.start.tolist(), res.length.tolist(), strict=True)) == want
    assert n_windows(n, wmin, wmax, circular) == len(want)
    tol = 1e-9 * max(1.0, float(np.abs(z).sum()))
    direct = np.array([z[(s + np.arange(w)) % n].sum() / math.sqrt(w) for s, w in want])
    if want:
        np.testing.assert_allclose(res.s, direct, rtol=0, atol=tol)
        assert abs(res.max_abs - float(np.abs(direct).max())) <= tol
    else:
        assert math.isnan(res.max_abs)


@SETTINGS
@given(scan_case(), st.integers(0, 10_000))
def test_circular_rotation_invariance(case: tuple[np.ndarray, int, int, bool], k: int) -> None:
    z, wmin, wmax, _ = case
    tol = 1e-9 * max(1.0, float(np.abs(z).sum()))
    a = scan_max(z, wmin, wmax, circular=True).max_abs
    b = scan_max(np.roll(z, k % z.size), wmin, wmax, circular=True).max_abs
    if math.isnan(a):
        assert math.isnan(b)
    else:
        assert abs(a - b) <= tol


@SETTINGS
@given(scan_case())
def test_reversal_and_negation_invariance(case: tuple[np.ndarray, int, int, bool]) -> None:
    z, wmin, wmax, circular = case
    tol = 1e-9 * max(1.0, float(np.abs(z).sum()))
    a = scan_max(z, wmin, wmax, circular).max_abs
    rev = scan_max(z[::-1], wmin, wmax, circular).max_abs
    neg = scan_max(-z, wmin, wmax, circular).max_abs
    if math.isnan(a):
        assert math.isnan(rev) and math.isnan(neg)
    else:
        assert abs(a - rev) <= tol
        assert abs(a - neg) <= tol


@SETTINGS
@given(
    st.integers(0, 12),
    st.integers(1, 40),
    st.integers(1, 6),
    st.integers(0, 20),
    st.booleans(),
    st.integers(0, 2**32 - 1),
)
def test_null_max_equals_rowwise(
    B: int, n: int, wmin: int, extra: int, circular: bool, seed: int
) -> None:
    wmax = wmin + extra
    Z = np.random.default_rng(seed).normal(size=(B, n)) * 3
    got = scan_null_max(Z, wmin, wmax, circular)
    assert got.shape == (B,)
    want = np.array([scan_max(row, wmin, wmax, circular).max_abs for row in Z], dtype=float)
    np.testing.assert_array_equal(np.isnan(got), np.isnan(want))
    tol = 1e-9 * max(1.0, float(np.abs(Z).sum(axis=1).max(initial=0.0)))
    fin = ~np.isnan(want)
    np.testing.assert_allclose(got[fin], want[fin], rtol=0, atol=tol)


@SETTINGS
@given(
    st.lists(st.floats(0.0, 20.0), min_size=0, max_size=30),
    st.lists(st.floats(0.0, 20.0), min_size=1, max_size=50),
)
def test_fwer_p_rule(s_abs: list[float], null: list[float]) -> None:
    p = fwer_p(s_abs, null)
    B = len(null)
    assert p.shape == (len(s_abs),)
    assert ((p >= 1 / (1 + B) - 1e-15) & (p <= 1 + 1e-15)).all()
    m = np.asarray(null)
    for v, pv in zip(s_abs, p, strict=True):
        exact = (1 + int((m >= v).sum())) / (1 + B)
        # tie tolerance may only add near-ties, never remove exceedances
        assert pv >= exact - 1e-15
        assert pv <= (1 + int((m >= v - 1e-9 * max(1.0, v)).sum())) / (1 + B) + 1e-15
    order = np.argsort(s_abs, kind="stable")
    assert (np.diff(p[order]) <= 1e-15).all()


@st.composite
def merge_case(draw: st.DrawFn) -> tuple[np.ndarray, int, int, int, bool, np.ndarray]:
    n = draw(st.integers(2, 30))
    z = np.array(draw(st.lists(st.integers(-3, 3), min_size=n, max_size=n)), dtype=float)
    wmin = draw(st.integers(1, 4))
    wmax = draw(st.integers(wmin, 8))
    circular = draw(st.booleans())
    k = n_windows(n, wmin, wmax, circular)
    sig = np.array(draw(st.lists(st.booleans(), min_size=k, max_size=k)), dtype=bool)
    return z, n, wmin, wmax, circular, sig


@SETTINGS
@given(merge_case())
def test_merge_invariants(case: tuple[np.ndarray, int, int, int, bool, np.ndarray]) -> None:
    z, n, wmin, wmax, circular, sig = case
    res = scan_max(z, wmin, wmax, circular)
    sig = sig & (res.s != 0)  # a zero-S window has no direction
    clusters = merge_windows(res, sig, n, circular)
    pos = [
        positions(int(s), int(w), n, circular)
        for s, w in zip(res.start.tolist(), res.length.tolist(), strict=True)
    ]
    seen: dict[int, int] = {}
    for ci, c in enumerate(clusters):
        assert list(c.members) == sorted(set(c.members))
        assert c.best in c.members
        sign = 1 if c.direction == "high" else -1
        for m in c.members:
            assert sig[m]
            assert np.sign(res.s[m]) == sign
            assert m not in seen
            seen[m] = ci
        best_abs = abs(res.s[c.best])
        assert all(abs(res.s[m]) <= best_abs for m in c.members)
        # connected: grow from the best window through shared TUs
        reach, frontier = {c.best}, [c.best]
        while frontier:
            a = frontier.pop()
            for m in c.members:
                if m not in reach and pos[a] & pos[m]:
                    reach.add(m)
                    frontier.append(m)
        assert reach == set(c.members)
    assert set(seen) == set(np.flatnonzero(sig).tolist())
    for a in seen:
        for b in seen:
            if a < b and np.sign(res.s[a]) == np.sign(res.s[b]) and pos[a] & pos[b]:
                assert seen[a] == seen[b]
    starts = [int(res.start[c.best]) for c in clusters]
    assert starts == sorted(starts)
