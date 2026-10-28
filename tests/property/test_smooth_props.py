"""Hypothesis properties of ``loess_detrend`` (issue M5-01; GIPPA_SPEC.md §7.3.5(c), (d) step 3).

- Rotation equivariance on circular input: rotating (y, pos) by any offset (pos taken mod period,
  re-sorted) rotates the residuals (<= 1e-9). Positions are integers (gene indices or TU ordinals,
  with gaps), so the rotated distances are exact.
- The linear fit reproduces any affine y exactly (<= 1e-9) once n >= min_points (below that the
  fit is the mean by definition).
- Circular and linear fits reproduce constants; residuals of a constant are 0.
- Linearity: fit(a y1 + b y2) = a fit(y1) + b fit(y2) (one smoother per pos, ADR-0020).
- ``loess_matrix(pos, ...) @ y == loess_detrend(y, pos, ...).fit`` (<= 1e-12; lead decision 6),
  linear and circular, also for an (n, B) batch column by column.
"""

from __future__ import annotations

import numpy as np
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.stats.smooth import loess_detrend, loess_matrix

SETTINGS = settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])

VAL = st.floats(-100.0, 100.0, allow_nan=False, allow_infinity=False)
FRAC = st.floats(0.05, 1.0)


@st.composite
def circular_case(draw: st.DrawFn, max_size: int = 80) -> tuple[np.ndarray, np.ndarray, int]:
    period = draw(st.integers(3, 150))
    idx = draw(st.sets(st.integers(0, period - 1), min_size=1, max_size=min(period, max_size)))
    pos = np.array(sorted(idx), dtype=float)
    y = np.array(draw(st.lists(VAL, min_size=pos.size, max_size=pos.size)))
    return y, pos, period


@SETTINGS
@given(circular_case(), FRAC, st.integers(0, 10_000))
def test_rotation_equivariance(
    case: tuple[np.ndarray, np.ndarray, int], frac: float, c: int
) -> None:
    y, pos, period = case
    offset = c % period
    base = loess_detrend(y, pos, frac, circular=True, period=float(period))
    rpos = (pos + offset) % period
    order = np.argsort(rpos, kind="stable")
    rot = loess_detrend(y[order], rpos[order], frac, circular=True, period=float(period))
    assert rot.k == base.k
    np.testing.assert_allclose(rot.resid, base.resid[order], rtol=0, atol=1e-9)


@st.composite
def linear_positions(draw: st.DrawFn) -> np.ndarray:
    gaps = draw(st.lists(st.floats(0.1, 10.0), min_size=3, max_size=120))
    start = draw(st.floats(-1e3, 1e3))
    return start + np.concatenate([[0.0], np.cumsum(gaps[:-1])])


@SETTINGS
@given(linear_positions(), st.floats(-10.0, 10.0), st.floats(-10.0, 10.0), FRAC)
def test_linear_reproduces_affine(pos: np.ndarray, a: float, b: float, frac: float) -> None:
    y = a + b * pos
    res = loess_detrend(y, pos, frac, circular=False)
    assert np.max(np.abs(res.resid)) <= 1e-9 * max(1.0, float(np.max(np.abs(y))))


@SETTINGS
@given(circular_case(), FRAC, VAL, st.booleans())
def test_constants_reproduced(
    case: tuple[np.ndarray, np.ndarray, int], frac: float, value: float, circular: bool
) -> None:
    _, pos, period = case
    y = np.full(pos.size, value)
    res = loess_detrend(y, pos, frac, circular=circular, period=float(period))
    assert np.max(np.abs(res.resid)) <= 1e-12 * max(1.0, abs(value))


@SETTINGS
@given(circular_case(), FRAC, st.booleans(), st.data())
def test_linearity(
    case: tuple[np.ndarray, np.ndarray, int], frac: float, circular: bool, data: st.DataObject
) -> None:
    y1, pos, period = case
    y2 = np.array(data.draw(st.lists(VAL, min_size=pos.size, max_size=pos.size)))
    a, b = data.draw(VAL), data.draw(VAL)
    p = float(period)
    f1 = loess_detrend(y1, pos, frac, circular=circular, period=p).fit
    f2 = loess_detrend(y2, pos, frac, circular=circular, period=p).fit
    f12 = loess_detrend(a * y1 + b * y2, pos, frac, circular=circular, period=p).fit
    scale = max(1.0, abs(a) * 100.0 + abs(b) * 100.0)
    np.testing.assert_allclose(f12, a * f1 + b * f2, rtol=0, atol=1e-10 * scale)


@SETTINGS
@given(circular_case(max_size=2), FRAC, st.booleans())
def test_mean_rule_below_min_points(
    case: tuple[np.ndarray, np.ndarray, int], frac: float, circular: bool
) -> None:
    y, pos, period = case
    res = loess_detrend(y, pos, frac, circular=circular, period=float(period))
    assert res.k == pos.size
    np.testing.assert_allclose(res.fit, np.mean(y), rtol=0, atol=1e-12)
    assert abs(float(res.resid.sum())) <= 1e-10


@SETTINGS
@given(circular_case(), FRAC, st.booleans(), st.data())
def test_loess_matrix_equals_detrend(
    case: tuple[np.ndarray, np.ndarray, int], frac: float, circular: bool, data: st.DataObject
) -> None:
    y, pos, period = case
    p = float(period)
    H = loess_matrix(pos, frac, circular=circular, period=p)
    assert H.shape == (pos.size, pos.size)
    tol = 1e-12 * max(1.0, float(np.max(np.abs(y))))
    np.testing.assert_allclose(
        H @ y, loess_detrend(y, pos, frac, circular=circular, period=p).fit, rtol=0, atol=tol
    )
    B = data.draw(st.integers(1, 4))
    Y = np.array(data.draw(st.lists(VAL, min_size=pos.size * B, max_size=pos.size * B)))
    Y = Y.reshape(pos.size, B)
    fits = H @ Y
    for b in range(B):
        ref = loess_detrend(Y[:, b], pos, frac, circular=circular, period=p).fit
        tol_b = 1e-12 * max(1.0, float(np.max(np.abs(Y[:, b]))))
        np.testing.assert_allclose(fits[:, b], ref, rtol=0, atol=tol_b)
