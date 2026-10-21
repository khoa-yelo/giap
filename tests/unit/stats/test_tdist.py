"""Student-t MLE and the two-sided p (issue M6-01; GIAP_SPEC.md §7.4.2, §7.4.3, §8.3).

Spec examples 1-3 of docs/specs/M6-01-stats-tdist-decoy-fdr.md and every row of its "Degenerate
cases" table that concerns the fit:

- ``fit_student_t(x, df_bounds=(2.0, 50.0)) -> StudentT(loc, scale, df, n, converged)`` is the
  maximum likelihood fit over (loc, scale, df), NaNs dropped, ``n`` the number of finite values;
- ``student_t_two_sided(x, fit) -> F64`` is ``2 min(F, 1 - F)`` with ``F = t.cdf((x - loc)/scale,
  df)``, exactly 1.0 at ``x = loc`` and NaN in / NaN out;
- fewer than 5 finite values gives ``StudentT(nan, nan, nan, n, False)`` and no exception;
- a df at a bound is a legitimate optimum (``converged`` True).

Examples 2 and 3 are hand-computed constants (the spec's 4-decimal values) and are also checked
against ``scipy.stats.t`` to full precision, which is the definition the spec writes. The
statistical claims (the fit recovers its parameters, the p is uniform under the fitted t) live in
tests/calibration/test_operon_tdist_calibration.py, not here.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy import stats as sps

from giap.stats.tdist import StudentT, fit_student_t, student_t_two_sided

FAST = settings(max_examples=50, deadline=None)
NAN = float("nan")


def _t_sample(rng: np.random.Generator, loc: float, scale: float, df: float, n: int) -> np.ndarray:
    return loc + scale * rng.standard_t(df, size=n)


def _two_sided(x: float, loc: float, scale: float, df: float) -> float:
    """The spec's definition, written out with scipy as the oracle for F_t."""
    f = float(sps.t.cdf((x - loc) / scale, df))
    return 2.0 * min(f, 1.0 - f)


# --------------------------------------------------------------------------- the fit


def test_example1_fit_recovers_the_drawing_parameters() -> None:
    """Example 1: 2,000 draws from t(0.3, 0.5, df 5), fixed seed."""
    x = _t_sample(np.random.default_rng(601), 0.3, 0.5, 5.0, 2000)
    fit = fit_student_t(x)
    assert fit.converged
    assert fit.n == 2000
    assert fit.loc == pytest.approx(0.3, abs=0.05)
    assert fit.scale == pytest.approx(0.5, abs=0.05)
    assert fit.df == pytest.approx(5.0, abs=2.0)


def test_fit_is_a_frozen_slots_dataclass_with_the_spec_fields() -> None:
    fit = fit_student_t(_t_sample(np.random.default_rng(602), 0.0, 1.0, 6.0, 200))
    assert [f.name for f in dataclasses.fields(fit)] == ["loc", "scale", "df", "n", "converged"]
    assert isinstance(fit.n, int) and isinstance(fit.converged, bool)
    with pytest.raises(dataclasses.FrozenInstanceError):
        fit.loc = 1.0  # type: ignore[misc]


def test_fit_drops_nan_and_counts_only_finite_values() -> None:
    x = _t_sample(np.random.default_rng(603), -1.0, 2.0, 7.0, 300)
    padded = np.concatenate([x[:100], [NAN], x[100:], [NAN, NAN]])
    bare, with_nan = fit_student_t(x), fit_student_t(padded)
    assert with_nan.n == 300
    assert with_nan.loc == pytest.approx(bare.loc, abs=1e-12)
    assert with_nan.scale == pytest.approx(bare.scale, abs=1e-12)
    assert with_nan.df == pytest.approx(bare.df, abs=1e-12)
    assert with_nan.converged == bare.converged


@pytest.mark.parametrize("n", [0, 1, 2, 3, 4])
def test_fewer_than_five_finite_values_is_unestimated(n: int) -> None:
    fit = fit_student_t(np.concatenate([np.arange(n, dtype=float), [NAN]]))
    assert fit.n == n
    assert not fit.converged
    assert math.isnan(fit.loc) and math.isnan(fit.scale) and math.isnan(fit.df)


def test_zero_mad_does_not_raise() -> None:
    """All values equal: the likelihood is unbounded, so the fit need not converge -- but it
    returns a ``StudentT`` instead of raising."""
    fit = fit_student_t(np.full(50, 2.5))
    assert fit.n == 50
    assert not fit.converged or fit.scale > 0.0


def test_df_at_the_upper_bound_is_converged() -> None:
    """Gaussian data: the ML df runs to the bound, which is a legitimate optimum."""
    x = np.random.default_rng(604).normal(0.0, 1.0, 4000)
    fit = fit_student_t(x)
    assert fit.df == pytest.approx(50.0)
    assert fit.converged
    assert fit.loc == pytest.approx(0.0, abs=0.1)


def test_df_bounds_are_respected_and_configurable() -> None:
    x = np.random.default_rng(605).normal(0.0, 1.0, 2000)
    fit = fit_student_t(x, df_bounds=(3.0, 10.0))
    assert fit.df == pytest.approx(10.0)
    assert fit.converged
    heavy = np.random.default_rng(606).standard_cauchy(2000)
    low = fit_student_t(heavy)
    assert low.df == pytest.approx(2.0)
    assert 2.0 <= low.df <= 50.0


def test_fit_stays_inside_the_bounds_on_heavy_and_light_tails() -> None:
    rng = np.random.default_rng(607)
    for df_true in (2.5, 5.0, 30.0):
        fit = fit_student_t(_t_sample(rng, 0.0, 1.0, df_true, 1500))
        assert 2.0 <= fit.df <= 50.0
        assert fit.scale > 0.0


# ------------------------------------------------------------------- the two-sided p


def test_example2_positive_delta() -> None:
    """Example 2: Delta = 1.2 under StudentT(0.05, 0.4, 6) -> z = 2.875, p = 0.0282."""
    fit = StudentT(loc=0.05, scale=0.4, df=6.0, n=100, converged=True)
    p = student_t_two_sided([1.2], fit)
    assert float(p[0]) == pytest.approx(0.0282, abs=5e-5)
    assert float(p[0]) == pytest.approx(_two_sided(1.2, 0.05, 0.4, 6.0), rel=1e-12)


def test_example3_negative_delta_is_asymmetric() -> None:
    """Example 3: Delta = -1.2 under the same fit -> p = 0.0205 (loc != 0)."""
    fit = StudentT(loc=0.05, scale=0.4, df=6.0, n=100, converged=True)
    p = student_t_two_sided([-1.2], fit)
    assert float(p[0]) == pytest.approx(0.0205, abs=5e-5)
    assert float(p[0]) == pytest.approx(_two_sided(-1.2, 0.05, 0.4, 6.0), rel=1e-12)
    assert float(p[0]) != pytest.approx(float(student_t_two_sided([1.2], fit)[0]), abs=1e-6)


def test_p_is_exactly_one_at_the_location() -> None:
    fit = StudentT(loc=-0.37, scale=1.3, df=4.5, n=80, converged=True)
    assert float(student_t_two_sided([-0.37], fit)[0]) == 1.0


def test_nan_in_nan_out_and_shape_is_preserved() -> None:
    fit = StudentT(loc=0.0, scale=1.0, df=6.0, n=50, converged=True)
    p = student_t_two_sided([0.0, NAN, 2.0], fit)
    assert p.shape == (3,)
    assert p.dtype == np.float64
    assert float(p[0]) == 1.0
    assert math.isnan(float(p[1]))
    assert 0.0 < float(p[2]) < 1.0


def test_an_unestimated_fit_gives_nan_p() -> None:
    fit = fit_student_t([1.0, 2.0, 3.0])
    p = student_t_two_sided([0.0, 1.0, 5.0], fit)
    assert np.isnan(p).all()


def test_p_is_symmetric_around_a_zero_location() -> None:
    fit = StudentT(loc=0.0, scale=0.8, df=9.0, n=120, converged=True)
    x = np.array([0.1, 0.5, 1.0, 3.0])
    assert student_t_two_sided(x, fit) == pytest.approx(student_t_two_sided(-x, fit), rel=1e-12)


def test_p_decreases_with_distance_from_the_location() -> None:
    fit = StudentT(loc=0.2, scale=0.5, df=5.0, n=200, converged=True)
    p = student_t_two_sided([0.2, 0.4, 0.9, 2.0, 6.0], fit)
    assert np.all(np.diff(p) < 0.0)
    assert np.all((p >= 0.0) & (p <= 1.0))


def test_empty_input_gives_an_empty_array() -> None:
    fit = StudentT(loc=0.0, scale=1.0, df=6.0, n=10, converged=True)
    p = student_t_two_sided([], fit)
    assert p.shape == (0,)
    assert p.dtype == np.float64


# --------------------------------------------------------------------- properties


@FAST
@given(
    x=st.floats(-50.0, 50.0, allow_nan=False),
    loc=st.floats(-5.0, 5.0),
    scale=st.floats(0.05, 5.0),
    df=st.floats(2.0, 50.0),
    a=st.floats(0.1, 10.0),
    b=st.floats(-20.0, 20.0),
)
def test_property_two_sided_p_is_affine_invariant(
    x: float, loc: float, scale: float, df: float, a: float, b: float
) -> None:
    """Acceptance: invariant to an affine shift applied to both x and the fit."""
    plain = StudentT(loc=loc, scale=scale, df=df, n=99, converged=True)
    moved = StudentT(loc=a * loc + b, scale=a * scale, df=df, n=99, converged=True)
    got = float(student_t_two_sided([a * x + b], moved)[0])
    assert got == pytest.approx(float(student_t_two_sided([x], plain)[0]), abs=1e-9)


@FAST
@given(
    xs=st.lists(st.floats(-100.0, 100.0), min_size=1, max_size=20),
    loc=st.floats(-5.0, 5.0),
    scale=st.floats(0.05, 5.0),
    df=st.floats(2.0, 50.0),
)
def test_property_two_sided_p_is_in_the_unit_interval(
    xs: list[float], loc: float, scale: float, df: float
) -> None:
    fit = StudentT(loc=loc, scale=scale, df=df, n=99, converged=True)
    p = student_t_two_sided(xs, fit)
    assert p.shape == (len(xs),)
    assert np.all((p >= 0.0) & (p <= 1.0))


@FAST
@given(
    u=st.floats(0.0, 20.0),
    v=st.floats(0.0, 20.0),
    loc=st.floats(-3.0, 3.0),
    scale=st.floats(0.1, 3.0),
    df=st.floats(2.0, 50.0),
)
def test_property_two_sided_p_is_monotone_in_the_absolute_z(
    u: float, v: float, loc: float, scale: float, df: float
) -> None:
    fit = StudentT(loc=loc, scale=scale, df=df, n=99, converged=True)
    p = student_t_two_sided([loc + scale * u, loc + scale * v], fit)
    if u < v:
        assert float(p[0]) >= float(p[1]) - 1e-12
    else:
        assert float(p[1]) >= float(p[0]) - 1e-12


@FAST
@given(
    n=st.integers(5, 60),
    loc=st.floats(-2.0, 2.0),
    scale=st.floats(0.2, 3.0),
    seed=st.integers(0, 2**32 - 1),
)
def test_property_fit_returns_a_usable_or_explicitly_unestimated_result(
    n: int, loc: float, scale: float, seed: int
) -> None:
    """No silent fallback: either the fields are finite and in range, or the fit is unestimated."""
    x = _t_sample(np.random.default_rng(seed), loc, scale, 6.0, n)
    fit = fit_student_t(x)
    assert fit.n == n
    if math.isnan(fit.loc):
        assert not fit.converged and math.isnan(fit.scale) and math.isnan(fit.df)
    else:
        assert math.isfinite(fit.scale) and fit.scale > 0.0
        assert 2.0 <= fit.df <= 50.0
