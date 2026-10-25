"""The abundance-dependent Student-t scale (issue M6-13, ADR-0055 option (a); GIAP_SPEC.md §7.4.2,
§7.4.3 as amended; docs/specs/M6-13-discordance-abundance-sigma.md).

Interface pinned here (the spec's "Public interface"):

- ``fit_student_t_scaled(x, covariate, *, df_bounds=(2, 50), gamma_bounds=(-2, 2), min_sd=0.1)
  -> ScaledStudentT``: Delta_i ~ mu + sigma_i t_nu with sigma_i = sigma0 exp(gamma z_i), z the
  covariate standardised by its mean and SD (ddof = 1) over the fit pairs; pairs with a
  non-finite x or covariate are dropped together;
- ``scale_at(covariate, fit)``: sigma0 exp(gamma clamp(z)), z clamped to [z_min, z_max];
- ``scaled_t_two_sided(x, covariate, fit)``: ``2 min(F, 1 - F)`` with
  ``F = F_t((x - mu) / sigma_i; nu)``.

Examples E1-E5 (a given fit, constants from ``scipy.stats.t``, 1e-6 relative) and E8-E11
(recovery, homoscedastic, affine/order invariance, orientation), the degenerate rows (n < 5,
flat abundance, the constant fit unestimated, collapse, gamma on a bound, ``not_converged`` and
``worse_than_constant`` through a patched optimiser) and the hypothesis properties of the
acceptance criteria. Through-the-step calibration is tests/calibration (M6-09 harness).
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable, Iterator
from typing import Any

import numpy as np
import pytest
import scipy.optimize
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy import stats as sps

import giap.stats.tdist as tdist
from giap.stats.tdist import (
    ScaledStudentT,
    fit_student_t,
    fit_student_t_scaled,
    scale_at,
    scaled_t_two_sided,
    student_t_two_sided,
)

NAN = float("nan")
FAST = settings(max_examples=25, deadline=None)
FALLBACKS = {"flat_abundance", "not_converged", "collapsed", "worse_than_constant"}
MODELS = {"loglinear", "constant", "none"}

#: the given fit of E1-E5
GIVEN = ScaledStudentT(
    loc=-0.3, scale0=0.5, gamma=-0.5, df=5.0, n=100, converged=True, model="loglinear",
    reason=None, center=20.0, sd=2.0, z_min=-2.5, z_max=2.0,
    nll=NAN, nll_constant=NAN, lr_stat=NAN, lr_p=NAN,
)  # fmt: skip


def _p(x: float, m: float, fit: ScaledStudentT = GIVEN) -> float:
    return float(scaled_t_two_sided(np.array([x]), np.array([m]), fit)[0])


def _sigma(m: float, fit: ScaledStudentT = GIVEN) -> float:
    return float(scale_at(np.array([m]), fit)[0])


def _oracle_p(x: float, loc: float, sigma: float, df: float) -> float:
    return float(2.0 * sps.t.sf(abs(x - loc) / sigma, df))


def e8() -> tuple[np.ndarray, np.ndarray]:
    """E8: n = 2,000, z ~ N(0, 1), Delta = -0.3 + 0.5 exp(-0.5 z) t_5, m = z."""
    rng = np.random.default_rng(11)
    z = rng.standard_normal(2000)
    t5 = rng.standard_t(5, 2000)
    return -0.3 + 0.5 * np.exp(-0.5 * z) * t5, z.copy()


def e9() -> tuple[np.ndarray, np.ndarray]:
    """E9: homoscedastic, Delta = 0.05 + 0.4 t_6, then m ~ N(20, 2)."""
    rng = np.random.default_rng(12)
    x = 0.05 + 0.4 * rng.standard_t(6, 2000)
    return x, rng.normal(20, 2, 2000)


def _nll(x: np.ndarray, m: np.ndarray, fit: ScaledStudentT) -> float:
    """The negative log-likelihood of the model at ``fit``, written out with scipy."""
    z = (m - fit.center) / fit.sd
    s = fit.scale0 * np.exp(fit.gamma * z)
    return -float(np.sum(sps.t.logpdf((x - fit.loc) / s, fit.df) - np.log(s)))


@pytest.fixture(scope="module")
def e8_fit() -> ScaledStudentT:
    return fit_student_t_scaled(*e8())


# --------------------------------------------------------------------------- interface


def test_scaled_fit_is_a_frozen_slots_dataclass_with_the_spec_fields() -> None:
    names = [f.name for f in dataclasses.fields(ScaledStudentT)]
    assert names == ["loc", "scale0", "gamma", "df", "n", "converged", "model", "reason",
                     "center", "sd", "z_min", "z_max", "nll", "nll_constant", "lr_stat",
                     "lr_p"]  # fmt: skip
    with pytest.raises(dataclasses.FrozenInstanceError):
        GIVEN.loc = 0.0  # type: ignore[misc]
    assert not hasattr(GIVEN, "__dict__")  # slots


def test_new_names_are_exported() -> None:
    for name in ("ScaledStudentT", "fit_student_t_scaled", "scaled_t_two_sided", "scale_at"):
        assert name in tdist.__all__, name


# --------------------------------------------------------------------------- E1-E5


@pytest.mark.parametrize(
    ("x", "m", "sigma", "p"),
    [
        pytest.param(0.7, 18.0, 0.824361, 0.279284, id="E1_low_abundance"),
        pytest.param(0.7, 22.0, 0.303265, 0.0215361, id="E2_high_abundance"),
        pytest.param(0.7, 20.0, 0.5, 0.101939, id="E3_centre"),
        pytest.param(-1.3, 18.0, 0.824361, 0.279284, id="E4_symmetric_about_mu"),
        pytest.param(0.7, 26.0, 0.183940, 0.00285671, id="E5_clamped_to_z_max"),
    ],
)
def test_examples_e1_to_e5(x: float, m: float, sigma: float, p: float) -> None:
    assert _sigma(m) == pytest.approx(sigma, rel=5e-6)  # the spec's printed digits
    assert _p(x, m) == pytest.approx(p, rel=1e-5)
    z = min(max((m - 20.0) / 2.0, -2.5), 2.0)
    assert _sigma(m) == pytest.approx(0.5 * math.exp(-0.5 * z), rel=1e-12)
    exact = _oracle_p(x, -0.3, 0.5 * math.exp(-0.5 * z), 5.0)
    assert _p(x, m) == pytest.approx(exact, rel=1e-6)


def test_e2_same_delta_is_13x_more_significant_at_high_abundance() -> None:
    assert _p(0.7, 18.0) / _p(0.7, 22.0) == pytest.approx(0.279284 / 0.0215361, rel=1e-5)
    assert _p(0.7, 18.0) / _p(0.7, 22.0) > 12.0


def test_e3_centre_equals_the_constant_model_at_sigma0() -> None:
    const = tdist.StudentT(loc=-0.3, scale=0.5, df=5.0, n=100, converged=True)
    assert _p(0.7, 20.0) == pytest.approx(float(student_t_two_sided([0.7], const)[0]), rel=1e-12)


def test_e5_clamp_and_the_unclamped_value() -> None:
    assert _p(0.7, 26.0) == pytest.approx(0.00285671, rel=1e-5)
    unclamped = _oracle_p(0.7, -0.3, 0.5 * math.exp(-0.5 * 3.0), 5.0)
    assert unclamped == pytest.approx(0.000288207, rel=1e-5)
    assert _p(0.7, 26.0) > 9 * unclamped
    # the lower end clamps as well: m = 10 is z = -5 -> -2.5
    assert _sigma(10.0) == pytest.approx(0.5 * math.exp(-0.5 * -2.5), rel=1e-12)


def test_p_is_one_at_loc_and_nan_in_gives_nan_out() -> None:
    p = scaled_t_two_sided(np.array([-0.3, NAN, 0.2]), np.array([19.0, 20.0, NAN]), GIVEN)
    assert p[0] == 1.0
    assert math.isnan(p[1]) and math.isnan(p[2])


def test_p_keeps_input_order_and_lies_in_unit_interval() -> None:
    x = np.linspace(-5, 5, 41)
    m = np.linspace(14, 26, 41)
    p = scaled_t_two_sided(x, m, GIVEN)
    assert p.shape == x.shape
    assert np.all((p >= 0) & (p <= 1))
    for i in (0, 7, 40):
        assert p[i] == pytest.approx(_p(float(x[i]), float(m[i])), rel=1e-12)


def test_unestimated_fit_gives_nan() -> None:
    none = dataclasses.replace(GIVEN, loc=NAN, scale0=NAN, gamma=NAN, df=NAN, model="none",
                               converged=False)  # fmt: skip
    assert np.isnan(scaled_t_two_sided(np.array([0.1, 0.2]), np.array([20.0, 21.0]), none)).all()
    assert np.isnan(scale_at(np.array([20.0]), none)).all()


# --------------------------------------------------------------------------- E8-E11


def test_e8_recovery(e8_fit: ScaledStudentT) -> None:
    f = e8_fit
    assert f.model == "loglinear" and f.reason is None and f.converged
    assert f.n == 2000
    assert -0.62 <= f.gamma <= -0.38
    assert -0.35 <= f.loc <= -0.25
    assert f.scale0 == pytest.approx(0.5, rel=0.10)
    assert 3.0 <= f.df <= 9.0
    assert f.lr_stat > 20.0


def test_e8_reported_statistics(e8_fit: ScaledStudentT) -> None:
    """The standardisation (S-1: mean and SD with ddof 1 over the fit pairs), the NLLs, the LR and
    its chi-square(1) p (S-7, reported only)."""
    x, m = e8()
    f = e8_fit
    assert f.center == pytest.approx(float(m.mean()), abs=1e-12)
    assert f.sd == pytest.approx(float(m.std(ddof=1)), rel=1e-12)
    assert f.z_min == pytest.approx(float((m.min() - m.mean()) / m.std(ddof=1)), rel=1e-9)
    assert f.z_max == pytest.approx(float((m.max() - m.mean()) / m.std(ddof=1)), rel=1e-9)
    assert f.nll == pytest.approx(_nll(x, m, f), rel=1e-9)
    c = fit_student_t(x)
    nll_c = -float(np.sum(sps.t.logpdf(x, c.df, loc=c.loc, scale=c.scale)))
    assert f.nll_constant == pytest.approx(nll_c, rel=1e-9)
    assert f.lr_stat == pytest.approx(2.0 * (f.nll_constant - f.nll), rel=1e-9)
    assert f.lr_p == pytest.approx(float(sps.chi2.sf(f.lr_stat, 1)), rel=1e-6, abs=1e-300)


def test_e8_scale_at_matches_the_model(e8_fit: ScaledStudentT) -> None:
    _, m = e8()
    z = (m - e8_fit.center) / e8_fit.sd
    np.testing.assert_allclose(scale_at(m, e8_fit), e8_fit.scale0 * np.exp(e8_fit.gamma * z),
                               rtol=1e-12)  # fmt: skip


def test_e9_homoscedastic_gamma_near_zero_and_p_near_constant() -> None:
    x, m = e9()
    f = fit_student_t_scaled(x, m)
    assert f.model == "loglinear"
    assert abs(f.gamma) <= 0.08
    p_c = student_t_two_sided(x, fit_student_t(x))
    assert float(np.max(np.abs(scaled_t_two_sided(x, m, f) - p_c))) <= 0.02


def test_e10_affine_and_order_invariance(e8_fit: ScaledStudentT) -> None:
    x, m = e8()
    p0 = scaled_t_two_sided(x, m, e8_fit)
    affine = fit_student_t_scaled(x, 3.0 * m + 7.0)
    assert affine.gamma == pytest.approx(e8_fit.gamma, abs=1e-6)
    np.testing.assert_allclose(scaled_t_two_sided(x, 3.0 * m + 7.0, affine), p0, atol=1e-6)
    rev = fit_student_t_scaled(x[::-1], m[::-1])
    assert rev.gamma == pytest.approx(e8_fit.gamma, abs=1e-6)
    np.testing.assert_allclose(scaled_t_two_sided(x[::-1], m[::-1], rev)[::-1], p0, atol=1e-6)


def test_e11_orientation(e8_fit: ScaledStudentT) -> None:
    """a and b swapped: Delta negated, m (symmetric) unchanged; to optimiser tolerance."""
    x, m = e8()
    flip = fit_student_t_scaled(-x, m)
    assert flip.loc == pytest.approx(-e8_fit.loc, abs=1e-3)
    for name in ("gamma", "scale0", "df"):
        assert getattr(flip, name) == pytest.approx(getattr(e8_fit, name), abs=1e-3), name
    np.testing.assert_allclose(scaled_t_two_sided(-x, m, flip), scaled_t_two_sided(x, m, e8_fit),
                               atol=1e-3)  # fmt: skip


# --------------------------------------------------------------------------- degenerate rows


def _assert_constant(f: ScaledStudentT, x: np.ndarray, reason: str) -> None:
    c = fit_student_t(x)
    assert f.model == "constant" and f.reason == reason
    assert f.gamma == 0.0
    assert (f.loc, f.scale0, f.df) == (c.loc, c.scale, c.df)  # bit for bit


def test_fewer_than_five_pairs_is_none() -> None:
    f = fit_student_t_scaled([0.1, 0.2, -0.3, 0.4], [10.0, 11.0, 12.0, 13.0])
    assert f.model == "none" and f.n == 4 and not f.converged
    assert all(math.isnan(v) for v in (f.loc, f.scale0, f.gamma, f.df))


def test_non_finite_pairs_are_dropped_together() -> None:
    rng = np.random.default_rng(3)
    x = rng.standard_t(5, 60)
    m = rng.normal(20, 2, 60)
    x2, m2 = x.copy(), m.copy()
    x2[[1, 4]] = NAN
    m2[[7, 9, 11]] = [NAN, np.inf, -np.inf]
    keep = np.isfinite(x2) & np.isfinite(m2)
    f = fit_student_t_scaled(x2, m2)
    assert f.n == int(keep.sum()) == 55
    assert f.center == pytest.approx(float(m[keep].mean()), abs=1e-12)
    g = fit_student_t_scaled(x[keep], m[keep])
    assert (f.gamma, f.loc) == pytest.approx((g.gamma, g.loc), abs=1e-9)


def test_constant_fit_unestimated_is_none() -> None:
    """Dominant ties collapse the constant fit (unchanged rule): no scaled fit is tried."""
    x = np.concatenate([np.zeros(40), [0.5, -0.4, 1.1, -0.9, 0.3]])
    assert not fit_student_t(x).converged
    f = fit_student_t_scaled(x, np.linspace(15, 25, x.size))
    assert f.model == "none" and not f.converged


@pytest.mark.parametrize("spread", [0.0, 0.04], ids=["single_value", "sd_below_0.1"])
def test_flat_abundance_falls_back_to_constant(spread: float) -> None:
    x, _ = e9()
    m = 20.0 + spread * np.random.default_rng(4).standard_normal(x.size)
    assert float(m.std(ddof=1)) < 0.1
    order = np.lexsort((m, x))
    _assert_constant(fit_student_t_scaled(x[order], m[order]), x[order], "flat_abundance")


def test_min_sd_is_the_flat_guard() -> None:
    x, _ = e9()
    m = 20.0 + 0.15 * np.random.default_rng(4).standard_normal(x.size)
    assert fit_student_t_scaled(x, m).model == "loglinear"  # s >= 0.1
    assert fit_student_t_scaled(x, m, min_sd=0.5).reason == "flat_abundance"


def test_twenty_pairs_fit_loglinear() -> None:
    """S-5: no separate pair-count gate; the four-parameter fit runs at n = 20."""
    rng = np.random.default_rng(21)
    z = rng.standard_normal(20)
    x = -0.3 + 0.5 * np.exp(-0.5 * z) * rng.standard_t(5, 20)
    f = fit_student_t_scaled(x, z)
    assert f.model in {"loglinear", "constant"}
    assert f.model != "constant" or f.reason in FALLBACKS - {"flat_abundance"}
    assert f.n == 20


def _planted_gamma(g: float, n: int = 400, seed: int = 5) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    return 0.5 * np.exp(g * z) * rng.standard_t(5, n), z


def test_gamma_on_a_bound_is_accepted() -> None:
    """A planted gamma = -3 is outside [-2, 2]: the fit is returned on the bound (snap), as a df on
    a bound is (M6-01)."""
    f = fit_student_t_scaled(*_planted_gamma(-3.0))
    assert f.model == "loglinear" and f.reason is None
    assert f.gamma == -2.0
    g = fit_student_t_scaled(*_planted_gamma(-3.0), gamma_bounds=(-1.0, 1.0))
    assert g.gamma == -1.0


def test_df_bounds_are_honoured() -> None:
    f = fit_student_t_scaled(*e8(), df_bounds=(6.0, 50.0))
    assert 6.0 <= f.df <= 50.0


# --------------------------------------------------------------------------- patched optimiser

_REAL_MINIMIZE = scipy.optimize.minimize


@pytest.fixture
def patch_scaled(monkeypatch: pytest.MonkeyPatch) -> Iterator[Callable[[Any], None]]:
    """Patch ``scipy.optimize.minimize`` for the four-parameter (mu, log sigma0, gamma, nu) fit
    only; the three-parameter constant fit runs unchanged."""

    def install(edit: Callable[[Any, np.ndarray], Any]) -> None:
        def fake(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
            res = _REAL_MINIMIZE(fun, x0, *a, **k)
            if np.asarray(x0).size != 4:
                return res
            return edit(res, np.asarray(x0, dtype=float))

        monkeypatch.setattr(scipy.optimize, "minimize", fake)
        if hasattr(tdist, "minimize"):
            monkeypatch.setattr(tdist, "minimize", fake)

    yield install


def _with(res: Any, **kw: Any) -> Any:
    out = scipy.optimize.OptimizeResult(dict(res))
    out.update(kw)
    return out


def test_not_converged_falls_back_to_constant(patch_scaled: Any) -> None:
    patch_scaled(lambda res, x0: _with(res, success=False, status=2))
    x, m = e8()
    order = np.lexsort((m, x))
    _assert_constant(fit_student_t_scaled(x[order], m[order]), x[order], "not_converged")


def test_collapsed_falls_back_to_constant(patch_scaled: Any) -> None:
    """sigma0 at e^-40 of the constant scale: min sigma_i under the existing collapse floor."""

    def collapse(res: Any, x0: np.ndarray) -> Any:
        x = x0.copy()
        x[1] -= 40.0
        return _with(res, x=x, success=True)

    patch_scaled(collapse)
    x, m = e8()
    order = np.lexsort((m, x))
    _assert_constant(fit_student_t_scaled(x[order], m[order]), x[order], "collapsed")


def test_worse_than_constant_falls_back(patch_scaled: Any) -> None:
    """The start is the constant optimum, so an optimiser result that fits worse is a failure."""

    def worse(res: Any, x0: np.ndarray) -> Any:
        x = x0.copy()
        x[0] += 5.0 * math.exp(x0[1])
        return _with(res, x=x, success=True)

    patch_scaled(worse)
    x, m = e8()
    order = np.lexsort((m, x))
    _assert_constant(fit_student_t_scaled(x[order], m[order]), x[order], "worse_than_constant")


def test_scaled_fit_starts_at_the_constant_optimum_with_an_explicit_simplex(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """S-6: start (mu_c, log sigma_c, 0, nu_c); an explicit initial simplex with a 0.1 step in
    gamma (scipy's default would leave gamma almost frozen at a zero start); maxiter 4,000 and
    maxfev 40,000."""
    seen: dict[str, Any] = {}

    def spy(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        if np.asarray(x0).size == 4:
            seen["x0"] = np.asarray(x0, dtype=float)
            seen["options"] = dict(k.get("options") or {})
        return _REAL_MINIMIZE(fun, x0, *a, **k)

    monkeypatch.setattr(scipy.optimize, "minimize", spy)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", spy)
    x, m = e8()
    fit_student_t_scaled(x, m)
    c = fit_student_t(np.sort(x))
    x0 = seen["x0"]
    assert x0[2] == 0.0
    assert x0[0] == pytest.approx(c.loc, abs=1e-6)
    assert x0[1] == pytest.approx(math.log(c.scale), abs=1e-6)
    assert x0[3] == pytest.approx(c.df, abs=1e-6)
    simplex = np.asarray(seen["options"]["initial_simplex"], dtype=float)
    assert simplex.shape == (5, 4)
    assert float(np.ptp(simplex[:, 2])) == pytest.approx(0.1, rel=1e-9)
    assert seen["options"].get("maxiter") == 4000
    assert seen["options"].get("maxfev") == 40000


# --------------------------------------------------------------------------- properties


@st.composite
def _data(draw: Any, gamma: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    seed = draw(st.integers(0, 2**31 - 1))
    n = draw(st.integers(20, 300))
    g = draw(st.sampled_from([-1.0, -0.5, 0.0, 0.5])) if gamma is None else gamma
    spread = draw(st.sampled_from([0.0, 0.05, 1.0, 3.0]))
    df = draw(st.sampled_from([3.0, 5.0, 30.0]))
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    x = draw(st.floats(-1, 1)) + 0.4 * np.exp(g * z) * rng.standard_t(df, n)
    return x, 18.0 + spread * z


@FAST
@given(_data(gamma=0.0))
def test_property_gamma_zero_gives_gamma_hat_near_zero_and_p_in_unit_interval(
    data: tuple[np.ndarray, np.ndarray],
) -> None:
    """gamma = 0 data: gamma-hat SD is ~0.25 at n = 20 and ~0.04 at n = 500 (S-5 table); the
    bound is the smallest-n SD times 4."""
    x, m = data
    f = fit_student_t_scaled(x, m)
    p = scaled_t_two_sided(x, m, f)
    if f.model == "none":
        assert np.isnan(p).all()
        return
    assert np.all((p >= 0) & (p <= 1))
    if f.model == "loglinear":
        assert abs(f.gamma) <= max(0.3, 4.0 * 1.15 / math.sqrt(x.size))


@FAST
@given(_data())
def test_property_never_worse_than_constant_and_lr_nonnegative(
    data: tuple[np.ndarray, np.ndarray],
) -> None:
    x, m = data
    f = fit_student_t_scaled(x, m)
    if f.model == "loglinear":
        assert f.nll <= f.nll_constant + 1e-9 * max(1.0, abs(f.nll_constant))
        assert f.lr_stat >= 0.0
        assert 0.0 <= f.lr_p <= 1.0
    elif math.isfinite(f.lr_stat):
        assert f.lr_stat >= 0.0


@FAST
@given(_data(), st.floats(0.5, 5.0), st.floats(-50.0, 50.0), st.integers(0, 2**31 - 1))
def test_property_affine_and_order_invariance(
    data: tuple[np.ndarray, np.ndarray], scale: float, shift: float, seed: int
) -> None:
    x, m = data
    f = fit_student_t_scaled(x, m)
    g = fit_student_t_scaled(x, scale * m + shift)
    perm = np.random.default_rng(seed).permutation(x.size)
    h = fit_student_t_scaled(x[perm], m[perm])
    assert f.model == h.model
    assert f.model == g.model or "flat_abundance" in {f.reason, g.reason}  # s moves with scale
    if f.model == "loglinear" and g.model == "loglinear":
        p = scaled_t_two_sided(x, m, f)
        assert g.gamma == pytest.approx(f.gamma, abs=1e-6)
        np.testing.assert_allclose(scaled_t_two_sided(x, scale * m + shift, g), p, atol=1e-6)
    if f.model == "loglinear":
        assert h.gamma == pytest.approx(f.gamma, abs=1e-6)
        np.testing.assert_allclose(scaled_t_two_sided(x[perm], m[perm], h),
                                   scaled_t_two_sided(x, m, f)[perm], atol=1e-6)  # fmt: skip


@FAST
@given(_data())
def test_property_fallback_reasons_exhaustive_and_exclusive(
    data: tuple[np.ndarray, np.ndarray],
) -> None:
    x, m = data
    f = fit_student_t_scaled(x, m)
    assert f.model in MODELS
    c = fit_student_t(x)
    if f.model == "loglinear":
        assert f.reason is None
        assert c.converged and float(np.std(m, ddof=1)) >= 0.1
    elif f.model == "constant":
        assert f.reason in FALLBACKS
        assert f.gamma == 0.0
        assert (f.reason == "flat_abundance") == (float(np.std(m, ddof=1)) < 0.1)
    else:
        assert not c.converged


# --------------------------------------------------------------------------- review round (PR #321)


def _gaussian(n: int = 400, seed: int = 7) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian Delta (the constant nu-hat is on its upper bound, 50) with a varied covariate."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n)
    return 0.5 * rng.standard_normal(n), z


def test_df_bounds_constant_is_the_default_of_both_fits() -> None:
    """CR S-3 / stat N2: one module constant holds the nu bounds (2, 50)."""
    import inspect

    assert tdist.DF_BOUNDS == (2.0, 50.0)
    for fn in (fit_student_t, fit_student_t_scaled):
        assert inspect.signature(fn).parameters["df_bounds"].default == tdist.DF_BOUNDS
    assert inspect.signature(fit_student_t_scaled).parameters["gamma_bounds"].default == (
        tdist.GAMMA_BOUNDS)  # fmt: skip


def test_inward_nu_step_at_the_upper_bound(monkeypatch: pytest.MonkeyPatch) -> None:
    """T-1 / S-6: with nu_c = 50 the nu step of the initial simplex points inward (every vertex
    nu <= 50) and is not zero (the nu column is not flat)."""
    x, m = _gaussian()
    assert fit_student_t(np.sort(x)).df == 50.0
    seen: dict[str, Any] = {}

    def spy(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        if np.asarray(x0).size == 4:
            seen["simplex"] = np.asarray((k.get("options") or {})["initial_simplex"], dtype=float)
        return _REAL_MINIMIZE(fun, x0, *a, **k)

    monkeypatch.setattr(scipy.optimize, "minimize", spy)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", spy)
    fit_student_t_scaled(x, m)
    nu = seen["simplex"][:, 3]
    assert np.all(nu <= 50.0), nu
    assert float(np.ptp(nu)) > 0.0
    assert float(np.ptp(nu)) == pytest.approx(5.0, rel=1e-9)  # 0.1 nu_c, inward


def test_nu_just_inside_the_bound_snaps_exactly(patch_scaled: Any) -> None:
    """T-1: an optimiser result at nu = 50 - 1e-9 on Gaussian Delta is snapped to exactly 50."""

    def near(res: Any, x0: np.ndarray) -> Any:
        x = np.asarray(res.x, dtype=float).copy()
        x[3] = 50.0 - 1e-9
        return _with(res, x=x, success=True)

    patch_scaled(near)
    f = fit_student_t_scaled(*_gaussian())
    assert f.model == "loglinear", f.reason
    assert f.df == 50.0


def test_gamma_just_inside_the_bound_snaps_exactly(patch_scaled: Any) -> None:
    """T-1: an optimiser result at gamma = -2 + 1e-9 on a planted gamma = -3 is snapped to -2."""

    def near(res: Any, x0: np.ndarray) -> Any:
        x = np.asarray(res.x, dtype=float).copy()
        x[2] = tdist.GAMMA_BOUNDS[0] + 1e-9
        return _with(res, x=x, success=True)

    patch_scaled(near)
    f = fit_student_t_scaled(*_planted_gamma(-3.0))
    assert f.model == "loglinear", f.reason
    assert f.gamma == -2.0


def test_sd_equal_to_min_sd_is_not_flat() -> None:
    """T-3: the flat guard is s >= min_sd (s == min_sd fits the scaled model)."""
    x, m = e9()
    m = 20.0 + 0.15 * np.random.default_rng(4).standard_normal(x.size)
    order = np.lexsort((m, x))
    x, m = x[order], m[order]  # the engine's own order, so s is computed identically
    s = float(np.std(m, ddof=1))
    at = fit_student_t_scaled(x, m, min_sd=s)
    assert at.sd == s
    assert at.reason != "flat_abundance"
    assert fit_student_t_scaled(x, m, min_sd=float(np.nextafter(s, np.inf))).reason == (
        "flat_abundance")  # fmt: skip


def test_lr_is_clamped_at_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    """T-3: an accepted fit whose NLL is above the constant one (within the 1e-9 relative
    tolerance) reports ``lr_stat = 0`` and ``lr_p = 1``, never a negative LR."""

    def fake(fun: Any, x0: Any, *a: Any, **k: Any) -> Any:
        res = _REAL_MINIMIZE(fun, x0, *a, **k)
        if np.asarray(x0).size != 4:
            return res
        args = k.get("args", ())
        start = np.asarray(x0, dtype=float)
        f0 = float(fun(start, *args))
        for e in range(4, 14):
            for sign in (1.0, -1.0):
                x = start.copy()
                x[0] += sign * 10.0 ** (-e) * math.exp(start[1])
                if 1e-9 < float(fun(x, *args)) - f0 < 1e-7:
                    return _with(res, x=x, success=True)
        raise AssertionError("no small worsening step found")

    monkeypatch.setattr(scipy.optimize, "minimize", fake)
    if hasattr(tdist, "minimize"):
        monkeypatch.setattr(tdist, "minimize", fake)
    f = fit_student_t_scaled(*e8())
    assert f.model == "loglinear", f.reason
    assert f.nll > f.nll_constant
    assert f.lr_stat == 0.0
    assert f.lr_p == 1.0
