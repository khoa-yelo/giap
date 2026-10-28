"""Code review of #140 (M6-01 decoy FDR): non-finite scores, the lfdr "larger wins" rule (lead
decision 3), the accuracy of the binned reflected KDE and a ``bin_index`` edge case."""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt
import pytest
from scipy import stats as sps

from gippa.stats.decoy_fdr import _reflected_density, bin_index, local_fdr, max_decoy_fdr

F64 = npt.NDArray[np.float64]


def _sample(seed: int, n: int, scale: float = 1.0) -> F64:
    return -np.abs(np.random.default_rng(seed).normal(scale=scale, size=n))


# ------------------------------------------------------------------- non-finite scores (B1)


@pytest.mark.parametrize("bad", [-math.inf, math.inf])
def test_an_infinite_target_is_excluded_by_max_decoy_fdr(bad: float) -> None:
    t_s = np.append(_sample(1, 60, 0.3), bad)
    t_b = np.zeros(t_s.size, dtype=np.int64)
    d_s = _sample(2, 80)
    decoys = {"convergent": (d_s, np.zeros(d_s.size, dtype=np.int64))}
    got = max_decoy_fdr(t_s, t_b, decoys, n_bins=1)
    plain = max_decoy_fdr(t_s[:-1], t_b[:-1], decoys, n_bins=1)
    assert math.isnan(got.q[-1]) and math.isnan(got.fdr[-1]) and math.isnan(got.lfdr[-1])
    assert got.which[-1] == ""
    np.testing.assert_allclose(got.q[:-1], plain.q, rtol=1e-12)
    np.testing.assert_allclose(got.fdr[:-1], plain.fdr, rtol=1e-12)
    np.testing.assert_allclose(got.lfdr[:-1], plain.lfdr, rtol=1e-12)
    assert np.isfinite(got.lfdr[:-1]).all()


@pytest.mark.parametrize("bad", [-math.inf, math.inf])
def test_an_infinite_target_is_excluded_by_local_fdr(bad: float) -> None:
    target = np.append(_sample(3, 60, 0.3), bad)
    decoy = _sample(4, 80)
    got = local_fdr(target, decoy, np.ones(decoy.size))
    assert math.isnan(got[-1])
    np.testing.assert_allclose(got[:-1], local_fdr(target[:-1], decoy, np.ones(80)), rtol=1e-12)
    assert np.isfinite(got[:-1]).all()


def test_a_minus_inf_decoy_is_dropped() -> None:
    target = _sample(5, 60, 0.3)
    decoy = np.append(_sample(6, 80), -math.inf)
    got = local_fdr(target, decoy, np.ones(decoy.size))
    np.testing.assert_allclose(got, local_fdr(target, decoy[:-1], np.ones(80)), rtol=1e-12)
    assert np.isfinite(got).all()
    t_b = np.zeros(target.size, dtype=np.int64)
    d_b = np.zeros(decoy.size, dtype=np.int64)
    with_inf = max_decoy_fdr(target, t_b, {"c": (decoy, d_b)}, n_bins=1)
    without = max_decoy_fdr(target, t_b, {"c": (decoy[:-1], d_b[:-1])}, n_bins=1)
    np.testing.assert_allclose(with_inf.fdr, without.fdr, rtol=1e-12)
    np.testing.assert_allclose(with_inf.lfdr, without.lfdr, rtol=1e-12)


def test_no_survivor_fields_are_distinct_arrays() -> None:
    got = max_decoy_fdr([-0.1, -0.2], [0, 0], {"c": ([math.nan], [0])}, n_bins=1)
    assert got.q is not got.lfdr and got.q is not got.fdr and got.fdr is not got.lfdr


# ------------------------------------------------------------ lead decision 3: lfdr = fmax


def _two_sets() -> tuple[F64, npt.NDArray[np.int64], F64, F64]:
    rng = np.random.default_rng(7)
    t_s = -np.abs(rng.normal(scale=1.0, size=300))
    t_b = rng.integers(0, 2, size=300).astype(np.int64)
    narrow = -np.abs(rng.normal(scale=0.4, size=400))  # null mass near |Delta| = 0
    wide = -np.abs(rng.normal(loc=1.5, scale=0.4, size=400))  # null mass near |Delta| = 1.5
    return t_s, t_b, narrow, wide


def test_lfdr_is_the_elementwise_fmax_over_the_decoy_sets() -> None:
    t_s, t_b, narrow, wide = _two_sets()
    rng = np.random.default_rng(8)
    nb = rng.integers(0, 2, size=narrow.size).astype(np.int64)
    wb = rng.integers(0, 2, size=wide.size).astype(np.int64)
    a = max_decoy_fdr(t_s, t_b, {"narrow": (narrow, nb)}, n_bins=2).lfdr
    b = max_decoy_fdr(t_s, t_b, {"wide": (wide, wb)}, n_bins=2).lfdr
    both = max_decoy_fdr(t_s, t_b, {"narrow": (narrow, nb), "wide": (wide, wb)}, n_bins=2).lfdr
    # the two sets must disagree in both directions, or fmax and fmin would coincide
    assert np.any(a > b + 1e-6) and np.any(b > a + 1e-6)
    np.testing.assert_allclose(both, np.fmax(a, b), rtol=1e-12)


def test_an_undefined_kde_does_not_mask_a_defined_one() -> None:
    t_s, t_b, narrow, _ = _two_sets()
    nb = np.zeros(narrow.size, dtype=np.int64)
    few = narrow[:3]  # fewer than MIN_KDE_N decoys: this set's lfdr is undefined
    defined = max_decoy_fdr(t_s, t_b, {"narrow": (narrow, nb)}, n_bins=2).lfdr
    assert np.isfinite(defined).all()
    for order in (("few", "narrow"), ("narrow", "few")):
        sets = {"few": (few, nb[:3]), "narrow": (narrow, nb)}
        got = max_decoy_fdr(t_s, t_b, {k: sets[k] for k in order}, n_bins=2).lfdr
        np.testing.assert_allclose(got, defined, rtol=1e-12)


# ---------------------------------------------------- accuracy of the binned reflected KDE


def _exact_reflected(sample: F64, w: F64, x: F64) -> F64:
    h = float(np.sqrt(sps.gaussian_kde(sample, bw_method="silverman").covariance[0, 0]))
    w = w / w.sum()
    x = np.abs(x)[:, None]
    phi = sps.norm.pdf(x - sample[None, :], scale=h) + sps.norm.pdf(x + sample[None, :], scale=h)
    out: F64 = phi @ w
    return out


@pytest.mark.parametrize("n", [50, 2000])
def test_reflected_density_matches_the_exact_kernel_sum(n: int) -> None:
    rng = np.random.default_rng(n)
    sample = np.abs(rng.normal(size=n))
    w = rng.uniform(0.5, 2.0, size=n)
    x = np.linspace(0.0, float(sample.max()), 997)
    np.testing.assert_allclose(
        _reflected_density(sample, w, x), _exact_reflected(sample, w, x), rtol=5e-3
    )


def test_reflected_density_integrates_to_one_on_the_half_line() -> None:
    sample = np.abs(np.random.default_rng(9).normal(size=500))
    x = np.linspace(0.0, 8.0, 4001)
    assert np.trapezoid(_reflected_density(sample, np.ones(500), x), x) == pytest.approx(
        1.0, abs=1e-3
    )


# ------------------------------------------------------------------------- bin_index


def test_a_distance_above_the_last_finite_edge_is_outside_every_bin() -> None:
    assert bin_index([150.0, 100.0, 99.0], [0.0, 50.0, 100.0]).tolist() == [-1, -1, 1]
