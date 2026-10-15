"""M4-13 detectability model (docs/specs/M4-13-detectability.md; GIAP_SPEC.md §7.3.3, §5.2, §8.2,
§8.4, G5, §6.6 / TE7).

Examples 1-9 of the spec, the edge cases and the interface rules, each against the independent
oracle in tests/helpers/detectability_oracle.py (own L2 logistic solver, own fold / C-search /
binning / recalibration code). Inputs from tests/helpers/detectability_data.py.

Resolved ambiguities (pinned here; the lead may overrule):
- Eligibility is carried by ``detected``: a pandas ``boolean`` Series indexed by locus_tag whose
  ``<NA>`` rows are not eligible (partial / non-coding). ``p_detect`` has ``detected``'s index and
  is NaN exactly there; the constant rule and diagnostics count eligible rows only.
- Folds are drawn over the TUs of eligible genes; the oracle tests use all-eligible inputs so the
  choice does not matter there.
- RNG consumption: ``tu_folds`` permutes sorted unique TU IDs with ``rng.permutation``; the outer
  folds are drawn first, then one inner ``tu_folds(training TUs, 3, rng)`` per outer fold in fold
  order 0..k-1, from the same generator.
- Inner metric: mean over the 3 inner folds of the held-out mean log loss (GridSearchCV
  convention); the C test allows a different choice only when the two scores are within 1e-3
  (sklearn lbfgs default tol 1e-4 vs the tight oracle).
- ``p_detect`` is compared with the oracle at atol 2e-3 (same reason).
- Equal-count bins: stable rank order split by ``np.array_split`` (first bins one larger); exact
  bin contents are only checked for n divisible by n_bins.
- ``transform_covariates`` returns the columns named as ``names`` in that order, input index kept;
  NaN is replaced by the median of the transformed column (odd non-NaN counts used here).
- Raw-``tm_segments`` / zero-variance decisions are made on all eligible genes (the inputs here
  make the training-fold reading identical).
- "Warning" / "logged" (``too_few_tus``, dropped covariate): accepted in the diagnostics, a Python
  warning or a log record.
"""

from __future__ import annotations

import logging
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detectability_data as dd
import detectability_oracle as O

from giap.config.model import DetectabilityConfig
from giap.rng import for_step

C_GRID = {0.01, 0.1, 1.0, 10.0}
SEED = 7


def _mod() -> Any:
    from giap.genetable import detectability

    return detectability


def _fit(gs: dd.GeneSet, cfg: DetectabilityConfig | None = None, seed: int = SEED) -> Any:
    return _mod().fit_detectability(
        gs.detected,
        gs.cov,
        gs.tu_of_gene,
        cfg or DetectabilityConfig(),
        for_step("detectability", seed),
    )


def _rows(diag: pd.DataFrame) -> pd.DataFrame:
    if "section" in diag.columns:
        assert set(diag["section"]) <= {"detectability"}
    return diag


def _metric(diag: pd.DataFrame, group_type: str, group: object, metric: str) -> float:
    d = _rows(diag)
    sel = d[
        (d["group_type"] == group_type)
        & (d["group"].astype(str) == str(group))
        & (d["metric"] == metric)
    ]
    assert len(sel) == 1, f"{group_type}/{group}/{metric}: {len(sel)} rows"
    return float(sel["value"].iloc[0])


def _mentions(token: str, res: Any, caplog: pytest.LogCaptureFixture, rec: list[Any]) -> bool:
    diag = res.diagnostics
    in_diag = any(token in str(v) for col in diag.columns for v in diag[col].tolist())
    in_warn = any(token in str(w.message) for w in rec)
    return in_diag or in_warn or token in caplog.text


# ----------------------------------------------------------------------------- examples 1, 2


def test_example1_constant_below_50_detected() -> None:
    gs = dd.with_counts(49, 951, seed=1)
    res = _fit(gs)
    assert res.constant is True
    assert np.allclose(res.p_detect.to_numpy(float), 0.049, rtol=0, atol=1e-15)


def test_example1_mirror_constant_below_50_undetected() -> None:
    gs = dd.with_counts(951, 49, seed=2)
    res = _fit(gs)
    assert res.constant is True
    assert np.allclose(res.p_detect.to_numpy(float), 0.951, rtol=0, atol=1e-15)


def test_constant_rate_counts_eligible_only() -> None:
    gs = dd.with_counts(49, 951, seed=3)
    extra = pd.Series(pd.array([pd.NA] * 5, dtype="boolean"), index=[f"X{i}" for i in range(5)])
    det = pd.concat([gs.detected, extra])
    cov = pd.concat([gs.cov, gs.cov.iloc[:5].set_axis(extra.index)])
    tu = pd.concat([gs.tu_of_gene, pd.Series(["TU_X"] * 5, index=extra.index)])
    res = _mod().fit_detectability(
        det, cov, tu, DetectabilityConfig(), for_step("detectability", 1)
    )
    assert res.constant is True
    assert res.p_detect.loc[gs.detected.index].eq(0.049).all()
    assert res.p_detect.loc[extra.index].isna().all()


def test_example2_fitted_at_exactly_50_each() -> None:
    gs = dd.with_counts(50, 50, seed=4)
    res = _fit(gs)
    assert res.constant is False
    p = res.p_detect.to_numpy(float)
    assert np.isfinite(p).all() and ((p > 0) & (p < 1)).all()
    assert np.unique(p).size > 1
    assert len(res.chosen_C) == 5 and set(res.chosen_C) <= C_GRID


@pytest.mark.parametrize("all_detected", [False, True])
def test_no_detected_or_all_detected_is_constant(all_detected: bool) -> None:
    gs = dd.with_counts(400 if all_detected else 0, 0 if all_detected else 400, seed=5)
    res = _fit(gs)
    assert res.constant is True
    assert np.all(res.p_detect.to_numpy(float) == (1.0 if all_detected else 0.0))


# ----------------------------------------------------------------------------- examples 3, 4


def test_example3_tu_folds_deterministic_and_balanced() -> None:
    tus = [f"TU_{i}" for i in range(1, 11)]
    a = _mod().tu_folds(tus, 5, np.random.default_rng(0))
    b = _mod().tu_folds(tus, 5, np.random.default_rng(0))
    assert a == b
    assert set(a) == set(tus)
    assert sorted(np.bincount(list(a.values()), minlength=5).tolist()) == [2] * 5
    assert a == O.tu_folds(tus, 5, np.random.default_rng(0))


def test_tu_folds_ignore_input_order_and_duplicates() -> None:
    tus = [f"TU_{i}" for i in range(1, 24)]
    a = _mod().tu_folds(tus, 5, np.random.default_rng(3))
    b = _mod().tu_folds(list(reversed(tus)) + tus[:4], 5, np.random.default_rng(3))
    assert a == b


def test_example4_genes_of_one_tu_share_a_fold() -> None:
    gs = dd.synthetic(150, seed=11)
    res = _fit(gs)
    assert res.fold_of_tu == O.tu_folds(list(gs.tu_of_gene), 5, for_step("detectability", SEED))
    assert set(res.fold_of_tu) == set(gs.tu_of_gene)
    assert set(res.fold_of_tu.values()) == set(range(5))


def test_example4_held_out_labels_never_reach_their_prediction() -> None:
    """Cross-fitting: flipping the labels of genes in fold j changes no p_detect of fold j."""
    gs = dd.synthetic(150, seed=12)
    base = _fit(gs)
    folds = base.fold_of_tu
    fold_of_gene = gs.tu_of_gene.map(folds)
    j = 2
    in_j = (fold_of_gene == j).to_numpy()
    flipped = gs.detected.copy()
    flipped[in_j] = ~flipped[in_j]
    res = _mod().fit_detectability(
        flipped, gs.cov, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", SEED)
    )
    assert res.constant is False
    np.testing.assert_array_equal(res.p_detect[in_j].to_numpy(), base.p_detect[in_j].to_numpy())
    assert res.chosen_C[j] == base.chosen_C[j]


# ----------------------------------------------------------------------------- model vs oracle


def test_cross_fitted_p_matches_oracle_given_chosen_c() -> None:
    gs = dd.synthetic(150, seed=21)
    res = _fit(gs)
    p, folds, _, _ = O.cross_fit(
        gs.detected, gs.cov, gs.tu_of_gene, dd.NAMES, 5, 5, for_step("detectability", SEED),
        fixed_c=res.chosen_C,
    )  # fmt: skip
    assert res.fold_of_tu == folds
    assert list(res.p_detect.index) == list(gs.detected.index)
    np.testing.assert_allclose(res.p_detect.to_numpy(float), p.to_numpy(), rtol=0, atol=2e-3)


def test_chosen_c_minimises_inner_log_loss() -> None:
    gs = dd.synthetic(150, seed=22)
    res = _fit(gs)
    _, _, chosen, scores = O.cross_fit(
        gs.detected, gs.cov, gs.tu_of_gene, dd.NAMES, 5, 5, for_step("detectability", SEED)
    )
    assert len(res.chosen_C) == 5
    for j, (c_impl, c_or, sc) in enumerate(zip(res.chosen_C, chosen, scores, strict=True)):
        assert c_impl in C_GRID, (j, c_impl)
        assert c_impl == c_or or abs(sc[c_impl] - sc[c_or]) < 1e-3, (j, c_impl, c_or, sc)


def test_fit_diagnostics_are_calibration_diagnostics_of_eligible_genes() -> None:
    gs = dd.synthetic(150, seed=23)
    res = _fit(gs)
    d = gs.detected.astype(bool).to_numpy()
    ref = O.diagnostics(d, res.p_detect.to_numpy(float))
    diag = res.diagnostics
    assert set(diag["section"]) == {"detectability"}
    overall = diag[(diag["group_type"] == "overall")]
    got = dict(zip(overall["metric"], overall["value"].astype(float), strict=False))
    assert got["brier"] == pytest.approx(ref["brier"], abs=1e-12)
    assert got["calibration_slope"] == pytest.approx(ref["calibration_slope"], abs=1e-5)
    assert got["calibration_intercept"] == pytest.approx(ref["calibration_intercept"], abs=1e-5)
    assert got["n_detected"] == ref["n_detected"] and got["n_eligible"] == ref["n_eligible"]
    bins_n = diag[(diag["group_type"] == "bin") & (diag["metric"] == "n")]
    assert sorted(bins_n["group"].astype(int)) == list(range(1, 11))
    assert bins_n["value"].astype(float).sum() == len(d)


# ----------------------------------------------------------------------------- example 5, 6


def test_example5_two_bins() -> None:
    diag = _mod().calibration_diagnostics([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], n_bins=2)
    assert _metric(diag, "bin", 1, "mean_p") == pytest.approx(0.15, abs=1e-12)
    assert _metric(diag, "bin", 1, "observed_rate") == 0.0
    assert _metric(diag, "bin", 2, "mean_p") == pytest.approx(0.85, abs=1e-12)
    assert _metric(diag, "bin", 2, "observed_rate") == 1.0
    assert _metric(diag, "bin", 1, "n") == 2 and _metric(diag, "bin", 2, "n") == 2
    overall = diag[diag["group_type"] == "overall"]
    got = dict(zip(overall["metric"], overall["value"].astype(float), strict=False))
    assert got["brier"] == pytest.approx(0.025, abs=1e-12)
    assert got["n_detected"] == 2 and got["n_eligible"] == 4
    assert {"calibration_intercept", "calibration_slope"} <= set(got)


def test_bins_ties_broken_by_input_order() -> None:
    diag = _mod().calibration_diagnostics([1, 1, 0, 0], [0.5, 0.5, 0.5, 0.5], n_bins=2)
    assert _metric(diag, "bin", 1, "observed_rate") == 1.0
    assert _metric(diag, "bin", 2, "observed_rate") == 0.0


def _calibrated(n: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    r = np.random.default_rng(seed)
    p = expit(r.normal(0.0, 1.5, n))
    return (r.random(n) < p).astype(int), p


def test_example6_perfect_calibration_slope_and_intercept() -> None:
    d, p = _calibrated(20_000, seed=6)
    diag = _mod().calibration_diagnostics(d, p)
    slope = _metric(diag, "overall", _overall_group(diag), "calibration_slope")
    icpt = _metric(diag, "overall", _overall_group(diag), "calibration_intercept")
    assert abs(slope - 1) <= 0.05 and abs(icpt) <= 0.05
    a, b = O.logistic_recalibration(d, p)
    assert icpt == pytest.approx(a, abs=1e-6) and slope == pytest.approx(b, abs=1e-6)


def _overall_group(diag: pd.DataFrame) -> object:
    groups = diag.loc[diag["group_type"] == "overall", "group"].unique()
    assert len(groups) == 1
    return groups[0]


def test_calibration_diagnostics_match_oracle() -> None:
    d, p = _calibrated(1_000, seed=8)
    p = np.round(p, 2)  # ties
    p[:3] = [0.0, 1.0, 1.0]  # clipped to [1e-6, 1 - 1e-6] before the logit
    diag = _mod().calibration_diagnostics(d, p)
    ref = O.diagnostics(d, p)
    g = _overall_group(diag)
    for row in ref["bins"]:  # type: ignore[attr-defined]
        for m in ("mean_p", "observed_rate", "n"):
            assert _metric(diag, "bin", row["group"], m) == pytest.approx(row[m], abs=1e-12)
    for m in ("brier", "n_detected", "n_eligible"):
        assert _metric(diag, "overall", g, m) == pytest.approx(ref[m], abs=1e-12)
    for m in ("calibration_intercept", "calibration_slope"):
        assert _metric(diag, "overall", g, m) == pytest.approx(ref[m], abs=1e-5)


# ----------------------------------------------------------------------------- example 7, 8, 9


@pytest.mark.slow
def test_example7_sim_covariates_only_spearman() -> None:
    gs = dd.sim_geneset(np.random.default_rng(4500), 4500)
    res = _fit(gs)
    elig = gs.extra["eligible"]
    rho = spearmanr(res.p_detect[elig], gs.p_true[elig])[0]
    assert res.constant is False
    assert rho >= 0.8, rho


def test_example8_non_eligible_rows_are_nan() -> None:
    gs = dd.synthetic(120, seed=31)
    det = gs.detected.copy()
    det.iloc[[0, 5, 9]] = pd.NA
    cov = gs.cov.copy()
    cov.iloc[5] = np.nan  # a partial gene may have no covariates at all
    res = _mod().fit_detectability(
        det, cov, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", 1)
    )
    assert list(res.p_detect.index) == list(det.index)
    na = det.isna().to_numpy()
    assert res.p_detect[na].isna().all()
    assert np.isfinite(res.p_detect[~na].to_numpy(float)).all()
    assert _metric(res.diagnostics, "overall", _overall_group(res.diagnostics), "n_eligible") == (
        len(det) - 3
    )


def test_example9_bitwise_deterministic() -> None:
    gs = dd.synthetic(150, seed=41)
    a, b = _fit(gs, seed=99), _fit(gs, seed=99)
    assert a.p_detect.equals(b.p_detect)
    assert np.array_equal(a.p_detect.to_numpy(float), b.p_detect.to_numpy(float), equal_nan=True)
    assert a.fold_of_tu == b.fold_of_tu and a.chosen_C == b.chosen_C
    c = _fit(gs, seed=100)
    assert c.fold_of_tu != a.fold_of_tu  # the seed reaches the folds


# ----------------------------------------------------------------------------- covariates


def test_transform_covariates_values_and_median_imputation() -> None:
    cov = pd.DataFrame(  # one NaN per column, five non-NaN values (unambiguous median)
        {
            "length_aa": [100, 200, np.nan, 50, 400, 30],
            "n_observable_peptides": [0, 3, 7, np.nan, 1, 2],
            "gravy": [-0.5, np.nan, 0.2, 0.1, 0.0, 0.3],
            "pi": [5.0, 6.0, 7.0, 8.0, np.nan, 9.0],
            "mw_da": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            "tm_segments": [0, 1, 0, 3, 2, np.nan],
        },
        index=list("abcdef"),
    )
    out = _mod().transform_covariates(cov, list(dd.NAMES))
    ref = O.transform_covariates(cov, dd.NAMES)
    assert list(out.columns) == list(dd.NAMES)
    assert list(out.index) == list(cov.index)
    pd.testing.assert_frame_equal(
        out.astype(float), ref.astype(float), check_exact=False, atol=1e-12
    )
    assert out.loc["a", "length_aa"] == pytest.approx(np.log(100))
    assert out.loc["b", "n_observable_peptides"] == pytest.approx(np.log1p(3))
    assert out.loc["a", "gravy"] == -0.5 and out.loc["a", "pi"] == 5.0
    assert not out.isna().any().any()
    assert out.loc["c", "length_aa"] == pytest.approx(np.log(100))  # median of log lengths
    assert out.loc["b", "gravy"] == pytest.approx(0.1)
    assert out.loc["f", "tm_segments"] == 1.0


def test_zero_variance_covariate_dropped(caplog: pytest.LogCaptureFixture) -> None:
    gs = dd.synthetic(150, seed=51, tm_values=(0,))
    caplog.set_level(logging.DEBUG)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        full = _fit(gs)
    four = _fit(gs, DetectabilityConfig(covariates=[n for n in dd.NAMES if n != "tm_segments"]))
    assert full.constant is False
    np.testing.assert_allclose(
        full.p_detect.to_numpy(float), four.p_detect.to_numpy(float), atol=1e-12
    )
    assert _mentions("tm_segments", full, caplog, rec)


def test_nan_covariates_imputed_by_median() -> None:
    gs = dd.synthetic(151, seed=52)
    cov = gs.cov.copy()
    cov.iloc[[3, 17], cov.columns.get_loc("gravy")] = np.nan
    cov.iloc[[8], cov.columns.get_loc("pi")] = np.nan
    res = _mod().fit_detectability(
        gs.detected, cov, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", SEED)
    )
    filled = cov.copy()
    for c in ("gravy", "pi"):
        filled[c] = filled[c].fillna(filled[c].median())
    ref = _mod().fit_detectability(
        gs.detected, filled, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", SEED)
    )
    assert np.isfinite(res.p_detect.to_numpy(float)).all()
    np.testing.assert_allclose(
        res.p_detect.to_numpy(float), ref.p_detect.to_numpy(float), atol=1e-12
    )


# ----------------------------------------------------------------------------- edge cases


def test_too_few_tus_is_constant_with_warning(caplog: pytest.LogCaptureFixture) -> None:
    gs = dd.with_counts(100, 100, seed=61, genes_per_tu=50)  # 4 TUs < 5 folds
    caplog.set_level(logging.DEBUG)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = _fit(gs)
    assert res.constant is True
    assert np.all(res.p_detect.to_numpy(float) == 0.5)
    assert _mentions("too_few_tus", res, caplog, rec)


def test_degenerate_training_fold_uses_constant_prediction() -> None:
    base = dd.synthetic(200, seed=71)
    folds = O.tu_folds(list(base.tu_of_gene), 5, for_step("detectability", SEED))
    fold_of_gene = base.tu_of_gene.map(folds).to_numpy()
    d = fold_of_gene == 0  # every detected gene is held out in fold 0
    assert d.sum() >= 50 and (~d).sum() >= 50
    gs = dd.synthetic(200, seed=71, detected=d)
    res = _fit(gs)
    assert res.constant is False
    assert res.fold_of_tu == folds
    p0 = res.p_detect[d].to_numpy(float)
    assert np.unique(p0).size == 1  # constant prediction for the degenerate fold
    assert (res.diagnostics["metric"] == "degenerate_fold").any()


def test_convergence_warnings_not_emitted(caplog: pytest.LogCaptureFixture) -> None:
    from sklearn.exceptions import ConvergenceWarning

    gs = dd.synthetic(150, seed=81)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = _fit(gs)
    assert not [w for w in rec if issubclass(w.category, ConvergenceWarning)]
    assert np.isfinite(res.p_detect.to_numpy(float)).all()


def test_result_types() -> None:
    gs = dd.synthetic(120, seed=91)
    res = _fit(gs)
    assert isinstance(res, _mod().DetectabilityResult)
    assert isinstance(res.p_detect, pd.Series) and res.p_detect.dtype.kind == "f"
    assert isinstance(res.chosen_C, tuple)
    with pytest.raises((AttributeError, TypeError)):
        res.constant = True  # frozen


# ----------------------------------------------------------------------------- M4-17 input, speed


@pytest.mark.slow
def test_attack2_short_prophage_genes_have_lower_p_detect() -> None:
    gs = dd.sim_geneset(np.random.default_rng(2002), 2000, short_genes_in="prophage")
    res = _fit(gs)
    elig = gs.extra["eligible"]
    pro = gs.extra["in_prophage"] & elig
    rest = ~gs.extra["in_prophage"] & elig
    assert pro.sum() >= 10
    assert res.p_detect[pro].mean() < res.p_detect[rest].mean()


@pytest.mark.slow
def test_fit_4500_genes_under_20_seconds() -> None:
    gs = dd.sim_geneset(np.random.default_rng(45), 4500)
    t0 = time.perf_counter()
    _fit(gs)
    assert time.perf_counter() - t0 < 20.0
