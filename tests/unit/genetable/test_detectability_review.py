"""M4-13 lead decisions 10-14 after the code review of PR #83 (docs/specs/M4-13-detectability.md).

Spec: GIPPA_SPEC.md §7.3.3 (detectability model, edge cases), §5.2, §8.2, §9.3 (detection_summary
columns), §6.6. Decisions pinned here:

- (10) Real-contract inputs. ``fit_detectability(detected, cov, tus, cfg, rng)``:
  ``detected`` is a pandas ``boolean`` Series indexed by locus_tag (``<NA>`` = not eligible; M4-11
  ``detected_2pep`` must be NA-masked for partial genes by the caller, as the docstring says);
  ``cov`` is accepted exactly as M4-10 ``covariate_table`` returns it (``locus_tag`` COLUMN on a
  RangeIndex); ``tus`` is accepted as a C6 ``layers`` table (``locus_tag`` and ``tu_id`` COLUMNS,
  other columns ignored). Rows are matched by locus_tag, never by position. The index-keyed forms
  of the earlier tests (``cov`` indexed by locus_tag without a ``locus_tag`` column, TUs as a Series
  indexed by locus_tag) keep working and give identical results. Duplicate locus_tags (in
  ``detected``, ``cov`` or ``tus``) and eligible genes missing from ``cov``/``tus`` (or with a
  missing ``tu_id``) raise ``ValueError`` naming the tags (no median fill, no singleton TUs).
- (11) Constant predicted logits (after clipping) → NaN slope, intercept and their CIs.
- (12) ``folds < 2`` or too few TUs for the nested search (``folds=2`` with 2 TUs) → ``ValueError``
  mentioning the folds (not a bare ``min()`` crash). ``n_tus < folds`` keeps the spec edge case
  (constant + ``too_few_tus``; tested in test_detectability.py). The < 50 constant fallback is
  logged at WARNING and recorded in diagnostics.
- (13) ``diagnostics`` columns ``value``/``ci_low``/``ci_high`` are float64.
- (14) One test per surviving review mutant
  (local/scratch/review-m4-13-code/mutants.json): tie_larger_C, tus_le, raw_le, degenerate_half,
  inner_leak, inner_fresh_rng, z_196, conv_not_filtered, no_conv_row, no_dropped_row,
  tol_default, C_grid_short, tu_na_shared.
"""

from __future__ import annotations

import functools
import logging
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit, logit
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detectability_data as dd
import detectability_oracle as O
import evidence_oracle as eo

from gippa.config.model import DetectabilityConfig
from gippa.rng import for_step

C_GRID = {0.01, 0.1, 1.0, 10.0}
SEED = 7


def _mod() -> Any:
    from gippa.genetable import detectability

    return detectability


def _metric_rows(diag: pd.DataFrame, group_type: str, group: object, metric: str) -> pd.DataFrame:
    return diag[
        (diag["group_type"] == group_type)
        & (diag["group"].astype(str) == str(group))
        & (diag["metric"] == metric)
    ]


def _overall(diag: pd.DataFrame, metric: str) -> pd.Series:
    sel = diag[(diag["group_type"] == "overall") & (diag["metric"] == metric)]
    assert len(sel) == 1, f"{metric}: {len(sel)} rows"
    return sel.iloc[0]


def _fit(gs: dd.GeneSet, cfg: DetectabilityConfig | None = None, seed: int = SEED) -> Any:
    return _mod().fit_detectability(
        gs.detected,
        gs.cov,
        gs.tu_of_gene,
        cfg or DetectabilityConfig(),
        for_step("detectability", seed),
    )


# ----------------------------------------------------------------------------- (10) real contracts


@dataclass
class Contract:
    detected: pd.Series  # boolean, index locus_tag, NA = partial / pseudogene
    cov: pd.DataFrame  # M4-10 covariate_table output (locus_tag column, RangeIndex)
    tus: pd.DataFrame  # C6 layers table (locus_tag, tu_id columns)
    cfg: Any  # GippaConfig


@functools.lru_cache(maxsize=1)
def _contract() -> Contract:
    """gippasim M4-07 layout + M4-08 ``covariates_only`` detection → perfect C9/C11 evidence →
    M4-11 ``gene_evidence`` (``detected_2pep``), M4-10 ``covariate_table``, C6 ``layers``."""
    import gippasim.expression as ex
    import gippasim.layout as la
    from gippa.config import GippaConfig
    from gippa.genetable.covariates import covariate_table
    from gippa.genetable.evidence import accepted_peptides, gene_evidence
    from gippasim.contract_tables import layers_frame

    rng = np.random.default_rng(1313)
    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=500), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(detection="covariates_only"), rng)
    genes = lay.genes.copy()
    partial_rows = [3, 40, 41, 200]
    genes.loc[genes.index[partial_rows], "partial"] = True
    truth = dict(
        zip(expr.genes["locus_tag"], expr.genes["detected_true"].astype(bool), strict=True)
    )
    peps, hits = eo.perfect_evidence(genes, lay.replicons, truth, rng, n_false=10)
    cfg = GippaConfig()
    acc = accepted_peptides(hits, peps, cfg)
    ev = gene_evidence(genes, hits, acc, lay.replicons)
    det = pd.Series(
        pd.array(ev["detected_2pep"].to_numpy(bool), dtype="boolean"),
        index=pd.Index(ev["locus_tag"].astype(str), name="locus_tag"),
        name="detected",
    )
    not_eligible = (
        genes["partial"].to_numpy(bool) | (genes["feature_type"] == "pseudogene").to_numpy()
    )
    det[not_eligible] = pd.NA
    cov = covariate_table(genes, cfg)
    tus = layers_frame(lay)
    assert "locus_tag" in cov.columns and isinstance(cov.index, pd.RangeIndex)
    assert {"locus_tag", "tu_id"} <= set(tus.columns)
    elig = det.notna()
    assert int(det[elig].sum()) >= 50 and int((~det[elig].astype(bool)).sum()) >= 50
    assert tus["tu_id"].nunique() < len(tus)  # multi-gene TUs exist (grouping is testable)
    return Contract(det, cov, tus, cfg)


def _fit_contract(det: pd.Series, cov: pd.DataFrame, tus: Any, cfg: Any, seed: int = SEED) -> Any:
    return _mod().fit_detectability(
        det, cov, tus, cfg.detectability, for_step("detectability", seed)
    )


def test_real_contract_inputs_fit_with_tu_grouping() -> None:
    c = _contract()
    res = _fit_contract(c.detected, c.cov, c.tus, c.cfg)
    assert res.constant is False
    assert list(res.p_detect.index) == list(c.detected.index)
    na = c.detected.isna().to_numpy()
    assert res.p_detect[na].isna().all()
    p = res.p_detect[~na].to_numpy(float)
    assert np.isfinite(p).all() and ((p > 0) & (p < 1)).all()
    # folds are over the C6 TUs of eligible genes (not one "TU" per gene)
    tu_of = c.tus.set_index("locus_tag")["tu_id"]
    elig_tags = c.detected.index[~na]
    assert set(res.fold_of_tu) == set(tu_of.loc[elig_tags].astype(str))
    # not the all-NaN-covariates path: every configured covariate survives on real data
    assert not (res.diagnostics["metric"] == "dropped_zero_variance").any()
    assert _overall(res.diagnostics, "n_eligible")["value"] == len(elig_tags)


def test_real_contract_equals_index_keyed_inputs() -> None:
    c = _contract()
    a = _fit_contract(c.detected, c.cov, c.tus, c.cfg)
    cov_idx = c.cov.set_index("locus_tag")
    tu_ser = c.tus.set_index("locus_tag")["tu_id"]
    b = _fit_contract(c.detected, cov_idx, tu_ser, c.cfg)
    assert np.array_equal(a.p_detect.to_numpy(float), b.p_detect.to_numpy(float), equal_nan=True)
    assert a.fold_of_tu == b.fold_of_tu and a.chosen_C == b.chosen_C


def test_real_contract_rows_matched_by_locus_tag_not_position() -> None:
    c = _contract()
    a = _fit_contract(c.detected, c.cov, c.tus, c.cfg)
    r = np.random.default_rng(5)
    cov = c.cov.iloc[r.permutation(len(c.cov))].reset_index(drop=True)
    tus = c.tus.iloc[r.permutation(len(c.tus))].reset_index(drop=True)
    b = _fit_contract(c.detected, cov, tus, c.cfg)
    assert list(b.p_detect.index) == list(c.detected.index)
    assert np.array_equal(a.p_detect.to_numpy(float), b.p_detect.to_numpy(float), equal_nan=True)


def _eligible_tag(c: Contract, k: int = 7) -> str:
    return str(c.detected.index[c.detected.notna().to_numpy()][k])


@pytest.mark.parametrize("where", ["detected", "cov", "tus"])
def test_duplicate_locus_tags_raise_naming_them(where: str) -> None:
    c = _contract()
    tag = _eligible_tag(c)
    det, cov, tus = c.detected, c.cov, c.tus
    if where == "detected":
        det = pd.concat([det, det.loc[[tag]]])
    elif where == "cov":
        cov = pd.concat([cov, cov[cov["locus_tag"] == tag]], ignore_index=True)
    else:
        tus = pd.concat([tus, tus[tus["locus_tag"] == tag]], ignore_index=True)
    with pytest.raises(ValueError, match=tag):
        _fit_contract(det, cov, tus, c.cfg)


@pytest.mark.parametrize("where", ["cov", "tus", "tu_id_na"])
def test_eligible_gene_missing_raises_naming_it(where: str) -> None:
    c = _contract()
    tag, tag2 = _eligible_tag(c, 7), _eligible_tag(c, 90)
    cov, tus = c.cov, c.tus
    if where == "cov":
        cov = cov[~cov["locus_tag"].isin([tag, tag2])].reset_index(drop=True)
    elif where == "tus":
        tus = tus[~tus["locus_tag"].isin([tag, tag2])].reset_index(drop=True)
    else:
        tus = tus.copy()
        tus["tu_id"] = tus["tu_id"].astype(object)
        tus.loc[tus["locus_tag"].isin([tag, tag2]), "tu_id"] = None
    with pytest.raises(ValueError) as exc:
        _fit_contract(c.detected, cov, tus, c.cfg)
    assert tag in str(exc.value) and tag2 in str(exc.value)


def test_index_keyed_missing_tu_raises_not_singleton() -> None:
    """Mutant tu_na_shared (and the old silent singleton-TU fallback)."""
    gs = dd.synthetic(150, seed=21)
    tag = str(gs.detected.index[11])
    tu = gs.tu_of_gene.drop(tag)
    with pytest.raises(ValueError, match=tag):
        _mod().fit_detectability(
            gs.detected, gs.cov, tu, DetectabilityConfig(), for_step("detectability", SEED)
        )


def test_index_keyed_missing_covariates_raise_not_median_fill() -> None:
    gs = dd.synthetic(150, seed=22)
    tag = str(gs.detected.index[12])
    with pytest.raises(ValueError, match=tag):
        _mod().fit_detectability(
            gs.detected,
            gs.cov.drop(tag),
            gs.tu_of_gene,
            DetectabilityConfig(),
            for_step("detectability", SEED),
        )


def test_docstring_mentions_detected_2pep_masking() -> None:
    doc = _mod().fit_detectability.__doc__ or ""
    assert "detected_2pep" in doc and "partial" in doc


# ----------------------------------------------------------------------------- (11) constant logits


@pytest.mark.parametrize(
    "p",
    [np.full(100, 0.3), np.r_[np.zeros(50), np.full(50, 1e-9)]],  # 2nd: equal after clipping
    ids=["constant", "constant_after_clip"],
)
def test_constant_logits_give_nan_slope_intercept_and_cis(p: np.ndarray) -> None:
    d = np.tile([0, 1], 50)
    diag = _mod().calibration_diagnostics(d, p)
    for m in ("calibration_slope", "calibration_intercept"):
        row = _overall(diag, m)
        assert np.isnan(float(row["value"])), m
        assert np.isnan(float(row["ci_low"])) and np.isnan(float(row["ci_high"])), m
    assert np.isfinite(float(_overall(diag, "brier")["value"]))


def test_constant_fallback_writes_nan_recalibration() -> None:
    res = _fit(dd.with_counts(49, 951, seed=1))
    assert res.constant is True
    for m in ("calibration_slope", "calibration_intercept"):
        row = _overall(res.diagnostics, m)
        assert np.isnan(float(row["value"])) and np.isnan(float(row["ci_low"]))
        assert np.isnan(float(row["ci_high"]))


# ----------------------------------------------------------------------------- (12) folds, fallback


@pytest.mark.parametrize("folds", [0, 1])
def test_folds_below_two_raise(folds: int) -> None:
    gs = dd.with_counts(100, 100, seed=3)
    with pytest.raises(ValueError, match=r"(?i)fold"):
        _fit(gs, DetectabilityConfig(folds=folds))


def test_folds_too_many_for_nested_search_raise_value_error() -> None:
    """folds=2 with 2 TUs: an outer training set of one TU cannot feed the inner search."""
    gs = dd.with_counts(100, 100, seed=4, genes_per_tu=100)
    assert gs.tu_of_gene.nunique() == 2
    with pytest.raises(ValueError, match=r"(?i)fold"):
        _fit(gs, DetectabilityConfig(folds=2))


def test_small_n_fallback_logged_at_warning_and_in_diagnostics(
    caplog: pytest.LogCaptureFixture,
) -> None:
    caplog.set_level(logging.DEBUG)
    res = _fit(dd.with_counts(49, 951, seed=5))
    assert res.constant is True
    recs = [r for r in caplog.records if r.levelno == logging.WARNING]
    assert any("constant" in r.getMessage() for r in recs), [r.getMessage() for r in recs]
    assert float(_overall(res.diagnostics, "detectability_constant")["value"]) == 1.0


# ----------------------------------------------------------------------------- (13) dtypes


@pytest.mark.parametrize("case", ["calibration", "constant", "fitted", "fitted_extra_rows"])
def test_diagnostics_numeric_dtypes(case: str) -> None:
    if case == "calibration":
        diag = _mod().calibration_diagnostics([0, 0, 1, 1], [0.1, 0.2, 0.8, 0.9], n_bins=2)
    elif case == "constant":
        diag = _fit(dd.with_counts(49, 951, seed=6)).diagnostics
    elif case == "fitted":
        diag = _fit(dd.synthetic(150, seed=41)).diagnostics  # no extra rows
    else:
        diag = _fit(dd.synthetic(150, seed=51, tm_values=(0,))).diagnostics
    for col in ("value", "ci_low", "ci_high"):
        assert diag[col].dtype == np.float64, (col, diag[col].dtype)


# ----------------------------------------------------------------------------- (14) mutants


def test_c_tie_goes_to_smaller_c(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mutant tie_larger_C: a C-independent model makes every C tie → C = 0.01 in every fold."""
    from sklearn.dummy import DummyClassifier

    det = _mod()
    monkeypatch.setattr(
        det, "LogisticRegression", lambda *a, **kw: DummyClassifier(strategy="prior")
    )
    res = _fit(dd.synthetic(150, seed=61))
    assert res.chosen_C == (0.01,) * 5


def test_logistic_settings_and_full_c_grid(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mutants tol_default and C_grid_short: G5 / decision 9 settings, all four C values tried."""
    det = _mod()
    real = det.LogisticRegression
    calls: list[dict[str, Any]] = []

    def spy(*a: Any, **kw: Any) -> Any:
        calls.append(kw)
        return real(*a, **kw)

    monkeypatch.setattr(det, "LogisticRegression", spy)
    _fit(dd.synthetic(150, seed=62))
    assert calls
    assert {float(kw["C"]) for kw in calls} == C_GRID
    for kw in calls:
        assert kw.get("tol") == 1e-6 and kw.get("max_iter") == 2000
        assert kw.get("solver", "lbfgs") == "lbfgs"
        assert kw.get("l1_ratio") == 0.0 or kw.get("penalty") == "l2"


def test_exactly_k_tus_is_fitted() -> None:
    """Mutant tus_le: n_tus == folds is enough (one TU per fold)."""
    gs = dd.with_counts(100, 100, seed=63, genes_per_tu=40)
    assert gs.tu_of_gene.nunique() == 5
    res = _fit(gs)
    assert res.constant is False
    assert sorted(res.fold_of_tu.values()) == [0, 1, 2, 3, 4]
    assert np.isfinite(res.p_detect.to_numpy(float)).all()


def test_tm_segments_with_as_many_values_as_knots_gets_a_spline() -> None:
    """Mutant raw_le: raw only for FEWER distinct values than knots (5 values, 5 knots → spline)."""
    gs = dd.synthetic(150, seed=64, tm_values=(0, 1, 2, 3, 4))
    assert gs.cov["tm_segments"].nunique() == 5
    res = _fit(gs)
    p, folds, _, _ = O.cross_fit(
        gs.detected, gs.cov, gs.tu_of_gene, dd.NAMES, 5, 5, for_step("detectability", SEED),
        fixed_c=res.chosen_C,
    )  # fmt: skip
    assert res.fold_of_tu == folds
    np.testing.assert_allclose(res.p_detect.to_numpy(float), p.to_numpy(), rtol=0, atol=2e-3)


@pytest.mark.parametrize("held_out_detected", [True, False])
def test_degenerate_fold_predicts_training_detected_share(held_out_detected: bool) -> None:
    """Mutant degenerate_half: one-class training fold → its detected share (0 or 1), not 0.5."""
    base = dd.synthetic(200, seed=71)
    folds = O.tu_folds(list(base.tu_of_gene), 5, for_step("detectability", SEED))
    in0 = base.tu_of_gene.map(folds).to_numpy() == 0
    d = in0 if held_out_detected else ~in0
    gs = dd.synthetic(200, seed=71, detected=d)
    res = _fit(gs)
    assert res.constant is False
    expected = 0.0 if held_out_detected else 1.0  # the training share of fold 0
    assert np.all(res.p_detect[in0].to_numpy(float) == expected)
    assert len(_metric_rows(res.diagnostics, "fold", 0, "degenerate_fold")) == 1


def test_inner_search_uses_training_tus_and_the_step_stream(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Mutants inner_leak and inner_fresh_rng: the inner 3-fold split of outer fold j is
    ``tu_folds(training TUs of j, 3, rng)`` with the caller's generator (spec: same stream)."""
    det = _mod()
    real = det.tu_folds
    calls: list[tuple[set[str], int, Any]] = []

    def spy(tu_ids: Any, k: int, rng: Any) -> Any:
        calls.append(({str(t) for t in tu_ids}, int(k), rng))
        return real(tu_ids, k, rng)

    monkeypatch.setattr(det, "tu_folds", spy)
    gs = dd.synthetic(150, seed=72)
    rng = for_step("detectability", SEED)
    res = det.fit_detectability(gs.detected, gs.cov, gs.tu_of_gene, DetectabilityConfig(), rng)
    outer = [c for c in calls if c[1] == 5]
    inner = [c for c in calls if c[1] == 3]
    assert len(outer) == 1 and outer[0][2] is rng
    assert outer[0][0] == set(gs.tu_of_gene.astype(str))
    assert len(inner) == 5
    for j, (tus, _, r) in enumerate(inner):
        assert tus == {t for t, f in res.fold_of_tu.items() if f != j}, j
        assert r is rng, j


def test_ci_half_width_uses_exact_normal_quantile() -> None:
    """Mutant z_196: half-width / SE = z_0.975 to 1e-9 (SE from X'WX at the reported fit)."""
    r = np.random.default_rng(73)
    p = expit(r.normal(0.0, 1.5, 2000))
    d = (r.random(2000) < p).astype(float)
    diag = _mod().calibration_diagnostics(d, p)
    a = float(_overall(diag, "calibration_intercept")["value"])
    b = float(_overall(diag, "calibration_slope")["value"])
    z = logit(np.clip(p, 1e-6, 1 - 1e-6))
    x = np.column_stack([np.ones_like(z), z])
    mu = expit(a + b * z)
    se = np.sqrt(np.diag(np.linalg.inv(x.T @ (x * (mu * (1 - mu))[:, None]))))
    for m, s in (("calibration_intercept", se[0]), ("calibration_slope", se[1])):
        row = _overall(diag, m)
        lo, hi = float(row["ci_low"]), float(row["ci_high"])
        assert (hi - lo) / 2 / s == pytest.approx(float(norm.ppf(0.975)), rel=1e-9, abs=0)
        assert (hi + lo) / 2 == pytest.approx(float(row["value"]), rel=1e-12, abs=1e-12)


def test_convergence_warnings_captured_into_diagnostics(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mutants conv_not_filtered and no_conv_row: forced non-convergence is counted, not printed."""
    from sklearn.exceptions import ConvergenceWarning

    det = _mod()
    real = det.LogisticRegression
    monkeypatch.setattr(
        det, "LogisticRegression", lambda *a, **kw: real(*a, **{**kw, "max_iter": 1})
    )
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        res = _fit(dd.synthetic(150, seed=74))
    assert not [w for w in rec if issubclass(w.category, ConvergenceWarning)]
    row = _overall(res.diagnostics, "convergence_warning")
    assert float(row["value"]) >= 1


def test_zero_variance_covariate_has_its_diagnostics_row(caplog: pytest.LogCaptureFixture) -> None:
    """Mutant no_dropped_row: machine-readable row (decision 2) plus a WARNING log record."""
    caplog.set_level(logging.DEBUG)
    res = _fit(dd.synthetic(150, seed=75, tm_values=(0,)))
    assert res.constant is False
    sel = _metric_rows(res.diagnostics, "covariate", "tm_segments", "dropped_zero_variance")
    assert len(sel) == 1 and float(sel["value"].iloc[0]) == 1.0
    assert any(
        r.levelno == logging.WARNING and "tm_segments" in r.getMessage() for r in caplog.records
    )
