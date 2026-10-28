"""The cross-fitted TU recalibration of q (issue M8-04b; docs/specs/M8-04 lead decision 12,
ADR-0035 decisions 1-2).

``tu_universe(genes, folds)`` keeps the raw TU quantities

- ``log_comp`` = Σ log1p(-p̂) over the TU's eligible genes, p̂ clipped to [1e-12, 1 - 1e-12]
  (``clipped`` marks TUs with a clipped gene), ``q_raw`` = -expm1(log_comp), ``n_genes``;

and fits d_TU ~ logistic(β0 + β1 logit q_raw + β2 log n_genes) on the TUs of the other folds,
predicting each fold's TUs (M4-13's ``fold_of_tu``, so a TU's own detection never enters its own
fitted value). ``q`` is the fitted probability and ``comp`` its complement computed directly
(never 1 - q); ``detection_enrichment(..., folds=)`` tilts (q, comp) with
``tilt_pair(q, N, comp=comp)`` and passes the tilted complements to the conditional law.

Resolutions (reported to the lead):

- ``folds=None`` (M4-13 constant model: ``fold_of_tu`` is empty) and universes too small to fit
  (T6: 6 TUs) skip the recalibration: ``q`` is ``q_raw`` exactly and ``flags`` holds
  ``tu_recal_skipped`` (also in the ORA ``run_flags``). The minimum size is the implementer's
  choice; only "T6 is too small" and "~1,300 TUs are enough" are pinned.
- A TU of the universe missing from a non-empty ``folds`` is a ``ValueError`` naming it (fail
  closed, as M4-16).
- When the training TUs all have one gene (log n constant), the size term is dropped rather than
  failing (the fit still runs).
- "Identity when calibrated" is checked on 6,000-gene draws with p̂ = p_true: mean |q - q_raw|
  < 0.03 and median |logit q - logit q_raw| < 0.2 (sampling noise of a 3-parameter fit on ~2,700
  TUs is ~0.01 / ~0.06; the TU-packed bias it corrects is ~0.19 on large TUs).
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from helpers.ora_calsim import collection, simulate, whole_tu_sets
from helpers.ora_data import ora, singletons, t6_frame
from scipy.special import logit

from gippa.genetable.detectability import tu_folds
from gippa.stats.poisson_binomial import conditional_poisson_binomial, tilt_pair


def _sim(seed: int, **kw):  # type: ignore[no-untyped-def]
    return simulate(np.random.default_rng(seed), **kw)


# ----------------------------------------------------------------------------- the skip path


def test_t6_with_folds_is_too_small_and_skips() -> None:
    g = t6_frame()
    folds = tu_folds([f"TU_{i}" for i in range(1, 7)], 5, np.random.default_rng(1))
    u = ora().tu_universe(g, folds)
    np.testing.assert_array_equal(u.q, u.q_raw)
    assert "tu_recal_skipped" in u.flags


def test_skipped_recalibration_is_a_run_flag() -> None:
    from helpers.ora_data import kegg

    out = ora().detection_enrichment(
        t6_frame(), [kegg({"x": {"g1", "g4", "g7"}})], (1, 500), 0.05, folds=None
    )
    assert "tu_recal_skipped" in out.attrs["run_flags"]


def test_a_tu_missing_from_the_folds_fails_closed() -> None:
    s = _sim(5)
    folds = dict(s.folds)
    folds.pop("TU_7")
    with pytest.raises(ValueError, match="TU_7"):
        ora().tu_universe(s.genes, folds)


# ------------------------------------------------------------------------------ the fit


def test_fitted_q_is_not_the_raw_q_and_is_open() -> None:
    s = _sim(11)
    u = ora().tu_universe(s.genes, s.folds)
    assert "tu_recal_skipped" not in u.flags
    assert not np.array_equal(u.q, u.q_raw)
    assert np.all(u.q > 0.0) and np.all(u.comp > 0.0)
    assert np.all(u.q <= 1.0) and np.all(u.comp < 1.0)
    assert np.max(np.abs(u.q + u.comp - 1.0)) <= 4 * np.finfo(np.float64).eps
    np.testing.assert_allclose(u.q_raw, -np.expm1(u.log_comp), rtol=1e-15)


def test_raw_fields_match_the_frame() -> None:
    s = _sim(12)
    u = ora().tu_universe(s.genes, s.folds)
    assert u.tu_ids == tuple(s.tu_ids)
    np.testing.assert_array_equal(u.n_genes, s.sizes)
    np.testing.assert_array_equal(u.detected, s.d_tu)
    lc = np.bincount(s.tu_of_gene, weights=np.log1p(-s.genes["p_detect"].to_numpy()))
    np.testing.assert_allclose(u.log_comp, lc, rtol=1e-12)


def test_cross_fitting_a_tus_own_detection_never_affects_its_q() -> None:
    """Flip the detection of every gene of one TU: the fitted q of every TU in that TU's fold is
    unchanged bit for bit (their model never saw it); TUs of other folds do move."""
    s = _sim(21)
    o = ora()
    u0 = o.tu_universe(s.genes, s.folds)
    t = int(np.flatnonzero(~s.d_tu & (s.sizes >= 2))[0])  # an undetected multi-gene TU
    g1 = s.genes.copy()
    g1.loc[g1["tu_id"] == f"TU_{t + 1}", "detected"] = True
    u1 = o.tu_universe(g1, s.folds)
    fold = np.array([s.folds[x] for x in s.tu_ids])
    same = fold == fold[t]
    np.testing.assert_array_equal(u1.q[same], u0.q[same])
    np.testing.assert_array_equal(u1.comp[same], u0.comp[same])
    assert not np.array_equal(u1.q[~same], u0.q[~same])
    assert u1.detected[t] and not u0.detected[t]


def test_fold_assignment_is_what_the_fit_uses() -> None:
    """Re-drawing the folds changes the fitted q (it is not an in-sample fit)."""
    s = _sim(22)
    other = tu_folds(s.tu_ids, 5, np.random.default_rng(999))
    a = ora().tu_universe(s.genes, s.folds)
    b = ora().tu_universe(s.genes, other)
    assert not np.array_equal(a.q, b.q)
    np.testing.assert_array_equal(a.q_raw, b.q_raw)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_recalibration_is_the_identity_when_p_hat_is_calibrated(seed: int) -> None:
    s = _sim(seed, n_genes=6000, oracle=True)
    u = ora().tu_universe(s.genes, s.folds)
    assert "tu_recal_skipped" not in u.flags
    dq = np.abs(u.q - u.q_raw)
    dl = np.abs(logit(u.q) - logit(u.q_raw))
    assert float(dq.mean()) < 0.03, float(dq.mean())
    assert float(np.median(dl)) < 0.2, float(np.median(dl))


def test_recalibration_removes_the_tu_packed_size_bias() -> None:
    """c = 2.5 with p̂ the exact marginal: raw q overstates large TUs (mean d - q ≈ -0.19); the
    fitted q is close to the observed TU detection rate by size class."""
    s = _sim(7, n_genes=6000, c=2.5)
    u = ora().tu_universe(s.genes, s.folds)
    d = u.detected.astype(float)
    big, one = u.n_genes >= 4, u.n_genes == 1
    assert float((d[big] - u.q_raw[big]).mean()) < -0.1  # the bias is there (teeth)
    assert abs(float((d[big] - u.q[big]).mean())) < 0.05
    assert abs(float((d[one] - u.q[one]).mean())) < 0.05


def test_all_singleton_tus_drop_the_size_term() -> None:
    rng = np.random.default_rng(31)
    p = rng.uniform(0.05, 0.9, 1500)
    det = np.flatnonzero(rng.random(1500) < p)
    g = singletons(p, det.tolist())
    folds = tu_folds([f"TU_{i}" for i in range(1, 1501)], 5, rng)
    u = ora().tu_universe(g, folds)
    assert "tu_recal_skipped" not in u.flags
    assert np.all(np.isfinite(u.q)) and np.all((u.q > 0) & (u.q < 1))
    assert float(np.abs(u.q - u.q_raw).mean()) < 0.03


def test_deterministic() -> None:
    s = _sim(41)
    a = ora().tu_universe(s.genes, s.folds)
    b = ora().tu_universe(s.genes.sample(frac=1.0, random_state=1).reset_index(drop=True), s.folds)
    np.testing.assert_array_equal(a.q, b.q)
    np.testing.assert_array_equal(a.comp, b.comp)


# ----------------------------------------------------------------- ORA uses (q, comp) tilted


def test_ora_p_is_the_conditional_law_of_the_fitted_q_with_its_complement() -> None:
    s = _sim(51)
    rng = np.random.default_rng(52)
    sets = {**whole_tu_sets(s, rng, "large", 3), **whole_tu_sets(s, rng, "single", 3)}
    out = ora().detection_enrichment(s.genes, [collection(sets)], (5, 500), 0.05, folds=s.folds)
    assert "tu_recal_skipped" not in out.attrs["run_flags"]
    u = ora().tu_universe(s.genes, s.folds)
    N = int(u.detected.sum())
    qt, ct = tilt_pair(u.q, N, comp=u.comp)
    for term, genes in sets.items():
        m = np.zeros(len(u.tu_ids), dtype=bool)
        m[list(ora().tu_members(u, genes))] = True
        law = conditional_poisson_binomial(
            qt[m], qt[~m], int(u.detected[m].sum()), N, comp_in=ct[m], comp_out=ct[~m]
        )
        r_over = out[(out["term"] == term) & (out["direction"] == "over")].iloc[0]
        r_under = out[(out["term"] == term) & (out["direction"] == "under")].iloc[0]
        assert r_over["p_value"] == pytest.approx(law.p_over, rel=1e-11)
        assert r_under["p_value"] == pytest.approx(law.p_under, rel=1e-11)
        assert r_over["effect"] == pytest.approx(r_over["k_detected_in_set"] - law.mean, rel=1e-9)
        raw = conditional_poisson_binomial(
            u.q_raw[m], u.q_raw[~m], int(u.detected[m].sum()), N
        ).p_over
        if term.startswith("large"):
            assert not math.isclose(r_over["p_value"], raw, rel_tol=1e-6)


def test_folds_none_in_ora_uses_the_raw_q() -> None:
    s = _sim(53)
    rng = np.random.default_rng(54)
    sets = whole_tu_sets(s, rng, "large", 2)
    out = ora().detection_enrichment(s.genes, [collection(sets)], (5, 500), 0.05, folds=None)
    assert "tu_recal_skipped" in out.attrs["run_flags"]
    u = ora().tu_universe(s.genes, folds=None)
    N = int(u.detected.sum())
    for term, genes in sets.items():
        m = np.zeros(len(u.tu_ids), dtype=bool)
        m[list(ora().tu_members(u, genes))] = True
        law = conditional_poisson_binomial(u.q_raw[m], u.q_raw[~m], int(u.detected[m].sum()), N)
        r = out[(out["term"] == term) & (out["direction"] == "under")].iloc[0]
        assert r["p_value"] == pytest.approx(law.p_under, rel=1e-11)


def test_saturated_tu_keeps_its_complement() -> None:
    """A 40-gene TU at p̂ 0.6 has q_raw = 1 - 1.2e-16 (rounds near 1): its comp stays positive
    through the fit and the ORA call does not raise when that TU is undetected."""
    s = _sim(61)
    g = s.genes
    big = pd.DataFrame(
        {
            "locus_tag": [f"big{i}" for i in range(40)],
            "tu_id": "TU_99999",
            "tu_position": np.arange(1, 41),
            "detected": False,
            "p_detect": 0.6,
            "partial": False,
            "abundance_log2": np.nan,
            "abundance_imputed": False,
            "abundance_scale": "molar",
        }
    )
    g2 = pd.concat([g, big], ignore_index=True)
    folds = {**s.folds, "TU_99999": 0}
    u = ora().tu_universe(g2, folds)
    i = u.tu_ids.index("TU_99999")
    assert u.comp[i] > 0.0 and u.log_comp[i] == pytest.approx(40 * math.log(0.4), rel=1e-13)
    sets = {"has_big": frozenset(big["locus_tag"]) | s.genes["locus_tag"].iloc[:30].pipe(frozenset)}
    out = ora().detection_enrichment(g2, [collection(sets)], (1, 500), 0.05, folds=folds)
    assert np.all(np.isfinite(out["p_value"]))
