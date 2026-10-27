"""M6-02 review (#157): a non-finite abundance scores NaN in both models (never raises), the
k_min = 2 conditional fit recovers the abundance slope, and a Newton step halved to nothing at
a non-stationary point is not reported as an MLE."""

from __future__ import annotations

import math

import numpy as np
import pytest
from helpers import detection_limit_data as D
from unit.operons.test_detection_limit_conditional import BETA, SEED, _dl, _rng

from gippa.config.model import GippaConfig

PEP_COEF = {"intercept": -10.0, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}


@pytest.mark.parametrize("a", [math.nan, math.inf, -math.inf])
def test_peptide_model_non_finite_abundance_is_nan(a: float) -> None:
    m = _dl().PMissModel("peptide", PEP_COEF, 1, 100, 0.8, True)
    assert math.isnan(_dl().p_detect(m, a, list(D.TOY_PEPS)))
    assert math.isnan(_dl().p_miss(m, 16.0, [0.0, a - 16.0], list(D.TOY_PEPS)))


@pytest.mark.parametrize("a", [math.nan, math.inf])
def test_replicate_model_non_finite_abundance_is_nan(a: float) -> None:
    m = _dl().PMissModel("replicate", {"intercept": -5.0, "abundance_log2": 0.4}, 1, 100, 0.8, True)
    assert math.isnan(_dl().p_detect(m, a, []))


def test_peptide_conditional_fit_k_min_2_recovers_abundance() -> None:
    """Decision 5 at k_min = 2 (two_peptide): the leave-one-out factor of the score matters;
    the k = 1 form gives beta_a ~0.44 here."""
    rng = np.random.default_rng(SEED + 2)
    train = D.simulate_peptides(rng, 3000, BETA, n_peps=(5, 20), k_min=2, tier="two_peptide")
    m = _dl().fit_peptide_model(train.genes, train.hits, GippaConfig(seed=SEED), _rng())
    assert m.usable is True
    assert m.k_min == 2
    assert m.coef["abundance_log2"] == pytest.approx(0.5, abs=0.04)


def test_newton_rejects_a_non_stationary_point() -> None:
    """f(b) = -(b^2 - 1)^2 + 0.1 b from b = 0: the Hessian there is positive, the Newton step
    points downhill and is halved to nothing; that is not an MLE."""

    def ll(eta: np.ndarray) -> tuple[float, np.ndarray]:
        b = float(eta[0])
        return -((b * b - 1.0) ** 2) + 0.1 * b, np.array([-4.0 * b * (b * b - 1.0) + 0.1])

    out = _dl()._newton(np.eye(1), ll)
    assert out is None or abs(-4.0 * out[0] * (out[0] ** 2 - 1.0) + 0.1) < 1e-6


def test_training_ignores_the_gene_level_call() -> None:
    """Statistician review B1: genes with >= k_min peptides seen train the model even when the
    gene-level FDR call failed (mayu/decoy_picked single-hit genes), so the fit equals the one
    with every such gene called."""
    rng = np.random.default_rng(SEED + 3)
    train = D.simulate_peptides(rng, 1500, BETA, n_peps=(5, 20))
    seen_genes = set(train.hits["gene"].astype(str))
    in_seen = train.genes.locus_tag.astype(str).isin(seen_genes).to_numpy()
    drop = in_seen & (np.random.default_rng(7).random(len(train.genes)) < 0.5)
    uncalled = train.genes.assign(detected=train.genes["detected"].to_numpy() & ~drop)
    assert drop.sum() > 100
    cfg = GippaConfig(seed=SEED)
    m_all = _dl().fit_peptide_model(train.genes, train.hits, cfg, _rng())
    m_sub = _dl().fit_peptide_model(uncalled, train.hits, cfg, _rng())
    assert m_sub.usable and m_sub.n_fit == m_all.n_fit
    assert m_sub.coef == pytest.approx(m_all.coef)
