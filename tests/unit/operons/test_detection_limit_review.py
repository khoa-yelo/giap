"""M6-02 review (#157): a non-finite abundance scores NaN in both models (never raises), the
k_min = 2 conditional fit recovers the abundance slope, and a Newton step halved to nothing at
a non-stationary point is not reported as an MLE."""

from __future__ import annotations

import math

import numpy as np
import pytest
from helpers import detection_limit_data as D
from unit.operons.test_detection_limit_conditional import BETA, SEED, _dl, _rng

from giap.config.model import GiapConfig

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
    m = _dl().fit_peptide_model(train.genes, train.hits, GiapConfig(seed=SEED), _rng())
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
