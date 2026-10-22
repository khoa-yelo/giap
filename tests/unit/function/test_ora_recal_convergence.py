"""The TU recalibration fit converges wherever the logistic MLE exists (statistician review of
PR #179, B4).

``_logistic_mle`` runs undamped Newton from beta = (0, 1, 0). When q_raw is badly scaled (the
4b world plus an M4-13 miscalibration: c = 2.5, intercept +0.5, slope 0.8), the iterates
oscillate and diverge (|beta| ~ 1e5 by iteration 5, then a singular Hessian). The fit returns
None and ``tu_universe`` silently falls back to the raw q, with ``tu_recal_skipped``. In the
review that happened in 2,898 of 3,000 replicates (``helpers.ora_calsim``: 38 of 40 seeds),
and the raw q then rejects ``under`` on large TUs at 0.49 for alpha 0.01. The MLE exists
(BFGS from 0 converges; seed 0 fold 0: beta = (-0.305, 0.55, -0.528)). The fix is a damped
Newton step (halve it until the log-likelihood increases).

Proposed repo path: tests/unit/function/test_ora_recal_convergence.py
"""

from __future__ import annotations

import numpy as np
import pytest
from helpers.ora_calsim import simulate
from helpers.ora_data import ora
from scipy.optimize import minimize
from scipy.special import expit

SEEDS = [1, 2, 3]


def _reference_q(u, folds) -> np.ndarray:  # type: ignore[no-untyped-def]
    """Cross-fitted q from the exact MLE of d ~ logistic(1, logit q_raw, log n).

    Trust-region Newton.
    """
    lq = np.log(u.q_raw) - u.log_comp
    x = np.column_stack([np.ones_like(lq), lq, np.log(u.n_genes.astype(float))])
    y = u.detected.astype(float)
    fold = np.array([folds[t] for t in u.tu_ids])
    eta = np.empty_like(lq)
    for j in np.unique(fold):
        tr = fold != j

        def nll(b, xt=x[tr], yt=y[tr]):  # type: ignore[no-untyped-def]
            e = xt @ b
            return float(np.sum(np.logaddexp(0.0, e) - yt * e))

        def grad(b, xt=x[tr], yt=y[tr]):  # type: ignore[no-untyped-def]
            return xt.T @ (expit(xt @ b) - yt)

        def hess(b, xt=x[tr]):  # type: ignore[no-untyped-def]
            mu = expit(xt @ b)
            return xt.T @ (xt * (mu * (1.0 - mu))[:, None])

        r = minimize(nll, np.zeros(3), jac=grad, hess=hess, method="trust-exact",
                     options={"gtol": 1e-9})  # fmt: skip
        b = r.x
        for _ in range(5):  # polish: Newton converges quadratically from the trust-region optimum
            b = b - np.linalg.solve(hess(b), grad(b))
        assert np.max(np.abs(grad(b))) < 1e-8 * tr.sum()
        eta[fold == j] = x[fold == j] @ b
    return expit(eta)


@pytest.mark.parametrize("seed", SEEDS)
def test_recalibration_converges_on_a_badly_scaled_q(seed: int) -> None:
    sim = simulate(np.random.default_rng(seed), c=2.5, a=0.5, b=0.8)
    u = ora().tu_universe(sim.genes, sim.folds)
    assert "tu_recal_skipped" not in u.flags
    np.testing.assert_allclose(u.q, _reference_q(u, sim.folds), rtol=1e-6, atol=1e-9)
