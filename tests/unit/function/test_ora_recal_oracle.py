"""Independent per-fold oracle for the TU recalibration (docs/specs/M8-04 lead decision 13, "an
independent per-fold oracle ... scipy.optimize, ≈ 1e-8, with a case where the size term
matters"; decision 12, ADR-0035).

For each fold j of M4-13's TU folds, d_TU ~ logistic(b0 + b1 logit q_raw + b2 log n_genes) is
fitted on the training TUs (fold != j) only, and the held-out TUs (fold == j) get
q = expit(eta), comp = expit(-eta). The design is rebuilt here from the gene frame (clip p̂ to
[1e-12, 1 - 1e-12], log comp = Σ log1p(-p̂), logit q_raw = log(-expm1(log comp)) - log comp),
not read back from ``TuUniverse``. The MLE is found by ``scipy.optimize.minimize`` (trust-exact on
the mean negative log-likelihood) and then refined with ``scipy.optimize.root`` on the score
equations: the trust region stalls at a gradient of ~1e-10 because the objective's float noise
hides further improvement, while the score root reaches ~1e-15. Both are scipy; neither is the
implementation's Newton loop. ``u.q`` and ``u.comp`` must match to a relative 1e-8.

Teeth (each checked, so the worlds discriminate):

- **The size term matters.** In the TU-packed world (c = 2.5) dropping log n, or using raw n
  instead of log n, moves q by > 1e-2 relative, far beyond the tolerance.
- **The fit must converge tightly.** A step-halving Newton from (0, 1, 0) stopped when its step
  is below 1e-2 moves q by > 1e-7 relative (observed 2e-6 to 2e-5), so a loose tolerance fails.
- **Badly scaled q_raw** (c = 2.5, intercept +0.5, slope 0.8): the MLE exists; the oracle must
  be met there too (HEAD's undamped Newton diverges on seed 1 and falls back to the raw q).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from functools import cache
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.ora_calsim import simulate
from helpers.ora_data import ora
from scipy.optimize import minimize, root
from scipy.special import expit

CLIP = 1e-12
RTOL = 1e-8
Fit = Callable[[np.ndarray, np.ndarray], np.ndarray]

WORLDS = {
    "tu_packed_s1": ({"c": 2.5}, 1),
    "tu_packed_s2": ({"c": 2.5}, 2),
    "badly_scaled_s0": ({"c": 2.5, "a": 0.5, "b": 0.8}, 0),
    "badly_scaled_s1": ({"c": 2.5, "a": 0.5, "b": 0.8}, 1),
}


def _design(genes: pd.DataFrame, ids: list[str], size: str | None) -> tuple[Any, Any]:
    el = genes[~genes["partial"].astype(bool)]
    p = np.clip(el["p_detect"].to_numpy(float), CLIP, 1.0 - CLIP)
    index = {t: i for i, t in enumerate(ids)}
    tu = np.array([index[t] for t in el["tu_id"]])
    log_comp = np.bincount(tu, weights=np.log1p(-p), minlength=len(ids))
    n = np.bincount(tu, minlength=len(ids)).astype(float)
    d = np.bincount(tu, weights=el["detected"].astype(float).to_numpy(), minlength=len(ids)) > 0
    lq = np.log(-np.expm1(log_comp)) - log_comp
    cols = [np.ones_like(lq), lq]
    if size == "log":
        cols.append(np.log(n))
    elif size == "raw":
        cols.append(n)
    return np.column_stack(cols), d.astype(float)


def _scipy_mle(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    m = y.size

    def nll(b: np.ndarray) -> float:
        e = x @ b
        return float(np.sum(np.logaddexp(0.0, e) - y * e)) / m

    def score(b: np.ndarray) -> np.ndarray:
        return x.T @ (expit(x @ b) - y) / m

    def fisher(b: np.ndarray) -> np.ndarray:
        e = x @ b
        return x.T @ (x * (expit(e) * expit(-e))[:, None]) / m

    r = minimize(nll, np.zeros(x.shape[1]), jac=score, hess=fisher, method="trust-exact",
                 options={"gtol": 1e-12})  # fmt: skip
    s = root(score, r.x, jac=fisher, method="hybr", options={"xtol": 1e-15})
    assert np.max(np.abs(score(s.x))) < 1e-13, np.max(np.abs(score(s.x)))
    return np.asarray(s.x)


def _loose_newton(x: np.ndarray, y: np.ndarray, tol: float = 1e-2) -> np.ndarray:
    """Step-halving Newton from (0, 1, 0), stopped once max |step| < ``tol``."""

    def nll(b: np.ndarray) -> float:
        e = x @ b
        return float(np.sum(np.logaddexp(0.0, e) - y * e))

    b = np.zeros(x.shape[1])
    b[1] = 1.0
    for _ in range(500):
        mu = expit(x @ b)
        step = np.linalg.solve(x.T @ (x * (mu * (1.0 - mu))[:, None]), x.T @ (y - mu))
        t = 1.0
        while nll(b + t * step) > nll(b) and t > 1e-12:
            t /= 2.0
        b = b + t * step
        if np.max(np.abs(t * step)) < tol:
            return b
    raise AssertionError("loose Newton did not stop")


def _cross_fit(
    genes: pd.DataFrame, ids: list[str], folds: Mapping[str, int], fit: Fit, size: str | None
) -> tuple[np.ndarray, np.ndarray]:
    x, y = _design(genes, ids, size)
    fold = np.array([int(folds[t]) for t in ids])
    eta = np.full(len(ids), np.nan)
    for j in np.unique(fold):
        train, test = fold != j, fold == j
        assert np.ptp(x[train, -1]) > 0  # sizes vary in these worlds: the size term is fitted
        eta[test] = x[test] @ fit(x[train], y[train])
    return expit(eta), expit(-eta)


def _rel(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(a - b) / b))


@cache
def _world(name: str) -> dict[str, Any]:
    kw, seed = WORLDS[name]
    sim = simulate(np.random.default_rng(seed), **kw)
    ids = sorted({str(t) for t in sim.genes["tu_id"]}, key=lambda t: int(t.split("_")[1]))
    q, comp = _cross_fit(sim.genes, ids, sim.folds, _scipy_mle, "log")
    return {"sim": sim, "ids": ids, "q": q, "comp": comp}


@pytest.fixture(params=list(WORLDS))
def world(request: pytest.FixtureRequest) -> dict[str, Any]:
    return _world(request.param)


def test_fitted_q_and_comp_match_the_per_fold_scipy_oracle(world: dict[str, Any]) -> None:
    sim = world["sim"]
    u = ora().tu_universe(sim.genes, folds=sim.folds)
    assert list(u.tu_ids) == world["ids"]
    assert "tu_recal_skipped" not in u.flags
    np.testing.assert_allclose(u.q, world["q"], rtol=RTOL, atol=0.0, err_msg="q")
    np.testing.assert_allclose(u.comp, world["comp"], rtol=RTOL, atol=0.0, err_msg="comp")


@pytest.mark.parametrize("size", [None, "raw"])
@pytest.mark.parametrize("name", ["tu_packed_s1", "tu_packed_s2"])
def test_the_size_term_and_its_log_scale_matter(name: str, size: str | None) -> None:
    """Teeth: the oracle's design is not interchangeable with a smaller or raw-n design here."""
    world = _world(name)
    sim = world["sim"]
    q, _ = _cross_fit(sim.genes, world["ids"], sim.folds, _scipy_mle, size)
    assert _rel(q, world["q"]) > 1e-2


def test_a_loose_newton_tolerance_misses_the_oracle(world: dict[str, Any]) -> None:
    """Teeth: stopping at a step of 1e-2 is detectably off at the 1e-8 tolerance."""
    sim = world["sim"]
    q, comp = _cross_fit(sim.genes, world["ids"], sim.folds, _loose_newton, "log")
    assert max(_rel(q, world["q"]), _rel(comp, world["comp"])) > 10 * RTOL
