"""M6-02 lead decision 5: the conditional (truncated) likelihood estimator
(docs/specs/M6-02-detection-limit-model.md lead decisions 5-6; docs/methods/detection-limit.md).

Pinned beyond test_detection_limit.py:

- ``PMissModel.reason``: why a model is not usable; ``"replicate_two_runs"`` for a replicate model
  at exactly 2 runs, decided from the run count with no fit attempted (``n_fit == 0``, empty
  ``coef`` and ``fold_coef``).
- Replicate model from 3 runs up: the doubly truncated binomial recovers the per-run abundance
  slope (truth 0.5; the plain logistic gives 0.149 at 3 runs, 0.245 at 4).
- Peptide model: the conditional fit recovers the abundance coefficient within 10 % (plain fit
  ~0.35) and is calibrated in the large (CITL within 0.1 of 0; plain fit ~-0.89) on an independent
  genome.
- Separated peptide-model data -> ``usable = False`` (unpenalised fit, no ridge).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import detection_limit_data as D

from gippa.config import GippaConfig
from gippa.rng import for_step

SEED = 20261021
BETA = {"intercept": -10.0, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}


def _dl() -> Any:
    from gippa.operons import detection_limit

    return detection_limit


def _rng(seed: int = SEED) -> np.random.Generator:
    return for_step(D.STEP, seed)


def _replicate_world(
    n: int, n_runs: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, np.ndarray]:
    """Per-run detection logit = -8 + 0.5 a, a ~ N(16, 2); one assigned peptide per gene."""
    rng = np.random.default_rng(seed)
    tags = [f"R{i:05d}" for i in range(n)]
    a = rng.normal(16.0, 2.0, n)
    peps = {t: D.make_peptides(rng, 2) for t in tags}
    q = 1.0 / (1.0 + np.exp(-(-8.0 + 0.5 * a)))
    seen = rng.random((n, n_runs)) < q[:, None]
    det = seen.any(axis=1)
    genes = D.genes_frame(
        tags, [f"T{i}" for i in range(n)], [D.protein(peps[t]) for t in tags], det.tolist(), a
    )
    hits = D.hits_frame([D.hit_row(peps[t][0], t) for t, d in zip(tags, det, strict=True) if d])
    runs = D.runs_frame(
        [(peps[t][0], f"run{r}") for i, t in enumerate(tags) for r in range(n_runs) if seen[i, r]]
    )
    return genes, hits, runs, seen


_NO_COV = {"detectability": {"covariates": []}}


# ---------------------------------------------------------------------------- replicate, R = 2


def test_pmissmodel_has_reason_field() -> None:
    m = _dl().PMissModel(
        kind="peptide", coef=dict(D.TOY_COEF), k_min=1, n_fit=1, auc=0.8, usable=True
    )
    assert hasattr(m, "reason")


def test_replicate_two_runs_not_usable_no_fit() -> None:
    """Lead decision 5: at exactly 2 runs the replicate model is not fitted (the mixed-gene
    training set has S = 1 of 2 for every gene: slope not identified)."""
    genes, hits, runs, seen = _replicate_world(1500, 2, seed=31)
    assert ((seen.sum(1) == 1).sum()) > 200  # plenty of 'mixed' genes: data are not the reason
    m = _dl().fit_replicate_model(genes, runs, GippaConfig(**_NO_COV), _rng(), peptide_hits=hits)
    assert m.kind == "replicate"
    assert m.usable is False
    assert m.reason == "replicate_two_runs"
    assert m.n_fit == 0
    assert dict(m.coef) == {}
    assert tuple(m.fold_coef) == ()
    assert math.isnan(_dl().p_detect(m, 16.0, [(10, 0.0)]))


def test_replicate_two_runs_with_default_covariates_same_reason() -> None:
    genes, hits, runs, _ = _replicate_world(600, 2, seed=32)
    rng = np.random.default_rng(32)
    genes = genes.assign(
        n_observable_peptides=2,
        length_aa=[len(x) for x in genes.aa],
        gravy=rng.normal(-0.3, 0.3, len(genes)),
        pi=rng.normal(7.0, 1.5, len(genes)),
        tm_segments=rng.integers(0, 3, len(genes)),
    )
    m = _dl().fit_replicate_model(genes, runs, GippaConfig(), _rng(), peptide_hits=hits)
    assert m.usable is False and m.reason == "replicate_two_runs" and m.n_fit == 0


def test_dispatcher_two_runs_replicate_reason_recorded() -> None:
    """With peptides the peptide model is in force; the replicate slot carries the reason."""
    s = D.simulate_peptides(np.random.default_rng(33), 300, BETA)
    rows = [
        (p, f"run{r}")
        for t in s.genes.locus_tag
        for p, o in zip(s.peptides[t], s.observed[t], strict=True)
        if o
        for r in range(2)
    ]
    fit = _dl().fit_detection_limit(s.genes, s.hits, D.runs_frame(rows), GippaConfig(), _rng())
    assert fit.model.kind == "peptide"
    assert fit.replicate is None or (
        fit.replicate.usable is False and fit.replicate.reason == "replicate_two_runs"
    )


# ---------------------------------------------------------------------------- replicate, R >= 3


@pytest.mark.parametrize("n_runs", [3, 4])
def test_replicate_truncated_binomial_recovers_slope(n_runs: int) -> None:
    """Doubly truncated binomial (1 <= S <= R - 1) on 4,000 genes: abundance slope within 0.12 of
    the truth 0.5 (oracle SE ~0.03 at R = 3). Teeth: the plain logistic gives 0.149 (R = 3) and
    0.245 (R = 4), both outside."""
    genes, hits, runs, seen = _replicate_world(4000, n_runs, seed=34 + n_runs)
    m = _dl().fit_replicate_model(genes, runs, GippaConfig(**_NO_COV), _rng(), peptide_hits=hits)
    s = seen.sum(1)
    assert m.usable is True and m.kind == "replicate"
    assert m.n_fit == int(((s > 0) & (s < n_runs)).sum()) * n_runs
    assert m.coef["abundance_log2"] == pytest.approx(0.5, abs=0.12)
    # intercept follows (truth -8; logit at a = 16 is 0)
    eta16 = m.coef["intercept"] + 16.0 * m.coef["abundance_log2"]
    assert eta16 == pytest.approx(0.0, abs=0.3)


# ---------------------------------------------------------------------------- peptide model


def _citl(d: np.ndarray, p: np.ndarray) -> float:
    """Calibration in the large: c in d ~ Bernoulli(expit(logit(p) + c)), by Newton."""
    z = np.log(np.clip(p, 1e-6, 1 - 1e-6)) - np.log1p(-np.clip(p, 1e-6, 1 - 1e-6))
    c = 0.0
    for _ in range(100):
        mu = 1.0 / (1.0 + np.exp(-(z + c)))
        step = float(np.sum(d - mu) / np.sum(mu * (1 - mu)))
        c += step
        if abs(step) < 1e-12:
            break
    return c


@pytest.mark.slow
def test_peptide_conditional_fit_recovers_abundance_and_citl() -> None:
    """Calibration design (truth -10 + 0.5 a - 0.03 L + 0.8 GRAVY, a ~ N(15, 2.5), 5-20
    peptides), trained on 6,000 genes, scored on an independent 4,000-gene genome. Oracle
    conditional MLE on these seeds: beta_a 0.496, CITL -0.013 (single-replicate SD ~0.02). Teeth:
    the plain logistic gives beta_a ~0.35 and CITL ~-0.89."""
    train = D.simulate_peptides(for_step("test.dl_citl.train", 1), 6000, BETA, n_peps=(5, 20))
    score = D.simulate_peptides(
        for_step("test.dl_citl.score", 0), 4000, BETA, n_peps=(5, 20), tu_size=1
    )
    dl = _dl()
    m = dl.fit_peptide_model(train.genes, train.hits, GippaConfig(seed=SEED), _rng())
    assert m.usable is True
    assert m.coef["abundance_log2"] == pytest.approx(0.5, rel=0.10)
    p = np.array(
        [
            dl.p_detect(m, a, score.peps(t))
            for t, a in zip(score.genes.locus_tag, score.genes.abundance_log2, strict=True)
        ]
    )
    c = _citl(score.genes.detected.to_numpy(float), p)
    assert abs(c) <= 0.1, c


def test_peptide_conditional_fit_recovers_abundance_small() -> None:
    """Fast version, 1,500 training genes: beta_a within 0.08 of 0.5 (oracle SD ~0.013;
    plain fit ~0.35)."""
    train = D.simulate_peptides(for_step("test.dl_citl.train", 2), 1500, BETA, n_peps=(5, 20))
    m = _dl().fit_peptide_model(train.genes, train.hits, GippaConfig(seed=SEED), _rng())
    assert m.usable is True
    assert m.coef["abundance_log2"] == pytest.approx(0.5, abs=0.08)


# ---------------------------------------------------------------------------- separation


def _separated(rule: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Every gene has peptides on both sides of the rule; observed iff the rule holds."""
    s = D.simulate_peptides(np.random.default_rng(41), 300, BETA, n_peps=(8, 20))
    rows = []
    for t in s.genes.locus_tag:
        for pep in s.peptides[t]:
            hit = D.gravy(pep) > 0 if rule == "gravy" else len(pep) <= 15
            if hit:
                rows.append(D.hit_row(pep, t))
    hits = D.hits_frame(rows)
    det = s.genes.locus_tag.isin(set(hits.gene))
    return s.genes.assign(detected=det.to_numpy()), hits


@pytest.mark.parametrize("rule", ["gravy", "length"])
def test_separated_peptide_model_not_usable(rule: str) -> None:
    """Lead decision 5: unpenalised fit; separation -> usable False (never a divergent iterate)."""
    genes, hits = _separated(rule)
    assert genes.detected.sum() >= 50
    m = _dl().fit_peptide_model(genes, hits, GippaConfig(seed=SEED), _rng())
    assert m.usable is False
    assert math.isnan(_dl().p_detect(m, 15.0, list(D.TOY_PEPS)))
