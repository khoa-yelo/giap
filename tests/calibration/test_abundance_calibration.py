"""Calibration of M4-14 gene-level abundance (docs/specs/M4-14-abundance.md, acceptance criterion 3;
GIAP_SPEC.md §7.3.4, §8.4, §11.3).

Statistical claims (for the statistician), each gated through the M4-05 harness
(``_criteria.calibrate``, ADR-0015 two-block confirmation), suite ``regions`` (SUITE_REGISTRY file
entry: the M4-08 AR(1) abundance generator of the region nulls, whose gene abundances these
estimates feed):

1. ``abundance.spearman.ibaq``: on simulated data, Spearman(``abundance_log2``,
   ``abundance_log2_true``) over detected genes with >= 3 quantified peptides is >= 0.7. Per
   replicate the value is the indicator 1{rho < 0.7}; criterion ``fwer`` at 0.05 (share of failing
   replicates <= 0.05 + k sqrt(0.05 * 0.95 / R)). The spec calls this a sanity check (not the M4
   gate); it is logged with the suite like every other claim.
2. ``abundance.spearman.top3``: the same with ``abundance.method = top3``.
3. Teeth: an estimator that ignores the data (``abundance_log2`` shuffled across genes) is a
   confirmed ``fail`` of claim 1 (private log, PR #69 decision 3).
4. ``abundance.ibaq_denominator`` (PR #82 lead decision 4): the iBAQ denominator matters. Over the
   same genes (detected, >= 3 quantified peptides; ``n_observable_peptides`` varies, sd of its log2
   >= 0.3), Spearman(iBAQ, truth) exceeds Spearman(raw sum, truth) beyond noise: per replicate the
   value is 1{the 5% paired-bootstrap lower bound (200 gene resamples) of rho_ibaq - rho_sum <= 0};
   criterion ``fwer`` at 0.05. Both estimates come from ``gene_abundance`` (method ibaq): the raw
   sum is the same call with every ``n_observable_peptides`` = 1, so only the denominator differs.
   Teeth: with ``n_observable_peptides`` permuted across genes the claim is a confirmed ``fail``.
   Re-run on the M4-09 simulator once #80 is merged (backlog).

Simulated data (the M4-09 evidence simulator is not merged; tests/helpers/abundance_oracle.py):
giapsim layout (M4-07, 400 chromosome genes) and expression (M4-08 defaults); each ``detected_true``
gene's observable tryptic peptides (0 missed cleavages, 7-30 aa) are quantified with probability
0.6 and intensity ln I = abundance_log2_true * ln 2 + N(0, 0.5^2) (the M4-09 intensity model),
as unique exact accepted peptides. ``n_observable_peptides`` is the oracle's digest count.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import abundance_oracle as ao

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "regions"
RHO_MIN = 0.7
LEVEL = 0.05
N_GENES = 400


def _simulate(rng: np.random.Generator) -> tuple[Any, ao.SimPeptides]:
    import giapsim.expression as ex
    import giapsim.layout as la

    lay = la.simulate_layout(la.LayoutConfig(n_genes_chromosome=N_GENES), rng)
    expr = ex.simulate_expression(lay, ex.ExpressionConfig(), rng)
    return lay, ao.simulate_peptides(lay, expr, rng)


def _estimate(
    lay: Any, sim: ao.SimPeptides, method: str, n_obs: dict[str, int] | None = None
) -> dict[str, float]:
    from giap.config import GiapConfig
    from giap.genetable.abundance import gene_abundance

    tags = lay.genes["locus_tag"].astype(str).tolist()
    res = gene_abundance(
        lay.genes,
        sim.accepted,
        ao.covariates_frame(sim.n_observable if n_obs is None else n_obs, tags),
        ao.detection(tags, sim.detected),
        None,
        ao.pmap_frame([]),
        GiapConfig(abundance={"method": method}),
    )
    t = res.table
    return dict(zip(t["locus_tag"].astype(str), t["abundance_log2"].astype(float), strict=True))


def _spearman(est: dict[str, float], sim: ao.SimPeptides) -> float:
    vals = ao.peptide_values(sim.accepted)
    keep = [g for g in sim.detected if len(vals.get(g, ())) >= 3]
    assert len(keep) >= 30, f"too few genes with >= 3 quantified peptides: {len(keep)}"
    x = np.array([est[g] for g in keep])
    assert np.isfinite(x).all(), "genes with >= 3 quantified peptides must have a value"
    y = np.array([sim.truth[g] for g in keep])
    return float(stats.spearmanr(x, y).statistic)


def _claim(method: str, shuffle: bool = False) -> Callable[[np.random.Generator], float]:
    def fn(rng: np.random.Generator) -> float:
        lay, sim = _simulate(rng)
        est = _estimate(lay, sim, method)
        if shuffle:
            keys = [k for k, v in est.items() if np.isfinite(v)]
            vals = np.array([est[k] for k in keys])
            est = dict(zip(keys, rng.permutation(vals), strict=True))
        return 1.0 if _spearman(est, sim) < RHO_MIN else 0.0

    return fn


@pytest.mark.parametrize("method", ["ibaq", "top3"])
@pytest.mark.parametrize("full", MODES)
def test_spearman_with_true_abundance(full: bool, method: str) -> None:
    R = C.replicates(SUITE, full)
    rec = C.calibrate(SUITE, "fwer", LEVEL, _claim(method), R,
                      stream=f"abundance.spearman.{method}")  # fmt: skip
    assert rec["status"] != "fail", rec


@pytest.mark.calibration
def test_teeth_shuffled_abundance_fails(teeth_log: Path) -> None:
    rec = C.calibrate(SUITE, "fwer", LEVEL, _claim("ibaq", shuffle=True), 20,
                      stream="abundance.spearman.teeth_shuffled", log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", rec


@pytest.mark.calibration
def test_simulated_data_sanity() -> None:
    """The helper simulator itself: enough detected genes with >= 3 quantified peptides, and
    intensities centred on 2^abundance (mean of ln I - a ln 2 near 0)."""
    _, sim = _simulate(C.replicate_rng("abundance.sim_sanity", 0))
    vals = ao.peptide_values(sim.accepted)
    assert sum(len(v) >= 3 for v in vals.values()) >= 30
    resid = [np.log(x) - sim.truth[g] * np.log(2.0) for g, v in vals.items() for x in v]
    assert abs(float(np.mean(resid))) < 0.05
    assert 0.45 < float(np.std(resid)) < 0.55


# ---------------------------------------------------------------- claim 4: the iBAQ denominator

N_BOOT = 200
LOG2_NOBS_SD_MIN = 0.3


def _rank_corr_rows(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Spearman per row of two (B, n) arrays (average ranks, Pearson of ranks)."""
    rx = stats.rankdata(x, axis=1)
    ry = stats.rankdata(y, axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    return np.asarray((rx * ry).sum(axis=1) / np.sqrt((rx**2).sum(axis=1) * (ry**2).sum(axis=1)))


def _denominator_claim(permute_n_obs: bool = False) -> Callable[[np.random.Generator], float]:
    def fn(rng: np.random.Generator) -> float:
        lay, sim = _simulate(rng)
        tags = list(sim.n_observable)
        n_obs = dict(sim.n_observable)
        if permute_n_obs:
            n_obs = dict(zip(tags, rng.permutation([n_obs[t] for t in tags]).tolist(), strict=True))
        ibaq = _estimate(lay, sim, "ibaq", n_obs)
        raw_sum = _estimate(lay, sim, "ibaq", dict.fromkeys(tags, 1))
        vals = ao.peptide_values(sim.accepted)
        keep = [g for g in sim.detected if len(vals.get(g, ())) >= 3]
        assert len(keep) >= 30, f"too few genes with >= 3 quantified peptides: {len(keep)}"
        log2_nobs = np.log2([max(1, sim.n_observable[g]) for g in keep])
        assert float(np.std(log2_nobs)) >= LOG2_NOBS_SD_MIN, "n_observable must vary"
        x = np.array([ibaq[g] for g in keep])
        s = np.array([raw_sum[g] for g in keep])
        y = np.array([sim.truth[g] for g in keep])
        assert np.isfinite(x).all() and np.isfinite(s).all()
        idx = rng.integers(0, len(keep), size=(N_BOOT, len(keep)))
        delta = _rank_corr_rows(x[idx], y[idx]) - _rank_corr_rows(s[idx], y[idx])
        lower = float(np.quantile(delta, 0.05))
        return 1.0 if lower <= 0.0 else 0.0

    return fn


@pytest.mark.parametrize("full", MODES)
def test_ibaq_denominator_beats_raw_sum(full: bool) -> None:
    R = C.replicates(SUITE, full)
    rec = C.calibrate(SUITE, "fwer", LEVEL, _denominator_claim(), R,
                      stream="abundance.ibaq_denominator")  # fmt: skip
    assert rec["status"] != "fail", rec


@pytest.mark.calibration
def test_teeth_permuted_denominator_fails(teeth_log: Path) -> None:
    claim = _denominator_claim(permute_n_obs=True)
    rec = C.calibrate(SUITE, "fwer", LEVEL, claim, 20,
                      stream="abundance.ibaq_denominator.teeth", log_dir=teeth_log)  # fmt: skip
    assert rec["status"] == "fail", rec
