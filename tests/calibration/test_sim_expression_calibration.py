"""Calibration of the giapsim expression model (issue M4-08; GIAP_SPEC.md §11.3, §8.4, §7.3.6).

Every distributional claim of docs/specs/M4-08-sim-expression.md becomes a per-replicate p-value
that is exactly U(0, 1) when the simulator follows the spec, gated through the M4-05 harness
(``_criteria.calibrate``: KS uniformity and the level bound at every ``levels`` alpha; ADR-0015
two-block confirmation). Replicate ``i`` draws its layout and expression from
``replicate_rng(stream, i)``.

Suite ``regions`` (generator of the region-test nulls; SUITE_REGISTRY file entry):
* AR(1) background (rho in {0, 0.5, 0.9}): the whitened TU series z_0 = x_0 sqrt(1 - rho^2)/sigma,
  z_t = (x_t - rho x_{t-1})/sigma iid N(0, 1), sigma = ``tu_sd``, mean ``mu_log2``: chi2 of sum z^2
  and the lag-1 autocorrelation of z;
* region effect (prophage -3), plasmid copy (log2 10) and island (+1.5 SD over 10 TUs) shifts of
  TU means at rho = 0 (z-tests with known sigma);
* gene level: abundance - tu_log2 regressed on TU position has slope -``position_decay_log2``
  and residual variance ``stoich_sd^2 + gene_sd^2``.
Suite ``detectability`` (function entries): ``detected_true ~ Bernoulli(p_detect_true)`` (total and
by p-decile), conditional independence of detection and abundance under ``covariates_only``
(score test), and the BGC detection shift of example 5 (power input of §7.3.6).
"""

from __future__ import annotations

import math
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import expression_oracle as eo

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
MU, SIGMA = 20.0, 2.0
_MEMO: dict[tuple[str, int], float] = {}


def _mods() -> tuple[Any, Any]:
    import giapsim.expression as ex
    import giapsim.layout as la

    return ex, la


def _layout(rng: np.random.Generator, **kw: Any) -> Any:
    _, la = _mods()
    return la.simulate_layout(la.LayoutConfig(**kw), rng)


def _express(lay: Any, rng: np.random.Generator, **kw: Any) -> Any:
    ex, _ = _mods()
    return ex.simulate_expression(lay, ex.ExpressionConfig(**kw), rng)


def _two_sided_z(z: float) -> float:
    return float(2 * stats.norm.sf(abs(z)))


def _two_sided_chi2(x: float, df: int) -> float:
    return float(min(1.0, 2 * min(stats.chi2.cdf(x, df), stats.chi2.sf(x, df))))


def _gate(suite: str, stream: str, full: bool, fn: Callable[[np.random.Generator], float]) -> None:
    """KS uniformity and the level bound at every alpha, each through ``calibrate``."""

    def memo(rng: np.random.Generator) -> float:
        key = (stream, int(rng.bit_generator.seed_seq.entropy))  # type: ignore[attr-defined]
        if key not in _MEMO:
            _MEMO[key] = float(fn(rng))
        return _MEMO[key]

    R = C.replicates(suite, full)
    statuses = [C.calibrate(suite, "uniform", 0.0, memo, R, stream=stream)["status"]]
    for a in C.load_tolerances().levels:
        statuses.append(C.calibrate(suite, "level", a, memo, R, stream=stream)["status"])
    assert "fail" not in statuses, (stream, statuses)


# ----------------------------------------------------------------------------- suite: regions


def _whitened(rng: np.random.Generator, rho: float) -> np.ndarray:
    lay = _layout(rng, n_genes_chromosome=400, plasmids=(30,))
    series = eo.tu_series(lay, _express(lay, rng, ar1_rho=rho).genes)
    return np.concatenate([eo.ar1_whiten(x - MU, rho, SIGMA) for x in series.values()])


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
@pytest.mark.parametrize("full", MODES)
def test_ar1_whitened_sum_of_squares(full: bool, rho: float) -> None:
    def fn(rng: np.random.Generator) -> float:
        z = _whitened(rng, rho)
        return _two_sided_chi2(float(np.sum(z * z)), len(z))

    _gate("regions", f"sim_expression.ar1_chi2.rho{rho}", full, fn)


@pytest.mark.parametrize("rho", [0.5, 0.9])
@pytest.mark.parametrize("full", MODES)
def test_ar1_whitened_mean_and_lag1(full: bool, rho: float) -> None:
    """Whitened series mean 0 (mu = ``mu_log2``); no residual lag-1 correlation (rho honoured)."""

    def mean_fn(rng: np.random.Generator) -> float:
        z = _whitened(rng, rho)
        return _two_sided_z(float(z.sum()) / math.sqrt(len(z)))

    def lag_fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=600)
        x = eo.tu_series(lay, _express(lay, rng, ar1_rho=rho).genes)["chr1"] - MU
        z = eo.ar1_whiten(x, rho, SIGMA)
        return _two_sided_z(float(np.sum(z[1:] * z[:-1])) / math.sqrt(len(z) - 1))

    _gate("regions", f"sim_expression.ar1_mean.rho{rho}", full, mean_fn)
    _gate("regions", f"sim_expression.ar1_lag1.rho{rho}", full, lag_fn)


def _tu_mean_z(lay: Any, genes: Any, tus: list[str], shift: float) -> float:
    tl = dict(zip(genes["locus_tag"], genes["tu_log2"], strict=True))
    mem = eo.tu_members(lay)
    x = np.array([float(tl[mem[t][0]]) for t in tus]) - MU - shift
    return float(x.sum() / (SIGMA * math.sqrt(len(x))))


@pytest.mark.parametrize("full", MODES)
def test_region_effect_shift(full: bool) -> None:
    """rho = 0: TUs entirely inside the prophage have mean mu - 3 (example 3)."""

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng)
        g = _express(lay, rng, ar1_rho=0.0, region_effect_log2={"prophage": -3.0}).genes
        status = eo.tu_region_status(lay, "prophage")
        tus = [t for t, s in status.items() if s == "all"]
        return _two_sided_z(_tu_mean_z(lay, g, tus, -3.0))

    _gate("regions", "sim_expression.region_effect", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_plasmid_copy_shift(full: bool) -> None:
    """rho = 0: plasmid TUs have mean mu + log2(10) (example 4)."""

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=300, plasmids=(60,))
        g = _express(lay, rng, ar1_rho=0.0, plasmid_copy_log2={"p1": math.log2(10)}).genes
        return _two_sided_z(_tu_mean_z(lay, g, eo.tu_order(lay)["p1"], math.log2(10)))

    _gate("regions", "sim_expression.plasmid_copy", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_island_shift(full: bool) -> None:
    """rho = 0: the 10 island TUs have mean mu + 1.5 sigma; the rest of chr1 mean mu (no leak)."""
    ex, _ = _mods()

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=300)
        first = int(rng.integers(0, 50))
        spec = ex.IslandSpec("chr1", first, 10, 1.5)
        g = _express(lay, rng, ar1_rho=0.0, islands=(spec,)).genes
        order = eo.tu_order(lay)["chr1"]
        z_in = _tu_mean_z(lay, g, order[first : first + 10], 1.5 * SIGMA)
        z_out = _tu_mean_z(lay, g, order[:first] + order[first + 10 :], 0.0)
        return float(stats.chi2.sf(z_in**2 + z_out**2, 2))

    _gate("regions", "sim_expression.island", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_gene_level_decay_and_noise(full: bool) -> None:
    """abundance - tu_log2 = c - 0.2*position + N(0, 0.7^2 + 0.5^2): slope and residual variance."""
    s2 = 0.7**2 + 0.5**2

    def parts(rng: np.random.Generator) -> tuple[float, float]:
        lay = _layout(rng, n_genes_chromosome=600)
        g = _express(lay, rng).genes
        pos = eo.gene_tu_position(lay)
        k = np.array([pos[t][1] for t in g["locus_tag"]], dtype=float)
        d = (g["abundance_log2_true"] - g["tu_log2"]).to_numpy(float)
        X = np.column_stack([np.ones_like(k), k])
        b, *_ = np.linalg.lstsq(X, d, rcond=None)
        rss = float(np.sum((d - X @ b) ** 2))
        se = math.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
        return _two_sided_z((b[1] + 0.2) / se), _two_sided_chi2(rss / s2, len(d) - 2)

    _gate("regions", "sim_expression.position_slope", full, lambda r: parts(r)[0])
    _gate("regions", "sim_expression.gene_noise", full, lambda r: parts(r)[1])


# ----------------------------------------------------------------------------- suite: detectability


@pytest.mark.parametrize("full", MODES)
def test_detection_bernoulli_total(full: bool) -> None:
    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=2000)
        g = _express(lay, rng).genes
        p = g["p_detect_true"].to_numpy(float)
        d = g["detected_true"].astype(float).to_numpy()
        return _two_sided_z(float((d - p).sum()) / math.sqrt(float((p * (1 - p)).sum())))

    _gate("detectability", "sim_expression.bernoulli_total", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_detection_bernoulli_by_decile(full: bool) -> None:
    """sum  over p-deciles of (O - E)^2/V ~ chi2_10 (reliability of the truth)."""

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=2000)
        g = _express(lay, rng).genes
        p = g["p_detect_true"].to_numpy(float)
        d = g["detected_true"].astype(float).to_numpy()
        order = np.argsort(p, kind="stable")
        stat = 0.0
        for idx in np.array_split(order, 10):
            stat += float((d[idx] - p[idx]).sum()) ** 2 / float((p[idx] * (1 - p[idx])).sum())
        return float(stats.chi2.sf(stat, 10))

    _gate("detectability", "sim_expression.bernoulli_decile", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_covariates_only_detection_independent_of_abundance(full: bool) -> None:
    """Example 6: score test of adding abundance to the true logit (exact given p and abundance)."""

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=2000)
        g = _express(lay, rng, detection="covariates_only").genes
        p = g["p_detect_true"].to_numpy(float)
        d = g["detected_true"].astype(float).to_numpy()
        a = g["abundance_log2_true"].to_numpy(float)
        a = a - a.mean()
        return _two_sided_z(float(((d - p) * a).sum()) / math.sqrt(float((p * (1 - p) * a * a)
                                                                        .sum())))  # fmt: skip

    _gate("detectability", "sim_expression.covariates_only_score", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_bgc_detection_shift(full: bool) -> None:
    """Example 5 (§7.3.6 power input): BGC genes detected at sigmoid(covariate logit + shift)."""
    ex, la = _mods()
    shift = eo.logit(0.7) - eo.logit(0.5)

    def fn(rng: np.random.Generator) -> float:
        lay = _layout(rng, n_genes_chromosome=800, regions=(la.RegionSpec("bgc", 300),))
        g = eo.merged(lay, _express(lay, rng, region_detect_shift={"bgc": shift}).genes)
        rt = eo.region_types_of_genes(lay)
        g = g[~g["aa"].map(eo.is_na) & g["locus_tag"].map(lambda t: "bgc" in rt.get(t, set()))]
        p = eo.sigmoid(eo.expected_logit(g, ex.DEFAULT_BETA, "covariates_only", {}, {}) + shift)
        d = g["detected_true"].astype(float).to_numpy()
        return _two_sided_z(float((d - p).sum()) / math.sqrt(float((p * (1 - p)).sum())))

    _gate("detectability", "sim_expression.bgc_shift", full, fn)


@pytest.mark.parametrize("full", MODES)
def test_default_detected_share_every_replicate(full: bool) -> None:
    """Example 1 holds across layouts, not only for seed 1: share within 0.55 +/- 0.15."""
    R = C.replicates("detectability", full)
    shares = C.run_replicates(
        "sim_expression.share",
        lambda rng: float(_express(_layout(rng), rng).genes["detected_true"].astype(bool).mean()),
        R,
    )
    assert all(abs(s - 0.55) <= 0.15 for s in shares), shares
