"""Calibration of the M4-13 detectability model (docs/specs/M4-13-detectability.md; GIAP_SPEC.md
§7.3.3, §8.4 row "Detectability | logistic truth on physicochemical covariates | calibration slope
1 ± 0.2 | 20 | 100"; M4-index lead decision 6: the gate is the MEAN slope).

Suite ``detectability`` (SUITE_REGISTRY file entry); R = ``replicates_{fast,full}.detectability``
(20 / 100). Replicate i of stream s draws from ``_criteria.replicate_rng(s, i)``.

Claims (each named so the statistician can review it):

D1  Cross-fitted slope (the §8.4 gate). Data: M4-07 layout (2,000 chromosome genes fast, 4,500
    full, default regions) + M4-08 ``detection = "covariates_only"`` (logistic truth in log length,
    log1p n_obs, GRAVY); eligible = non-partial cds/sORF; covariates from M4-10. Fit
    ``fit_detectability`` (rng ``giap.rng.for_step("detectability", seed_i)``) and take the
    ``calibration_slope`` from its diagnostics (checked per replicate against the oracle's
    unpenalised recalibration of the returned ``p_detect``). Criterion: mean over the R replicates
    ∈ ``tolerances.detectability_slope`` = [0.8, 1.2]. Per-replicate slopes, intercepts and Brier
    scores are appended to ``<log_dir>/detectability.jsonl`` (criterion
    ``detectability_slope_mean``; no ``values`` key, so ``merge_shards`` ignores it). Not routed
    through ``calibrate``: a mean-in-interval criterion is none of level/uniform/fdr/fwer. Under
    ``--shard`` it runs in the shard holding replicate 0 only (whole R there); under ``--merge``
    confirmation it is skipped.
D2  Equal-count bins are honest under perfect calibration. p_i = expit(N(0, 1.5²)), d_i ~
    Bernoulli(p_i), n = 20,000 (example 6 design). From the 10 ``calibration_diagnostics`` bins,
    H = Σ_b n_b² (observed_rate_b - mean_p_b)² / Σ_{i∈b} p_i(1 - p_i) (variances from the oracle's
    binning of the same p; bins depend on p only) ~ χ²₁₀; replicate value = χ²₁₀ upper tail of H.
    ``calibrate``: KS uniformity and the level bound at every alpha in ``levels``.
D3  ``calibration_slope`` is the recalibration MLE: same data, z = (slope - 1)/SE₁ with SE from
    the Fisher information at (intercept 0, slope 1); two-sided normal p. Uniform + level.
D4  Same for ``calibration_intercept``: z = intercept/SE₀. Uniform + level.
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.special import expit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detectability_data as dd
import detectability_oracle as O

from giap.config.model import DetectabilityConfig
from giap.rng import for_step

SUITE = "detectability"
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
N_PERFECT = 20_000
_MEMO: dict[tuple[str, int], Any] = {}


def _mod() -> Any:
    from giap.genetable import detectability

    return detectability


def _overall(diag: pd.DataFrame, metric: str) -> float:
    sel = diag[(diag["group_type"] == "overall") & (diag["metric"] == metric)]
    assert len(sel) == 1, metric
    return float(sel["value"].iloc[0])


def _bins(diag: pd.DataFrame) -> pd.DataFrame:
    b = diag[diag["group_type"] == "bin"]
    wide = b.pivot(index="group", columns="metric", values="value").astype(float)
    wide.index = wide.index.astype(int)
    return wide.sort_index()


def _memo(
    stream: str, fn: Callable[[np.random.Generator], Any]
) -> Callable[[np.random.Generator], Any]:
    def wrapped(rng: np.random.Generator) -> Any:
        key = (stream, int(rng.bit_generator.seed_seq.entropy))  # type: ignore[attr-defined]
        if key not in _MEMO:
            _MEMO[key] = fn(rng)
        return _MEMO[key]

    return wrapped


def _gate(stream: str, full: bool, fn: Callable[[np.random.Generator], float]) -> None:
    R = C.replicates(SUITE, full)
    statuses = [C.calibrate(SUITE, "uniform", 0.0, fn, R, stream=stream)["status"]]
    for a in C.load_tolerances().levels:
        statuses.append(C.calibrate(SUITE, "level", a, fn, R, stream=stream)["status"])
    assert "fail" not in statuses, (stream, statuses)


# ----------------------------------------------------------------------------- D1


def _cross_fit_replicate(rng: np.random.Generator, n_genes: int) -> tuple[float, float, float]:
    gs = dd.sim_geneset(rng, n_genes)
    seed = int(rng.integers(0, 2**31 - 1))
    res = _mod().fit_detectability(
        gs.detected, gs.cov, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", seed)
    )
    assert res.constant is False
    elig = gs.extra["eligible"]
    d = gs.detected[elig].astype(bool).to_numpy()
    p = res.p_detect[elig].to_numpy(float)
    assert np.isfinite(p).all() and res.p_detect[~elig].isna().all()
    slope = _overall(res.diagnostics, "calibration_slope")
    icpt = _overall(res.diagnostics, "calibration_intercept")
    brier = _overall(res.diagnostics, "brier")
    a, b = O.logistic_recalibration(d, p)
    assert slope == pytest.approx(b, abs=1e-4) and icpt == pytest.approx(a, abs=1e-4)
    assert brier == pytest.approx(float(np.mean((p - d) ** 2)), abs=1e-12)
    return slope, icpt, brier


@pytest.mark.parametrize("full", MODES)
def test_d1_mean_cross_fitted_slope_in_tolerance(full: bool, tolerances: dict[str, object]) -> None:
    if os.environ.get(C.CONFIRM_ENV) is not None:
        pytest.skip("D1 is not a calibrate() check; nothing to confirm")
    R = C.replicates(SUITE, full)
    if C.replicate_range(R).start != 0:
        pytest.skip("D1 runs whole in the shard holding replicate 0")
    lo, hi = (float(x) for x in tolerances["detectability_slope"])  # type: ignore[attr-defined]
    n_genes = 4500 if full else 2000
    stream = "detectability.cross_fit_slope"
    t0 = time.perf_counter()
    out = C.run_replicates(stream, lambda rng: _cross_fit_replicate(rng, n_genes), R)
    slopes, icpts, briers = (np.array(v, float) for v in zip(*out, strict=True))
    mean = float(slopes.mean())
    ok = lo <= mean <= hi
    record = {
        "suite": SUITE,
        "stream": stream,
        "criterion": "detectability_slope_mean",
        "level": 0.0,
        "observed": mean,
        "bound": hi if mean > 1 else lo,
        "R": R,
        "seed": C.suite_seed(stream, 0),
        "duration": time.perf_counter() - t0,
        "ok": ok,
        "block": 1,
        "status": "pass" if ok else "fail",
        "n_genes": n_genes,
        "slopes": slopes.tolist(),
        "intercepts": icpts.tolist(),
        "briers": briers.tolist(),
    }
    directory = C.log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / f"{SUITE}.jsonl").open("a") as fh:
        fh.write(json.dumps(record) + "\n")
    assert ok, f"mean slope {mean:.4f} outside [{lo}, {hi}]; slopes {np.round(slopes, 3)}"


# ----------------------------------------------------------------------------- D2-D4


def _perfect(rng: np.random.Generator) -> dict[str, float]:
    p = expit(rng.normal(0.0, 1.5, N_PERFECT))
    d = (rng.random(N_PERFECT) < p).astype(int)
    diag = _mod().calibration_diagnostics(d, p)
    bins = _bins(diag)
    assert list(bins.index) == list(range(1, 11))
    var = np.array([float(np.sum(p[b] * (1 - p[b]))) for b in O.bins(p, 10)])
    n_b = bins["n"].to_numpy()
    h = float(np.sum((n_b * (bins["observed_rate"] - bins["mean_p"]).to_numpy()) ** 2 / var))
    se0, se1 = O.recalibration_se_at_truth(p)
    slope = _overall(diag, "calibration_slope")
    icpt = _overall(diag, "calibration_intercept")
    return {
        "hl": float(stats.chi2.sf(h, 10)),
        "slope": float(2 * stats.norm.sf(abs(slope - 1.0) / se1)),
        "intercept": float(2 * stats.norm.sf(abs(icpt) / se0)),
    }


def _check(stream: str, full: bool, key: str) -> None:
    fn = _memo(stream, _perfect)
    _gate(stream, full, lambda rng: float(fn(rng)[key]))


@pytest.mark.parametrize("full", MODES)
def test_d2_bins_chi2_under_perfect_calibration(full: bool) -> None:
    _check("detectability.perfect_bins", full, "hl")


@pytest.mark.parametrize("full", MODES)
def test_d3_slope_wald_under_perfect_calibration(full: bool) -> None:
    _check("detectability.perfect_slope", full, "slope")


@pytest.mark.parametrize("full", MODES)
def test_d4_intercept_wald_under_perfect_calibration(full: bool) -> None:
    _check("detectability.perfect_intercept", full, "intercept")
