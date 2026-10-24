"""Attack (GIAP_SPEC.md §8.5, statistician): a TU-structured null for the M8-08 determinants
models (§7.8.2 steps 2-3, §8.2 "resample TUs, not genes"; docs/specs/M8-08-determinants-models.md).

Generator (replicate ``i`` of stream ``determinants.attack_tu.<name>``): 1,500 genes (600 for the
CV arm) in TUs of geometric size (p 0.45, capped at 8). Features that are **constant within a TU**
(``promotech_tu`` uniform, ``tx_rate_sum_tu`` normal: the transcription block) and gene-level noise
(``cai``, ``gravy``, ``pi``); the response is a TU random effect plus gene noise (SD 1 / 0.7);
detection is a TU-shared latent above its 30 % quantile. Everything is independent of the
features, but every feature and the response share the TU structure, which is what fools a
gene-level resampler: a gene bootstrap or gene-split folds treat a TU's genes as independent.

Gated (a failing attack blocks the merge):

- ``test_univariate_ci_level``: the 95 % TU-block bootstrap CIs of Spearman rho (true 0) and
  AUROC (true 0.5) of ``promotech_tu`` miss the truth in at most ``level_bound(0.05, R)`` of
  replicates, **unless** every univariate row carries ``ci_unestimated`` (then the rates are
  recorded and a warning is raised: an unflagged interval must be at level, R4). The gene
  bootstrap on the same data is recorded (it should fail: the attack bites). Statistician probe
  (Slurm 1866223, R = 1,000, B = 1,000): TU-block miss rates rho 0.051, AUROC 0.072; binary
  TU-constant feature rho 0.065, AUROC 0.071, against a bound of 0.064 -- the bootstrap SE is
  5-9 % low for AUROC, and the percentile, normal and basic intervals all miss (Slurm 1866311).
- ``test_cv_r2_no_tu_leak``: with ``tu_folds`` and ``tu=``, the mean CV R² over replicates is
  <= 0 for both models and the mean ΔR²(transcription) is <= 0.02 (one-sided: a fold that splits a
  TU lets LightGBM read the TU's response off a TU-constant feature, giving R² > 0). The same
  data on gene-split folds is recorded.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d

from giap.stats import bootstrap

SUITE = "determinants"
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite(SUITE)]
MODES = [pytest.param(False, marks=pytest.mark.calibration_nightly, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip
B = 1000
BLOCKS = {"promotech_tu": "transcription", "tx_rate_sum_tu": "transcription", "cai": "codon",
          "gravy": "protein", "pi": "protein"}  # fmt: skip
CV_LEAK_TOL = 0.02


def _sizes(rng: np.random.Generator, n: int) -> list[int]:
    s: list[int] = []
    tot = 0
    while tot < n:
        s.append(int(min(8, rng.geometric(0.45))))
        tot += s[-1]
    s[-1] -= tot - n
    return [x for x in s if x]


def _world(
    rng: np.random.Generator, n: int
) -> tuple[list[str], pd.DataFrame, np.ndarray, np.ndarray]:
    sizes = _sizes(rng, n)
    k, codes = len(sizes), np.repeat(np.arange(len(sizes)), sizes)
    X = pd.DataFrame(
        {"promotech_tu": rng.random(k)[codes], "tx_rate_sum_tu": rng.normal(size=k)[codes],
         "cai": rng.normal(size=n), "gravy": rng.normal(size=n), "pi": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)]),
    )  # fmt: skip
    y = rng.normal(size=k)[codes] + 0.7 * rng.normal(size=n)
    lat = rng.normal(size=k)[codes] + 0.5 * rng.normal(size=n)
    return d.tu_ids(sizes), X, y, lat > np.quantile(lat, 0.3)


def _miss(lo: float, hi: float, truth: float) -> bool:
    return not (lo <= truth <= hi)


def _log(record: dict[str, object]) -> None:
    path = C.log_dir() / f"{SUITE}_attack_tu.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")


@pytest.mark.parametrize("full", MODES)
def test_univariate_ci_level(full: bool) -> None:
    m = d.mod()
    R = 5 * C.replicates(SUITE, full=full)  # univariate replicates are cheap (~1 s)
    miss = {"rho": [], "auroc": [], "rho_gene": [], "auroc_gene": []}
    flagged = True
    for i in range(R):
        rng = C.replicate_rng(f"{SUITE}.attack_tu.uni", i)
        tu, X, y, det = _world(rng, 1500)
        f = d.features({"promotech_tu": X["promotech_tu"].tolist()}, tu=tu, detected=det.tolist(),
                       abundance=np.where(det, y, np.nan).tolist())  # fmt: skip
        uni = m.univariate(f, B, rng)
        rho = d.one(uni, metric="spearman_rho")
        auc = d.one(uni, metric="auroc")
        flagged &= all("ci_unestimated" in d.flags_of(v) for v in uni["flags"])
        miss["rho"].append(_miss(rho.ci_low, rho.ci_high, 0.0))
        miss["auroc"].append(_miss(auc.ci_low, auc.ci_high, 0.5))
        x = X["promotech_tu"].to_numpy()
        gr = pd.DataFrame({"tu": np.arange(det.sum()).astype(str), "x": x[det], "y": y[det]})
        ga = pd.DataFrame({"tu": np.arange(det.size).astype(str), "x": x, "y": det})
        cr = bootstrap.tu_block_bootstrap(gr, "tu", m._rho, B, rng, columns=("x", "y"))
        ca = bootstrap.tu_block_bootstrap(ga, "tu", m._auroc, B, rng, columns=("x", "y"))
        miss["rho_gene"].append(_miss(cr.ci_low, cr.ci_high, 0.0))
        miss["auroc_gene"].append(_miss(ca.ci_low, ca.ci_high, 0.5))
    rates = {k: float(np.mean(v)) for k, v in miss.items()}
    bound = C.level_bound(0.05, R)
    _log({"claim": "univariate_ci_level", "R": R, "rates": rates, "bound": bound,
          "flagged": flagged})  # fmt: skip
    print(f"TU-null univariate CI miss rates, R={R}: {rates} (bound {bound:.3f})")
    over = {k: v for k, v in rates.items() if k in ("rho", "auroc") and v > bound}
    if flagged:
        if over:
            warnings.warn(f"univariate CIs (flagged ci_unestimated) miss at {over}", stacklevel=1)
        return
    assert not over, f"unflagged univariate CIs miss the truth at {over} (bound {bound:.3f})"


@pytest.mark.parametrize("full", MODES)
def test_cv_r2_no_tu_leak(full: bool) -> None:
    m = d.mod()
    R = C.replicates(SUITE, full=full)
    cfg = d.config(fast=False)
    rec: dict[str, dict[str, list[float]]] = {}
    for i in range(R):
        rng = C.replicate_rng(f"{SUITE}.attack_tu.cv", i)
        tu, X, y, _ = _world(rng, 600)
        gene = [f"g{j}" for j in range(len(tu))]
        for arm, groups in (("tu", tu), ("gene", gene)):
            folds = m.tu_folds(groups, cfg.outer_folds, cfg.repeats, rng)
            seed = C.load_tolerances().calibration_base + i
            cv = m.cross_validate(X, y, folds, BLOCKS, cfg, seed, 1, tu=groups)
            for name, res in cv.items():
                r = rec.setdefault(f"{arm}.{name}", {"r2": [], "delta_tx": []})
                r["r2"].append(float(np.mean(res.r2)))
                r["delta_tx"].append(float(np.mean(res.delta_r2["transcription"])))
    means = {k: {s: float(np.mean(v)) for s, v in r.items()} for k, r in rec.items()}
    _log({"claim": "cv_r2_no_tu_leak", "R": R, "means": means})
    print(f"TU-null CV, R={R}: {means}")
    for name in d.MODELS:
        assert means[f"tu.{name}"]["r2"] <= 0.0, (name, means)
        assert means[f"tu.{name}"]["delta_tx"] <= CV_LEAK_TOL, (name, means)
