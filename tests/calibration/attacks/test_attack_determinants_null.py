"""Attack (GIPPA_SPEC.md §8.5, statistician): structure the M8-09 TU-size-preserving null does not
keep (§7.8.2 step 4, §8.1 row "Determinant model"; docs/specs/M8-09-determinants-attribution-null.md
rules 5-6). In every arm the response is independent of every feature in the population; the
method passes when the null p stays at level (§8.4).

Arms (replicate ``i`` of stream ``determinants.attack_null.<arm>.<model>``, 400 genes in TUs of
1-5, six N(0, 1) gene-level columns, TU effect + gene noise as in the M8-09 calibration):

- ``autocorr``: §8.5 attack 1 for this method. Abundance has an inter-TU AR(1) component
  (rho 0.9 along the TU order, about 37 % of the variance; no islands, no trend) and
  ``gc_local_dev`` is an independent gene-level AR(1) (rho 0.9, the 10-kb window). Two independent
  autocorrelated series correlate spuriously in one genome; the TU permutation destroys the
  response's autocorrelation, so its null is too narrow.
- ``autocorr_pos``: the same plus ``replichore_position`` (a §7.8.1 context feature). With TU
  folds drawn at random along the genome LightGBM interpolates the response from its
  neighbours, so the CV R² itself is positive under the null.
- ``detect``: §8.5 attack 2 for this method (detectability confounding). 700 genes; detection
  logit = 0.4 + 1.5 z(y) + 1.5 z(``length_aa``); model rows = detected genes (§7.8.2 step 1).
  Selection on a collider makes ``length_aa`` predict abundance among the detected genes.

Statistician probe (Slurm 1867304; R = 100, B = 99, 50 trees; rates at alpha 0.01 / 0.05 / 0.10,
bounds 0.030 / 0.094 / 0.160): autocorr EN 0.05 / 0.12 / 0.13, LightGBM 0.02 / 0.12 / 0.23;
autocorr_pos EN 0.33 / 0.41 / 0.43, LightGBM 0.80 / 0.95 / 0.96; detect EN 0.64 / 0.83 / 0.86,
LightGBM 0.21 / 0.53 / 0.69. The clustered and size worlds at R = 200 stay at level.

Gate (ADR-0032 pattern): every arm's rejection rate <= ``level_bound`` at every alpha **unless**
every NullResult of the arm carries ``null_uncalibrated`` (then the rates are logged and a
warning raised: an unflagged permutation p must be at level, R4).
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
import determinants_attrib_data as a
import determinants_models_data as d

from gippa.config.model import DeterminantsConfig, LightgbmConfig

SUITE = "determinants"
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite(SUITE)]
MODES = [pytest.param(False, marks=pytest.mark.calibration_nightly, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip
MODELS = ("lightgbm", "elastic_net")
ARMS = ("autocorr", "autocorr_pos", "detect")
UNCALIBRATED = "null_uncalibrated"


def _ar1(rng: np.random.Generator, n: int, rho: float) -> np.ndarray:
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0]
    for i in range(1, n):
        x[i] = rho * x[i - 1] + np.sqrt(1 - rho**2) * e[i]
    return x


def _world(rng: np.random.Generator, arm: str) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    n = 700 if arm == "detect" else 400
    sizes = d.mixed_sizes(rng, n, 1, 5)
    tu = d.tu_ids(sizes)
    k = len(sizes)
    u = np.repeat(rng.normal(size=k), sizes)
    X, _, _, _ = d.planted(rng, n=n, signal=0.0)
    keep = np.ones(n, dtype=bool)
    if arm in ("autocorr", "autocorr_pos"):
        y = 0.6 * np.repeat(_ar1(rng, k, 0.9), sizes) + 0.5 * u + rng.normal(0.0, 0.6, n)
        X["gc_local_dev"] = _ar1(rng, n, 0.9)
        if arm == "autocorr_pos":
            X["replichore_position"] = np.arange(n) / n
    elif arm == "detect":
        X["length_aa"] = rng.normal(5.5, 0.6, n)
        y = 0.7 * u + rng.normal(0.0, 0.7, n)
        ln = X["length_aa"].to_numpy()
        z = 0.4 + 1.5 * (y - y.mean()) / y.std() + 1.5 * (ln - ln.mean()) / ln.std()
        keep = rng.random(n) < 1.0 / (1.0 + np.exp(-z))
    else:
        raise ValueError(arm)
    return (
        X.loc[keep],
        np.asarray(y, dtype=float)[keep],
        [t for t, kp in zip(tu, keep, strict=True) if kp],
    )


def _cfg(full: bool) -> DeterminantsConfig:
    return (
        DeterminantsConfig()
        if full
        else DeterminantsConfig(lightgbm=LightgbmConfig(n_estimators=50))
    )


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("arm", ARMS)
@pytest.mark.parametrize("full", MODES)
def test_null_level_under_attack(full: bool, arm: str, model: str) -> None:
    m, nl = d.mod(), a.null_mod()
    R = C.replicates(SUITE, full=full)
    B = 199 if full else 99
    cfg = _cfg(full)
    ps: list[float] = []
    flagged = True
    for i in range(R):
        rng = C.replicate_rng(f"{SUITE}.attack_null.{arm}.{model}", i)
        X, y, tu = _world(rng, arm)
        folds = np.asarray(m.tu_folds(tu, cfg.outer_folds, cfg.repeats, rng))
        seed = C.load_tolerances().calibration_base + i
        cv = m.cross_validate(X, y, folds, {}, cfg, seed, 1, tu=tu)
        obs = {model: float(np.mean(cv[model].r2))}
        res = nl.null_p(X, y, tu, folds[0], obs, cfg, B, seed, 1)[model]
        flagged &= UNCALIBRATED in d.flags_of(res.flags)
        ps.append(float(res.p))
    levels = C.load_tolerances().levels
    p = np.asarray(ps)
    rates = {float(al): float(np.mean(p <= al + 1e-12)) for al in levels}
    over = {al: r for al, r in rates.items() if r > C.level_bound(al, R)}
    path = C.log_dir() / f"{SUITE}_attack_null.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps({"arm": arm, "model": model, "R": R, "B": B, "rates": rates,
                             "flagged": flagged, "p": ps}) + "\n")  # fmt: skip
    if flagged:
        if over:
            warnings.warn(f"{arm}/{model}: flagged {UNCALIBRATED} null rejects at {over}",
                          stacklevel=1)  # fmt: skip
        return
    assert not over, f"{arm}/{model}: unflagged null p rejects above level at {over}"
