"""Calibration of the M8-09 TU-size-preserving permutation null (GIPPA_SPEC.md §7.8.2 step 4,
§8.1 row "Determinant model", §8.4 level criterion; suite ``determinants``;
docs/specs/M8-09-determinants-attribution-null.md rule 6, decision 16; M8-index decisions 23-24).

Replicate ``i`` of stream ``determinants.m809.<arm>.<model>`` draws one world from
``_criteria.replicate_rng``: 400 genes in TUs of 1-5, the observed mean CV R² from
``cross_validate`` on ``tu_folds`` (5 x GroupKFold(5)), and ``null_p`` on the repeat-1 folds with
B = 99 (fast) / 199 (full), seed ``calibration_base + i``. Fast tier: 50 trees (index decision 23,
"reduced model size"); full tier: the §5.2 defaults. R = ``replicates_fast/full.determinants``.
Each model is its own check (the two p of one world are not independent, so no pooling).

Worlds (features never carry information about y beyond what the null keeps):

- ``clustered``: y = TU effect + gene noise (ICC 0.5); six N(0, 1) gene-level columns and one
  TU-level column, all independent of y. The TU-size-preserving permutation is exchangeable
  here, so p is valid (and conservative: one null repeat against a 5-repeat observed mean).
- ``size``: y = 0.8 * TU size + TU effect + noise; columns ``tx_rate_sum_tu`` (= TU size) and
  ``tu_position`` carry the TU-size structure, five noise columns. The null keeps y paired with
  same-size TUs, so the size signal is part of the null and p stays at level.

Gated:

- ``test_null_level`` (both worlds, both models): §8.4 level at every alpha in ``levels``
  (rate of p <= alpha <= alpha + k sqrt(alpha (1 - alpha) / R); index decision 24: level, not KS).
- ``test_teeth_gene_level_null_fails`` (``size`` world, both models): the same refits with a
  gene-level permutation (``null_p`` given one TU per gene, which ignores the TU structure) reject
  far above level — a confirmed ``fail`` (ADR-0015, both blocks), so the level gate can fail.
- ``test_power_planted`` (y = cai + N(0, 1), true R² 0.5): p = 1/(1 + B) in every replicate for
  both models — the null is not trivially conservative (an identity permutation or p = 1 passes
  the level gate; it fails here).

Test-author pre-run (a prototype of the spec's rules, Slurm job 1865936, one core per test, fast
tier R = 20): every fast test passed; about 16.5 min per level / power test (about 49 s per
replicate per model) and 25 min per teeth test (both blocks run, as a confirmed fail must).
Rejection rates at alpha 0.01 / 0.05 / 0.10: clustered 0 / 0 / 0.05 (LightGBM) and 0 / 0 / 0
(elastic net); size 0 / 0.05 / 0.10 and 0 / 0 / 0.05; mean p 0.49-0.53 (conservative, as the spec
states); gene-level teeth p = 0.01 in every replicate; planted p = 0.01 in every replicate. The
full tier (R = 100, B = 199, 500 trees) is not measured; it is ``calibration_full`` only.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import determinants_attrib_data as a
import determinants_models_data as d

from gippa.config.model import DeterminantsConfig, LightgbmConfig

SUITE = "determinants"
N = 400
ALPHA = 0.05
MODELS = ("lightgbm", "elastic_net")

pytestmark = pytest.mark.calsuite(SUITE)
MODES = [
    pytest.param(False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly],
                 id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]  # fmt: skip


def _B(full: bool) -> int:
    return 199 if full else 99


def _cfg(full: bool) -> DeterminantsConfig:
    return (
        DeterminantsConfig()
        if full
        else DeterminantsConfig(lightgbm=LightgbmConfig(n_estimators=50))
    )


def _threads(full: bool) -> int:
    """Results do not depend on ``threads`` (rule 6, example 8); the full tier uses the job's
    CPUs (``SLURM_CPUS_PER_TASK``, at most 8)."""
    return min(8, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))) if full else 1


def _world(rng: np.random.Generator, kind: str) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    sizes = d.mixed_sizes(rng, N, 1, 5)
    tu = d.tu_ids(sizes)
    size = np.repeat(sizes, sizes).astype(float)
    pos = np.concatenate([np.arange(1, s + 1) for s in sizes]).astype(float)
    u = np.repeat(rng.normal(0.0, 1.0, len(sizes)), sizes)
    X, _, _, _ = d.planted(rng, n=N, signal=0.0)
    X["gc_local_dev"] = np.repeat(rng.normal(0.0, 1.0, len(sizes)), sizes)  # TU-level column
    if kind == "clustered":
        y = 0.7 * u + rng.normal(0.0, 0.7, N)
    elif kind == "size":
        X["tu_position"] = pos
        X["tx_rate_sum_tu"] = size
        y = 0.8 * size + 0.7 * u + rng.normal(0.0, 0.7, N)
    elif kind == "planted":
        y = X["cai"].to_numpy() + rng.normal(0.0, 1.0, N)
    else:
        raise ValueError(kind)
    return X, np.asarray(y, dtype=float), tu


def _p(rng: np.random.Generator, kind: str, model: str, full: bool, gene_level: bool = False,
       seed_offset: int = 0) -> float:  # fmt: skip
    X, y, tu = _world(rng, kind)
    cfg = _cfg(full)
    m = d.mod()
    folds = np.asarray(m.tu_folds(tu, cfg.outer_folds, cfg.repeats, rng))
    seed = C.load_tolerances().calibration_base + int(rng.integers(0, 2**31 - 1)) + seed_offset
    cv = m.cross_validate(X, y, folds, {}, cfg, seed, _threads(full), tu=tu)
    obs = {model: float(np.mean(cv[model].r2))}
    null_tu = [f"gene_{j:05d}" for j in range(N)] if gene_level else tu
    res = a.null_mod().null_p(X, y, null_tu, folds[0], obs, cfg, _B(full), seed, _threads(full))
    p = float(res[model].p)
    assert res[model].n_perm == _B(full)
    assert 0.0 < p <= 1.0, p
    return p


def _log(record: dict[str, Any]) -> None:
    path = C.log_dir() / f"{SUITE}_m809.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("kind", ["clustered", "size"])
@pytest.mark.parametrize("full", MODES)
def test_null_level(full: bool, kind: str, model: str) -> None:
    R = C.replicates(SUITE, full=full)
    levels = C.load_tolerances().levels
    t0 = time.perf_counter()
    gate = calibrate_levels(
        SUITE, "level", levels, lambda g: _p(g, kind, model, full), R,
        stream=f"{SUITE}.m809.{kind}.{model}{'.full' if full else ''}",
    )  # fmt: skip
    _log({"suite": SUITE, "claim": f"null_level_{kind}", "model": model, "R": R, "B": _B(full),
          "mean_p": None if gate.mean() is None else float(gate.mean()[0]),
          "block_seconds": gate.block_seconds, "seconds": time.perf_counter() - t0})  # fmt: skip
    gate.assert_ok()


@pytest.mark.calibration
@pytest.mark.calibration_nightly
@pytest.mark.parametrize("model", MODELS)
def test_teeth_gene_level_null_fails(model: str, teeth_log: Path) -> None:
    rec = C.calibrate(
        SUITE, "level", ALPHA, lambda g: _p(g, "size", model, False, gene_level=True),
        C.replicates(SUITE, full=False), stream=f"{SUITE}.m809.teeth_gene_level.{model}",
        log_dir=teeth_log,
    )  # fmt: skip
    assert rec["status"] == "fail", rec


@pytest.mark.calibration
@pytest.mark.calibration_nightly
@pytest.mark.parametrize("model", MODELS)
def test_power_planted(model: str) -> None:
    R = C.replicates(SUITE, full=False)
    t0 = time.perf_counter()
    ps = C.run_replicates(f"{SUITE}.m809.planted.{model}",
                          lambda g: _p(g, "planted", model, False), R)  # fmt: skip
    _log({"suite": SUITE, "claim": "null_power_planted", "model": model, "R": R, "p": ps,
          "seconds": time.perf_counter() - t0})  # fmt: skip
    floor = 1.0 / (1 + _B(False))
    assert all(p == pytest.approx(floor) for p in ps), ps
