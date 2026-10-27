"""Calibration of the M8-08 drop-one-block ΔR² (GIPPA_SPEC.md §7.8.2 step 5, §8.4; suite
``determinants``; docs/specs/M8-08-determinants-models.md example 1 over R replicates).

Replicate ``i`` of stream ``determinants.m808.<name>`` draws the spec's example-1 world from
``_criteria.replicate_rng``: 400 genes in TUs of 1-5, ``cai`` ~ N(0, 1) (block ``codon``) and five
N(0, 1) noise columns in the other four blocks, y = s * cai + N(0, 1); folds from ``tu_folds`` on
the same generator; §5.2 defaults (5 x GroupKFold(5), 500 trees). R is
``replicates_fast/full.determinants`` (20 / 100). Every replicate refits both models six times
(~25 s on one core), so the fast mode is ``calibration_nightly``.

Gated:

- ``test_planted_block_recovered`` (s = 2, true R² 0.8): the codon CI excludes 0 and the codon ΔR²
  is the largest block in every replicate (both models); the mean elastic-net ΔR²(codon) is
  within ``determinants_delta_r2_max_abs_bias`` of 0.8 and every noise block's mean ΔR² within it
  of 0 (both models). LightGBM's ΔR²(codon) is **recorded**, not gated against 0.8: its drop
  model overfits the noise columns (CV R² ≈ -0.25), so ΔR² ≈ 0.98 on these data
  (test-author pre-run, see the M8-08 report).
- ``test_null_delta_r2_near_zero`` (s = 0): every block's mean ΔR² within the same bound of 0.

Recorded, not gated (spec acceptance: "5 repeats share data, so the interval is optimistic"):
the share of null replicates whose 95 % t-interval misses 0, pooled over the five blocks, per
model, against ``level_bound(0.05, 5 R)``; a warning flags it above the bound. The pre-run
measured about 0.4 for LightGBM and 0.05 for the elastic net (R = 40).

``determinants_delta_r2_max_abs_bias`` is not in the protected tolerances.yaml yet: the gated
tests fail with that message until the lead adds it (proposed 0.05, the spec's example-1 margin).
"""

from __future__ import annotations

import json
import math
import sys
import time
import warnings
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import determinants_models_data as d

SUITE = "determinants"
BIAS_KEY = "determinants_delta_r2_max_abs_bias"
N = 400
TRUE_R2 = 0.8
BLOCKS = ("codon", "transcription", "translation", "protein", "context")
NOISE = BLOCKS[1:]

pytestmark = pytest.mark.calsuite(SUITE)
MODES = [
    pytest.param(False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly],
                 id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]  # fmt: skip


def _bias_tol(tolerances: dict[str, object]) -> float:
    if BIAS_KEY not in tolerances:
        pytest.fail(
            f"tests/calibration/tolerances.yaml has no {BIAS_KEY!r} (protected; lead to add)"
        )
    return float(tolerances[BIAS_KEY])  # type: ignore[arg-type]


def _replicate(rng: np.random.Generator, signal: float, i: int) -> dict[str, dict[str, Any]]:
    """One world -> per model: mean CV R² and per block (ΔR² mean, ci_low, ci_high)."""
    X, y, tu, blocks = d.planted(rng, n=N, signal=signal)
    cfg = d.config(fast=False)
    m = d.mod()
    folds = m.tu_folds(tu, cfg.outer_folds, cfg.repeats, rng)
    cv = m.cross_validate(X, y, folds, blocks, cfg, C.load_tolerances().calibration_base + i, 1)
    t = stats.t.ppf(0.975, cfg.repeats - 1)
    out: dict[str, dict[str, Any]] = {}
    for name in d.MODELS:
        res = cv[name]
        delta = {}
        for b in BLOCKS:
            v = np.asarray(res.delta_r2[b], dtype=float)
            half = t * v.std(ddof=1) / math.sqrt(v.size)
            delta[b] = (float(v.mean()), float(v.mean() - half), float(v.mean() + half))
        out[name] = {"r2": float(np.mean(res.r2)), "delta": delta}
    return out


def _run(stream: str, signal: float, R: int) -> list[dict[str, dict[str, Any]]]:
    return [_replicate(C.replicate_rng(stream, i), signal, i) for i in range(R)]


def _log(record: dict[str, Any]) -> None:
    path = C.log_dir() / f"{SUITE}_m808.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")


@pytest.mark.parametrize("full", MODES)
def test_planted_block_recovered(full: bool, tolerances: dict[str, object]) -> None:
    tol = _bias_tol(tolerances)
    R = C.replicates(SUITE, full=full)
    t0 = time.perf_counter()
    reps = _run(f"{SUITE}.m808.planted", 2.0, R)
    summary: dict[str, Any] = {}
    for name in d.MODELS:
        codon = np.array([r[name]["delta"]["codon"] for r in reps])
        assert np.all(codon[:, 1] > 0), f"{name}: codon ΔR² CI covers 0 in some replicate"
        for r in reps:
            dl = r[name]["delta"]
            assert max(BLOCKS, key=lambda b: dl[b][0]) == "codon", name
        for b in NOISE:
            mean = float(np.mean([r[name]["delta"][b][0] for r in reps]))
            assert abs(mean) <= tol, f"{name}: noise block {b} mean ΔR² {mean:+.4f}"
        summary[name] = {
            "cv_r2_mean": float(np.mean([r[name]["r2"] for r in reps])),
            "delta_codon_mean": float(codon[:, 0].mean()),
        }
    en = summary["elastic_net"]["delta_codon_mean"]
    _log({"suite": SUITE, "claim": "planted_codon_delta_r2", "R": R, "true": TRUE_R2,
          "summary": summary, "seconds": time.perf_counter() - t0})  # fmt: skip
    print(f"planted ΔR²(codon), R={R}: {summary}")
    assert abs(en - TRUE_R2) <= tol, f"elastic-net mean ΔR²(codon) {en:.4f} vs {TRUE_R2}"


@pytest.mark.parametrize("full", MODES)
def test_null_delta_r2_near_zero(full: bool, tolerances: dict[str, object]) -> None:
    tol = _bias_tol(tolerances)
    R = C.replicates(SUITE, full=full)
    t0 = time.perf_counter()
    reps = _run(f"{SUITE}.m808.null", 0.0, R)
    bound = C.level_bound(0.05, R * len(BLOCKS))
    recorded: dict[str, Any] = {}
    for name in d.MODELS:
        ci = np.array([[r[name]["delta"][b] for b in BLOCKS] for r in reps])  # (R, blocks, 3)
        miss = float(np.mean(~((ci[..., 1] <= 0.0) & (ci[..., 2] >= 0.0))))
        means = ci[..., 0].mean(axis=0)
        recorded[name] = {"miss_rate": miss, "bound": bound,
                          "mean_delta": dict(zip(BLOCKS, means.tolist(), strict=True))}  # fmt: skip
        if miss > bound:
            warnings.warn(
                f"M8-08 recorded: {name} ΔR² t-interval misses 0 in {miss:.3f} of null "
                f"block-replicates (bound {bound:.3f}); the interval is optimistic (spec)",
                stacklevel=1,
            )
        for b, mean in zip(BLOCKS, means, strict=True):
            assert abs(mean) <= tol, f"{name}: null block {b} mean ΔR² {mean:+.4f}"
    _log({"suite": SUITE, "claim": "null_delta_r2", "R": R, "recorded": recorded,
          "seconds": time.perf_counter() - t0})  # fmt: skip
    print(f"null ΔR², R={R}: {recorded}")
