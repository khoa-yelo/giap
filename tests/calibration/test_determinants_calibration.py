"""M8 gate, half 2: the ``determinants`` null at level and recovery (docs/specs/M8-15-calibration-
gate.md § Suites and criteria, Behaviour 2, examples 5-6, lead decisions 1-2; GIAP_SPEC.md
§7.8.3, §8.4 row "Determinants null", §11.3; M8-14 scenarios; ADR-0015, ADR-0034).

**Null (gate).** ``determinants_null(400)`` fast / ``(1500)`` full: the response is the TU effect
plus noise, independent of every feature, with TU sizes and the detection pattern kept. Model
rows are the detected genes (``detected_true``), y = their ``abundance_log2_true``; X = the
simulator's planted features (``promoter``, ``rbs``, ``cai``; M8-14 lead decision 2) plus the
TU-structure columns ``tu_position`` and ``tu_size``. The observed mean CV R² from
``cross_validate`` on ``tu_folds``, then ``null_p`` (TU-size-preserving permutations, repeat-1
folds). Level of ``null_p`` per model at every alpha of ``levels``, one p per model per replicate
(R = ``replicates_fast/full.determinants`` = 20 / 100; example 5: at most 2 rejections at 0.05 for
R = 20). KS recorded (lead decision 2).

* Fast tier (lead decision 1): 400 genes, B = 99, LightGBM ``n_estimators = 50``, elastic net as
  configured, **1 CV repeat** for the observed value (generator settings, not tolerances: R3).
* Full tier: the §5.2 parameters (``DeterminantsConfig()``: 5 repeats, 500 trees), B = 999,
  1,500 genes; Slurm array task ``determinants``, sharded (``--shard K/N``).

**Recovery (recorded; spec table).** ``determinants_planted(1500, coefs, 1.0)`` at true R² 0.2 and
0.5 (the M8-14 coefficient pairs), the §5.2 parameters, R = ``replicates_fast.determinants``:
LightGBM ``cv_r2_mean`` against the recorded truth (band +-0.1, example 6), and the SHAP dominant
block of the top-decile genes (``attribution_table``) against the planted truth (>= 70 %, example
6): (a) the per-gene truth of ``truth/determinants_genes.tsv`` (largest |z beta|), (b) the signed
argmax of z beta (GIAP's rule takes the largest **positive** SHAP sum; M8-14 backlog "signed-argmax
comparison"), (c) the block of the largest beta. All recorded in
``function_gate_diagnostics.jsonl`` with a warning outside the band; nothing gated.

Pre-run (test-author, iron-01, 1-2 threads): per null replicate at 400 genes (B = 99, 50 trees)
~20 s for LightGBM and ~50 s for the elastic net, so the fast tier runs ``null_p`` on up to 4
threads; recovery ~25 s per replicate with 4 threads. Two recovery replicates at true R² 0.5
gave LightGBM CV R² 0.32 and 0.41 (elastic net 0.45 and 0.52).
"""

from __future__ import annotations

import os
import time
import warnings
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers import function_gate_chain as F

from giap.config.model import DeterminantsConfig, LightgbmConfig

SUITE = "determinants"
MODELS = ("lightgbm", "elastic_net")
N_FAST, N_FULL, N_RECOVERY = 400, 1500, 1500
B_FAST, B_FULL = 99, 999
R2_BAND = 0.1  # example 6
DOMINANT_MIN = 0.70  # example 6
RECOVERY = {
    "r2_0.2": ({"promoter": 0.5, "rbs": 0.4, "cai": 0.3}, 0.2),
    "r2_0.5": ({"promoter": 1.0, "rbs": 0.8, "cai": 0.6}, 0.5),
}

pytestmark = pytest.mark.calsuite(SUITE)
MODES = [
    pytest.param(False, marks=[pytest.mark.calibration, pytest.mark.calibration_nightly],
                 id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]  # fmt: skip


def _fs() -> Any:
    import giapsim.function_sim as fs

    return fs


def _sc() -> Any:
    import giapsim.scenarios as sc

    return sc


def _cfg(full: bool) -> DeterminantsConfig:
    if full:
        return DeterminantsConfig()
    return DeterminantsConfig(repeats=1, lightgbm=LightgbmConfig(n_estimators=50))


def _threads(full: bool) -> int:
    """Results do not depend on ``threads`` (M8-09 rule 6); capped per tier."""
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))
    return max(1, min(8, cpus))  # lead decision 8 gap 3: the fast tier gets 8 threads too


def model_inputs(w: Any) -> tuple[pd.DataFrame, np.ndarray, list[str], list[str]]:
    """Model rows (detected genes, layout order): X, y, TU per row, locus tags."""
    fs = _fs()
    tu, pos = F.tu_maps(w)
    size = {str(t): len(gs) for t, gs in zip(w.layout.tus["tu_id"], w.layout.tus["genes"],
                                             strict=True)}  # fmt: skip
    f = w.features.set_index(w.features["locus_tag"].astype(str))
    e = w.expression.genes.set_index(w.expression.genes["locus_tag"].astype(str)).loc[f.index]
    det = e["detected_true"].to_numpy(bool)
    tags = [str(t) for t in f.index[det]]
    X = f.loc[tags, list(fs.FEATURES)].astype(float).reset_index(drop=True)
    X["tu_position"] = [float(pos[t]) for t in tags]
    X["tu_size"] = [float(size[tu[t]]) for t in tags]
    y = e["abundance_log2_true"].to_numpy(float)[det]
    return X, y, [tu[t] for t in tags], tags


def _null_p(i: int, model: str, full: bool) -> float:
    from giap.determinants import models as m
    from giap.determinants import null as nl

    stream = f"{SUITE}.m815.null.sim{'.full' if full else ''}"
    rng = C.replicate_rng(stream, i)
    w = _fs().draw(_sc().determinants_null(N_FULL if full else N_FAST), rng)
    X, y, tus, _ = model_inputs(w)
    cfg = _cfg(full)
    folds = np.asarray(m.tu_folds(tus, cfg.outer_folds, cfg.repeats, rng))
    seed = C.load_tolerances().calibration_base + int(rng.integers(0, 2**31 - 1))
    cv = m.cross_validate(X, y, folds, {}, cfg, seed, _threads(full), tu=tus)
    obs = {model: float(np.mean(cv[model].r2))}
    B = B_FULL if full else B_FAST
    res = nl.null_p(X, y, tus, folds[0], obs, cfg, B, seed, _threads(full))[model]
    assert res.n_perm == B and 0.0 < res.p <= 1.0, res
    assert "null_degenerate" not in res.flags, res.flags
    return float(res.p)


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("full", MODES)
def test_determinants_null_level_per_model(full: bool, model: str) -> None:
    R = C.replicates(SUITE, full)
    levels = list(C.load_tolerances().levels)
    stream = f"{SUITE}.m815.null.{model}{'.full' if full else ''}"
    got: dict[int, float] = {}

    def one(i: int) -> float:
        got[i] = _null_p(i, model, full)
        return got[i]

    gate = calibrate_levels(SUITE, "level", levels, F.by_index(stream, R, one), R, stream=stream)
    if got:
        p = [got[i] for i in sorted(got)]
        F.diagnostics({"check": f"determinants.null.{model}", "full": full, "R_run": len(p),
                       "B": B_FULL if full else B_FAST, "n_genes": N_FULL if full else N_FAST,
                       "rates": F.rate_record(p, levels), "ks": F.ks_record(p),
                       "mean_p": float(np.mean(p)),
                       "block_seconds": gate.block_seconds})  # fmt: skip
    gate.assert_ok()


# --------------------------------------------------------------------------- recovery


def _recovery(i: int, coefs: dict[str, float], r2: float) -> dict[str, Any]:
    from giap.determinants import attribution as at
    from giap.determinants import models as m

    fs = _fs()
    rng = C.replicate_rng(f"{SUITE}.m815.recovery.r2_{r2:g}", i)
    w = fs.draw(_sc().determinants_planted(N_RECOVERY, coefs, 1.0), rng)
    true = float(w.truth["determinants"]["true_r2"].iloc[0])
    assert abs(true - r2) < 1e-12, true
    X, y, tus, tags = model_inputs(w)
    blocks = {f: str(fs.BLOCKS[f]) for f in fs.FEATURES}
    blocks.update({"tu_position": "transcription", "tu_size": "transcription"})
    cfg = DeterminantsConfig()
    folds = np.asarray(m.tu_folds(tus, cfg.outer_folds, cfg.repeats, rng))
    seed = C.load_tolerances().calibration_base + int(rng.integers(0, 2**31 - 1))
    cv = m.cross_validate(X, y, folds, blocks, cfg, seed, _threads(True), tu=tus)
    Xi = X.set_axis(pd.Index(tags, name="locus_tag"))
    shap = at.shap_blocks(cv["lightgbm"], Xi, blocks)
    feats = pd.DataFrame({"locus_tag": tags, "abundance_log2": y, "ribosomal_protein": False,
                          "detected": True, "abundance_imputed": False})  # fmt: skip
    tab = at.attribution_table(feats, np.ones(len(tags), dtype=bool), cv["lightgbm"], shap)
    top = tab[tab["top_decile"].astype(bool)].set_index("locus_tag")
    # no positive SHAP sum -> NA dominant block: counted as a disagreement
    got = np.array([None if pd.isna(v) else str(v) for v in top["dominant_block"]], dtype=object)
    truth_abs = (w.truth["determinants_genes"].set_index("locus_tag")
                 .loc[top.index, "dominant_block"].to_numpy(object))  # fmt: skip
    beta = np.array([coefs.get(f, 0.0) for f in fs.FEATURES])
    feat = w.features.set_index(w.features["locus_tag"].astype(str))
    z = np.column_stack([(feat[f] - feat[f].mean()) / feat[f].std(ddof=0) for f in fs.FEATURES])
    contrib = pd.DataFrame(z * beta, index=feat.index).loc[top.index].to_numpy(float)
    names = np.array([fs.BLOCKS[f] for f in fs.FEATURES], dtype=object)
    truth_signed = names[contrib.argmax(axis=1)]
    largest = names[int(np.argmax(beta))]
    return {
        "true_r2": true,
        "cv_r2": {k: float(np.mean(v.r2)) for k, v in cv.items()},
        "n_top": len(top),
        "dominant_vs_truth_abs": float(np.mean(got == truth_abs)),
        "dominant_vs_truth_signed": float(np.mean(got == truth_signed)),
        "dominant_vs_largest_beta": float(np.mean(got == largest)),
    }


@pytest.mark.calibration
@pytest.mark.calibration_nightly
@pytest.mark.parametrize("arm", sorted(RECOVERY))
def test_determinants_recovery_recorded(arm: str) -> None:
    coefs, r2 = RECOVERY[arm]
    R = C.replicates(SUITE, False)
    t0 = time.perf_counter()
    reps = [_recovery(i, coefs, r2) for i in range(R)]
    lgb = float(np.mean([r["cv_r2"]["lightgbm"] for r in reps]))
    en = float(np.mean([r["cv_r2"]["elastic_net"] for r in reps]))
    dom = {
        k: float(np.mean([r[k] for r in reps]))
        for k in ("dominant_vs_truth_abs", "dominant_vs_truth_signed", "dominant_vs_largest_beta")
    }
    rec = {"check": f"determinants.recovery.{arm}", "R": R, "true_r2": r2,
           "lightgbm_cv_r2_mean": lgb, "elastic_net_cv_r2_mean": en,
           "lightgbm_in_band": abs(lgb - r2) <= R2_BAND, "band": R2_BAND,
           "dominant": dom, "dominant_min": DOMINANT_MIN,
           "dominant_meets": dom["dominant_vs_truth_signed"] >= DOMINANT_MIN,
           "per_replicate": reps, "seconds": time.perf_counter() - t0}  # fmt: skip
    F.diagnostics(rec)
    if not rec["lightgbm_in_band"]:
        warnings.warn(f"recovery {arm}: LightGBM CV R² {lgb:.3f} outside {r2} +- {R2_BAND} "
                      "(recorded, spec table)", stacklevel=1)  # fmt: skip
    if not rec["dominant_meets"]:
        warnings.warn(f"recovery {arm}: dominant block agrees in {dom} of top-decile genes "
                      f"(< {DOMINANT_MIN}; recorded)", stacklevel=1)  # fmt: skip
    assert all(np.isfinite(r["cv_r2"]["lightgbm"]) and r["n_top"] > 0 for r in reps), reps
