"""Calibration of the M6-02 detection-limit model (docs/specs/M6-02-detection-limit-model.md
acceptance: "Calibration (``detectability``-style, R = 20 fast / 100 full): with a logistic truth,
the calibration slope of ``p_detect`` is 1 ± 0.2 — the same criterion M4-13 already meets").

Suite ``detectability`` (``calsuite`` marker; no new ``tolerances.yaml`` key: R =
``replicates_{fast,full}.detectability``, bounds = ``detectability_slope`` [0.8, 1.2]). Replicate i
of stream s draws from ``_criteria.replicate_rng(s, i)``; the model's folds come from
``gippa.rng.for_step("operons.detection_limit", seed_i)``.

Claims:

DL1 Gene-level calibration slope of the cross-fitted ``p_detect`` (the gate). Truth: per-peptide
    logit = -10 + 0.5 a - 0.03 L + 0.8 GRAVY on hand-built tryptic proteins (5-20 observable
    peptides per gene, a ~ N(15, 2.5), TUs of 2 genes; 600 genes fast, 1,500 full); detected =
    >= 1 observed peptide (tier ``mayu``, k_min = 1); C11 = the observed peptides. Fit
    ``fit_peptide_model``; score every gene, detected or not, with ``model_for_tu(model, tu)`` (the
    model that did not see its TU); slope of the unpenalised recalibration d ~ a + b logit(p).
    Criterion: mean slope over R in [0.8, 1.2]. Intercepts and Brier scores are recorded in
    ``<log_dir>/detectability.jsonl`` (criterion ``detection_limit_slope_mean``). Recorded, not
    gated: the spec trains on detected genes only, which selects on the response; an oracle
    (unpenalised MLE on the same design) gives mean slope ~1.07 but mean intercept ~-0.97, i.e.
    p_detect too high and P(miss) too small (see the test-author report).
DL2 Teeth (no fitting): with the TRUE coefficients p_detect is calibrated (slope within 0.2 of 1
    on 4,000 genes); with every slope coefficient halved (mean logit kept) the slope is > 1.4
    (oracle ~1.84), so the DL1 criterion can fail.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import detectability_oracle as O
from helpers import detection_limit_data as D

from gippa.config import GippaConfig
from gippa.rng import for_step

SUITE = "detectability"
BETA = {"intercept": -10.0, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _dl() -> Any:
    from gippa.operons import detection_limit

    return detection_limit


def _replicate(rng: np.random.Generator, n_genes: int) -> tuple[float, float, float]:
    s = D.simulate_peptides(rng, n_genes, BETA, a_mu=15.0, a_sd=2.5, n_peps=(5, 20), tu_size=2)
    seed = int(rng.integers(0, 2**31 - 1))
    dl = _dl()
    cfg = GippaConfig(seed=seed)
    model = dl.fit_peptide_model(s.genes, s.hits, cfg, for_step(D.STEP, seed))
    assert model.usable and model.k_min == 1
    p = np.array(
        [
            dl.p_detect(dl.model_for_tu(model, tu), a, s.peps(t))
            for t, tu, a in zip(
                s.genes.locus_tag, s.genes.tu_id, s.genes.abundance_log2, strict=True
            )
        ]
    )
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    d = s.genes.detected.to_numpy(float)
    icpt, slope = O.logistic_recalibration(d, p)
    return slope, icpt, float(np.mean((p - d) ** 2))


@pytest.mark.calsuite(SUITE)
@pytest.mark.parametrize("full", MODES)
def test_dl1_mean_gene_level_slope_in_tolerance(full: bool, tolerances: dict[str, object]) -> None:
    if os.environ.get(C.CONFIRM_ENV) is not None:
        pytest.skip("DL1 is not a calibrate() check; nothing to confirm")
    R = C.replicates(SUITE, full)
    if C.replicate_range(R).start != 0:
        pytest.skip("DL1 runs whole in the shard holding replicate 0")
    lo, hi = (float(x) for x in tolerances["detectability_slope"])  # type: ignore[attr-defined]
    n_genes = 1500 if full else 600
    stream = "detection_limit.gene_slope"
    t0 = time.perf_counter()
    out = C.run_replicates(stream, lambda rng: _replicate(rng, n_genes), R)
    slopes, icpts, briers = (np.array(v, float) for v in zip(*out, strict=True))
    mean = float(slopes.mean())
    ok = lo <= mean <= hi
    record = {
        "suite": SUITE,
        "stream": stream,
        "criterion": "detection_limit_slope_mean",
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
        "intercept_mean": float(icpts.mean()),
        "briers": briers.tolist(),
    }
    directory = C.log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / f"{SUITE}.jsonl").open("a") as fh:
        fh.write(json.dumps(record) + "\n")
    assert ok, f"mean slope {mean:.4f} outside [{lo}, {hi}]; slopes {np.round(slopes, 3)}"


@pytest.mark.calsuite(SUITE)
@pytest.mark.calibration
def test_dl2_teeth_true_model_calibrated_flattened_model_fails() -> None:
    rng = C.replicate_rng("detection_limit.teeth", 0)
    s = D.simulate_peptides(rng, 4_000, BETA, n_peps=(5, 20), tu_size=1)
    dl = _dl()
    d = s.genes.detected.to_numpy(float)

    def slope(coef: dict[str, float]) -> float:
        m = dl.PMissModel(kind="peptide", coef=coef, k_min=1, n_fit=1, auc=math.nan, usable=True)
        p = [
            dl.p_detect(m, a, s.peps(t))
            for t, a in zip(s.genes.locus_tag, s.genes.abundance_log2, strict=True)
        ]
        return O.logistic_recalibration(d, p)[1]

    assert abs(slope(dict(BETA)) - 1.0) < 0.2
    # every slope coefficient halved, the intercept moved so the mean logit is unchanged
    peps = [x for t in s.genes.locus_tag for x in s.peps(t)]
    centre = (
        BETA["abundance_log2"] * float(s.genes.abundance_log2.mean())
        + BETA["length"] * float(np.mean([n for n, _ in peps]))
        + BETA["gravy"] * float(np.mean([g for _, g in peps]))
    )
    flat = {k: 0.5 * v for k, v in BETA.items()}
    flat["intercept"] = BETA["intercept"] + 0.5 * centre
    assert slope(flat) > 1.4
