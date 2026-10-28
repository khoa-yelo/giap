"""DL1 calibration in the large, recorded (docs/specs/M6-02-detection-limit-model.md lead decision
6; docs/methods/detection-limit.md "Calibration criteria").

Same design and streams as DL1 in ``test_detection_limit_calibration.py`` (stream
``detection_limit.gene_slope``, R = ``replicates_{fast,full}.detectability``); per replicate the
CITL = c in d ~ Bernoulli(expit(logit(p_detect) + c)) over every gene, scored with
``model_for_tu``. Appended to ``<log_dir>/detectability.jsonl`` as criterion
``detection_limit_citl_mean``. **Record only, no gate**: a CITL tolerance (proposed +-0.2) needs a
``tolerances.yaml`` key, Khoa's decision. The test asserts only that the record is written and the
values are finite.
"""

from __future__ import annotations

import json
import os
import time
from typing import Any

import _criteria as C
import numpy as np
import pytest
from helpers import detection_limit_data as D

from gippa.config import GippaConfig
from gippa.rng import for_step

SUITE = "detectability"
STREAM = "detection_limit.gene_slope"
BETA = {"intercept": -10.0, "abundance_log2": 0.5, "length": -0.03, "gravy": 0.8}
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _citl(d: np.ndarray, p: np.ndarray) -> float:
    pc = np.clip(p, 1e-6, 1 - 1e-6)
    z = np.log(pc) - np.log1p(-pc)
    c = 0.0
    for _ in range(100):
        mu = 1.0 / (1.0 + np.exp(-(z + c)))
        step = float(np.sum(d - mu) / np.sum(mu * (1 - mu)))
        c += step
        if abs(step) < 1e-12:
            break
    return c


def _replicate(rng: np.random.Generator, n_genes: int) -> tuple[float, float]:
    from gippa.operons import detection_limit as dl

    s = D.simulate_peptides(rng, n_genes, BETA, a_mu=15.0, a_sd=2.5, n_peps=(5, 20), tu_size=2)
    seed = int(rng.integers(0, 2**31 - 1))
    model: Any = dl.fit_peptide_model(
        s.genes, s.hits, GippaConfig(seed=seed), for_step(D.STEP, seed)
    )
    assert model.usable
    p = np.array(
        [
            dl.p_detect(dl.model_for_tu(model, tu), a, s.peps(t))
            for t, tu, a in zip(
                s.genes.locus_tag, s.genes.tu_id, s.genes.abundance_log2, strict=True
            )
        ]
    )
    d = s.genes.detected.to_numpy(float)
    return _citl(d, p), float(p.mean() - d.mean())


@pytest.mark.calsuite(SUITE)
@pytest.mark.parametrize("full", MODES)
def test_dl1_citl_recorded(full: bool) -> None:
    if os.environ.get(C.CONFIRM_ENV) is not None:
        pytest.skip("record only; nothing to confirm")
    R = C.replicates(SUITE, full)
    if C.replicate_range(R).start != 0:
        pytest.skip("runs whole in the shard holding replicate 0")
    n_genes = 1500 if full else 600
    t0 = time.perf_counter()
    out = C.run_replicates(STREAM, lambda rng: _replicate(rng, n_genes), R)
    citl, gap = (np.array(v, float) for v in zip(*out, strict=True))
    assert np.isfinite(citl).all() and np.isfinite(gap).all()
    record = {
        "suite": SUITE,
        "stream": STREAM,
        "criterion": "detection_limit_citl_mean",
        "level": 0.0,
        "observed": float(citl.mean()),
        "R": R,
        "seed": C.suite_seed(STREAM, 0),
        "duration": time.perf_counter() - t0,
        "status": "recorded",
        "n_genes": n_genes,
        "citl": citl.tolist(),
        "mean_p_minus_mean_d": gap.tolist(),
    }
    directory = C.log_dir()
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{SUITE}.jsonl"
    with path.open("a") as fh:
        fh.write(json.dumps(record) + "\n")
    last = path.read_text().strip().splitlines()[-1]
    assert json.loads(last)["criterion"] == "detection_limit_citl_mean"
