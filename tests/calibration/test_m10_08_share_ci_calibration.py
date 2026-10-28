"""Calibration of the A3 ``share_ci95`` as the evaluator reports it (M10-08b, R4, spec §8.4).

Claim (docs/specs/M10-08-a3-feasibility.md, step 14): ``share_ci95`` is the Wilson 95 % interval
of ``share = n_match / n_top``. ``gippa.stats.cutoffs.wilson_ci`` is calibrated in
``tests/calibration/test_cutoffs_calibration.py``; this test checks the interval **as written by
``a3_feasibility.py evaluate``** (after any rounding), end to end on synthetic inputs: for
``n_top = 40`` the evaluator is run once for every ``k = 0..40`` matches, and the exact binomial
coverage of the reported intervals, averaged over p in [0.05, 0.95], must be 0.95 within 0.01
(Wilson's mean coverage at n = 40 is 0.95 to about 0.003; the margin absorbs 4-decimal rounding).
Deterministic (exact coverage, no Monte Carlo); fast (41 small evaluator runs).
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np
import pytest
from helpers.a3_eval_kit import (
    Genome,
    Row,
    ladder,
    make_case,
    planted_contig,
    random_peptides,
)
from scipy.stats import binom

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("regions")]

N_TOP = 40


def test_reported_share_ci_has_nominal_exact_coverage(tmp_path: Path) -> None:
    rng = random.Random(20261024)
    peps = random_peptides(N_TOP, rng)
    rows = [Row(p, s) for p, s in zip(peps, ladder(N_TOP), strict=True)]
    lo = np.empty(N_TOP + 1)
    hi = np.empty(N_TOP + 1)
    for k in range(N_TOP + 1):
        genome = Genome([planted_contig("rep1", peps[:k], random.Random(k))])
        res = make_case(tmp_path / str(k), rows, genome).evaluate()
        assert res.j is not None, res.stderr
        j = res.j
        assert j["n_top"] == N_TOP and j["n_match"] == k
        lo[k], hi[k] = j["share_ci95"]
        assert lo[k] <= k / N_TOP <= hi[k]
    ks = np.arange(N_TOP + 1)
    ps = np.linspace(0.05, 0.95, 181)
    cov = np.array([binom.pmf(ks, N_TOP, p)[(lo <= p) & (p <= hi)].sum() for p in ps])
    print(f"A3 share_ci95 n_top={N_TOP}: mean exact coverage {cov.mean():.4f}, min {cov.min():.4f}")
    (tmp_path / "coverage.json").write_text(json.dumps({"mean": float(cov.mean())}))
    assert abs(float(cov.mean()) - 0.95) <= 0.01
