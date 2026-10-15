"""M4-13 lead decision 7 after the statistician's review (docs/specs/M4-13-detectability.md):
the detectability diagnostics (``detection_summary`` rows, section ``detectability``; column
contract of M4-16 / GIAP_SPEC.md §9.3: ``value``, ``ci_low``, ``ci_high``) carry Wald 95 % CIs for
the recalibration slope and intercept.

Oracle: unpenalised MLE (a, b) of d ~ a + b·logit(clip(p)) (tests/helpers/detectability_oracle.py),
SE from the inverse observed Fisher information X'WX at the MLE (W = mu(1 - mu)); CI = estimate ±
z_{0.975}·SE. Rows other than slope / intercept carry no CI (NaN).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit, logit
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import detectability_data as dd
import detectability_oracle as O

from giap.config.model import DetectabilityConfig
from giap.rng import for_step

Z975 = float(norm.ppf(0.975))
CI_METRICS = ("calibration_intercept", "calibration_slope")


def _mod() -> Any:
    from giap.genetable import detectability

    return detectability


def _wald(d: np.ndarray, p: np.ndarray) -> dict[str, tuple[float, float, float]]:
    a, b = O.logistic_recalibration(d, p)
    z = logit(np.clip(p, O.CLIP, 1 - O.CLIP))
    x = np.column_stack([np.ones_like(z), z])
    mu = expit(a + b * z)
    cov = np.linalg.inv(x.T @ (x * (mu * (1 - mu))[:, None]))
    se = np.sqrt(np.diag(cov))
    return {
        "calibration_intercept": (a, a - Z975 * se[0], a + Z975 * se[0]),
        "calibration_slope": (b, b - Z975 * se[1], b + Z975 * se[1]),
    }


def _overall(diag: pd.DataFrame, metric: str) -> pd.Series:
    sel = diag[(diag["group_type"] == "overall") & (diag["metric"] == metric)]
    assert len(sel) == 1, metric
    return sel.iloc[0]


def _check(diag: pd.DataFrame, d: np.ndarray, p: np.ndarray) -> None:
    assert {"ci_low", "ci_high"} <= set(diag.columns)
    ref = _wald(d, p)
    for m in CI_METRICS:
        row = _overall(diag, m)
        est, lo, hi = ref[m]
        assert float(row["value"]) == pytest.approx(est, abs=1e-5)
        assert float(row["ci_low"]) == pytest.approx(lo, abs=1e-4), m
        assert float(row["ci_high"]) == pytest.approx(hi, abs=1e-4), m
        assert float(row["ci_low"]) < float(row["value"]) < float(row["ci_high"])
    other = diag[~diag["metric"].isin(CI_METRICS)]
    assert other["ci_low"].isna().all() and other["ci_high"].isna().all()


def test_calibration_diagnostics_wald_ci_example6_design() -> None:
    r = np.random.default_rng(6)
    p = expit(r.normal(0.0, 1.5, 20_000))
    d = (r.random(p.size) < p).astype(int)
    diag = _mod().calibration_diagnostics(d, p)
    _check(diag, d, p)
    # perfect calibration: the CIs are narrow and (for this seed) cover the truth (0, 1)
    s, i = _overall(diag, "calibration_slope"), _overall(diag, "calibration_intercept")
    assert float(s["ci_high"]) - float(s["ci_low"]) < 0.1
    assert float(s["ci_low"]) <= 1.0 <= float(s["ci_high"])
    assert float(i["ci_low"]) <= 0.0 <= float(i["ci_high"])


def test_calibration_diagnostics_wald_ci_miscalibrated_small_n() -> None:
    r = np.random.default_rng(11)
    p = expit(r.normal(0.0, 1.0, 300))
    d = (r.random(p.size) < expit(0.4 + 0.6 * logit(p))).astype(int)
    _check(_mod().calibration_diagnostics(d, p), d, p)


def test_fit_detectability_diagnostics_carry_wald_ci() -> None:
    gs = dd.synthetic(150, seed=23)
    res = _mod().fit_detectability(
        gs.detected, gs.cov, gs.tu_of_gene, DetectabilityConfig(), for_step("detectability", 7)
    )
    d = gs.detected.astype(bool).to_numpy().astype(int)
    _check(res.diagnostics, d, res.p_detect.to_numpy(float))
