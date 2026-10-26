"""The 14s2 operating-point driver (validation/operon_fallback_operating_point.py; ADR-0022
addendum) on a synthetic labelled set, and the recorded keys of ``calibration.yaml``.
Statistician-owned; no RegulonDB data."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.operon_fallback_operating_point import operating_point

from giap.data.calibration import load


def _L(n: int = 400) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    y = rng.random(n) < 0.6
    d = np.where(y, rng.integers(-20, 140, n), rng.integers(40, 400, n))
    u = np.where(y, rng.beta(8, 2, n), rng.beta(2, 8, n))
    return pd.DataFrame({"y": y.astype(int), "d": d, "uniop_p": u,
                         "directon_id": np.arange(n) // 4})  # fmt: skip


def test_point_estimates_match_the_rule_d_le_d0() -> None:
    L = _L()
    res = operating_point(L, 102.68667744326773, 47.16037265004493, 0.5, B=50)
    called = L["d"].to_numpy() <= 102
    y = L["y"].to_numpy() == 1
    dm = res["distance_model"]
    assert dm["tp"] == int((called & y).sum()) and dm["fp"] == int((called & ~y).sum())
    assert math.isclose(dm["recall"], (called & y).sum() / y.sum())
    assert res["d_max_called_bp"] <= 102
    for key in ("precision", "recall"):
        lo, hi = dm[f"{key}_ci"]
        assert lo <= dm[key] <= hi
        diff = res["uniop_minus_distance"]
        assert math.isclose(diff[key], res["uniop"][key] - dm[key])


def test_deterministic_for_a_seed() -> None:
    a = operating_point(_L(), 102.7, 47.2, 0.5, B=30, seed=5)
    b = operating_point(_L(), 102.7, 47.2, 0.5, B=30, seed=5)
    assert a == b


def test_recorded_operating_point_keys() -> None:
    od = load()["operon_distance"]
    for key in ("precision_at_pair_threshold", "recall_at_pair_threshold"):
        assert od[f"{key}_ci_low"] <= od[key] <= od[f"{key}_ci_high"]
        assert 0.0 < od[key] < 1.0
    assert od["precision_at_pair_threshold"] == 0.8284518828451883
    assert od["recall_at_pair_threshold"] == 0.915606936416185
