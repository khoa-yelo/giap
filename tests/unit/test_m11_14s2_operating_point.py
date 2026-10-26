"""The 14s2 operating-point driver (validation/operon_fallback_operating_point.py; ADR-0022
addendum) on a synthetic labelled set, and the recorded keys of ``calibration.yaml``.
Synthetic inputs only; no RegulonDB or pilot data."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.operon_fallback_operating_point import _point, operating_point

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


def test_recorded_flat_ci_key_names() -> None:
    """14d / 14d2 read these exact flat names; renaming them is a contract break."""
    od = load()["operon_distance"]
    for key in (
        "precision_at_pair_threshold_ci_low",
        "precision_at_pair_threshold_ci_high",
        "recall_at_pair_threshold_ci_low",
        "recall_at_pair_threshold_ci_high",
    ):
        assert key in od
        assert isinstance(od[key], float) and 0.0 < od[key] < 1.0


def test_point_on_a_tiny_example() -> None:
    score = np.array([0.9, 0.7, 0.4, 0.6, 0.1])
    y = np.array([1, 0, 1, 1, 0])
    res = _point(score, y, 0.5)
    # called: 0.9 (tp), 0.7 (fp), 0.6 (tp); p = 3
    assert (res["tp"], res["fp"], res["n_called"]) == (2, 1, 3)
    assert math.isclose(res["precision"], 2 / 3)
    assert math.isclose(res["recall"], 2 / 3)
    assert math.isclose(res["f1"], 2 * 2 / (2 + 1 + 3))


def test_point_threshold_is_inclusive() -> None:
    res = _point(np.array([0.5, 0.49]), np.array([1, 1]), 0.5)
    assert res["tp"] == 1 and res["n_called"] == 1


def test_point_nothing_called_gives_nan_precision() -> None:
    res = _point(np.array([0.1, 0.2]), np.array([1, 0]), 0.5)
    assert math.isnan(res["precision"])
    assert res["recall"] == 0.0 and res["f1"] == 0.0 and res["n_called"] == 0


def test_point_no_positives_is_nan_recall() -> None:
    res = _point(np.array([0.1, 0.2]), np.array([0, 0]), 0.5)
    assert math.isnan(res["recall"])
    assert math.isnan(res["f1"])


def test_point_no_positives_with_calls_is_nan_recall() -> None:
    res = _point(np.array([0.9, 0.8]), np.array([0, 0]), 0.5)
    assert res["precision"] == 0.0
    assert math.isnan(res["recall"])


def test_operating_point_no_distance_calls() -> None:
    L = _L(40)
    L["d"] = 5000  # p_dist ~ 0 everywhere: nothing called
    res = operating_point(L, 102.7, 47.2, 0.5, B=20)
    assert res["distance_model"]["n_called"] == 0
    assert math.isnan(res["distance_model"]["precision"])
    assert res["d_max_called_bp"] is None or math.isnan(res["d_max_called_bp"])


def test_operating_point_output_keys() -> None:
    res = operating_point(_L(80), 102.7, 47.2, 0.5, B=20)
    for name in ("distance_model", "uniop"):
        for key in ("precision_ci", "recall_ci", "precision_ci_iid", "recall_ci_iid"):
            lo, hi = res[name][key]
            assert lo <= hi
    for key in (
        "precision",
        "recall",
        "precision_ci",
        "recall_ci",
        "precision_ci_iid",
        "recall_ci_iid",
    ):
        assert key in res["uniop_minus_distance"]
    assert res["n_labelled"] == 80
    assert res["n_positive"] + res["n_negative"] == 80
    assert res["n_directons"] == 20
    assert (res["pair_threshold"], res["d0"], res["s"], res["B"]) == (0.5, 102.7, 47.2, 20)
    assert res["step"] == "m11_14s2_operon_fallback"
