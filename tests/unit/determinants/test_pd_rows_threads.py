"""Partial-dependence support rows on a thread pool (perf; ADR-0063).

``pd_rows(..., threads)`` predicts booster by booster on ``threads`` threads (each predict
single-threaded, LightGBM releases the GIL) and collects in booster order: the rows are
bit-identical to ``threads`` = 1, including log-scale features and a feature with no value.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d

from giap.determinants.models import cross_validate
from giap.determinants.support import pd_rows

SEED = 20261026


@pytest.fixture(scope="module")
def fitted() -> tuple[object, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(11)
    X, y, tu, blocks = d.planted(rng, n=150)
    cfg = d.config(fast=True, n_estimators=15)
    cv = cross_validate(X, y, d.folds_for(tu, cfg, SEED), blocks, cfg, SEED, 1, tu=tu)
    feats = X.reset_index()
    feats["tai"] = np.abs(rng.normal(0.4, 0.1, len(feats)))  # absent from the model columns
    feats["tir_ostir"] = np.nan  # no value: no rows
    feats["promotech_max"] = np.exp(X["promotech_max"].to_numpy())  # positive, own scale
    return cv["lightgbm"], X, feats


@pytest.mark.parametrize("threads", [2, 4])
def test_pd_rows_identical_across_threads(fitted: tuple[object, pd.DataFrame, pd.DataFrame],
                                          threads: int) -> None:  # fmt: skip
    cv, X, feats = fitted
    one = pd_rows(cv, X, feats)  # type: ignore[arg-type]
    many = pd_rows(cv, X, feats, threads)  # type: ignore[arg-type]
    assert len(cv.boosters) > 1  # type: ignore[attr-defined]
    assert not one.empty and set(one["feature"]) == {"cai", "tai", "promotech_max"}
    pd.testing.assert_frame_equal(one, many, check_exact=True)
    for c in ("grid_value", "pd_mean", "pd_low", "pd_high"):
        assert one[c].to_numpy(float).tobytes() == many[c].to_numpy(float).tobytes(), c
