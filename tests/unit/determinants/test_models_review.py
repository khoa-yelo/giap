"""M8-08 reviews (PR #166): lead decision 7 and the code reviewer's test gaps T1-T3.

- every univariate row is ``ci_unestimated`` (the TU-block bootstrap CI of a rank statistic
  under-covers under a TU-structured null); a NA value is also ``unestimated``;
- ``tu_folds`` repeats differ even when every TU size is distinct (randomised greedy);
- T1: the elastic net standardises features (a x1000 rescale changes nothing);
- T2: ``oof`` is the mean over repeats of the per-repeat out-of-fold predictions;
- T3: with ``tu=``, permuting the input rows changes nothing."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d
from determinants_models_data import flags_of, mod
from test_models import _example_6

from gippa.rng import for_step


def _planted(n: int = 120, seed: int = 4) -> tuple[pd.DataFrame, np.ndarray, list[str], Any, Any]:
    X, y, tu, blocks = d.planted(np.random.default_rng(seed), n=n)
    cfg = d.config(n_estimators=20)
    folds = np.asarray(mod().tu_folds(tu, cfg.outer_folds, cfg.repeats,
                                      for_step("determinants.cv", 1)))  # fmt: skip
    return X, y, list(tu), blocks, (cfg, folds)


def test_univariate_rows_are_ci_unestimated() -> None:
    f, _ = _example_6()
    full = mod().univariate(f, 50, for_step("determinants.bootstrap", 1))
    few = mod().univariate(f.iloc[:8], 50, for_step("determinants.bootstrap", 1))  # < 10 rows
    assert full["value"].notna().all() and few["value"].isna().all()
    uni = pd.concat([full, few])
    for _, r in uni.iterrows():
        fl = set(flags_of(r["flags"]))
        assert "ci_unestimated" in fl, r
        assert ("unestimated" in fl) == bool(pd.isna(r["value"])), r


def test_tu_fold_repeats_differ_with_distinct_sizes() -> None:
    tu = d.tu_ids(list(range(1, 31)))
    folds = np.asarray(mod().tu_folds(tu, 5, 5, for_step("determinants.cv", 3)))
    assert len({d.partition(folds[r]) for r in range(5)}) > 1, "all repeats identical"
    tu_arr = np.asarray(tu)
    for r in range(5):
        for t in np.unique(tu_arr):
            assert len(set(folds[r, tu_arr == t].tolist())) == 1


def test_t1_elastic_net_standardises_features() -> None:
    X, y, tu, blocks, (cfg, folds) = _planted()
    col = X.columns[0]
    big = X.copy()
    big[col] = big[col] * 1000.0
    a = mod().cross_validate(X, y, folds, blocks, cfg, 5, 1, tu=tu)["elastic_net"]
    b = mod().cross_validate(big, y, folds, blocks, cfg, 5, 1, tu=tu)["elastic_net"]
    np.testing.assert_allclose(np.asarray(b.oof), np.asarray(a.oof), rtol=1e-6, atol=1e-8)


def test_t2_oof_is_the_mean_over_repeats() -> None:
    X, y, tu, blocks, (cfg, folds) = _planted()
    assert folds.shape[0] >= 2 and len({d.partition(f) for f in folds}) > 1
    full = mod().cross_validate(X, y, folds, blocks, cfg, 5, 1, tu=tu)["elastic_net"]
    per = [
        np.asarray(
            mod()
            .cross_validate(X, y, folds[r : r + 1], blocks, cfg, 5, 1, tu=tu)["elastic_net"]
            .oof
        )
        for r in range(folds.shape[0])
    ]
    np.testing.assert_allclose(np.asarray(full.oof), np.mean(per, axis=0), rtol=1e-9, atol=1e-12)


def test_t3_row_order_with_tu() -> None:
    X, y, tu, blocks, (cfg, folds) = _planted()
    perm = np.random.default_rng(9).permutation(len(y))
    a = mod().cross_validate(X, y, folds, blocks, cfg, 5, 1, tu=tu)["elastic_net"]
    b = mod().cross_validate(X.iloc[perm], y[perm], folds[:, perm], blocks, cfg, 5, 1,
                             tu=[tu[i] for i in perm])["elastic_net"]  # fmt: skip
    np.testing.assert_allclose(b.r2, a.r2, rtol=0, atol=1e-10)
    np.testing.assert_allclose(np.asarray(b.oof), np.asarray(a.oof)[perm], atol=1e-10)
