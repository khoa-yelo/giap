"""Hypothesis properties of M8-08 (docs/specs/M8-08-determinants-models.md, acceptance list).

* **folds** -- ``tu_folds`` never splits a TU, uses every fold, keeps fold gene counts within the
  largest TU of each other (largest-first greedy), reproduces from the same seed, differs between
  seeds and between repeats, and does not depend on the input row order;
* **oracles** -- univariate rho equals ``scipy.stats.spearmanr`` and AUROC the Mann-Whitney pair
  count with ties ½, on random data with ties and NaNs;
* **row order** -- permuting the rows leaves every metric unchanged (univariate table, CV R²,
  ΔR², out-of-fold predictions per row).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d
from determinants_models_data import finite, is_na, mod, one

from gippa.rng import for_step

SETTINGS = settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)
SLOW = settings(SETTINGS, max_examples=3)


def _fold_map(tu: list[str], folds: np.ndarray) -> list[dict[str, int]]:
    out = []
    for r in range(folds.shape[0]):
        m: dict[str, int] = {}
        for t, k in zip(tu, folds[r].tolist(), strict=True):
            assert m.setdefault(t, k) == k, f"TU {t} split across folds in repeat {r}"
        out.append(m)
    return out


@SETTINGS
@given(
    sizes=st.lists(st.integers(1, 9), min_size=5, max_size=60),
    n_folds=st.integers(2, 5),
    repeats=st.integers(1, 5),
    seed=st.integers(0, 2**32 - 1),
)
def test_tu_folds_never_split_a_tu_and_balance(
    sizes: list[int], n_folds: int, repeats: int, seed: int
) -> None:
    tu = d.tu_ids(sizes)
    folds = np.asarray(mod().tu_folds(tu, n_folds, repeats, np.random.default_rng(seed)))
    assert folds.shape == (repeats, len(tu)) and folds.dtype == np.int64
    assert folds.min() >= 0 and folds.max() < n_folds
    _fold_map(tu, folds)
    if len(sizes) >= n_folds:
        for r in range(repeats):
            counts = np.bincount(folds[r], minlength=n_folds)
            assert (counts > 0).all()
            assert counts.max() - counts.min() <= max(sizes), counts


@SETTINGS
@given(
    sizes=st.lists(st.integers(1, 3), min_size=25, max_size=60),
    seed=st.integers(0, 2**32 - 2),
)
def test_tu_folds_seeded(sizes: list[int], seed: int) -> None:
    tu = d.tu_ids(sizes)
    f = mod().tu_folds
    a = np.asarray(f(tu, 5, 5, np.random.default_rng(seed)))
    b = np.asarray(f(tu, 5, 5, np.random.default_rng(seed)))
    c = np.asarray(f(tu, 5, 5, np.random.default_rng(seed + 1)))
    assert np.array_equal(a, b)
    assert d.partition(a[0]) != d.partition(c[0])
    assert len({d.partition(a[r]) for r in range(5)}) == 5


@SETTINGS
@given(
    sizes=st.lists(st.integers(1, 4), min_size=5, max_size=40),
    seed=st.integers(0, 2**32 - 1),
    perm_seed=st.integers(0, 2**32 - 1),
)
def test_tu_folds_ignore_input_order(sizes: list[int], seed: int, perm_seed: int) -> None:
    tu = d.tu_ids(sizes)
    perm = np.random.default_rng(perm_seed).permutation(len(tu))
    tu_p = [tu[i] for i in perm]
    a = np.asarray(mod().tu_folds(tu, 3, 2, np.random.default_rng(seed)))
    b = np.asarray(mod().tu_folds(tu_p, 3, 2, np.random.default_rng(seed)))
    assert _fold_map(tu, a) == _fold_map(tu_p, b)


def _auroc(score: np.ndarray, label: np.ndarray) -> float:
    pos, neg = score[label], score[~label]
    gt = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(gt / (pos.size * neg.size))


@st.composite
def worlds(draw: Any, n_min: int = 12, n_max: int = 50) -> dict[str, Any]:
    n = draw(st.integers(n_min, n_max))
    seed = draw(st.integers(0, 2**32 - 1))
    rng = np.random.default_rng(seed)
    x = rng.integers(0, draw(st.integers(2, 8)), n).astype(float)  # ties
    x[rng.random(n) < draw(st.sampled_from([0.0, 0.1]))] = math.nan
    det = rng.random(n) < 0.6
    det[:2] = [True, False]
    ab = np.where(det, rng.normal(size=n) + 0.3 * np.nan_to_num(x), math.nan)
    tu = d.tu_ids(d.mixed_sizes(rng, n, 1, 3))
    tags = [f"G{i:05d}" for i in range(n)]
    f = d.features({"gc3": x}, tu=tu, detected=det.tolist(), abundance=ab.tolist(), tags=tags)
    return {"f": f, "x": x, "det": det, "ab": ab, "seed": seed}


@SETTINGS
@given(worlds())
def test_univariate_matches_the_oracles(w: dict[str, Any]) -> None:
    uni = mod().univariate(w["f"], 20, np.random.default_rng(w["seed"]))
    x, det, ab = w["x"], w["det"], w["ab"]
    ok = np.isfinite(x)
    m = ok & det
    rho = one(uni, feature_or_block="gc3", metric="spearman_rho")
    if m.sum() >= 10 and np.unique(x[m]).size > 1:
        want = stats.spearmanr(x[m], ab[m]).statistic
        assert math.isclose(rho["value"], want, abs_tol=1e-12)
        assert int(rho["n"]) == int(m.sum())
    elif m.sum() >= 10:
        assert is_na(rho["value"])
    auc = one(uni, feature_or_block="gc3", metric="auroc")
    if ok.sum() >= 10 and det[ok].any() and (~det[ok]).any():
        assert math.isclose(auc["value"], _auroc(x[ok], det[ok]), abs_tol=1e-12)
        assert int(auc["n"]) == int(ok.sum())


@SETTINGS
@given(worlds(n_min=20), st.integers(0, 2**32 - 1))
def test_univariate_ignores_row_order(w: dict[str, Any], perm_seed: int) -> None:
    f = w["f"]
    perm = np.random.default_rng(perm_seed).permutation(len(f))
    fp = f.iloc[perm].reset_index(drop=True)
    a = mod().univariate(f, 30, np.random.default_rng(w["seed"]))
    b = mod().univariate(fp, 30, np.random.default_rng(w["seed"]))
    pd.testing.assert_frame_equal(d.canonical(a), d.canonical(b))


@SLOW
@given(st.integers(0, 2**32 - 1), st.integers(0, 2**32 - 1))
def test_cross_validate_ignores_row_order(seed: int, perm_seed: int) -> None:
    rng = np.random.default_rng(seed)
    X, y, tu, blocks = d.planted(rng, n=120)
    cfg = d.config(n_estimators=20)
    folds = np.asarray(
        mod().tu_folds(tu, cfg.outer_folds, cfg.repeats, for_step("determinants.cv", 1))
    )
    perm = np.random.default_rng(perm_seed).permutation(len(y))
    a = mod().cross_validate(X, y, folds, blocks, cfg, 5, 1)
    b = mod().cross_validate(X.iloc[perm], y[perm], folds[:, perm], blocks, cfg, 5, 1)
    for m in d.MODELS:
        np.testing.assert_allclose(b[m].r2, a[m].r2, rtol=0, atol=1e-10)
        np.testing.assert_allclose(np.asarray(b[m].oof), np.asarray(a[m].oof)[perm], atol=1e-10)
        assert set(a[m].delta_r2) == set(b[m].delta_r2)
        for k in a[m].delta_r2:
            np.testing.assert_allclose(b[m].delta_r2[k], a[m].delta_r2[k], rtol=0, atol=1e-10)
        assert all(finite(v) for v in a[m].r2)
