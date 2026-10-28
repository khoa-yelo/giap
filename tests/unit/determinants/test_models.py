"""M8-08 determinants models: examples 1-8 and every degenerate row of
docs/specs/M8-08-determinants-models.md (GIPPA_SPEC.md §7.8.2 steps 1-3 and 5, §8.2, §9.3).

Builders and the contract choices the spec leaves open are in
``tests/helpers/determinants_models_data.py``. Behavioural CV tests use a fast config (2 repeats,
60 trees); examples 1 and 2 use the §5.2 defaults.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d
from determinants_models_data import finite, flags_of, is_na, mod, one, rows

from gippa.rng import for_step

SEED = 20261014


def _uni(features: pd.DataFrame, B: int = 200, seed: int = SEED) -> pd.DataFrame:
    return mod().univariate(features, B, for_step("determinants.bootstrap", seed))


def _cv(
    X: pd.DataFrame,
    y: np.ndarray,
    tu: list[str],
    blocks: dict[str, str],
    cfg: Any,
    seed: int = SEED,
    threads: int = 1,
    folds: np.ndarray | None = None,
) -> dict[str, Any]:
    f = d.folds_for(tu, cfg) if folds is None else folds
    return dict(mod().cross_validate(X, y, f, blocks, cfg, seed, threads))


def _auroc(score: np.ndarray, label: np.ndarray) -> float:
    pos, neg = score[label], score[~label]
    gt = (pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()
    return float(gt / (pos.size * neg.size))


# ---------------------------------------------------------------------------------------------
# rule 1: model rows


def test_model_rows_rule() -> None:
    """Eligible (every TSV row), detected, finite non-imputed abundance, not ribosomal."""
    tu = d.tu_ids([1] * 7)
    f = d.features(
        {"gravy": [0.1] * 7},
        tu=tu,
        detected=[True, False, True, True, True, True, True],
        abundance=[1.0, 1.0, math.nan, 2.0, 3.0, math.inf, 4.0],
        imputed=[False, False, False, True, False, False, False],
        ribosomal=[False, False, False, False, True, False, False],
    )
    got = np.asarray(mod().model_rows(f))
    assert got.dtype == bool and got.shape == (7,)
    assert got.tolist() == [True, False, False, False, False, False, True]


# ---------------------------------------------------------------------------------------------
# rule 2: univariate rows (examples 4-6)


def test_example_4_spearman_by_hand_and_scipy() -> None:
    """x = 1..12, y swaps each adjacent pair: Σd² = 12, rho = 1 - 6·12/(12·143) = 0.958042.

    (The spec's 6-row example is doubled: univariate rows need ≥ 10 rows.)"""
    x = np.arange(1.0, 13.0)
    y = np.array([2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11], dtype=float)
    f = d.features({"gravy": x}, tu=d.tu_ids([1] * 12), abundance=y)
    r = one(_uni(f), model="univariate", feature_or_block="gravy", metric="spearman_rho")
    assert r["value"] == pytest.approx(1 - 72 / 1716, abs=1e-12)
    assert r["value"] == pytest.approx(stats.spearmanr(x, y).statistic, abs=1e-12)
    assert int(r["n"]) == 12


def test_spearman_with_ties_matches_scipy() -> None:
    rng = np.random.default_rng(4)
    x = rng.integers(0, 4, 40).astype(float)
    y = rng.normal(size=40) + x
    f = d.features({"gc3": x}, tu=d.tu_ids(d.mixed_sizes(rng, 40)), abundance=y)
    r = one(_uni(f), feature_or_block="gc3", metric="spearman_rho")
    assert r["value"] == pytest.approx(stats.spearmanr(x, y).statistic, abs=1e-12)


def test_example_5_auroc_ties_count_half() -> None:
    """Detected scores (3, 2, 1), undetected (2, 0), each twice: 18/24 = 0.75 (spec: 4.5/6)."""
    score = np.array([3, 2, 1, 3, 2, 1, 2, 0, 2, 0], dtype=float)
    det = np.array([True] * 6 + [False] * 4)
    f = d.features(
        {"gravy": score},
        tu=d.tu_ids([1] * 10),
        detected=det.tolist(),
        abundance=[1.0] * 6 + [math.nan] * 4,
    )
    r = one(_uni(f), feature_or_block="gravy", metric="auroc")
    assert r["value"] == pytest.approx(0.75, abs=1e-12)
    assert _auroc(score, det) == pytest.approx(0.75)
    assert int(r["n"]) == 10


def _example_6() -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    rng = np.random.default_rng(6)
    n, n_rib = 60, 8
    rib = np.zeros(n, dtype=bool)
    rib[:n_rib] = True
    det = rng.random(n) < 0.7
    det[:n_rib] = True
    cai = rng.uniform(0.2, 0.6, n)
    cai[:n_rib] = 0.95
    gravy = rng.normal(size=n)
    ab = rng.normal(size=n)
    ab[:n_rib] = 12.0
    ab[~det] = math.nan
    f = d.features(
        {"cai": cai, "gravy": gravy},
        tu=d.tu_ids(d.mixed_sizes(rng, n)),
        detected=det.tolist(),
        abundance=ab.tolist(),
        ribosomal=rib.tolist(),
    )
    return f, {"rib": rib, "det": det, "cai": cai, "gravy": gravy, "ab": ab}


def test_example_6_ribosomal_genes_excluded_from_cai_rows() -> None:
    f, v = _example_6()
    uni = _uni(f)
    model = v["det"] & ~v["rib"]
    rho = one(uni, feature_or_block="cai", metric="spearman_rho")
    assert rho["value"] == pytest.approx(
        stats.spearmanr(v["cai"][model], v["ab"][model]).statistic, abs=1e-12
    )
    assert int(rho["n"]) == int(model.sum())
    auc = one(uni, feature_or_block="cai", metric="auroc")
    assert auc["value"] == pytest.approx(_auroc(v["cai"][~v["rib"]], v["det"][~v["rib"]]))
    assert int(auc["n"]) == int((~v["rib"]).sum())


def test_example_6_ribosomal_genes_included_in_gravy_auroc() -> None:
    f, v = _example_6()
    uni = _uni(f)
    auc = one(uni, feature_or_block="gravy", metric="auroc")
    assert auc["value"] == pytest.approx(_auroc(v["gravy"], v["det"]))
    assert int(auc["n"]) == len(f)
    rho = one(uni, feature_or_block="gravy", metric="spearman_rho")  # model rows: no ribosomal
    model = v["det"] & ~v["rib"]
    assert int(rho["n"]) == int(model.sum())


def test_univariate_rows_have_the_table_shape() -> None:
    f, _ = _example_6()
    uni = _uni(f)
    assert list(uni.columns) == list(d.TABLE_COLUMNS)
    assert set(uni["model"]) == {"univariate"}
    assert set(uni["metric"]) == {"spearman_rho", "auroc"}
    assert set(uni["feature_or_block"]) == {"cai", "gravy"}  # model_matrix drops all-NA features
    assert len(uni) == 4
    assert (uni["n_perm"].astype(int) == 0).all()
    for _, r in uni.iterrows():
        assert finite(r["ci_low"]) and finite(r["ci_high"])
        assert r["ci_low"] <= r["ci_high"]
        lo, hi = (-1.0, 1.0) if r["metric"] == "spearman_rho" else (0.0, 1.0)
        assert lo <= r["ci_low"] and r["ci_high"] <= hi


def test_univariate_ci_brackets_a_strong_effect() -> None:
    rng = np.random.default_rng(11)
    n = 80
    x = rng.normal(size=n)
    f = d.features(
        {"tai": x}, tu=d.tu_ids(d.mixed_sizes(rng, n)), abundance=3 * x + rng.normal(size=n)
    )
    r = one(_uni(f, B=500), feature_or_block="tai", metric="spearman_rho")
    assert r["ci_low"] <= r["value"] <= r["ci_high"]
    assert r["ci_low"] > 0.5


def test_univariate_one_hot_columns_get_rows() -> None:
    """Rule 2: per *encoded* column: start_codon -> start_codon_ATG ... start_codon_other."""
    rng = np.random.default_rng(12)
    n = 40
    sc = rng.choice(["ATG", "GTG", "TTG", "CTG"], n).tolist()
    f = d.features({"start_codon": sc}, tu=d.tu_ids([1] * n), abundance=rng.normal(size=n))
    names = set(_uni(f)["feature_or_block"])
    assert names == {f"start_codon_{lv}" for lv in ("ATG", "GTG", "TTG", "other")}


def test_univariate_uses_the_tu_block_bootstrap(monkeypatch: pytest.MonkeyPatch) -> None:
    """CIs are TU-block bootstrap percentile CIs with B = the argument (§8.2 resample TUs)."""
    import gippa.stats.bootstrap as boot

    real = boot.tu_block_bootstrap
    calls: list[dict[str, Any]] = []

    def spy(
        df: pd.DataFrame, tu_col: str, stat_fn: Any, B: int, rng: Any, *a: Any, **k: Any
    ) -> Any:
        calls.append({"tu": df[tu_col].tolist(), "B": B})
        return real(df, tu_col, stat_fn, B, rng, *a, **k)

    m = mod()
    monkeypatch.setattr(boot, "tu_block_bootstrap", spy)
    if hasattr(m, "tu_block_bootstrap"):
        monkeypatch.setattr(m, "tu_block_bootstrap", spy)
    f, _ = _example_6()
    _uni(f, B=37)
    assert calls, "univariate never called tu_block_bootstrap"
    assert {c["B"] for c in calls} == {37}
    tus = set(f["tu_id"])
    assert all(set(c["tu"]) <= tus for c in calls)


def test_univariate_drops_na_feature_rows_per_feature() -> None:
    rng = np.random.default_rng(13)
    n = 30
    x = rng.normal(size=n)
    x[:5] = math.nan
    ab = rng.normal(size=n)
    f = d.features({"mrna_dg": x}, tu=d.tu_ids([1] * n), abundance=ab)
    r = one(_uni(f), feature_or_block="mrna_dg", metric="spearman_rho")
    assert int(r["n"]) == n - 5
    assert r["value"] == pytest.approx(stats.spearmanr(x[5:], ab[5:]).statistic, abs=1e-12)


def test_degenerate_constant_feature_rho_na_auroc_half() -> None:
    n = 20
    det = [True] * 14 + [False] * 6
    f = d.features(
        {"fop": [0.4] * n, "gravy": np.linspace(-1, 1, n)},
        tu=d.tu_ids([1] * n),
        detected=det,
        abundance=list(np.linspace(0, 5, 14)) + [math.nan] * 6,
    )
    uni = _uni(f)
    assert is_na(one(uni, feature_or_block="fop", metric="spearman_rho")["value"])
    assert one(uni, feature_or_block="fop", metric="auroc")["value"] == pytest.approx(0.5)


@pytest.mark.parametrize("all_detected", [True, False])
def test_degenerate_one_class_auroc_na(all_detected: bool) -> None:
    n = 20
    ab = list(np.linspace(0, 5, n)) if all_detected else [math.nan] * n
    f = d.features(
        {"gravy": np.linspace(-1, 1, n)},
        tu=d.tu_ids([1] * n),
        detected=[all_detected] * n,
        abundance=ab,
    )
    uni = _uni(f)
    auc = rows(uni, feature_or_block="gravy", metric="auroc")
    assert not any(finite(v) for v in auc["value"])
    if all_detected:
        assert finite(one(uni, feature_or_block="gravy", metric="spearman_rho")["value"])


def test_fewer_than_ten_model_rows_give_no_rho() -> None:
    n = 30
    det = [True] * 8 + [False] * 22
    f = d.features(
        {"gravy": np.linspace(-1, 1, n)},
        tu=d.tu_ids([1] * n),
        detected=det,
        abundance=list(np.linspace(0, 5, 8)) + [math.nan] * 22,
    )
    uni = _uni(f)
    assert not any(finite(v) for v in rows(uni, metric="spearman_rho")["value"])
    assert finite(one(uni, feature_or_block="gravy", metric="auroc")["value"])  # 30 AUROC rows


def test_univariate_is_seeded_and_leaves_global_numpy_state_alone() -> None:
    f, _ = _example_6()
    state = np.random.get_state()  # noqa: NPY002
    a, b = _uni(f, seed=3), _uni(f, seed=3)
    after = np.random.get_state()  # noqa: NPY002
    assert state[0] == after[0] and np.array_equal(state[1], after[1])
    pd.testing.assert_frame_equal(a, b)
    c = _uni(f, seed=4)
    assert not np.allclose(a["ci_low"].to_numpy(float), c["ci_low"].to_numpy(float))
    pd.testing.assert_series_equal(a["value"], c["value"])  # estimates do not depend on the seed


# ---------------------------------------------------------------------------------------------
# rule 3: folds (example 3)


def _example_3() -> list[str]:
    rng = np.random.default_rng(3)
    sizes = rng.choice([1, 1, 1, 2, 2, 3, 4, 6, 8], 80).tolist()
    return d.tu_ids(sizes)


def test_example_3_every_tu_in_one_fold_and_balanced() -> None:
    tu = _example_3()
    folds = np.asarray(mod().tu_folds(tu, 5, 5, for_step("determinants.cv", SEED)))
    assert folds.shape == (5, len(tu)) and folds.dtype == np.int64
    tu_arr = np.asarray(tu)
    for r in range(5):
        assert set(folds[r].tolist()) == set(range(5))
        for t in np.unique(tu_arr):
            assert len(set(folds[r, tu_arr == t].tolist())) == 1, f"TU {t} split in repeat {r}"
        counts = np.bincount(folds[r], minlength=5)
        assert np.all(np.abs(counts - counts.mean()) <= 0.15 * counts.mean()), counts
    assert len({d.partition(folds[r]) for r in range(5)}) == 5, "repeats must differ"


def test_tu_folds_reproducible_from_the_cv_stream() -> None:
    tu = _example_3()
    a = mod().tu_folds(tu, 5, 5, for_step("determinants.cv", SEED))
    b = mod().tu_folds(tu, 5, 5, for_step("determinants.cv", SEED))
    c = mod().tu_folds(tu, 5, 5, for_step("determinants.cv", SEED + 1))
    assert np.array_equal(a, b)
    assert d.partition(np.asarray(a)[0]) != d.partition(np.asarray(c)[0])


# ---------------------------------------------------------------------------------------------
# rules 4-7: cross_validate (examples 1, 2, 8)


def _example_1() -> dict[str, Any]:
    X, y, tu, blocks = d.planted(np.random.default_rng(1), n=400, signal=2.0)
    cfg = d.config(fast=False)
    cv = _cv(X, y, tu, blocks, cfg)
    return {"cv": cv, "X": X, "y": y, "tu": tu, "blocks": blocks, "cfg": cfg}


def _check_example_1_result_shapes(example_1: dict[str, Any]) -> None:
    cv, cfg, n = example_1["cv"], example_1["cfg"], len(example_1["y"])
    assert set(cv) == set(d.MODELS)
    for name, res in cv.items():
        assert res.model == name
        assert np.asarray(res.r2).shape == (cfg.repeats,)
        assert np.asarray(res.oof).shape == (n,)
        assert np.asarray(res.folds).shape == (cfg.repeats, n)
        assert set(res.delta_r2) == set(example_1["blocks"].values())
        assert all(np.asarray(v).shape == (cfg.repeats,) for v in res.delta_r2.values())
    boosters = cv["lightgbm"].boosters
    assert len(boosters) == cfg.repeats * cfg.outer_folds
    assert all(hasattr(b, "predict") for b in boosters)


def _check_example_1_cv_r2(example_1: dict[str, Any]) -> None:
    """Elastic net recovers the true R² 0.8 within 0.05 (spec). LightGBM's step-function fit of
    a linear signal on 320 training rows is looser: within 0.15 (test-author tolerance, see the
    M8-08 report: the spec's 0.05 is not met by a correct 500-tree fit)."""
    r2 = {m: float(np.mean(example_1["cv"][m].r2)) for m in d.MODELS}
    assert r2["elastic_net"] == pytest.approx(0.8, abs=0.05)
    assert r2["lightgbm"] == pytest.approx(0.8, abs=0.15)


def _check_example_1_delta_r2(example_1: dict[str, Any]) -> None:
    """ΔR²(codon) ≈ 0.8 (elastic net; LightGBM's drop model overfits noise, so ΔR² is larger);
    every noise block |ΔR²| ≤ 0.05."""
    en, lg = example_1["cv"]["elastic_net"].delta_r2, example_1["cv"]["lightgbm"].delta_r2
    assert float(np.mean(en["codon"])) == pytest.approx(0.8, abs=0.05)
    assert float(np.mean(lg["codon"])) > 0.7
    for res in (en, lg):
        for block in ("transcription", "translation", "protein", "context"):
            assert abs(float(np.mean(res[block]))) <= 0.05, block


def _check_example_1_pooled_oof_r2_and_table(example_1: dict[str, Any]) -> None:
    cv, cfg = example_1["cv"], example_1["cfg"]
    tbl = mod().model_rows_table(cv, _uni_empty())
    t = stats.t.ppf(0.975, cfg.repeats - 1)
    assert t == pytest.approx(2.776, abs=1e-3)
    for m in d.MODELS:
        res = cv[m]
        r2 = np.asarray(res.r2, dtype=float)
        mean = one(tbl, model=m, metric="cv_r2_mean")
        assert mean["value"] == pytest.approx(r2.mean(), abs=1e-12)
        assert int(mean["n"]) == len(example_1["y"])
        sd = one(tbl, model=m, metric="cv_r2_sd")
        assert sd["value"] == pytest.approx(r2.std(ddof=1), abs=1e-12)
        for block, v in res.delta_r2.items():
            v = np.asarray(v, dtype=float)
            r = one(tbl, model=m, metric="delta_r2", feature_or_block=block)
            half = t * v.std(ddof=1) / math.sqrt(cfg.repeats)
            assert r["value"] == pytest.approx(v.mean(), abs=1e-12)
            assert r["ci_low"] == pytest.approx(v.mean() - half, abs=1e-12)
            assert r["ci_high"] == pytest.approx(v.mean() + half, abs=1e-12)


def test_example_1() -> None:
    """Spec example 1 with the §5.2 defaults (one fit shared by four checks)."""
    ex = _example_1()
    _check_example_1_result_shapes(ex)
    _check_example_1_cv_r2(ex)
    _check_example_1_delta_r2(ex)
    _check_example_1_pooled_oof_r2_and_table(ex)


def _uni_empty() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([], dtype=object) for c in d.TABLE_COLUMNS})


@pytest.mark.slow
def test_example_2_same_seed_and_threads_give_identical_tables() -> None:
    X, y, tu, blocks = d.planted(np.random.default_rng(2), n=300)
    cfg = d.config(fast=False)
    folds = d.folds_for(tu, cfg)
    t1 = mod().model_rows_table(_cv(X, y, tu, blocks, cfg, threads=1, folds=folds), _uni_empty())
    t1b = mod().model_rows_table(_cv(X, y, tu, blocks, cfg, threads=1, folds=folds), _uni_empty())
    t4 = mod().model_rows_table(_cv(X, y, tu, blocks, cfg, threads=4, folds=folds), _uni_empty())
    pd.testing.assert_frame_equal(t1, t1b, check_exact=True)
    pd.testing.assert_frame_equal(t1, t4, check_exact=True)


@pytest.fixture(scope="module")
def small() -> dict[str, Any]:
    X, y, tu, blocks = d.planted(np.random.default_rng(21), n=200)
    return {"X": X, "y": y, "tu": tu, "blocks": blocks, "cfg": d.config()}


def test_threads_do_not_change_results_fast(small: dict[str, Any]) -> None:
    s = small
    folds = d.folds_for(s["tu"], s["cfg"])
    a = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], threads=1, folds=folds)
    b = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], threads=3, folds=folds)
    for m in d.MODELS:
        assert np.array_equal(a[m].oof, b[m].oof) and np.array_equal(a[m].r2, b[m].r2)


def test_seed_drives_lightgbm_only_and_global_state_is_untouched(small: dict[str, Any]) -> None:
    s = small
    folds = d.folds_for(s["tu"], s["cfg"])
    state = np.random.get_state()  # noqa: NPY002
    a = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], seed=1, folds=folds)
    a2 = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], seed=1, folds=folds)
    b = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], seed=2, folds=folds)
    after = np.random.get_state()  # noqa: NPY002
    assert np.array_equal(state[1], after[1])
    for m in d.MODELS:
        assert np.array_equal(a[m].oof, a2[m].oof)
        assert np.array_equal(a[m].folds, folds)
    assert not np.array_equal(a["lightgbm"].oof, b["lightgbm"].oof)  # subsample 0.8 is seeded
    np.testing.assert_allclose(a["elastic_net"].oof, b["elastic_net"].oof, rtol=0, atol=1e-12)


def test_cv_r2_is_the_pooled_oof_r2(small: dict[str, Any]) -> None:
    """Decision 15: per repeat 1 - SSE/SST over all rows. With the same folds in both repeats the
    (deterministic) elastic net gives identical repeats, so oof is the per-repeat prediction."""
    s = small
    f = d.folds_for(s["tu"], s["cfg"])
    same = np.vstack([f[0], f[0]])
    res = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], folds=same)["elastic_net"]
    y, oof = s["y"], np.asarray(res.oof)
    pooled = 1 - np.sum((y - oof) ** 2) / np.sum((y - y.mean()) ** 2)
    np.testing.assert_allclose(res.r2, [pooled, pooled], rtol=0, atol=1e-12)
    fold_mean = np.mean(
        [
            1
            - np.sum((y[f[0] == k] - oof[f[0] == k]) ** 2)
            / np.sum((y[f[0] == k] - y[f[0] == k].mean()) ** 2)
            for k in range(5)
        ]
    )
    assert abs(fold_mean - pooled) > 1e-9  # the two definitions differ on these data


def test_oof_is_out_of_fold(small: dict[str, Any]) -> None:
    """A row's prediction never sees its own response: a huge y on one test row does not pull
    that row's OOF prediction (it would under leakage)."""
    s = small
    folds = d.folds_for(s["tu"], s["cfg"])
    y = s["y"].copy()
    base = _cv(s["X"], y, s["tu"], s["blocks"], s["cfg"], folds=folds)
    y[0] = 1e4
    hot = _cv(s["X"], y, s["tu"], s["blocks"], s["cfg"], folds=folds)
    for m in d.MODELS:
        assert abs(float(hot[m].oof[0]) - float(base[m].oof[0])) < 1.0, m


def test_example_8_imputer_uses_the_training_median() -> None:
    """``cai`` NA on half of fold 0's test rows, the other half at +10: the training median is
    ≈ 0 (two clusters ±[1, 2]), the all-row median ≈ 1.1. y = 2·cai, so the NA rows' elastic-net
    prediction is ≈ 2·median_train (no leak) and far from 2·median_all (leak)."""
    rng = np.random.default_rng(8)
    n = 200
    tu = d.tu_ids([2] * (n // 2))
    tu_arr = np.asarray(tu)
    fold = np.repeat(np.arange(5), n // 5)  # TUs are contiguous pairs: no TU split
    folds = np.vstack([fold, fold])
    cai = np.where(rng.random(n) < 0.5, rng.uniform(-2, -1, n), rng.uniform(1, 2, n))
    test = fold == 0
    na_rows = np.flatnonzero(test)[::2]
    hot_rows = np.flatnonzero(test)[1::2]
    cai[hot_rows] = 10.0
    y = 2 * cai + rng.normal(0, 0.01, n)
    cai[na_rows] = math.nan
    X = pd.DataFrame(
        {"cai": cai, "gravy": rng.normal(size=n)}, index=pd.Index([f"G{i:05d}" for i in range(n)])
    )
    assert len({tuple(np.unique(fold[tu_arr == t])) for t in np.unique(tu_arr)}) <= 5
    cfg = d.config()
    res = mod().cross_validate(X, y, folds, {"cai": "codon", "gravy": "protein"}, cfg, SEED, 1)
    pred = np.asarray(res["elastic_net"].oof)[na_rows]
    med_train = float(np.nanmedian(cai[~test]))
    med_all = float(np.nanmedian(cai))
    assert abs(med_all - med_train) > 0.5
    np.testing.assert_allclose(pred, 2 * med_train, atol=0.25)
    assert np.all(np.abs(pred - 2 * med_train) < np.abs(pred - 2 * med_all))


def test_lightgbm_passes_nan_through(small: dict[str, Any]) -> None:
    s = small
    X = s["X"].copy()
    X.iloc[::7, 0] = math.nan
    res = _cv(X, s["y"], s["tu"], s["blocks"], s["cfg"])
    for m in d.MODELS:
        assert np.all(np.isfinite(res[m].oof)) and np.all(np.isfinite(res[m].r2))


# ---------------------------------------------------------------------------------------------
# lead decision 2: zero-variance and all-NaN columns dropped per training fold


def test_constant_and_all_nan_columns_change_nothing(small: dict[str, Any]) -> None:
    """Degenerate row "a constant feature": dropped before fitting, so both models give exactly
    the results without it (LightGBM's colsample would otherwise see a different column count)."""
    s = small
    folds = d.folds_for(s["tu"], s["cfg"])
    base = _cv(s["X"], s["y"], s["tu"], s["blocks"], s["cfg"], folds=folds)
    X = s["X"].copy()
    X["start_codon_other"] = 0.0
    X["nterm_residue_W"] = math.nan
    blocks = {**s["blocks"], "start_codon_other": "translation", "nterm_residue_W": "protein"}
    got = _cv(X, s["y"], s["tu"], blocks, s["cfg"], folds=folds)
    for m in d.MODELS:
        np.testing.assert_allclose(got[m].oof, base[m].oof, rtol=0, atol=1e-12)
        np.testing.assert_allclose(got[m].r2, base[m].r2, rtol=0, atol=1e-12)
        for b in base[m].delta_r2:
            np.testing.assert_allclose(got[m].delta_r2[b], base[m].delta_r2[b], rtol=0, atol=1e-12)
        assert {"start_codon_other", "nterm_residue_W"} <= set(got[m].dropped)
        assert not set(got[m].dropped) & set(s["X"].columns)
        assert list(got[m].dropped) == sorted(got[m].dropped)


def test_zero_variance_is_judged_on_the_training_fold_only(small: dict[str, Any]) -> None:
    """``varies_in_test``: constant on every row outside fold 0 (so constant on fold 0's
    training rows) -> dropped. ``const_in_test``: constant on fold 0's test rows only, varying in
    every training fold -> kept (judging it on a test fold would leak)."""
    s = small
    folds = d.folds_for(s["tu"], s["cfg"])
    rng = np.random.default_rng(22)
    test0 = folds[0] == 0
    varies = np.where(test0, rng.normal(size=test0.size), 1.0)
    const = np.where(test0, 5.0, rng.normal(size=test0.size))
    X = s["X"].assign(varies_in_test=varies, const_in_test=const)
    blocks = {**s["blocks"], "varies_in_test": "context", "const_in_test": "context"}
    got = _cv(X, s["y"], s["tu"], blocks, s["cfg"], folds=folds)
    for m in d.MODELS:
        assert "varies_in_test" in got[m].dropped
        assert "const_in_test" not in got[m].dropped


def test_block_with_no_columns_has_no_delta_row(small: dict[str, Any]) -> None:
    """Degenerate row "a block with every column dropped (layer missing)": model_matrix drops the
    layer's columns, so the block is absent from ``blocks`` and gets no ΔR²."""
    s = small
    X = s["X"].drop(columns=["promotech_max"])
    blocks = {c: b for c, b in s["blocks"].items() if c != "promotech_max"}
    cv = _cv(X, s["y"], s["tu"], blocks, s["cfg"])
    tbl = mod().model_rows_table(cv, _uni_empty())
    for m in d.MODELS:
        assert "transcription" not in cv[m].delta_r2
    assert set(rows(tbl, metric="delta_r2")["feature_or_block"]) == {
        "codon", "translation", "protein", "context",
    }  # fmt: skip


def test_degenerate_constant_y_is_unestimated(small: dict[str, Any]) -> None:
    s = small
    y = np.full(len(s["y"]), 3.0)
    cv = _cv(s["X"], y, s["tu"], s["blocks"], s["cfg"])
    for m in d.MODELS:
        assert not np.any(np.isfinite(np.asarray(cv[m].r2, dtype=float)))
    tbl = mod().model_rows_table(cv, _uni_empty())
    for m in d.MODELS:
        for metric in ("cv_r2_mean", "cv_r2_sd"):
            r = one(tbl, model=m, metric=metric)
            assert is_na(r["value"])
            assert "unestimated" in flags_of(r["flags"])
            assert "too_few_rows" not in flags_of(r["flags"])


# ---------------------------------------------------------------------------------------------
# model_rows_table on hand-made results (§9.3 columns, ΔR² t-interval)


def _fake(model: str, r2: list[float], delta: dict[str, list[float]], n: int = 150) -> Any:
    R = len(r2)
    return mod().CvResult(
        model=model,
        r2=np.asarray(r2, dtype=float),
        oof=np.zeros(n),
        folds=np.zeros((R, n), dtype=np.int64),
        delta_r2={b: np.asarray(v, dtype=float) for b, v in delta.items()},
        boosters=(),
    )


def test_table_columns_values_and_t_interval() -> None:
    r2 = [0.30, 0.32, 0.28, 0.31, 0.29]
    dl = {"codon": [0.10, 0.12, 0.08, 0.11, 0.09], "protein": [0.0, 0.01, -0.01, 0.0, 0.0]}
    cv = {"lightgbm": _fake("lightgbm", r2, dl), "elastic_net": _fake("elastic_net", r2, dl)}
    f, _ = _example_6()
    uni = _uni(f)
    tbl = mod().model_rows_table(cv, uni)
    assert list(tbl.columns) == list(d.TABLE_COLUMNS)
    assert set(tbl["model"]) == {"lightgbm", "elastic_net", "univariate"}
    assert set(tbl["metric"]) <= d.METRICS
    assert len(rows(tbl, model="univariate")) == len(uni)
    pd.testing.assert_frame_equal(
        d.canonical(rows(tbl, model="univariate")), d.canonical(uni), check_dtype=False
    )
    for m in d.MODELS:
        assert len(rows(tbl, model=m, metric="cv_r2_mean")) == 1
        assert len(rows(tbl, model=m, metric="cv_r2_sd")) == 1
        assert len(rows(tbl, model=m, metric="delta_r2")) == 2
        r = one(tbl, model=m, metric="cv_r2_mean")
        assert r["value"] == pytest.approx(0.30, abs=1e-12)
        assert int(r["n"]) == 150 and int(r["n_perm"]) == 0
        assert one(tbl, model=m, metric="cv_r2_sd")["value"] == pytest.approx(
            float(np.std(r2, ddof=1)), abs=1e-12
        )
        c = one(tbl, model=m, metric="delta_r2", feature_or_block="codon")
        half = 2.7764451051977987 * float(np.std(dl["codon"], ddof=1)) / math.sqrt(5)
        assert c["value"] == pytest.approx(0.10, abs=1e-12)
        assert (c["ci_low"], c["ci_high"]) == (
            pytest.approx(0.10 - half, abs=1e-12),
            pytest.approx(0.10 + half, abs=1e-12),
        )
        for _, row in rows(tbl, model=m).iterrows():
            # lead decision 5: ΔR² rows carry ci_unestimated (the t interval is not calibrated)
            want = {"ci_unestimated"} if row["metric"] == "delta_r2" else set()
            assert set(flags_of(row["flags"])) == want


# ---------------------------------------------------------------------------------------------
# rule 8: minimum data (example 7, decision 18)


@pytest.mark.parametrize(
    ("n_rows", "n_tus", "ok"),
    [(100, 25, True), (99, 25, False), (100, 24, False), (60, 30, False), (400, 80, True)],
)
def test_minimum_data_rule(n_rows: int, n_tus: int, ok: bool) -> None:
    m = mod()
    assert (m.MIN_ROWS, m.MIN_TUS) == (100, 25)
    sizes = [1] * (n_tus - 1) + [n_rows - (n_tus - 1)]
    assert m.enough_data(d.tu_ids(sizes)) is ok


def test_example_7_sixty_model_rows() -> None:
    rng = np.random.default_rng(7)
    n = 60
    x = rng.normal(size=n)
    f = d.features(
        {"cai": x, "gravy": rng.normal(size=n)},
        tu=d.tu_ids(d.mixed_sizes(rng, n)),
        abundance=x + rng.normal(size=n),
    )
    tu = f.loc[np.asarray(mod().model_rows(f)), "tu_id"].tolist()
    assert mod().enough_data(tu) is False
    uni = _uni(f)
    tbl = mod().model_rows_table({}, uni)
    assert list(tbl.columns) == list(d.TABLE_COLUMNS)
    for m in d.MODELS:
        sel = rows(tbl, model=m)
        assert len(sel) == 1
        r = sel.iloc[0]
        assert r["metric"] == "cv_r2_mean" and is_na(r["value"])
        assert flags_of(r["flags"]) == {"too_few_rows"}
    assert len(rows(tbl, model="univariate")) == 4
    assert finite(one(tbl, feature_or_block="cai", metric="spearman_rho")["value"])
