"""M8-09 SHAP block attribution, top decile, dominant block and residual flag: examples 5-8, the
degenerate rows and the §7.8.3 unit "SHAP block sums equal the prediction minus the expected
value" (docs/specs/M8-09-determinants-attribution-null.md rules 1-4, decisions Q2-Q3;
GIAP_SPEC.md §7.8.2 steps 5-6, §9.3, §9.4).

Contract choices the interface leaves open are in ``tests/helpers/determinants_attrib_data.py``.
SHAP tests fit real LightGBM models (3 repeats, 40 trees); the table tests use a synthetic
``CvResult`` and SHAP frame so every number is set by hand.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_attrib_data as a
import determinants_models_data as d
from determinants_models_data import flags_of, is_na

from giap.config.model import DeterminantsConfig, LightgbmConfig
from giap.rng import for_step

SEED = 20261021
TOL = 1e-9


def _cfg() -> DeterminantsConfig:
    return DeterminantsConfig(repeats=3, lightgbm=LightgbmConfig(n_estimators=40))


def _fit(threads: int = 1) -> tuple[pd.DataFrame, np.ndarray, list[str], dict[str, str], Any]:
    X, y, tu, blocks = a.toy_world(np.random.default_rng(5))
    cfg = _cfg()
    folds = d.mod().tu_folds(tu, cfg.outer_folds, cfg.repeats, for_step("determinants.cv", SEED))
    cv = d.mod().cross_validate(X, y, folds, blocks, cfg, SEED, threads, tu=tu)
    return X, y, tu, blocks, cv


@pytest.fixture(scope="module")
def fitted() -> tuple[pd.DataFrame, np.ndarray, list[str], dict[str, str], Any]:
    return _fit()


def _ribosomal_rows(X: pd.DataFrame, k: int = 4) -> pd.DataFrame:
    rng = np.random.default_rng(9)
    return pd.DataFrame(
        {c: rng.normal(size=k) for c in X.columns},
        index=pd.Index([f"R{i:05d}" for i in range(k)], name="locus_tag"),
    )


# ---------------------------------------------------------------------------------------------
# rule 1 / example 5: the SHAP identity


def test_example_5_shap_identity_every_row(fitted: Any) -> None:
    X, _, _, blocks, cv = fitted
    shap = a.attr_mod().shap_blocks(cv["lightgbm"], X, blocks)
    assert list(shap.index) == list(X.index)
    for c in (*a.SHAP_COLUMNS, "expected_value", "predicted_log2"):
        assert c in shap.columns, c
    total = shap[list(a.SHAP_COLUMNS)].astype(float).fillna(0.0).sum(axis=1).to_numpy()
    expected = shap["expected_value"].to_numpy(float)
    oof = np.asarray(cv["lightgbm"].oof, dtype=float)
    np.testing.assert_allclose(total + expected, oof, rtol=0, atol=TOL)
    np.testing.assert_allclose(shap["predicted_log2"].to_numpy(float), oof, rtol=0, atol=TOL)


def test_rows_out_of_fitted_order_raise(fitted: Any) -> None:
    """Review B1: fold labels and OOF are positional, so shuffled rows must not pass silently."""
    X, _, _, blocks, cv = fitted
    shuffled = X.iloc[np.random.default_rng(5).permutation(len(X))]
    with pytest.raises(ValueError, match="fitted order"):
        a.attr_mod().shap_blocks(cv["lightgbm"], shuffled, blocks)


def test_blocks_without_columns_are_na(fitted: Any) -> None:
    X, _, _, blocks, cv = fitted
    shap = a.attr_mod().shap_blocks(cv["lightgbm"], X, blocks)
    for b in ("translation", "context"):
        assert shap[f"shap_{b}"].isna().all(), b
    for b in ("transcription", "codon", "protein"):
        v = shap[f"shap_{b}"].to_numpy(float)
        assert np.all(np.isfinite(v)), b
    # the planted signal: codon (y = cai - 0.5 gravy + 0.3 promotech) carries most attribution
    mag = {b: float(np.abs(shap[f"shap_{b}"].to_numpy(float)).mean()) for b in
           ("transcription", "codon", "protein")}  # fmt: skip
    assert max(mag, key=mag.__getitem__) == "codon"


def test_shap_is_mean_over_repeats_of_the_held_out_fold_model(fitted: Any) -> None:
    """Rule 1 by hand: booster (r, k) is ``boosters[r * K + k]`` (repeat x fold, M8-08); per gene,
    the ``pred_contrib`` vectors of the model that held it out in each repeat, averaged, then
    summed by block."""
    X, _, _, blocks, cv = fitted
    res = cv["lightgbm"]
    folds = np.asarray(res.folds)
    R, n = folds.shape
    K = int(folds.max()) + 1
    assert len(res.boosters) == R * K
    contrib = {b: np.zeros(n) for b in ("transcription", "codon", "protein")}
    expected = np.zeros(n)
    for r in range(R):
        for k in range(K):
            bst = res.boosters[r * K + k]
            test = folds[r] == k
            names = list(bst.feature_name_)
            c = np.asarray(bst.predict(X.loc[test, names], pred_contrib=True), dtype=float)
            expected[test] += c[:, -1] / R
            for j, col in enumerate(names):
                contrib[blocks[col]][test] += c[:, j] / R
    shap = a.attr_mod().shap_blocks(res, X, blocks)
    np.testing.assert_allclose(shap["expected_value"].to_numpy(float), expected, atol=TOL, rtol=0)
    for b, v in contrib.items():
        np.testing.assert_allclose(shap[f"shap_{b}"].to_numpy(float), v, atol=TOL, rtol=0)


def test_ribosomal_extra_rows_use_every_fold_model(fitted: Any) -> None:
    """Decision Q3: rows the models never saw get the mean prediction (and SHAP) of every fold
    model of every repeat; the identity holds for them too."""
    X, _, _, blocks, cv = fitted
    res = cv["lightgbm"]
    extra = _ribosomal_rows(X)
    shap = a.attr_mod().shap_blocks(res, X, blocks, extra=extra)
    assert list(shap.index) == [*X.index, *extra.index]
    ex = shap.loc[extra.index]
    preds = np.mean(
        [np.asarray(b.predict(extra[list(b.feature_name_)]), dtype=float) for b in res.boosters],
        axis=0,
    )
    np.testing.assert_allclose(ex["predicted_log2"].to_numpy(float), preds, atol=TOL, rtol=0)
    total = ex[list(a.SHAP_COLUMNS)].astype(float).fillna(0.0).sum(axis=1).to_numpy()
    np.testing.assert_allclose(total + ex["expected_value"].to_numpy(float), preds, atol=TOL,
                               rtol=0)  # fmt: skip
    # the CV rows are unchanged by the extra rows
    base = a.attr_mod().shap_blocks(res, X, blocks)
    pd.testing.assert_frame_equal(shap.loc[X.index], base)


def test_elastic_net_has_no_shap(fitted: Any) -> None:
    X, _, _, blocks, cv = fitted
    with pytest.raises(ValueError):
        a.attr_mod().shap_blocks(cv["elastic_net"], X, blocks)


def test_example_8_shap_and_flags_identical_across_threads(fitted: Any) -> None:
    X, y, tu, blocks, cv1 = fitted
    *_, cv4 = _fit(threads=4)
    s1 = a.attr_mod().shap_blocks(cv1["lightgbm"], X, blocks)
    s4 = a.attr_mod().shap_blocks(cv4["lightgbm"], X, blocks)
    pd.testing.assert_frame_equal(s1, s4, check_exact=True)
    f = d.features({c: X[c].tolist() for c in X.columns}, tu=tu, abundance=y.tolist(),
                   tags=list(X.index))  # fmt: skip
    rows = d.mod().model_rows(f)
    t1 = a.attr_mod().attribution_table(f, rows, cv1["lightgbm"], s1)
    t4 = a.attr_mod().attribution_table(f, rows, cv4["lightgbm"], s4)
    pd.testing.assert_frame_equal(t1, t4, check_exact=True)


def test_attribution_from_real_fit_with_ribosomal_rows(fitted: Any) -> None:
    X, y, tu, blocks, cv = fitted
    extra = _ribosomal_rows(X)
    shap = a.attr_mod().shap_blocks(cv["lightgbm"], X, blocks, extra=extra)
    n, k = len(X), len(extra)
    f = d.features(
        {c: [*X[c], *extra[c]] for c in X.columns},
        tu=[*tu, *[f"TR_{i}" for i in range(k)]],
        abundance=[*y, *[30.0] * k],
        ribosomal=[False] * n + [True] * k,
        tags=[*X.index, *extra.index],
    )
    rows = d.mod().model_rows(f)
    assert rows.sum() == n
    tbl = a.attr_mod().attribution_table(f, rows, cv["lightgbm"], shap).set_index("locus_tag")
    assert set(tbl.index) == set(X.index) | set(extra.index)
    np.testing.assert_allclose(tbl.loc[X.index, "predicted_log2"].to_numpy(float),
                               np.asarray(cv["lightgbm"].oof), atol=TOL, rtol=0)  # fmt: skip
    np.testing.assert_allclose(tbl.loc[extra.index, "predicted_log2"].to_numpy(float),
                               shap.loc[extra.index, "predicted_log2"].to_numpy(float), atol=TOL,
                               rtol=0)  # fmt: skip
    for g in extra.index:
        assert "cai_reference" in flags_of(tbl.at[g, "flags"])
        assert not bool(tbl.at[g, "top_decile"])  # abundance 30 would top every model row
    for g in X.index:
        assert "cai_reference" not in flags_of(tbl.at[g, "flags"])


# ---------------------------------------------------------------------------------------------
# rules 2-4 on a hand-built table


N = 100
RIBO = ["R00000", "R00001"]
UNDETECTED = ["U00000", "U00001", "U00002"]


def _resid_sd1() -> np.ndarray:
    """Example 6: 98 residuals at ±c, one +2.5, one -2.3, sample SD exactly 1."""
    m = 0.2 / N
    c = math.sqrt((N - 1 - (2.5**2 + 2.3**2) + N * m**2) / 98)
    r = np.array([c if i % 2 == 0 else -c for i in range(98)] + [2.5, -2.3])
    assert abs(r.std(ddof=1) - 1.0) < 1e-12
    return r


def _table(
    resid: np.ndarray, shap_rows: dict[int, list[float]] | None = None,
    ribo_ab: float = 1000.0, ribo_pred: float = 950.0,
) -> tuple[pd.DataFrame, dict[str, Any]]:  # fmt: skip
    """100 model rows with abundance 0..99 (gene ``G000ii`` has abundance ii) and residual
    ``resid[(ii + 2) % 100]`` (example 6's +2.5 lands on gene 96, its -2.3 on gene 97); 2 ribosomal
    genes (residual ``ribo_ab - ribo_pred``); 3 undetected genes. ``shap_rows`` maps abundance ->
    five block sums (NaN = NA block); every other gene gets small random sums."""
    model_tags = [f"G{i:05d}" for i in range(N)]
    ab = np.arange(N, dtype=float)
    res = np.roll(resid, -2)  # gene i gets resid[(i + 2) % N]
    pred = ab - res
    tags = [*model_tags, *RIBO, *UNDETECTED]
    m = len(tags)
    f = d.features(
        {"gravy": [0.0] * m},
        tu=[f"TU_{i:04d}" for i in range(m)],
        detected=[True] * (N + 2) + [False] * 3,
        abundance=[*ab, ribo_ab, ribo_ab, math.nan, math.nan, math.nan],
        ribosomal=[False] * N + [True, True] + [False] * 3,
        tags=tags,
    )
    rng = np.random.default_rng(11)
    sums = rng.normal(0.0, 0.05, (N + 2, 5))
    for i, v in (shap_rows or {}).items():
        sums[i] = v
    preds = np.array([*pred, ribo_pred, ribo_pred])
    shap = pd.DataFrame(sums, columns=list(a.SHAP_COLUMNS),
                        index=pd.Index([*model_tags, *RIBO], name="locus_tag"))  # fmt: skip
    shap["expected_value"] = preds - np.nansum(sums, axis=1)
    shap["predicted_log2"] = preds
    cv = d.mod().CvResult("lightgbm", np.array([0.5]), pred.copy(), np.zeros((1, N), np.int64),
                          {}, ())  # fmt: skip
    rows = d.mod().model_rows(f)
    assert int(rows.sum()) == N
    tbl = a.attr_mod().attribution_table(f, rows, cv, shap)
    return tbl, {"ab": ab, "pred": pred, "res": res, "tags": model_tags}


def test_columns_and_rows() -> None:
    tbl, info = _table(_resid_sd1())
    assert tuple(tbl.columns) == a.ATTRIBUTION_COLUMNS
    assert set(tbl["locus_tag"]) == set(info["tags"]) | set(RIBO)  # undetected genes: no row
    t = tbl.set_index("locus_tag")
    np.testing.assert_allclose(t.loc[info["tags"], "predicted_log2"].to_numpy(float),
                               info["pred"], atol=TOL, rtol=0)  # fmt: skip
    np.testing.assert_allclose(t.loc[info["tags"], "residual"].to_numpy(float), info["res"],
                               atol=TOL, rtol=0)  # fmt: skip
    np.testing.assert_allclose(t.loc[info["tags"], "abundance_log2"].to_numpy(float), info["ab"])
    for b in a.BLOCKS:
        assert t[f"shap_{b}"].notna().all(), b


def test_top_decile_is_90th_percentile_of_model_rows() -> None:
    """Abundance 0..99: the 90th percentile (linear) is 89.1 -> genes 90..99 (10 genes); the
    ribosomal genes (abundance 1000) are never in the top decile and do not move the cut."""
    tbl, _ = _table(_resid_sd1())
    t = tbl.set_index("locus_tag")
    top = {g for g in t.index if bool(t.at[g, "top_decile"])}
    assert top == {f"G{i:05d}" for i in range(90, 100)}
    assert not any(bool(t.at[g, "top_decile"]) for g in RIBO)


def test_example_7_dominant_block() -> None:
    shap_rows = {
        99: [0.4, -0.1, 0.6, 0.0, -0.2],  # example 7 -> codon
        98: [-0.1, 0.0, -0.3, -0.2, 0.0],  # none positive -> NA
        97: [0.5, 0.1, 0.5, 0.2, 0.0],  # tie transcription / codon -> BLOCKS order: transcription
        96: [0.1, 0.05, math.nan, 0.2, math.nan],  # NA blocks ignored -> protein
        10: [0.0, 0.0, 9.0, 0.0, 0.0],  # not in the top decile -> NA
    }
    tbl, _ = _table(_resid_sd1(), shap_rows)
    t = tbl.set_index("locus_tag")
    assert t.at["G00099", "dominant_block"] == "codon"
    assert is_na(t.at["G00098", "dominant_block"])
    assert t.at["G00097", "dominant_block"] == "transcription"
    assert t.at["G00096", "dominant_block"] == "protein"
    assert is_na(t.at["G00010", "dominant_block"])
    assert is_na(t.at["G00096", "shap_codon"]) and is_na(t.at["G00096", "shap_context"])
    for g in RIBO:
        assert is_na(t.at[g, "dominant_block"])
    for i in range(90):
        assert is_na(t.at[f"G{i:05d}", "dominant_block"])


def test_example_6_two_sided_residual_flag() -> None:
    """SD 1 over the model rows; +2.5 and -2.3 flagged, the ±0.94 rest not. The ribosomal rows'
    residual of 50 is not in the SD (with it the SD would be about 7 and nothing flagged)."""
    tbl, info = _table(_resid_sd1())
    t = tbl.set_index("locus_tag")
    flagged = {g for g in info["tags"] if a.FLAG in flags_of(t.at[g, "flags"])}
    big = {g for g, r in zip(info["tags"], info["res"], strict=True) if abs(r) > 2.0}
    assert len(big) == 2
    assert flagged == big
    signs = sorted(float(t.at[g, "residual"]) for g in flagged)
    assert signs[0] == pytest.approx(-2.3) and signs[1] == pytest.approx(2.5)
    for g in RIBO:
        assert "cai_reference" in flags_of(t.at[g, "flags"])
        assert a.FLAG not in flags_of(t.at[g, "flags"]), g  # decision Q3 (residual 50)


def test_residual_flag_strictly_beyond_two_sd() -> None:
    """Generic residuals: flagged iff |r| > 2 SD (ddof 1, model rows only)."""
    r = np.random.default_rng(2).normal(size=N)
    r = (r - r.mean()) / r.std(ddof=1)
    r[r > 1.9] = 1.0
    r[r < -1.9] = -1.0
    tbl, info = _table(r)
    sd = float(info["res"].std(ddof=1))
    lim = 2.0 * sd
    t = tbl.set_index("locus_tag")
    for g, x in zip(info["tags"], info["res"], strict=True):
        assert (a.FLAG in flags_of(t.at[g, "flags"])) == (abs(x) > lim), (g, x, lim)


def test_residual_sd_zero_sets_no_flag() -> None:
    tbl, info = _table(np.full(N, 0.5))  # dyadic: residuals exactly equal
    assert float(np.std(info["res"])) == 0.0
    for v in tbl["flags"]:
        assert a.FLAG not in flags_of(v)


def test_fewer_than_100_model_rows() -> None:
    """Degenerate row: no SHAP, no null; the model rows get NA predictions and ``too_few_rows``."""
    n = 60
    f = d.features({"gravy": [0.0] * n}, tu=[f"TU_{i:04d}" for i in range(n)],
                   abundance=list(np.arange(n, dtype=float)))  # fmt: skip
    rows = d.mod().model_rows(f)
    tbl = a.attr_mod().attribution_table(f, rows, None, None)
    assert tuple(tbl.columns) == a.ATTRIBUTION_COLUMNS
    assert len(tbl) == n
    assert tbl["predicted_log2"].isna().all() and tbl["residual"].isna().all()
    for c in a.SHAP_COLUMNS:
        assert tbl[c].isna().all(), c
    assert tbl["dominant_block"].isna().all()
    for v in tbl["flags"]:
        assert "too_few_rows" in flags_of(v)
        assert a.FLAG not in flags_of(v)
