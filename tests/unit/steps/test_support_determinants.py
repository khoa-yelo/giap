"""M9-17: ``determinants_shap``, ``determinants_pd`` and ``null_determinants`` from the
``determinants`` step (V15(b,c), V18; docs/specs/M9-17-plot-support-analysis.md "Support files",
"Definitions", examples 2, 3 and 5, the degenerate cases ``too_few_rows`` / no abundance /
``layer_missing:ostir``, the tables-unchanged acceptance; lead decisions 1 and 3).

Fixture: the M8-10 one (``helpers.determinants_step_data``: sim1000 through the M4 chain, B = 49,
2 repeats x 3 folds, 30 rounds), with the four PD features filled in the gene table (the sim1000
gene table has ``cai``, ``tai``, ``tir_ostir`` and ``promoter_score`` NA everywhere, which would
make every PD check vacuous). The test re-runs the step's model chain (``feature_frame`` ->
``model_matrix`` -> ``tu_folds`` on ``determinants.cv`` -> ``cross_validate``), which is
deterministic, to get the boosters the step fitted.

Interface and readings pinned here (test author; contract in tests/helpers/
support_analysis_data.py):

* ``giap.determinants.attribution.shap_features(cv, X, blocks, *, extra=None)``: the per-column
  matrix ``shap_blocks`` sums by block (same rows, index ``X.index`` then ``extra.index``;
  columns the encoded model columns the boosters use). Example 2: its block sums equal
  ``shap_blocks`` to 1e-9.
* ``determinants_shap``: model rows only; ``feature`` = an encoded model column; the top
  ``min(20, n columns)`` columns by mean |SHAP| over the model rows (ties by column name);
  ``feature_value`` the model-matrix value (NA when NaN).
* ``determinants_pd``: rows in ``cai, tai, tir_ostir, promotech_max`` order, ``grid_value``
  ascending, on the **feature's own scale** (``determinants_features.tsv``; the step applies the
  model's ln to ``tir_ostir`` before predicting); grid = numpy's default (linear) quantiles at
  5 %, ..., 95 % plus the min and max over the model rows with a value; per booster the mean
  prediction over **all** model rows with that column set; ``pd_mean`` the mean over boosters
  and ``pd_low`` / ``pd_high`` their min / max.
* ``null_determinants``: per model (lightgbm, then elastic_net) ``permutation`` 1..B in b order,
  ``cv_r2`` the null R² (NA when NaN) and ``observed_cv_r2_mean`` the model table's
  ``cv_r2_mean`` cell.
"""

from __future__ import annotations

import dataclasses
import math
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import support_analysis_data as sa

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_attrib_data as da
import determinants_models_data as dm
import determinants_step_data as s

SEED = 20261022
KEY = sa.SUPPORT_KEYS["determinants"]
NAMES = ("determinants_pd", "determinants_shap", "null_determinants")
TOL = 1e-9


@pytest.fixture(scope="module")
def step_mod() -> ModuleType:
    return s.mod()


@pytest.fixture(scope="module")
def work(step_mod: ModuleType, tmp_path_factory: pytest.TempPathFactory) -> Path:
    return s.build_work(tmp_path_factory.mktemp("m917_det_work"), SEED)


def _with_pd_features(work: Path, dest: Path, *, ostir: bool = True) -> Path:
    """The M4 gene table with synthetic CAI, tAI, OSTIR TIR and promoter score (``ostir=False``:
    ``tir_ostir`` left NA, the ``layer_missing:ostir`` case)."""
    g = pd.read_parquet(work / "genetable" / "genes.parquet")
    rng = np.random.default_rng(917)
    n = len(g)
    g["cai"] = rng.uniform(0.2, 0.9, n)
    g["tai"] = rng.uniform(0.1, 0.6, n)
    g["tir_ostir"] = np.exp(rng.normal(6.0, 1.5, n)) if ostir else math.nan
    g["promoter_score"] = rng.uniform(0.0, 1.0, n)
    dest.parent.mkdir(parents=True, exist_ok=True)
    g.to_parquet(dest, index=False)
    return dest


def _run(cfg: Any, ins: dict[str, Path], out: Path, *, support: bool = True) -> dict[str, Path]:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    outs = s.outputs(out)
    if support:
        outs[KEY] = out / KEY
    reg.run_step(s.STEP, cfg, dict(ins), outs)
    return outs


@pytest.fixture(scope="module")
def ins(work: Path, tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    gt = _with_pd_features(work, tmp_path_factory.mktemp("m917_gt") / "genes.parquet")
    return s.inputs(work, gene_table=gt)


@pytest.fixture(scope="module")
def full(ins: dict[str, Path], tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    return _run(s.config(SEED), ins, tmp_path_factory.mktemp("m917_det_full"))


@pytest.fixture(scope="module")
def chain(ins: dict[str, Path]) -> dict[str, Any]:
    """The step's model chain on the same inputs (deterministic: same seed, threads 1)."""
    from giap.determinants.features import model_matrix
    from giap.determinants.models import cross_validate, model_rows, tu_folds
    from giap.rng import for_step

    cfg = s.config(SEED)
    f = s.expected_features(ins)
    rows = model_rows(f)
    X_all, blocks = model_matrix(f)
    X = X_all[rows]
    y = pd.to_numeric(f.loc[rows, "abundance_log2"]).to_numpy(float)
    tu = f.loc[rows, "tu_id"].astype(str).tolist()
    d = cfg.determinants
    folds = tu_folds(tu, d.outer_folds, d.repeats, for_step("determinants.cv", SEED))
    cv = cross_validate(X, y, folds, blocks, d, SEED, 1, tu=tu)
    tags = f.loc[rows, "locus_tag"].astype(str).tolist()
    return {"features": f, "rows": rows, "X": X, "blocks": blocks, "cv": cv, "tags": tags}


def _table(paths: dict[str, Path], name: str) -> pd.DataFrame:
    return sa.read(sa.need(paths[KEY], name))


# ----------------------------------------------------------------------------- wiring


def test_registry_declares_the_optional_support_output(step_mod: ModuleType) -> None:
    spec = sa.step_spec("determinants")
    assert KEY in spec.outputs and KEY in spec.optional_outputs
    for key in s.OUTPUTS:
        assert key in spec.outputs and key not in spec.optional_outputs


def test_tables_bytes_unchanged_by_the_support_writer(
    ins: dict[str, Path], full: dict[str, Path], tmp_path: Path
) -> None:
    bare = _run(s.config(SEED), ins, tmp_path / "bare", support=False)
    assert sa.tables_bytes({k: full[k] for k in s.OUTPUTS}) == sa.tables_bytes(
        {k: bare[k] for k in s.OUTPUTS}
    )


def test_published_under_plots_support(ins: dict[str, Path], tmp_path: Path) -> None:
    work = tmp_path / "work"
    sa.run_with_spec_outputs("determinants", s.config(SEED, B=9), ins, work)
    view = sa.publish_and_view(work, tmp_path / "OUT")
    for name in NAMES:
        assert sa.find(tmp_path / "OUT" / "plots" / "support", name) is not None, name
        got = view.support(name)
        assert got is not None and list(got.columns) == list(sa.COLUMNS[name]), name


# ----------------------------------------------------------------------------- schema


def test_the_three_files_and_their_headers(full: dict[str, Path]) -> None:
    assert sa.names(full[KEY]) == sorted(NAMES)
    for name in NAMES:
        path = sa.need(full[KEY], name)
        assert sa.header(path) == list(sa.COLUMNS[name]), name
        sa.check_written_form(path)


@pytest.mark.parametrize("name", NAMES)
def test_validates_against_its_schema(full: dict[str, Path], name: str) -> None:
    import json

    frictionless = pytest.importorskip("frictionless")
    import giap.schemas

    path = sa.need(full[KEY], name)
    d = json.loads((Path(giap.schemas.__file__).parent / "plot_support" / f"{name}.json")
                   .read_text(encoding="utf-8"))  # fmt: skip
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ----------------------------------------------------------------------------- example 2: SHAP


@pytest.fixture(scope="module")
def toy() -> tuple[pd.DataFrame, dict[str, str], Any, pd.DataFrame]:
    from giap.config.model import DeterminantsConfig, LightgbmConfig
    from giap.rng import for_step

    X, y, tu, blocks = da.toy_world(np.random.default_rng(5))
    cfg = DeterminantsConfig(repeats=3, lightgbm=LightgbmConfig(n_estimators=40))
    folds = dm.mod().tu_folds(tu, cfg.outer_folds, cfg.repeats, for_step("determinants.cv", 3))
    cv = dm.mod().cross_validate(X, y, folds, blocks, cfg, 3, 1, tu=tu)
    rng = np.random.default_rng(9)
    extra = pd.DataFrame({c: rng.normal(size=4) for c in X.columns},
                         index=pd.Index([f"R{i}" for i in range(4)], name="locus_tag"))  # fmt: skip
    return X, blocks, cv["lightgbm"], extra


def test_example_2_shap_features_sum_by_block_to_shap_blocks(toy: Any) -> None:
    from giap.determinants.attribution import BLOCKS, shap_blocks, shap_features

    X, blocks, cv, extra = toy
    per = shap_features(cv, X, blocks, extra=extra)
    sums = shap_blocks(cv, X, blocks, extra=extra)
    assert list(per.index) == list(sums.index)
    assert set(per.columns) <= set(X.columns) and len(per.columns) > 0
    for b in BLOCKS:
        cols = [c for c in per.columns if blocks[c] == b]
        want = sums[f"shap_{b}"].to_numpy(float)
        if not cols:
            assert np.isnan(want).all(), b
            continue
        got = per[cols].to_numpy(float).sum(axis=1)
        np.testing.assert_allclose(got, want, rtol=0, atol=TOL, err_msg=b)


def test_shap_features_plus_expected_is_the_oof_prediction(toy: Any) -> None:
    from giap.determinants.attribution import shap_blocks, shap_features

    X, blocks, cv, _ = toy
    per = shap_features(cv, X, blocks)
    expected = shap_blocks(cv, X, blocks)["expected_value"].to_numpy(float)
    total = per.to_numpy(float).sum(axis=1) + expected
    np.testing.assert_allclose(total, np.asarray(cv.oof, float), rtol=0, atol=1e-6)


def test_shap_features_is_a_lightgbm_output(toy: Any) -> None:
    from giap.determinants.attribution import shap_features

    X, blocks, cv, _ = toy
    with pytest.raises(ValueError):
        shap_features(dataclasses.replace(cv, model="elastic_net"), X, blocks)


def _top(per: pd.DataFrame, k: int) -> list[str]:
    m = per.abs().mean(axis=0)
    return sorted(m.index, key=lambda c: (-float(m[c]), str(c)))[:k]


def test_shap_rows_are_the_model_rows_times_the_top_features(
    full: dict[str, Path], chain: dict[str, Any]
) -> None:
    from giap.determinants.attribution import shap_features

    shap = _table(full, "determinants_shap")
    per = shap_features(chain["cv"]["lightgbm"], chain["X"], chain["blocks"])
    per.index = pd.Index(chain["tags"])
    k = min(sa.SHAP_TOP, per.shape[1])
    want = _top(per, k)
    assert set(shap["locus_tag"]) == set(chain["tags"])
    assert len(shap) == len(chain["tags"]) * k
    assert not shap.duplicated(["locus_tag", "feature"]).any()
    assert sorted(set(shap["feature"])) == sorted(want)
    X = chain["X"].set_axis(chain["tags"])
    for tag, feat, value, phi in shap.itertuples(index=False, name=None):
        assert sa.close(sa.num(phi), float(per.at[tag, feat]), abs_=1e-6), (tag, feat)
        x = float(X.at[tag, feat])
        assert sa.close(sa.num(value), x), (tag, feat, value, x)


def test_example_2_on_the_step_full_blocks_match_attribution(
    full: dict[str, Path], chain: dict[str, Any]
) -> None:
    """Block sums of the per-column SHAP the file is cut from equal the attribution table's
    ``shap_<block>`` cells (%.6g) on every model row."""
    from giap.determinants.attribution import BLOCKS, shap_features

    per = shap_features(chain["cv"]["lightgbm"], chain["X"], chain["blocks"])
    per.index = pd.Index(chain["tags"])
    att = s.read(full["determinants_attribution"]).set_index("locus_tag")
    for b in BLOCKS:
        cols = [c for c in per.columns if chain["blocks"][c] == b]
        for tag in chain["tags"]:
            cell = sa.num(att.at[tag, f"shap_{b}"])
            if not cols:
                assert math.isnan(cell)
                continue
            assert sa.close(float(per.loc[tag, cols].sum()), cell, abs_=1e-5), (tag, b)


# ----------------------------------------------------------------------------- example 3: PD


def _grid(values: np.ndarray) -> np.ndarray:
    v = values[np.isfinite(values)]
    q = np.quantile(v, np.arange(5, 100, 5) / 100.0)
    return np.sort(np.concatenate([[v.min()], q, [v.max()]]))


def test_example_3_pd_rows_grid_and_order(full: dict[str, Path], chain: dict[str, Any]) -> None:
    pdt = _table(full, "determinants_pd")
    feats = list(dict.fromkeys(pdt["feature"]))
    assert feats == list(sa.PD_FEATURES), "all four PD features present, in spec order"
    f = chain["features"].loc[chain["rows"]]
    for feat, grp in pdt.groupby("feature", sort=False):
        assert len(grp) == sa.PD_GRID_POINTS, (feat, len(grp))
        g = sa.nums(grp["grid_value"])
        assert np.all(np.diff(g) >= 0), feat
        want = _grid(pd.to_numeric(f[feat]).to_numpy(float))
        np.testing.assert_allclose(g, want, rtol=sa.REL, err_msg=feat)


def test_example_3_spread_brackets_the_mean(full: dict[str, Path]) -> None:
    pdt = _table(full, "determinants_pd")
    lo, mid, hi = (sa.nums(pdt[c]) for c in ("pd_low", "pd_mean", "pd_high"))
    assert np.isfinite(mid).all()
    assert np.all(lo <= mid * (1 + np.sign(mid) * 1e-6) + 1e-9)
    assert np.all(mid <= hi * (1 + np.sign(hi) * 1e-6) + 1e-9)
    assert np.any(hi - lo > 0), "boosters disagree somewhere: a zero spread everywhere is vacuous"


def test_example_3_pd_is_the_booster_mean_prediction(
    full: dict[str, Path], chain: dict[str, Any]
) -> None:
    """For each booster, the mean prediction over the model rows with the column set to the grid
    value (``tir_ostir`` on its ln model scale); then the mean, min and max over boosters."""
    pdt = _table(full, "determinants_pd")
    X = chain["X"]
    boosters = chain["cv"]["lightgbm"].boosters
    assert len(boosters) == 6  # 2 repeats x 3 folds (example 3 has 25 with the defaults)
    for feat, grp in pdt.groupby("feature", sort=False):
        for g, m, lo, hi in zip(sa.nums(grp["grid_value"]), sa.nums(grp["pd_mean"]),
                                sa.nums(grp["pd_low"]), sa.nums(grp["pd_high"]),
                                strict=True):  # fmt: skip
            xg = math.log(g) if feat == "tir_ostir" else g
            per = []
            for bst in boosters:
                names = [str(c) for c in bst.feature_name_]
                Z = X[names].astype(float).copy()
                if feat in Z.columns:
                    Z[feat] = xg
                per.append(float(np.mean(bst.predict(Z))))
            assert sa.close(m, float(np.mean(per)), abs_=1e-6), (feat, g)
            assert sa.close(lo, min(per), abs_=1e-6) and sa.close(hi, max(per), abs_=1e-6)


# ----------------------------------------------------------------------------- the null


def test_null_determinants_rows(full: dict[str, Path]) -> None:
    null = _table(full, "null_determinants")
    model = s.read(full["determinants_model"])
    B = 49
    assert list(dict.fromkeys(null["model"])) == ["lightgbm", "elastic_net"]
    for m, grp in null.groupby("model", sort=False):
        assert grp["permutation"].astype(int).tolist() == list(range(1, B + 1)), m
        obs_cell = s.rows(model, model=m, metric="cv_r2_mean").iloc[0]["value"]
        assert set(grp["observed_cv_r2_mean"]) == {obs_cell}, m
        row = s.rows(model, model=m, metric="null_p").iloc[0]
        r2 = sa.nums(grp["cv_r2"])
        n = int(np.count_nonzero(~np.isnan(r2)))
        assert n == int(row["n_perm"]), m
        obs = sa.num(obs_cell)
        b = int(np.count_nonzero(r2[~np.isnan(r2)] >= obs))
        p = (1 + b) / (1 + n)
        assert abs(p - sa.num(row["value"])) <= 2 / (1 + n) + 1e-9, (m, p, row["value"])


# ----------------------------------------------------------------------------- determinism


def test_example_5_same_run_twice_byte_identical_threads_too(
    ins: dict[str, Path], full: dict[str, Path], tmp_path: Path
) -> None:
    """Same seed, run again with 2 threads: byte-identical support files (§6.6)."""
    again = _run(s.config(SEED, threads=2), ins, tmp_path / "again")
    assert sa.tree_bytes(full[KEY]) == sa.tree_bytes(again[KEY])


# ----------------------------------------------------------------------------- degenerate


def test_too_few_rows_writes_no_support(work: Path, ins: dict[str, Path], tmp_path: Path) -> None:
    src = tmp_path / "src"
    s_gt = s.gene_table_variant(work, src / "few.parquet", keep_abundance=60)
    g = pd.read_parquet(s_gt)
    pdg = pd.read_parquet(ins["gene_table"])
    for c in ("cai", "tai", "tir_ostir", "promoter_score"):
        g[c] = pdg[c].to_numpy()
    g.to_parquet(s_gt, index=False)
    out = _run(s.config(SEED), s.inputs(work, gene_table=s_gt), tmp_path / "few")
    assert out[KEY].is_dir(), "the declared support directory is always created"
    assert sa.names(out[KEY]) == []


def test_no_abundance_writes_no_support(work: Path, tmp_path: Path) -> None:
    gt = s.gene_table_variant(work, tmp_path / "noab.parquet", keep_abundance=None)
    out = _run(s.config(SEED), s.inputs(work, gene_table=gt), tmp_path / "noab")
    assert out[KEY].is_dir()
    assert sa.names(out[KEY]) == []


def test_ostir_missing_drops_its_pd_rows(work: Path, tmp_path: Path) -> None:
    gt = _with_pd_features(work, tmp_path / "gt" / "genes.parquet", ostir=False)
    out = _run(s.config(SEED, B=9), s.inputs(work, gene_table=gt), tmp_path / "noostir")
    pdt = _table(out, "determinants_pd")
    assert list(dict.fromkeys(pdt["feature"])) == ["cai", "tai", "promotech_max"]
    shap = _table(out, "determinants_shap")
    assert "tir_ostir" not in set(shap["feature"])
    assert len(_table(out, "null_determinants")) == 2 * 9
