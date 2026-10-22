"""Determinants models: univariate rho / AUROC, repeated grouped CV, drop-one-block ΔR²
(GIAP_SPEC.md §7.8.2 steps 1-3 and 5, §5.2, §6.6, §9.3; docs/specs/M8-08-determinants-models.md).

Pure functions, no file I/O. Row order is canonicalised (``univariate`` on ``locus_tag``,
``cross_validate`` on ``X.index``, ``tu_folds`` on the TU ID), so permuting rows changes nothing.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Literal

import lightgbm as lgb
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from giap.config.model import DeterminantsConfig
from giap.determinants.features import model_matrix
from giap.stats import bootstrap

__all__ = ["MIN_ROWS", "MIN_TUS", "TABLE_COLUMNS", "CvResult", "cross_validate", "enough_data",
           "model_rows", "model_rows_table", "tu_folds", "univariate"]  # fmt: skip

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
B8 = npt.NDArray[np.bool_]
ModelName = Literal["lightgbm", "elastic_net"]

MIN_ROWS = 100
MIN_TUS = 25
MIN_UNIVARIATE_ROWS = 10
MODELS: tuple[ModelName, ...] = ("lightgbm", "elastic_net")
TABLE_COLUMNS: tuple[str, ...] = (
    "model", "feature_or_block", "metric", "value", "ci_low", "ci_high", "n", "n_perm", "flags",
)  # fmt: skip
#: ribosomal genes are excluded from the univariate statistics of these features (§7.8.2 step 1)
RIBOSOMAL_EXCLUDED: frozenset[str] = frozenset({"cai"})


@dataclass(frozen=True, slots=True)
class CvResult:
    """Repeated grouped CV of one model (spec rules 4-7)."""

    model: ModelName
    r2: F64  # per repeat (pooled out-of-fold R²)
    oof: F64  # out-of-fold prediction averaged over repeats, per row (input order)
    folds: I64  # (repeats, n_rows) fold index (input order)
    delta_r2: Mapping[str, F64]  # block -> per-repeat ΔR²
    boosters: tuple[object, ...]  # repeat x fold fitted LightGBM models (M8-09 SHAP)
    dropped: tuple[str, ...] = ()  # sorted columns zero-variance / all-NaN in a training fold


def _flag(s: pd.Series) -> B8:
    return np.asarray(s.astype("boolean").fillna(False).to_numpy(dtype=bool))


def model_rows(features: pd.DataFrame) -> B8:
    """Rule 1: detected, finite non-imputed ``abundance_log2``, not ribosomal (every TSV row is
    eligible)."""
    ab = pd.to_numeric(features["abundance_log2"], errors="coerce").to_numpy(float, na_value=np.nan)
    keep = (
        _flag(features["detected"])
        & np.isfinite(ab)
        & ~_flag(features["abundance_imputed"])
        & ~_flag(features["ribosomal_protein"])
    )
    return np.asarray(keep, dtype=bool)


def enough_data(tu: Sequence[str]) -> bool:
    """Rule 8 (decision 18): at least ``MIN_ROWS`` model rows in at least ``MIN_TUS`` TUs."""
    return len(tu) >= MIN_ROWS and len(set(tu)) >= MIN_TUS


# ---------------------------------------------------------------------------- univariate


def _rho(df: pd.DataFrame) -> float:
    x, y = df["x"].to_numpy(float), df["y"].to_numpy(float)
    if x.size < MIN_UNIVARIATE_ROWS or np.ptp(x) == 0 or np.ptp(y) == 0:
        return math.nan
    return float(stats.spearmanr(x, y).statistic)


def _auroc(df: pd.DataFrame) -> float:
    """Mann-Whitney AUROC of ``x`` for ``y`` (bool), ties ½."""
    x, y = df["x"].to_numpy(float), df["y"].to_numpy(bool)
    n1 = int(y.sum())
    n0 = y.size - n1
    if y.size < MIN_UNIVARIATE_ROWS or n1 == 0 or n0 == 0:
        return math.nan
    ranks = stats.rankdata(x)
    return float((ranks[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _uni_row(
    feature: str, metric: str, df: pd.DataFrame, fn: Callable[[pd.DataFrame], float], B: int,
    rng: np.random.Generator,
) -> dict[str, Any]:  # fmt: skip
    if len(df) < MIN_UNIVARIATE_ROWS:
        est, lo, hi = math.nan, math.nan, math.nan
    else:
        ci = bootstrap.tu_block_bootstrap(df, "tu", fn, B, rng)
        est, lo, hi = ci.estimate, ci.ci_low, ci.ci_high
    # lead decision 7: the TU-block bootstrap CI of a rank statistic under-covers (SE 5-9 % low
    # under a TU-structured null), so every univariate CI is ci_unestimated; NA value: unestimated
    flags = "ci_unestimated" if math.isfinite(est) else "unestimated;ci_unestimated"
    return {"model": "univariate", "feature_or_block": feature, "metric": metric, "value": est,
            "ci_low": lo, "ci_high": hi, "n": len(df), "n_perm": 0, "flags": flags}  # fmt: skip


def univariate(features: pd.DataFrame, B: int, rng: np.random.Generator) -> pd.DataFrame:
    """Rule 2: per encoded column, Spearman rho with ``abundance_log2`` on model rows and AUROC
    for ``detected`` on all rows (ribosomal genes excluded for ``RIBOSOMAL_EXCLUDED``), each with
    a 95 % percentile TU-block bootstrap CI over ``B`` draws; < 10 usable rows -> NA."""
    f = features.sort_values("locus_tag", kind="stable").reset_index(drop=True)
    X, _ = model_matrix(f)
    on_model = model_rows(f)
    ribosomal = _flag(f["ribosomal_protein"])
    tu = f["tu_id"].astype(str).to_numpy()
    ab = pd.to_numeric(f["abundance_log2"], errors="coerce").to_numpy(float, na_value=np.nan)
    det = _flag(f["detected"])
    out: list[dict[str, Any]] = []
    for col in X.columns:
        x = X[col].to_numpy(float)
        ok = np.isfinite(x)
        auc_rows = ok & ~ribosomal if col in RIBOSOMAL_EXCLUDED else ok
        m, a = ok & on_model, auc_rows
        rho_df = pd.DataFrame({"tu": tu[m], "x": x[m], "y": ab[m]})
        auc_df = pd.DataFrame({"tu": tu[a], "x": x[a], "y": det[a]})
        out.append(_uni_row(col, "spearman_rho", rho_df, _rho, B, rng))
        out.append(_uni_row(col, "auroc", auc_df, _auroc, B, rng))
    return _frame(out)


def _frame(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(TABLE_COLUMNS))
    df["value"] = df["value"].astype(float)
    df[["ci_low", "ci_high"]] = df[["ci_low", "ci_high"]].astype(float)
    df[["n", "n_perm"]] = df[["n", "n_perm"]].astype("Int64")
    return df


# ---------------------------------------------------------------------------- folds


def tu_folds(tu: Sequence[str], n_folds: int, repeats: int, rng: np.random.Generator) -> I64:
    """Rule 3 (lead decision 7): per repeat, shuffle the (sorted) TU IDs and assign each, in the
    shuffled order, to the fold with the fewest genes so far (ties between folds drawn from
    ``rng``), so repeats differ even when every TU size is distinct. Returns (repeats, n_rows)."""
    ids, codes, sizes = np.unique(
        np.asarray(tu, dtype=str), return_inverse=True, return_counts=True
    )
    out = np.empty((repeats, codes.size), dtype=np.int64)
    for r in range(repeats):
        order = rng.permutation(ids.size)
        load = np.zeros(n_folds, dtype=np.int64)
        fold_of = np.empty(ids.size, dtype=np.int64)
        for t in order:
            least = np.flatnonzero(load == load.min())
            k = int(least[0] if least.size == 1 else rng.choice(least))
            fold_of[t] = k
            load[k] += sizes[t]
        out[r] = fold_of[codes]
    return out


# ---------------------------------------------------------------------------- cross-validation


def _varying(Xtr: F64) -> B8:
    """Columns with two distinct finite values on the training rows (decision 2)."""
    has = ~np.isnan(Xtr)
    lo = np.where(has, Xtr, np.inf).min(axis=0, initial=np.inf)
    hi = np.where(has, Xtr, -np.inf).max(axis=0, initial=-np.inf)
    return np.asarray(hi > lo, dtype=bool)


def _fit_predict(
    model: ModelName, Xtr: F64, ytr: F64, Xte: F64, groups: npt.NDArray[Any], names: list[str],
    cfg: DeterminantsConfig, seed: int, threads: int,
) -> tuple[F64, object | None]:  # fmt: skip
    if Xtr.shape[1] == 0:
        return np.full(Xte.shape[0], float(ytr.mean())), None
    if model == "lightgbm":
        p = cfg.lightgbm
        booster = lgb.LGBMRegressor(
            n_estimators=p.n_estimators, learning_rate=p.learning_rate, num_leaves=p.num_leaves,
            min_child_samples=p.min_child_samples, subsample=p.subsample,
            subsample_freq=p.subsample_freq, colsample_bytree=p.colsample_bytree,
            deterministic=p.deterministic, force_row_wise=True, n_jobs=threads,
            random_state=seed, verbose=-1,
        )  # fmt: skip
        booster.fit(pd.DataFrame(Xtr, columns=names), ytr)
        return np.asarray(booster.predict(pd.DataFrame(Xte, columns=names)), dtype=float), booster
    inner = list(GroupKFold(n_splits=cfg.elastic_net.inner_folds).split(Xtr, ytr, groups))
    # random_state is unused by selection="cyclic"; fixed so no global state is ever touched
    enet = ElasticNetCV(
        l1_ratio=list(cfg.elastic_net.l1_ratio), cv=inner, alphas=100, max_iter=10000,
        selection="cyclic", random_state=0,
    )  # fmt: skip
    pipe = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), enet)
    pipe.fit(Xtr, ytr)
    return np.asarray(pipe.predict(Xte), dtype=float), None


def _r2(y: F64, pred: F64) -> float:
    sst = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - float(np.sum((y - pred) ** 2)) / sst if sst > 0 else math.nan


def _repeated_oof(
    model: ModelName, X: F64, y: F64, folds: I64, groups: npt.NDArray[Any], cols: B8,
    names: list[str], cfg: DeterminantsConfig, seed: int, threads: int,
) -> tuple[F64, list[object], set[str]]:  # fmt: skip
    """Out-of-fold predictions (repeats, n) on the columns ``cols``; boosters; dropped columns."""
    pred = np.empty(folds.shape, dtype=float)
    boosters: list[object] = []
    dropped: set[str] = set()
    for r in range(folds.shape[0]):
        for k in np.unique(folds[r]).tolist():
            test = folds[r] == k
            train = ~test
            varying = _varying(X[train])
            dropped |= {n for n, c, v in zip(names, cols, varying, strict=True) if c and not v}
            use = np.flatnonzero(cols & varying)
            s = int(np.random.SeedSequence([seed, r, k]).generate_state(1)[0] >> 1)
            pred[r, test], b = _fit_predict(
                model, X[np.ix_(train, use)], y[train], X[np.ix_(test, use)], groups[train],
                [names[i] for i in use], cfg, s, threads,
            )  # fmt: skip
            if b is not None:
                boosters.append(b)
    return pred, boosters, dropped


def cross_validate(
    X: pd.DataFrame, y: F64, folds: I64, blocks: Mapping[str, str], cfg: DeterminantsConfig,
    seed: int, threads: int, *, tu: Sequence[str] | None = None,
) -> dict[str, CvResult]:  # fmt: skip
    """Rules 4-7: LightGBM and elastic net on the given outer folds (repeats x folds), pooled
    out-of-fold R² per repeat, and drop-one-block ΔR² on the same folds.

    ``tu`` (the TU ID per row) groups the elastic net's inner ``GroupKFold``; without it every row
    is its own group. Constant y -> R² and ΔR² NaN (unestimated)."""
    order = np.asarray(X.index.argsort(kind="stable"))
    Xs = X.to_numpy(dtype=float)[order]
    ys = np.asarray(y, dtype=float)[order]
    folds = np.asarray(folds, dtype=np.int64)
    fs = folds[:, order]
    groups = np.arange(len(ys)) if tu is None else np.asarray(tu, dtype=str)[order]
    names = [str(c) for c in X.columns]
    col_block = [blocks.get(c) for c in names]
    present = list(dict.fromkeys(b for b in col_block if b is not None))
    everything = np.ones(len(names), dtype=bool)
    R = fs.shape[0]
    out: dict[str, CvResult] = {}
    for model in MODELS:
        if ys.size == 0 or np.ptp(ys) == 0:
            nan = np.full(R, math.nan)
            out[model] = CvResult(model, nan, np.full(len(ys), math.nan), folds,
                                  {b: nan.copy() for b in present}, ())  # fmt: skip
            continue
        pred, boosters, dropped = _repeated_oof(
            model, Xs, ys, fs, groups, everything, names, cfg, seed, threads
        )
        r2 = np.array([_r2(ys, pred[r]) for r in range(R)])
        delta: dict[str, F64] = {}
        for b in present:
            keep = np.array([cb != b for cb in col_block], dtype=bool)
            p_drop, _, _ = _repeated_oof(model, Xs, ys, fs, groups, keep, names, cfg, seed, threads)
            delta[b] = r2 - np.array([_r2(ys, p_drop[r]) for r in range(R)])
        oof = np.empty(len(ys), dtype=float)
        oof[order] = pred.mean(axis=0)
        out[model] = CvResult(model, r2, oof, folds, delta, tuple(boosters), tuple(sorted(dropped)))
    return out


# ---------------------------------------------------------------------------- table


def _row(model: str, what: str, metric: str, value: float, n: Any, flags: list[str],
         lo: float = math.nan, hi: float = math.nan) -> dict[str, Any]:  # fmt: skip
    return {"model": model, "feature_or_block": what, "metric": metric, "value": value,
            "ci_low": lo, "ci_high": hi, "n": n, "n_perm": 0, "flags": ";".join(flags)}  # fmt: skip


def model_rows_table(cv: Mapping[str, CvResult], uni: pd.DataFrame) -> pd.DataFrame:
    """``determinants_model.tsv`` rows (§9.3 columns plus ``flags``): per model ``cv_r2_mean``,
    ``cv_r2_sd`` (ddof 1) and per block ``delta_r2`` with ``mean ± t_{0.975,R-1}·sd/√R``, then the
    univariate rows. Empty ``cv`` (too few rows) -> one NA ``cv_r2_mean`` row per model. Every
    ΔR² row is flagged ``ci_unestimated`` (lead decision 5: the interval under-covers)."""
    rows: list[dict[str, Any]] = []
    if not cv:
        rows += [_row(m, "all", "cv_r2_mean", math.nan, pd.NA, ["too_few_rows"]) for m in MODELS]
    for name, res in cv.items():
        r2 = np.asarray(res.r2, dtype=float)
        n = len(res.oof)
        bad = [] if np.all(np.isfinite(r2)) else ["unestimated"]
        sd = float(r2.std(ddof=1)) if r2.size > 1 else math.nan
        rows.append(
            _row(name, "all", "cv_r2_mean", float(r2.mean()) if not bad else math.nan, n, bad)
        )
        rows.append(_row(name, "all", "cv_r2_sd", sd if not bad else math.nan, n, bad))
        for block, v in res.delta_r2.items():
            v = np.asarray(v, dtype=float)
            mean = float(v.mean())
            half = math.nan
            if v.size > 1:
                half = (
                    float(stats.t.ppf(0.975, v.size - 1)) * float(v.std(ddof=1)) / math.sqrt(v.size)
                )
            fl = ([] if math.isfinite(mean) else ["unestimated"]) + ["ci_unestimated"]
            rows.append(_row(name, block, "delta_r2", mean, n, fl, mean - half, mean + half))
    tbl = _frame(rows)
    if len(uni):
        tbl = pd.concat([tbl, uni[list(TABLE_COLUMNS)]], ignore_index=True) if len(tbl) else uni
    return tbl.reset_index(drop=True)
