"""Cross-fitted detectability model and its calibration diagnostics (GIAP_SPEC.md §7.3.3, G5).

``p_detect`` for every eligible gene comes from a spline + L2-logistic model fitted on the other
TUs (``GroupKFold``-style TU folds from a seeded TU order); C is chosen from {0.01, 0.1, 1, 10} by
an inner 3-fold TU-grouped search on mean held-out log loss (ties → smaller C). Fewer than 50
detected or undetected eligible genes, or fewer TUs than folds, gives the constant overall rate.
Expression determinants and tested categories are never covariates (spec docs/specs/M4-13).
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.special import expit, logit
from sklearn.compose import ColumnTransformer
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import SplineTransformer
from threadpoolctl import threadpool_limits  # type: ignore[import-untyped]

from giap.config.model import DetectabilityConfig

__all__ = [
    "DetectabilityResult",
    "calibration_diagnostics",
    "fit_detectability",
    "transform_covariates",
    "tu_folds",
]

log = logging.getLogger(__name__)

C_GRID: tuple[float, ...] = (0.01, 0.1, 1.0, 10.0)
INNER_K = 3
MIN_CLASS = 50
CLIP = 1e-6
SECTION = "detectability"
OVERALL = "all"

Array = npt.NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class DetectabilityResult:
    p_detect: pd.Series  # index locus_tag; NaN for non-eligible genes
    constant: bool  # detectability_constant
    fold_of_tu: Mapping[str, int]
    chosen_C: tuple[float, ...]  # one per outer fold
    diagnostics: pd.DataFrame  # detection_summary rows, section "detectability"


def _row(group_type: str, group: object, metric: str, value: float) -> dict[str, object]:
    return {
        "section": SECTION,
        "group_type": group_type,
        "group": str(group),
        "metric": metric,
        "value": float(value),
    }


def transform_covariates(cov: pd.DataFrame, names: Sequence[str]) -> pd.DataFrame:
    """log1p(n_observable_peptides), log(length_aa), others raw; NaN/inf → transformed median."""
    out: dict[str, Array] = {}
    for name in names:
        x = pd.to_numeric(cov[name], errors="coerce").to_numpy(dtype=float, na_value=np.nan)
        with np.errstate(divide="ignore", invalid="ignore"):
            if name == "n_observable_peptides":
                x = np.log1p(x)
            elif name == "length_aa":
                x = np.log(x)
        x = np.where(np.isfinite(x), x, np.nan)
        med = float(np.nanmedian(x)) if np.isfinite(x).any() else 0.0
        out[name] = np.where(np.isnan(x), med, x)
    return pd.DataFrame(out, index=cov.index, columns=list(names))


def tu_folds(tu_ids: Sequence[str], k: int, rng: np.random.Generator) -> dict[str, int]:
    """Sorted unique TU IDs permuted by ``rng``; fold = position mod k."""
    tus = sorted({str(t) for t in tu_ids})
    order = rng.permutation(len(tus))
    return {tus[int(i)]: pos % k for pos, i in enumerate(order)}


def _recalibration(d: Array, p: Array) -> tuple[float, float]:
    """Unpenalised MLE (intercept, slope) of d ~ a + b·logit(clip(p)); NaN if not identifiable."""
    if not 0 < d.sum() < d.size:
        return float("nan"), float("nan")
    z = logit(np.clip(p, CLIP, 1 - CLIP))
    x = np.column_stack([np.ones_like(z), z])
    beta = np.array([0.0, 1.0])
    for _ in range(100):
        mu = expit(x @ beta)
        w = mu * (1 - mu)
        try:  # separated data: the MLE diverges; keep the last finite iterate
            step = np.linalg.solve(x.T @ (x * w[:, None]), x.T @ (d - mu))
        except np.linalg.LinAlgError:
            break
        if not np.all(np.isfinite(beta + step)):
            break
        beta = beta + step
        if np.max(np.abs(step)) < 1e-12:
            break
    return float(beta[0]), float(beta[1])


def calibration_diagnostics(d: npt.ArrayLike, p: npt.ArrayLike, n_bins: int = 10) -> pd.DataFrame:
    """Equal-count bins of p (stable rank order) plus Brier score and recalibration fit."""
    dd = np.asarray(d, dtype=float)
    pp = np.asarray(p, dtype=float)
    rows = []
    order = np.argsort(pp, kind="stable")
    for i, b in enumerate(np.array_split(order, n_bins), start=1):
        rows.append(_row("bin", i, "mean_p", float(pp[b].mean()) if b.size else np.nan))
        rows.append(_row("bin", i, "observed_rate", float(dd[b].mean()) if b.size else np.nan))
        rows.append(_row("bin", i, "n", b.size))
    a, s = _recalibration(dd, pp)
    rows.append(_row("overall", OVERALL, "brier", float(np.mean((pp - dd) ** 2))))
    rows.append(_row("overall", OVERALL, "calibration_intercept", a))
    rows.append(_row("overall", OVERALL, "calibration_slope", s))
    rows.append(_row("overall", OVERALL, "n_detected", dd.sum()))
    rows.append(_row("overall", OVERALL, "n_eligible", dd.size))
    return pd.DataFrame(rows)


def _pipeline(cols: Sequence[str], raw: Sequence[str], n_knots: int, c: float) -> Pipeline:
    parts = [
        (
            col,
            "passthrough"
            if col in raw
            else SplineTransformer(
                n_knots=n_knots, degree=3, knots="quantile", extrapolation="linear"
            ),
            [col],
        )
        for col in cols
    ]
    # l1_ratio=0 is penalty="l2" (the ``penalty`` argument is deprecated since sklearn 1.8).
    lr = LogisticRegression(C=c, l1_ratio=0.0, solver="lbfgs", max_iter=2000, tol=1e-6)
    return Pipeline([("features", ColumnTransformer(parts)), ("logistic", lr)])


def _fit_predict(
    x: pd.DataFrame,
    y: npt.NDArray[np.bool_],
    train: npt.NDArray[np.bool_],
    cols: Sequence[str],
    raw: Sequence[str],
    n_knots: int,
    c: float,
) -> tuple[Array, bool]:
    """Predictions for the rows not in ``train``; second value True if the fit was degenerate."""
    ytr = y[train]
    if ytr.all() or not ytr.any():
        return np.full(int((~train).sum()), float(ytr.mean())), True
    model = _pipeline(cols, raw, n_knots, c).fit(x[train], ytr)
    return np.asarray(model.predict_proba(x[~train])[:, 1], dtype=float), False


def _log_loss(y: npt.NDArray[np.bool_], p: Array) -> float:
    q = np.clip(p, 1e-15, 1 - 1e-15)
    return float(-np.mean(np.where(y, np.log(q), np.log1p(-q))))


def _choose_c(
    x: pd.DataFrame,
    y: npt.NDArray[np.bool_],
    tu: npt.NDArray[np.str_],
    cols: Sequence[str],
    raw: Sequence[str],
    n_knots: int,
    rng: np.random.Generator,
) -> float:
    folds = tu_folds(list(tu), INNER_K, rng)
    f = np.array([folds[t] for t in tu])
    scores = {}
    for c in C_GRID:
        losses = []
        for j in range(INNER_K):
            te = f == j
            if te.any():
                p, _ = _fit_predict(x, y, ~te, cols, raw, n_knots, c)
                losses.append(_log_loss(y[te], p))
        scores[c] = float(np.mean(losses))
    best = min(scores.values())
    return min(c for c, s in scores.items() if s == best)


def _constant(
    detected: pd.Series,
    elig: npt.NDArray[np.bool_],
    y: npt.NDArray[np.bool_],
    extra: list[dict[str, object]],
) -> DetectabilityResult:
    rate = float(y.mean()) if y.size else float("nan")
    p = pd.Series(np.where(elig, rate, np.nan), index=detected.index, dtype=float, name="p_detect")
    extra.append(_row("overall", OVERALL, "detectability_constant", 1.0))
    diag = pd.concat(
        [calibration_diagnostics(y, np.full(y.size, rate)), pd.DataFrame(extra)], ignore_index=True
    )
    return DetectabilityResult(p, True, {}, (), diag)


def fit_detectability(
    detected: pd.Series,
    cov: pd.DataFrame,
    tu_of_gene: pd.Series,
    cfg: DetectabilityConfig,
    rng: np.random.Generator,
) -> DetectabilityResult:
    """Cross-fitted ``p_detect`` (index of ``detected``; ``<NA>`` rows are not eligible → NaN)."""
    elig = detected.notna().to_numpy()
    tags = detected.index[elig]
    # canonical row order: results do not depend on the input gene order
    tags = tags[np.argsort(tags.astype(str).to_numpy(), kind="stable")]
    y = detected.loc[tags].astype(bool).to_numpy()
    extra: list[dict[str, object]] = []
    n_det, n_und = int(y.sum()), int((~y).sum())
    if n_det < MIN_CLASS or n_und < MIN_CLASS:
        return _constant(detected, elig, detected[elig].astype(bool).to_numpy(), extra)
    tu_raw = tu_of_gene.reindex(tags)
    tu = np.where(tu_raw.isna(), tags.astype(str), tu_raw.astype(str)).astype(str)
    k = int(cfg.folds)
    n_tus = len(set(tu))
    if n_tus < k:
        log.warning("detectability: too_few_tus (%d TUs < %d folds); constant p_detect", n_tus, k)
        extra.append(_row("overall", OVERALL, "too_few_tus", n_tus))
        return _constant(detected, elig, detected[elig].astype(bool).to_numpy(), extra)

    x_all = transform_covariates(cov.reindex(tags), list(cfg.covariates))
    cols = [c for c in x_all.columns if float(np.std(x_all[c].to_numpy())) > 0]
    for c in x_all.columns:
        if c not in cols:
            log.warning("detectability: covariate %s has zero variance; dropped", c)
            extra.append(_row("covariate", c, "dropped_zero_variance", 1.0))
    x = x_all[cols]
    raw = [c for c in cols if c == "tm_segments" and x[c].nunique() < cfg.spline_knots]

    folds = tu_folds(list(tu), k, rng)
    f = np.array([folds[t] for t in tu])
    p = np.full(len(tags), np.nan)
    chosen: list[float] = []
    # single-threaded BLAS: the small L-BFGS problems are ~10x slower with thread hand-offs
    with warnings.catch_warnings(record=True) as caught, threadpool_limits(1):
        warnings.simplefilter("always", ConvergenceWarning)
        for j in range(k):
            te = f == j
            c = _choose_c(x[~te], y[~te], tu[~te], cols, raw, cfg.spline_knots, rng)
            chosen.append(c)
            p[te], degenerate = _fit_predict(x, y, ~te, cols, raw, cfg.spline_knots, c)
            if degenerate:
                extra.append(_row("fold", j, "degenerate_fold", 1.0))
    n_conv = sum(issubclass(w.category, ConvergenceWarning) for w in caught)
    for w in caught:
        if not issubclass(w.category, ConvergenceWarning):
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)
    if n_conv:
        extra.append(_row("overall", OVERALL, "convergence_warning", n_conv))

    p_ser = pd.Series(p, index=tags, dtype=float)
    out = p_ser.reindex(detected.index).rename("p_detect")
    order = detected.index[elig]
    diag = pd.concat(
        [
            calibration_diagnostics(detected[elig].astype(bool).to_numpy(), p_ser.loc[order]),
            pd.DataFrame(extra, columns=["section", "group_type", "group", "metric", "value"]),
        ],
        ignore_index=True,
    )
    return DetectabilityResult(out, False, folds, tuple(chosen), diag)
