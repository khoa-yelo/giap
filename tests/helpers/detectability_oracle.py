"""Independent oracle for the M4-13 detectability model (docs/specs/M4-13-detectability.md).

Written from the spec text, not from ``gippa.genetable.detectability``:

- ``tu_folds``: sorted unique TU IDs, permuted by ``rng.permutation``, fold = position mod k;
- covariate transforms: log1p(n_observable_peptides), log(length_aa), the rest raw; NaN → median
  of the transformed column;
- features: per covariate one ``SplineTransformer(n_knots, degree=3, knots="quantile",
  extrapolation="linear")`` fitted on the training rows (the transformer is mandated by G5, so the
  oracle reuses it); ``tm_segments`` with fewer distinct values than knots enters raw; zero-variance
  covariates are dropped;
- the L2 logistic fit is NOT sklearn's: it minimises C·Σ logloss + ½‖w‖² (unpenalised intercept,
  sklearn's documented objective) with scipy L-BFGS-B at tight tolerance;
- C from {0.01, 0.1, 1, 10} by inner 3-fold TU-grouped CV on the training TUs (``tu_folds`` rule,
  same stream, outer folds in order 0..k-1), metric = mean over inner folds of the held-out mean
  log loss, ties → smaller C;
- ``calibration_diagnostics``: stable rank order (ties by input order), ``np.array_split`` into
  equal-count bins; unpenalised logistic regression of d on logit(clip(p, 1e-6, 1 - 1e-6)) by
  Newton-Raphson.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, log_expit, logit

C_GRID: tuple[float, ...] = (0.01, 0.1, 1.0, 10.0)
INNER_K = 3
MIN_CLASS = 50
CLIP = 1e-6


# ----------------------------------------------------------------------------- folds


def tu_folds(tu_ids: Sequence[str], k: int, rng: np.random.Generator) -> dict[str, int]:
    tus = sorted(set(map(str, tu_ids)))
    order = [tus[i] for i in rng.permutation(len(tus))]
    return {t: i % k for i, t in enumerate(order)}


# ----------------------------------------------------------------------------- covariates


def transform_covariates(cov: pd.DataFrame, names: Sequence[str]) -> pd.DataFrame:
    out = {}
    for name in names:
        x = pd.to_numeric(cov[name], errors="coerce").astype(float).to_numpy()
        if name == "n_observable_peptides":
            x = np.log1p(x)
        elif name == "length_aa":
            with np.errstate(divide="ignore", invalid="ignore"):
                x = np.log(x)
        x = np.where(np.isfinite(x), x, np.nan)
        med = np.nanmedian(x) if np.isfinite(x).any() else 0.0
        out[name] = np.where(np.isnan(x), med, x)
    return pd.DataFrame(out, index=cov.index)


def _kept(x: pd.DataFrame) -> list[str]:
    return [c for c in x.columns if np.nanstd(x[c].to_numpy(float)) > 0]


class _Features:
    """Per-column spline (or raw for few-valued ``tm_segments``) fitted on training rows."""

    def __init__(self, cols: Sequence[str], raw: Sequence[str], n_knots: int) -> None:
        self.cols, self.raw, self.n_knots = list(cols), set(raw), n_knots
        self.splines: dict[str, object] = {}

    def fit(self, x: pd.DataFrame) -> _Features:
        from sklearn.preprocessing import SplineTransformer

        for c in self.cols:
            if c not in self.raw:
                st = SplineTransformer(
                    n_knots=self.n_knots, degree=3, knots="quantile", extrapolation="linear"
                )
                st.fit(x[[c]].to_numpy(float))
                self.splines[c] = st
        return self

    def transform(self, x: pd.DataFrame) -> np.ndarray:
        blocks = []
        for c in self.cols:
            v = x[[c]].to_numpy(float)
            blocks.append(v if c in self.raw else self.splines[c].transform(v))  # type: ignore[attr-defined]
        return np.hstack(blocks)


def l2_logistic(x: np.ndarray, y: np.ndarray, c: float) -> tuple[np.ndarray, float]:
    """argmin C·Σ logloss + ½‖w‖² (intercept unpenalised)."""
    s = np.where(y, 1.0, -1.0)

    def f(t: np.ndarray) -> tuple[float, np.ndarray]:
        w, b = t[:-1], t[-1]
        z = x @ w + b
        g = -s * expit(-s * z)
        val = c * float(-log_expit(s * z).sum()) + 0.5 * float(w @ w)
        return val, np.r_[c * (x.T @ g) + w, c * g.sum()]

    o = minimize(
        f,
        np.zeros(x.shape[1] + 1),
        jac=True,
        method="L-BFGS-B",
        options={"ftol": 1e-15, "gtol": 1e-10, "maxiter": 20000},
    )
    return o.x[:-1], float(o.x[-1])


def _predict(feat: _Features, w: np.ndarray, b: float, x: pd.DataFrame) -> np.ndarray:
    out: np.ndarray = expit(feat.transform(x) @ w + b)
    return out


def _fit_predict(
    feat_cols: Sequence[str],
    raw: Sequence[str],
    n_knots: int,
    xtr: pd.DataFrame,
    ytr: np.ndarray,
    xte: pd.DataFrame,
    c: float,
) -> np.ndarray:
    feat = _Features(feat_cols, raw, n_knots).fit(xtr)
    w, b = l2_logistic(feat.transform(xtr), ytr, c)
    return _predict(feat, w, b, xte)


def mean_log_loss(y: np.ndarray, p: np.ndarray) -> float:
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return float(-np.mean(np.where(y, np.log(p), np.log1p(-p))))


def design(
    cov: pd.DataFrame, names: Sequence[str], n_knots: int
) -> tuple[pd.DataFrame, list[str], list[str]]:
    x = transform_covariates(cov, names)
    cols = _kept(x)
    raw = [c for c in cols if c == "tm_segments" and x[c].nunique() < n_knots]
    return x[cols], cols, raw


def inner_scores(
    x: pd.DataFrame,
    y: np.ndarray,
    tu: np.ndarray,
    cols: Sequence[str],
    raw: Sequence[str],
    n_knots: int,
    rng: np.random.Generator,
) -> dict[float, float]:
    folds = tu_folds(list(tu), INNER_K, rng)
    f = np.array([folds[t] for t in tu])
    scores: dict[float, float] = {}
    for c in C_GRID:
        losses = []
        for j in range(INNER_K):
            te = f == j
            p = _fit_predict(cols, raw, n_knots, x[~te], y[~te], x[te], c)
            losses.append(mean_log_loss(y[te], p))
        scores[c] = float(np.mean(losses))
    return scores


def choose_c(scores: Mapping[float, float]) -> float:
    best = min(scores.values())
    return min(c for c, s in scores.items() if s == best)


def cross_fit(
    detected: pd.Series,
    cov: pd.DataFrame,
    tu_of_gene: pd.Series,
    names: Sequence[str],
    n_knots: int,
    k: int,
    rng: np.random.Generator,
    fixed_c: Sequence[float] | None = None,
) -> tuple[pd.Series, dict[str, int], list[float], list[dict[float, float]]]:
    """Cross-fitted p_detect for an all-eligible input (both classes ≥ 50, ≥ k TUs).

    Returns (p_detect, fold_of_tu, chosen C per fold, inner scores per fold). With ``fixed_c``
    the inner search is skipped (its stream draws are still made so later folds agree).
    """
    tags = list(detected.index)
    y = detected.astype(bool).to_numpy()
    tu = tu_of_gene.loc[tags].astype(str).to_numpy()
    x, cols, raw = design(cov.loc[tags], names, n_knots)
    folds = tu_folds(list(tu), k, rng)
    f = np.array([folds[t] for t in tu])
    p = np.full(len(tags), np.nan)
    chosen: list[float] = []
    all_scores: list[dict[float, float]] = []
    for j in range(k):
        te = f == j
        xtr, ytr = x[~te], y[~te]
        if fixed_c is None:
            sc = inner_scores(xtr, ytr, tu[~te], cols, raw, n_knots, rng)
            c = choose_c(sc)
        else:
            tu_folds(list(tu[~te]), INNER_K, rng)
            sc, c = {}, float(fixed_c[j])
        all_scores.append(sc)
        chosen.append(c)
        p[te] = _fit_predict(cols, raw, n_knots, xtr, ytr, x[te], c)
    return pd.Series(p, index=detected.index), folds, chosen, all_scores


# ----------------------------------------------------------------------------- diagnostics


def bins(p: Sequence[float], n_bins: int) -> list[np.ndarray]:
    order = np.argsort(np.asarray(p, float), kind="stable")
    return [np.asarray(b) for b in np.array_split(order, n_bins)]


def logistic_recalibration(d: Sequence[float], p: Sequence[float]) -> tuple[float, float]:
    """Unpenalised MLE (intercept, slope) of d ~ a + b·logit(clip(p)) by Newton-Raphson."""
    y = np.asarray(d, float)
    z = logit(np.clip(np.asarray(p, float), CLIP, 1 - CLIP))
    x = np.column_stack([np.ones_like(z), z])
    beta = np.array([0.0, 1.0])
    for _ in range(100):
        mu = expit(x @ beta)
        w = mu * (1 - mu)
        step = np.linalg.solve(x.T @ (x * w[:, None]), x.T @ (y - mu))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-13:
            break
    return float(beta[0]), float(beta[1])


def recalibration_se_at_truth(p: Sequence[float]) -> tuple[float, float]:
    """SE of (intercept, slope) from the Fisher information at (0, 1) for fixed p (n large)."""
    q = np.clip(np.asarray(p, float), CLIP, 1 - CLIP)
    z = logit(q)
    w = q * (1 - q)
    info = np.array([[w.sum(), (w * z).sum()], [(w * z).sum(), (w * z * z).sum()]])
    cov = np.linalg.inv(info)
    return math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])


def diagnostics(d: Sequence[float], p: Sequence[float], n_bins: int = 10) -> dict[str, object]:
    dd, pp = np.asarray(d, float), np.asarray(p, float)
    rows = [
        {
            "group": i + 1,
            "mean_p": float(pp[b].mean()),
            "observed_rate": float(dd[b].mean()),
            "n": len(b),
        }
        for i, b in enumerate(bins(pp, n_bins))
    ]
    out: dict[str, object] = {
        "bins": rows,
        "brier": float(np.mean((pp - dd) ** 2)),
        "n_detected": int(dd.sum()),
        "n_eligible": len(dd),
    }
    if 0 < dd.sum() < len(dd):
        a, s = logistic_recalibration(dd, pp)
        out["calibration_intercept"], out["calibration_slope"] = a, s
    return out
