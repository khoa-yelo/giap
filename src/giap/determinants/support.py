"""Plot support of the determinants step: SHAP beeswarm rows, partial dependence and the
permutation null (V15(b,c), V18; docs/specs/M9-17-plot-support-analysis.md, ADR-0043).

Pure. ``shap_rows`` cuts the per-column SHAP matrix (``attribution.shap_features``) to the
``SHAP_TOP`` columns by mean |SHAP| over the model rows (ties by name). ``pd_rows``: per PD
feature with a value in some model row, the grid of 5 %, ..., 95 % (numpy linear) quantiles plus
min and max on the feature's own scale, rounded as written (%.6g, §9.2: a value near a tree split
gives the curve at the plotted x); per booster the mean prediction over the model rows with
the column set to the grid value (ln for ``LOG_FEATURES``); ``pd_mean`` the mean over boosters,
``pd_low``/``pd_high`` their min/max, a spread, not a CI. LightGBM only (lead decision 3).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from giap.determinants.features import LOG_FEATURES
from giap.determinants.models import MODELS, CvResult
from giap.determinants.null import NullResult

__all__ = ["PD_FEATURES", "SHAP_TOP", "null_rows", "pd_rows", "shap_rows"]

SHAP_TOP = 20
PD_FEATURES: tuple[str, ...] = ("cai", "tai", "tir_ostir", "promotech_max")
PD_QUANTILES = np.arange(5, 100, 5) / 100.0


def shap_rows(per: pd.DataFrame, X: pd.DataFrame, tags: Sequence[str]) -> pd.DataFrame:
    """``determinants_shap``: ``per`` and ``X`` hold the model rows (same order) of ``tags``."""
    mean = per.abs().mean(axis=0)
    top = sorted(map(str, per.columns), key=lambda c: (-float(mean[c]), c))[:SHAP_TOP]
    n = len(tags)
    return pd.DataFrame({
        "locus_tag": np.tile(np.asarray(tags, dtype=object), len(top)),
        "feature": np.repeat(np.asarray(top, dtype=object), n),
        "feature_value": np.concatenate([X[c].to_numpy(float) for c in top]) if top else [],
        "shap": np.concatenate([per[c].to_numpy(float) for c in top]) if top else [],
    })  # fmt: skip


def pd_rows(cv: CvResult, X: pd.DataFrame, features: pd.DataFrame) -> pd.DataFrame:
    """``determinants_pd`` of the LightGBM ``cv`` over the model rows ``X`` (``features``: the
    same rows of the features frame, own scale)."""
    base = [(bst, X[[str(c) for c in bst.feature_name_]].astype(float))  # type: ignore[attr-defined]
            for bst in cv.boosters]  # fmt: skip
    rows: list[tuple[str, float, float, float, float]] = []
    for feat in PD_FEATURES:
        v = pd.to_numeric(features[feat], errors="coerce").to_numpy(float, na_value=np.nan)
        v = v[np.isfinite(v)]
        if not v.size:
            continue  # e.g. layer_missing:ostir
        grid = np.sort(np.concatenate([[v.min()], np.quantile(v, PD_QUANTILES), [v.max()]]))
        for g in (float(f"{q:.6g}") for q in grid):  # the written (§9.2 %.6g) grid value
            x = (np.log(g) if g > 0 else np.nan) if feat in LOG_FEATURES else g
            per = []
            for bst, Z in base:
                if feat in Z.columns:
                    Z = Z.assign(**{feat: x})
                per.append(float(np.mean(bst.predict(Z))))  # type: ignore[attr-defined]
            rows.append((feat, g, float(np.mean(per)), min(per), max(per)))
    return pd.DataFrame(rows, columns=["feature", "grid_value", "pd_mean", "pd_low", "pd_high"])


def null_rows(nulls: Mapping[str, NullResult], observed: Mapping[str, float]) -> pd.DataFrame:
    """``null_determinants``: per model (``MODELS`` order) the B null R² in b order and the
    model table's ``cv_r2_mean`` (``observed``)."""
    parts = [
        pd.DataFrame({"model": m, "permutation": np.arange(1, nulls[m].null_r2.size + 1),
                      "cv_r2": nulls[m].null_r2, "observed_cv_r2_mean": observed[m]})
        for m in MODELS if m in nulls
    ]  # fmt: skip
    return pd.concat(parts, ignore_index=True)
