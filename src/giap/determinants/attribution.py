"""LightGBM SHAP block attribution, top decile, dominant block and the residual flag
(GIAP_SPEC.md §7.8.2 steps 5-6, §7.8.3, §9.3 ``determinants_attribution.tsv``, §9.4;
docs/specs/M8-09-determinants-attribution-null.md rules 1-4, decisions Q2-Q3).

SHAP values are LightGBM's own ``pred_contrib=True`` (no ``shap`` package); per gene they are
averaged over the repeats' held-out fold models and summed by block, so the block sums plus the
expected value equal the out-of-fold prediction. ADR-0033: the elastic-net ΔR² is the primary
attribution; SHAP stays a LightGBM output. Pure functions, no file I/O.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.determinants.features import BLOCKS as BLOCK_COLUMNS
from giap.determinants.models import CvResult

__all__ = ["ATTRIBUTION_COLUMNS", "BLOCKS", "SHAP_COLUMNS", "attribution_table", "shap_blocks"]

F64 = npt.NDArray[np.float64]
B8 = npt.NDArray[np.bool_]

BLOCKS: tuple[str, ...] = tuple(BLOCK_COLUMNS)
SHAP_COLUMNS: tuple[str, ...] = tuple(f"shap_{b}" for b in BLOCKS)
ATTRIBUTION_COLUMNS: tuple[str, ...] = (
    "locus_tag", "abundance_log2", "predicted_log2", "residual", *SHAP_COLUMNS,
    "dominant_block", "top_decile", "flags",
)  # fmt: skip
UNEXPLAINED = "abundance_not_explained_by_sequence"
TOP_QUANTILE = 0.9
RESIDUAL_SDS = 2.0


# ---------------------------------------------------------------------------- SHAP


def _contrib(booster: Any, X: pd.DataFrame, block_of: Mapping[str, int]) -> tuple[F64, F64]:
    """(rows, len(BLOCKS)) block sums of ``pred_contrib`` and the expected value per row."""
    names = [str(c) for c in booster.feature_name_]
    c = np.asarray(booster.predict(X[names].astype(float), pred_contrib=True), dtype=float)
    sums = np.zeros((len(X), len(BLOCKS)))
    for j, name in enumerate(names):
        sums[:, block_of[name]] += c[:, j]
    return sums, c[:, -1]


def _block_index(boosters: Sequence[Any], blocks: Mapping[str, str]) -> tuple[dict[str, int], B8]:
    """Column -> block position for every column a booster uses; which blocks have a column."""
    block_of: dict[str, int] = {}
    for bst in boosters:
        for name in map(str, bst.feature_name_):
            b = blocks.get(name)
            if b not in BLOCKS:
                raise ValueError(f"column {name!r} has no attribution block (got {b!r})")
            block_of[name] = BLOCKS.index(b)
    present = np.zeros(len(BLOCKS), dtype=bool)
    present[list(set(block_of.values()))] = True
    return block_of, present


def shap_blocks(
    cv: CvResult, X: pd.DataFrame, blocks: Mapping[str, str], *, extra: pd.DataFrame | None = None
) -> pd.DataFrame:
    """Rule 1: per-gene LightGBM SHAP sums by block.

    ``X`` is the matrix ``cross_validate`` saw (same rows and order); booster ``(r, k)`` is
    ``cv.boosters[r * K + k]``. A CV row gets the mean over repeats of the model that held it
    out; an ``extra`` row (genes the models never saw: ribosomal, decision Q3) the mean over every
    fold model of every repeat. Columns ``shap_<block>`` (``BLOCKS`` order, NaN for a block with
    no column), ``expected_value`` and ``predicted_log2`` (``cv.oof`` on the CV rows); index
    ``X.index`` then ``extra.index``.

    Raises:
        ValueError: for an elastic-net (or unfitted) result, a booster count other than
            repeats x folds, or a model column without an attribution block.
    """
    if cv.model != "lightgbm":
        raise ValueError(f"SHAP is a LightGBM output; got a {cv.model} result")
    folds = np.asarray(cv.folds, dtype=np.int64)
    R, n = folds.shape
    if n != len(X):
        raise ValueError(f"X has {len(X)} rows, the CV result {n}")
    labels = [np.unique(folds[r]).tolist() for r in range(R)]
    if not cv.boosters or len(cv.boosters) != sum(len(k) for k in labels):
        raise ValueError("the CV result has no booster for every repeat x fold")
    block_of, present = _block_index(cv.boosters, blocks)
    S = np.zeros((n, len(BLOCKS)))
    E = np.zeros(n)
    i = 0
    for r in range(R):
        for k in labels[r]:
            test = folds[r] == k
            s, e = _contrib(cv.boosters[i], X.loc[test], block_of)
            S[test] += s / R
            E[test] += e / R
            i += 1
    out = _shap_frame(S, E, np.asarray(cv.oof, dtype=float), X.index, present)
    if extra is None or len(extra) == 0:
        return out
    boosters: Sequence[Any] = cv.boosters
    m = len(boosters)
    S2 = np.zeros((len(extra), len(BLOCKS)))
    E2 = np.zeros(len(extra))
    P2 = np.zeros(len(extra))
    for bst in boosters:
        s, e = _contrib(bst, extra, block_of)
        S2 += s / m
        E2 += e / m
        P2 += np.asarray(bst.predict(extra[list(bst.feature_name_)].astype(float))) / m
    return pd.concat([out, _shap_frame(S2, E2, P2, extra.index, present)])


def _shap_frame(S: F64, E: F64, pred: F64, index: pd.Index, present: B8) -> pd.DataFrame:
    S = S.copy()
    S[:, ~present] = np.nan
    df = pd.DataFrame(S, columns=list(SHAP_COLUMNS), index=index)
    df["expected_value"] = E
    df["predicted_log2"] = pred
    return df


# ---------------------------------------------------------------------------- table


def _flag(s: pd.Series) -> B8:
    return np.asarray(s.astype("boolean").fillna(False).to_numpy(dtype=bool))


def _dominant(sums: F64) -> str | None:
    """Rule 3: the block with the largest positive SHAP sum (NaN ignored, ties by ``BLOCKS``)."""
    v = np.where(np.isnan(sums), -np.inf, sums)
    j = int(np.argmax(v))
    return BLOCKS[j] if v[j] > 0 else None


def attribution_table(
    features: pd.DataFrame, rows: B8, cv: CvResult | None, shap: pd.DataFrame | None
) -> pd.DataFrame:
    """Rules 2-4: ``determinants_attribution.tsv`` (§9.3 columns, in order).

    One row per model row (``rows``; ``cv`` holds them in ``features`` order) plus one per
    detected, non-imputed ribosomal gene present in ``shap`` (flag ``cai_reference``; never in
    the top decile, the residual SD or the residual flag). ``predicted_log2`` = ``cv.oof``
    (``shap["predicted_log2"]`` for ribosomal rows); ``residual`` = observed - predicted.
    ``top_decile``: model rows with ``abundance_log2`` ≥ the linear 90th percentile;
    ``dominant_block`` on top-decile rows only. ``|residual|`` > 2 SD (ddof 1, model rows; SD
    0 -> none) flags ``abundance_not_explained_by_sequence``. ``cv=None, shap=None``: fewer than
    100 model rows -> NA predictions and ``too_few_rows``."""
    rows = np.asarray(rows, dtype=bool)
    tags = features["locus_tag"].astype(str).to_numpy()
    ab = pd.to_numeric(features["abundance_log2"], errors="coerce").to_numpy(float, na_value=np.nan)
    ribo = np.zeros(len(features), dtype=bool)
    if shap is not None:
        ribo = (
            _flag(features["ribosomal_protein"]) & _flag(features["detected"]) & np.isfinite(ab)
            & ~_flag(features["abundance_imputed"]) & np.isin(tags, shap.index.astype(str))
        )  # fmt: skip
    keep = rows | ribo
    is_model = rows[keep]
    t_tags, t_ab = tags[keep], ab[keep]
    n = int(keep.sum())
    pred = np.full(n, math.nan)
    S = np.full((n, len(BLOCKS)), math.nan)
    if shap is not None:
        s = shap.set_axis(shap.index.astype(str)).reindex(t_tags)
        S = np.array(s[list(SHAP_COLUMNS)], dtype=float)
        pred = np.array(s["predicted_log2"], dtype=float)
    if cv is not None:
        oof = np.asarray(cv.oof, dtype=float)
        if oof.size != int(is_model.sum()):
            raise ValueError(f"cv has {oof.size} rows, features {int(is_model.sum())} model rows")
        pred[is_model] = oof
    resid = t_ab - pred
    top = np.zeros(n, dtype=bool)
    if is_model.any():
        top = is_model & (t_ab >= float(np.quantile(t_ab[is_model], TOP_QUANTILE)))
    r_model = resid[is_model]
    sd = float(np.std(r_model, ddof=1)) if r_model.size > 1 else math.nan
    unexplained = is_model & (np.abs(resid) > RESIDUAL_SDS * sd) if sd > 0 else np.zeros(n, bool)
    dominant: list[str | None] = [None] * n
    flags: list[str] = []
    for i in range(n):
        if top[i] and cv is not None:
            dominant[i] = _dominant(S[i])
        f: list[str] = []
        if cv is None and is_model[i]:
            f.append("too_few_rows")
        elif not math.isfinite(pred[i]):
            f.append("unestimated")
        if unexplained[i]:
            f.append(UNEXPLAINED)
        if not is_model[i]:
            f.append("cai_reference")
        flags.append(";".join(f))
    out = pd.DataFrame({"locus_tag": t_tags, "abundance_log2": t_ab, "predicted_log2": pred,
                        "residual": resid})  # fmt: skip
    for j, c in enumerate(SHAP_COLUMNS):
        out[c] = S[:, j]
    out["dominant_block"] = pd.array(dominant, dtype="string")
    out["top_decile"] = top
    out["flags"] = flags
    return out[list(ATTRIBUTION_COLUMNS)]
