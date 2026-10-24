"""The ``determinants`` sensitivity family (issue M8-16; GIAP_SPEC.md §8.6 "the determinants CV
R² (without its permutation null)", "plus ``agreement`` for the dominant determinant block";
docs/specs/M8-16-sensitivity-function-determinants.md item 2, lead decisions 2 and 4).

Per run: the M8-10 chain without its null on the M8-07 features table, with ``detected``,
``abundance_log2`` and ``abundance_imputed`` taken from the run's gene table: model rows, outer
TU folds from ``rng``, LightGBM / elastic-net CV (``determinants.repeats``), LightGBM SHAP and the
attribution table. The result is a ``DeterminantsSensitivity``: the genes flagged
``abundance_not_explained_by_sequence`` (the set behind ``n_*`` / ``jaccard``), the dominant block
of every explained model row, i.e. not flagged (behind ``agreement``; M8-16c, lead decision 6),
and the mean LightGBM CV R² (an event of the step). The dominant block is ``_dominant`` on the
row's LightGBM SHAP block sums; ``determinants_attribution.tsv`` still reports it for the top
decile only. Pure: randomness only from ``rng``; ``cfg.seed`` seeds the models as in the step.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from giap.config.model import GiapConfig
from giap.determinants.attribution import (
    SHAP_COLUMNS,
    UNEXPLAINED,
    _dominant,
    attribution_table,
    shap_blocks,
)
from giap.determinants.features import model_matrix
from giap.determinants.models import cross_validate, enough_data, model_rows, tu_folds

__all__ = ["DeterminantsSensitivity", "determinants_family", "with_gene_calls"]

_CALLS = ("detected", "abundance_log2", "abundance_imputed")


@dataclass(frozen=True)
class DeterminantsSensitivity:
    """One determinants rerun: ``flagged`` locus tags, ``dominant`` block per explained model
    row's locus tag (not flagged; ``None`` without a positive SHAP block sum) and the mean
    LightGBM CV R²."""

    flagged: frozenset[str]
    dominant: dict[str, str | None] = field(default_factory=dict)
    cv_r2: float = float("nan")


def with_gene_calls(features: pd.DataFrame, genes: pd.DataFrame) -> pd.DataFrame:
    """``features`` with ``detected`` / ``abundance_log2`` / ``abundance_imputed`` from ``genes``
    (by locus tag; a feature row without a gene is not detected and has no abundance)."""
    g = genes.drop_duplicates("locus_tag")
    g = g.set_axis(pd.Index(g["locus_tag"].astype(str)))
    tags = features["locus_tag"].astype(str)
    out = features.copy()
    det = g["detected"].astype("boolean").reindex(tags).fillna(False)
    imp = g["abundance_imputed"].astype("boolean").reindex(tags).fillna(False)
    out["detected"] = pd.array(det.to_numpy(bool), dtype="boolean")
    out["abundance_imputed"] = pd.array(imp.to_numpy(bool), dtype="boolean")
    ab = pd.to_numeric(g["abundance_log2"], errors="coerce").astype("float64")
    out["abundance_log2"] = ab.reindex(tags).to_numpy("float64")
    return out


def determinants_family(
    genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator, *, features: pd.DataFrame
) -> DeterminantsSensitivity | None:
    """The rerun on ``genes`` and the M8-07 ``features`` table, or ``None`` below the 100-row /
    25-TU minimum (no abundance included)."""
    f = with_gene_calls(features, genes)
    rows = model_rows(f)
    tu = f.loc[rows, "tu_id"].astype(str).tolist()
    if not enough_data(tu):
        return None
    X_all, blocks = model_matrix(f)
    X = X_all[rows]
    y = pd.to_numeric(f.loc[rows, "abundance_log2"]).to_numpy(float)
    d = cfg.determinants
    cpus = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
    threads = max(1, min(cfg.threads, cpus or 1))
    folds = tu_folds(tu, d.outer_folds, d.repeats, rng)
    lgbm = cross_validate(X, y, folds, blocks, d, cfg.seed, threads, tu=tu)["lightgbm"]
    shap = shap_blocks(lgbm, X, blocks)
    shap.index = pd.Index(f.loc[shap.index, "locus_tag"].astype(str))
    att = attribution_table(f, rows, lgbm, shap)
    flagged = frozenset(
        str(t) for t, fl in zip(att["locus_tag"], att["flags"], strict=True)
        if UNEXPLAINED in str(fl).split(";")
    )  # fmt: skip
    model_tags = set(f.loc[rows, "locus_tag"].astype(str))
    sums = shap.loc[:, list(SHAP_COLUMNS)]
    dominant = {t: _dominant(np.asarray(sums.loc[t], dtype=float))
                for t in att["locus_tag"].astype(str)
                if t in model_tags and t not in flagged}  # fmt: skip
    return DeterminantsSensitivity(flagged, dominant, float(np.mean(lgbm.r2)))
