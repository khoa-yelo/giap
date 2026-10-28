"""Inputs for the M8-08 determinants-model tests (docs/specs/M8-08-determinants-models.md).

Builders:

* ``features(cols, ...)`` -- a ``determinants_features.tsv`` frame (M8-07 ``FEATURE_COLUMNS``,
  ``KINDS`` dtypes): every feature NA except the ones given, so ``model_matrix`` keeps only those;
* ``tu_ids(sizes)`` -- TU IDs ``TU_0000`` ... repeated by gene count;
* ``planted(rng, n, signal)`` -- spec example 1: ``cai`` ~ N(0, 1) in block ``codon``, five
  N(0, 1) noise columns in the other four blocks, y = signal * cai + N(0, 1), TUs of 1-5 genes;
* ``config(fast)`` -- ``DeterminantsConfig``: the §5.2 defaults, or (``fast``) 2 repeats and
  60 trees so the behavioural unit tests stay quick (each test builds folds with ``cfg.repeats``
  and ``cfg.outer_folds``, so the two always agree).

Contract choices resolved by the test author where the spec's interface leaves them open:

* ``univariate`` takes the TSV frame and encodes it with ``model_matrix`` itself (one rho and one
  AUROC row per encoded column); it bootstraps with ``tu_id`` as the TU column; a statistic with
  fewer than 10 usable rows has no finite value;
* the minimum-data rule (decision 18) is ``enough_data(tu) -> bool`` with ``MIN_ROWS = 100`` and
  ``MIN_TUS = 25`` (``tu`` = the TU ID of each model row); when it is false the caller runs no CV
  and ``model_rows_table({}, uni)`` writes one ``cv_r2_mean`` row per model, value NA,
  ``flags`` ``too_few_rows``;
* ``CvResult`` gains ``dropped: tuple[str, ...] = ()`` (lead decision 2): the sorted columns
  dropped as zero-variance or all-NaN in at least one outer training fold;
* ``model_rows_table`` columns are §9.3 plus ``flags`` (index decision 12): ``model``,
  ``feature_or_block``, ``metric``, ``value``, ``ci_low``, ``ci_high``, ``n``, ``n_perm``,
  ``flags``; ``flags`` is a ``;``-joined string, empty or NA when there are none;
* row order is canonicalised on the frame index (``cross_validate``: ``X.index``;
  ``univariate``: ``locus_tag``) and on the TU ID (``tu_folds``).
"""

from __future__ import annotations

import importlib
import math
from collections.abc import Mapping, Sequence
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd

MODULE = "gippa.determinants.models"
TABLE_COLUMNS = (
    "model", "feature_or_block", "metric", "value", "ci_low", "ci_high", "n", "n_perm", "flags",
)  # fmt: skip
MODELS = ("lightgbm", "elastic_net")
METRICS = {"cv_r2_mean", "cv_r2_sd", "null_p", "delta_r2", "spearman_rho", "auroc"}

#: spec example 1: the signal column and the noise columns, with their blocks
PLANTED_BLOCKS: Mapping[str, str] = {
    "cai": "codon",
    "promotech_max": "transcription",
    "sd_dg": "translation",
    "gravy": "protein",
    "pi": "protein",
    "gc_local_dev": "context",
}


def mod() -> ModuleType:
    """``gippa.determinants.models`` (imported lazily so collection never fails)."""
    return importlib.import_module(MODULE)


def feats() -> ModuleType:
    return importlib.import_module("gippa.determinants.features")


def is_na(x: Any) -> bool:
    if x is None or x is pd.NA:
        return True
    try:
        return bool(pd.isna(x))
    except (TypeError, ValueError):
        return False


def flags_of(v: Any) -> set[str]:
    if is_na(v):
        return set()
    return {f for f in str(v).split(";") if f}


def tu_ids(sizes: Sequence[int], prefix: str = "TU_") -> list[str]:
    return [f"{prefix}{i:04d}" for i, s in enumerate(sizes) for _ in range(int(s))]


def mixed_sizes(rng: np.random.Generator, n_genes: int, lo: int = 1, hi: int = 5) -> list[int]:
    """TU sizes in [lo, hi] summing to ``n_genes``."""
    sizes: list[int] = []
    while sum(sizes) < n_genes:
        sizes.append(int(rng.integers(lo, hi + 1)))
    sizes[-1] -= sum(sizes) - n_genes
    if sizes[-1] == 0:
        sizes.pop()
    return sizes


def features(
    cols: Mapping[str, Sequence[Any]],
    *,
    tu: Sequence[str],
    detected: Sequence[bool] | None = None,
    abundance: Sequence[float] | None = None,
    imputed: Sequence[bool] | None = None,
    ribosomal: Sequence[bool] | None = None,
    tags: Sequence[str] | None = None,
) -> pd.DataFrame:
    """A TSV frame (M8-07 ``FEATURE_COLUMNS``, typed) with only ``cols`` non-NA."""
    f = feats()
    n = len(tu)
    ident: dict[str, Any] = {
        "locus_tag": list(tags) if tags is not None else [f"G{i:05d}" for i in range(n)],
        "tu_id": list(tu),
        "detected": [True] * n if detected is None else list(detected),
        "abundance_log2": [0.0] * n if abundance is None else list(abundance),
        "abundance_imputed": [False] * n if imputed is None else list(imputed),
        "ribosomal_protein": [False] * n if ribosomal is None else list(ribosomal),
    }
    data: dict[str, pd.Series] = {}
    for c in f.FEATURE_COLUMNS:
        vals = ident[c] if c in ident else cols.get(c, [None] * n)
        kind = f.KINDS[c]
        s = pd.Series(list(vals), dtype=object)
        if kind == "float":
            data[c] = pd.to_numeric(s, errors="coerce").astype("float64")
        elif kind == "int":
            data[c] = pd.to_numeric(s, errors="coerce").astype("Int64")
        elif kind == "bool":
            data[c] = s.astype("boolean")
        else:
            data[c] = pd.Series([v if isinstance(v, str) else None for v in vals], dtype=object)
    unknown = set(cols) - set(f.FEATURE_COLUMNS)
    if unknown:
        raise KeyError(f"not feature columns: {sorted(unknown)}")
    return pd.DataFrame(data)


def planted(
    rng: np.random.Generator, n: int = 400, signal: float = 2.0, hi: int = 5
) -> tuple[pd.DataFrame, np.ndarray, list[str], dict[str, str]]:
    """Spec example 1 (signal 2 -> true R² 0.8; signal 0 -> the null): X, y, tu, column->block."""
    tu = tu_ids(mixed_sizes(rng, n, 1, hi))
    X = pd.DataFrame(
        {c: rng.normal(0.0, 1.0, n) for c in PLANTED_BLOCKS},
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = signal * X["cai"].to_numpy() + rng.normal(0.0, 1.0, n)
    return X, np.asarray(y, dtype=float), tu, dict(PLANTED_BLOCKS)


def config(fast: bool = True, **lgb: Any) -> Any:
    from gippa.config.model import DeterminantsConfig, LightgbmConfig

    if not fast:
        return DeterminantsConfig()
    params = {"n_estimators": 60, **lgb}
    return DeterminantsConfig(repeats=2, lightgbm=LightgbmConfig(**params))


def folds_for(tu: Sequence[str], cfg: Any, seed: int = 7) -> np.ndarray:
    from gippa.rng import for_step

    return np.asarray(
        mod().tu_folds(list(tu), cfg.outer_folds, cfg.repeats, for_step("determinants.cv", seed))
    )


def partition(fold_row: np.ndarray) -> frozenset[frozenset[int]]:
    """A fold assignment as a set of row sets (fold labels do not matter)."""
    return frozenset(frozenset(np.flatnonzero(fold_row == k).tolist()) for k in np.unique(fold_row))


def rows(tbl: pd.DataFrame, **eq: Any) -> pd.DataFrame:
    m = np.ones(len(tbl), dtype=bool)
    for k, v in eq.items():
        m &= (tbl[k] == v).to_numpy(dtype=bool, na_value=False)
    return tbl.loc[m]


def one(tbl: pd.DataFrame, **eq: Any) -> pd.Series:
    sel = rows(tbl, **eq)
    assert len(sel) == 1, f"expected one row for {eq}, got {len(sel)}"
    return sel.iloc[0]


def canonical(tbl: pd.DataFrame) -> pd.DataFrame:
    """The table sorted on its key columns (for order-free comparison)."""
    return tbl.sort_values(["model", "feature_or_block", "metric"], kind="stable").reset_index(
        drop=True
    )


def finite(x: Any) -> bool:
    return not is_na(x) and math.isfinite(float(x))
