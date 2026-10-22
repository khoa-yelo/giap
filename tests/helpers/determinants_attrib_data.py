"""Inputs for the M8-09 tests (docs/specs/M8-09-determinants-attribution-null.md).

Contract choices resolved by the test author where the spec's interface leaves them open:

* ``tu_size_blocks(tu)`` returns one 2-D ``int64`` array per distinct TU size, ascending by size;
  array ``(k, s)`` holds the **row indices** of the ``k`` TUs of size ``s`` (one TU per row, TUs
  sorted by TU ID, the TU's rows in input order). Groups of one TU are included (they are the
  fixed blocks). Every row index appears exactly once over all arrays.
* ``permute_response(y, tu, rng)`` moves whole TU response vectors among the TUs of one size
  group (a uniform random permutation of the group), gene order inside a block kept; TUs whose
  size is unique are fixed. ``y`` is not modified.
* ``null_p(...)`` runs the models named by the keys of ``observed_mean_r2`` (``lightgbm``,
  ``elastic_net``) and returns ``{model: NullResult}``. ``NullResult`` has at least ``p``,
  ``b_count``, ``n_perm`` (= B), ``observed`` (the observed mean CV R² it was given), ``null_r2``
  (the B null R², in b order) and ``flags`` (contains ``null_degenerate`` when no TU can move).
  Permutation b (1..B) draws from ``default_rng(child[b - 1])`` with
  ``child = SeedSequence([seed, crc32("determinants.null")]).spawn(B)`` (= the seed sequence of
  ``for_step("determinants.null", seed)``), and its R² is the pooled-OOF R² of one CV repeat on
  ``folds_rep1`` — the same refit ``cross_validate`` makes with ``folds_rep1[None, :]``,
  ``tu=tu`` (spec rule 6, decision 16).
* ``shap_blocks(cv, X, blocks, *, extra=None)``: ``X`` is the model matrix ``cross_validate`` saw
  (same rows, same order); ``extra`` (optional) holds rows the models never saw (ribosomal
  genes, decision Q3), predicted by every fold model of every repeat, averaged. The frame is
  indexed by ``X.index`` then ``extra.index`` and has the columns ``shap_<block>`` (``BLOCKS``
  order), ``expected_value`` (mean over the contributing models of the last ``pred_contrib``
  column) and ``predicted_log2`` (the matching mean prediction; on the CV rows it is
  ``cv.oof``). A block with no column is NA. Elastic-net ``CvResult`` -> ``ValueError`` (no SHAP).
* ``attribution_table(features, rows, cv, shap)``: §9.3 columns in order; one row per model row
  plus one per ribosomal gene present in ``shap`` (flag ``cai_reference``); ``cv`` rows are the
  model rows in ``features`` order; ``cv=None, shap=None`` is the < 100-row case (NA predictions,
  ``too_few_rows``). ``flags`` is a ``;``-joined string; ``top_decile`` is boolean;
  ``dominant_block`` a block name or NA.
"""

from __future__ import annotations

import importlib
import zlib
from collections.abc import Sequence
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd

BLOCKS = ("transcription", "translation", "codon", "protein", "context")
SHAP_COLUMNS = tuple(f"shap_{b}" for b in BLOCKS)
ATTRIBUTION_COLUMNS = (
    "locus_tag", "abundance_log2", "predicted_log2", "residual", *SHAP_COLUMNS,
    "dominant_block", "top_decile", "flags",
)  # fmt: skip
NULL_STREAM = "determinants.null"
FLAG = "abundance_not_explained_by_sequence"


def null_mod() -> ModuleType:
    return importlib.import_module("giap.determinants.null")


def attr_mod() -> ModuleType:
    return importlib.import_module("giap.determinants.attribution")


def tu_from_sizes(sizes: Sequence[int], names: Sequence[str] | None = None) -> list[str]:
    names = names if names is not None else [f"TU_{i:04d}" for i in range(len(sizes))]
    return [n for n, s in zip(names, sizes, strict=True) for _ in range(int(s))]


def null_children(seed: int, B: int) -> list[np.random.SeedSequence]:
    return np.random.SeedSequence([seed, zlib.crc32(NULL_STREAM.encode())]).spawn(B)


def as_str(v: Any) -> str | None:
    if v is None or v is pd.NA:
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    return str(v)


def toy_world(
    rng: np.random.Generator, n: int = 200, signal: float = 1.0
) -> tuple[pd.DataFrame, np.ndarray, list[str], dict[str, str]]:
    """Three features in three blocks (spec example 5), TUs of 1-4 genes, y on all three."""
    import determinants_models_data as d

    tu = d.tu_ids(d.mixed_sizes(rng, n, 1, 4))
    X = pd.DataFrame(
        {
            "cai": rng.normal(0.0, 1.0, n),
            "gravy": rng.normal(0.0, 1.0, n),
            "promotech_max": rng.normal(0.0, 1.0, n),
        },
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = signal * (X["cai"] - 0.5 * X["gravy"] + 0.3 * X["promotech_max"]).to_numpy()
    y = y + rng.normal(0.0, 0.5, n)
    blocks = {"cai": "codon", "gravy": "protein", "promotech_max": "transcription"}
    return X, np.asarray(y, dtype=float), tu, blocks
