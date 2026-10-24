"""The ``enriched_categories`` sensitivity family (issue M8-16; GIAP_SPEC.md §8.6 "enriched
categories (q ≤ 0.05)", "Sensitivity reruns use B = 199"; docs/specs/M8-16-sensitivity-function-
determinants.md item 1, lead decisions 1 and 4).

A family of the M4-19 framework (``giap.context.sensitivity.FAMILIES``): the flagged IDs are
``{collection}:{term}:{kind}`` with q ≤ ``stats.alpha.enrichment``, kind ``gsea``, ``ora_over``
or ``ora_under``. Each call reruns the M8-04 detection ORA (analytical; it only changes when the
detection calls do, i.e. under ``single_peptide_detected``) and GSEA at ``n_perm_simple = 199``
and ``sample_size = 101``, whatever ``function.gsea`` says. Pure: randomness only from ``rng``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from giap.config.model import GiapConfig
from giap.function.collections import GeneSets
from giap.function.gsea import SIMPLE_ENGINE, GseaEngine, gsea_table
from giap.function.ora import detection_enrichment

__all__ = ["N_PERM", "SAMPLE_SIZE", "enriched_categories_family"]

#: §8.6 "Sensitivity reruns use B = 199" (``giap.context.sensitivity.SENSITIVITY_B``)
N_PERM = 199
#: spec item 1: the multilevel stage at ``sampleSize = 101``
SAMPLE_SIZE = 101


def _significant(res: pd.DataFrame, alpha: float) -> pd.DataFrame:
    q = res["q_value"].astype("float64")
    return res[q.notna() & (q <= alpha)]


def enriched_categories_family(
    genes: pd.DataFrame,
    cfg: GiapConfig,
    rng: np.random.Generator,
    *,
    collections: Sequence[GeneSets],
    tus: pd.DataFrame,
    folds: Mapping[str, int] | None = None,
    engine: GseaEngine = SIMPLE_ENGINE,
) -> frozenset[str]:
    """The ORA and GSEA terms with q ≤ ``stats.alpha.enrichment`` on ``genes``.

    ``collections`` are the step's ``build_collections`` output on the primary gene table;
    ``folds`` the M4-13 TU folds (``None``: no recalibration, as the function step without
    folds); ``engine`` runs GSEA (the step binds R fgsea, ADR-0038). ``tus`` (C5) is bound for
    the framework; the TU of each gene is read from ``genes["tu_id"]``, as the function step does.
    """
    del tus  # the gene table carries the TU membership
    if not collections:
        return frozenset()
    fn, alpha = cfg.function, cfg.stats.alpha.enrichment
    colls = list(collections)
    ora = _significant(
        detection_enrichment(genes, colls, (fn.ora_size[0], fn.ora_size[1]), alpha,
                             folds=dict(folds) if folds else None),
        alpha,
    )  # fmt: skip
    params = fn.gsea.model_copy(update={"n_perm_simple": N_PERM, "sample_size": SAMPLE_SIZE})
    gsea = _significant(
        gsea_table(genes, colls, (fn.gsea_size[0], fn.gsea_size[1]), params, alpha, rng,
                   engine=engine),
        alpha,
    )  # fmt: skip
    ids = {f"{c}:{t}:ora_{d}" for c, t, d in
           zip(ora["collection"], ora["term"], ora["direction"], strict=True)}  # fmt: skip
    ids |= {f"{c}:{t}:gsea" for c, t in zip(gsea["collection"], gsea["term"], strict=True)}
    return frozenset(ids)
