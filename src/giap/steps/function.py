"""Step ``function`` (GIAP_SPEC.md §6.2, §6.3, §7.7; docs/specs/M8-06-function-step-and-run.md):
the gene table's function columns and the GO/KEGG/COG setup resources →
``tables/category_summary.tsv``, ``enrichment_detection.tsv`` and ``enrichment_gsea.tsv``.

``run_function`` is the pure chain of spec item 2: eligibility, ``build_collections``, the
category summary (stream ``function.bootstrap``), the TU-level ORA (deterministic; M4-13's TU
folds from ``tu_folds.parquet`` passed as ``folds=``, ``None`` when that file is header-only,
lead decision 3) and GSEA (stream ``function.gsea``, C1 topology for ``genomic_cluster``). The
streams are sub-streams of ``giap.rng.for_step(<name>, seed)``, not the step's own, so a draw
added to one part never moves the other (the ``context_stats.<part>`` precedent). With no
eligible gene the three tables are header-only and ``no_eligible_genes`` is a run flag.

Each run flag is logged once as the WARNING ``function: run flag <flag>``; then the §6.7 events
``function.run_flags`` ``{"flags": [...]}`` and ``function.collections`` (per built collection:
``n_terms``, ``coverage``, ``n_unknown_ids``). A resource that does not parse raises (exit 4).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from giap.config import GiapConfig
from giap.function import gsea as _gsea
from giap.function import ora as _ora
from giap.function import summary as _summary
from giap.function.collections import GeneSets, Resources, build_collections
from giap.function.gsea import gsea_table
from giap.function.ora import detection_enrichment
from giap.function.summary import category_summary
from giap.genetable.eligibility import eligible
from giap.io.function_resources import read_resources
from giap.io.function_tables import (
    write_category_summary_tsv,
    write_enrichment_detection_tsv,
    write_enrichment_gsea_tsv,
)
from giap.io.genes_tsv import read_genes_parquet
from giap.rng import for_step
from giap.setupdb.installed import FUNCTION_RESOURCES
from giap.steps.common import one_path
from giap.steps.registry import StepContext, step

__all__ = ["read_folds", "run_function"]

NO_ELIGIBLE = "no_eligible_genes"
_FLAGS_EVENT = "function.run_flags"
_COLL_EVENT = "function.collections"


def _empty(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([], dtype=object) for c in columns})


def _union(*parts: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(f for p in parts for f in p))


def run_function(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    folds: Mapping[str, int],
    res: Resources,
    cfg: GiapConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """(category_summary, enrichment_detection, enrichment_gsea, run_flags). Pure.

    ``folds`` empty means no fitted folds (``folds=None`` to the ORA, ``tu_recal_skipped``).

    Raises:
        ValueError: an eligible gene without a finite ``p_detect`` or a ``tu_id``.
    """
    summary, detection, gsea, flags, _ = _run(genes, replicons, folds, res, cfg)
    return summary, detection, gsea, flags


def _run(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    folds: Mapping[str, int],
    res: Resources,
    cfg: GiapConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str], list[GeneSets]]:
    el = eligible(genes)
    colls, flags = build_collections(genes, el, res)
    if not bool(el.any()):
        return (_empty(_summary.COLUMNS), _empty(_ora.COLUMNS), _empty(_gsea.COLUMNS),
                _union(flags, [NO_ELIGIBLE]), colls)  # fmt: skip
    fn, alpha = cfg.function, cfg.stats.alpha.enrichment
    summary = category_summary(genes, colls, cfg.stats.bootstrap,
                               for_step("function.bootstrap", cfg.seed))  # fmt: skip
    detection = detection_enrichment(genes, colls, (fn.ora_size[0], fn.ora_size[1]), alpha,
                                     folds=dict(folds) if folds else None)  # fmt: skip
    gsea = gsea_table(genes, colls, (fn.gsea_size[0], fn.gsea_size[1]), fn.gsea, alpha,
                      for_step("function.gsea", cfg.seed), replicons)  # fmt: skip
    run_flags = _union(flags, detection.attrs.get("run_flags", []),
                       gsea.attrs.get("run_flags", []))  # fmt: skip
    return summary, detection, gsea, run_flags, colls


def read_folds(path: Path) -> dict[str, int]:
    """``genetable/tu_folds.parquet`` → ``{tu_id: fold}`` (``{}`` when header-only)."""
    df = pd.read_parquet(path)
    return {str(t): int(f) for t, f in zip(df["tu_id"], df["fold"], strict=True)}


@step(
    "function",
    module="function",
    inputs=["gene_table", "replicons", "tu_folds"],
    optional_inputs=list(FUNCTION_RESOURCES),
    outputs=["category_summary", "enrichment_detection", "enrichment_gsea"],
)
def function_step(ctx: StepContext) -> None:
    """Category summary, detection ORA and GSEA tables (Module 7)."""
    genes = read_genes_parquet(one_path(ctx.inputs["gene_table"]))
    replicons = pd.read_parquet(one_path(ctx.inputs["replicons"]))
    folds = read_folds(one_path(ctx.inputs["tu_folds"]))
    paths = {k: one_path(v) for k in FUNCTION_RESOURCES if (v := ctx.inputs.get(k)) is not None}
    res = read_resources(paths)
    summary, detection, gsea, flags, colls = _run(genes, replicons, folds, res, ctx.cfg)
    for flag in flags:
        ctx.log.warning("function: run flag %s", flag)
    ctx.log.info("function: %d run flags", len(flags),
                 extra={"event": _FLAGS_EVENT, "data": {"flags": list(flags)}})  # fmt: skip
    data = {c.collection: {"n_terms": len(c.sets), "coverage": float(c.coverage),
                           "n_unknown_ids": int(c.n_unknown_ids)} for c in colls}  # fmt: skip
    ctx.log.info("function: %d collections built", len(data),
                 extra={"event": _COLL_EVENT, "data": data})  # fmt: skip
    write_category_summary_tsv(ctx.outputs["category_summary"], summary)
    write_enrichment_detection_tsv(ctx.outputs["enrichment_detection"], detection)
    write_enrichment_gsea_tsv(ctx.outputs["enrichment_gsea"], gsea)
