"""Step ``function`` (GIPPA_SPEC.md §6.2, §6.3, §7.7; docs/specs/M8-06-function-step-and-run.md):
the gene table's function columns and the GO/KEGG/COG setup resources →
``tables/category_summary.tsv``, ``enrichment_detection.tsv`` and ``enrichment_gsea.tsv``.

``run_function`` is the pure chain of spec item 2: eligibility, ``build_collections``, the
category summary (stream ``function.bootstrap``), the TU-level ORA (deterministic; M4-13's TU
folds from ``tu_folds.parquet`` passed as ``folds=``, ``None`` when that file is header-only,
lead decision 3) and GSEA (stream ``function.gsea``, C1 topology for ``genomic_cluster``). The
streams are sub-streams of ``gippa.rng.for_step(<name>, seed)``, not the step's own, so a draw
added to one part never moves the other (the ``context_stats.<part>`` precedent). With no
eligible gene the three tables are header-only and ``no_eligible_genes`` is a run flag.

Each run flag is logged once as the WARNING ``function: run flag <flag>``; then the §6.7 events
``function.run_flags`` ``{"flags": [...]}`` and ``function.collections`` (per built collection:
``n_terms``, ``coverage``, ``n_unknown_ids``). A resource that does not parse raises (exit 4).

The step's GSEA engine is R fgsea (``gippa.io.fgsea_r``, ADR-0038), built lazily on the first
engine call, so a run without abundance or without a testable set never needs R. Each call logs
the event ``function.gsea_engine``; an engine error fails the step (exit 4, no table written).
"""

from __future__ import annotations

import logging
import math
import os
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd

from gippa.config import GippaConfig
from gippa.function import gsea as _gsea
from gippa.function import ora as _ora
from gippa.function import summary as _summary
from gippa.function.collections import GeneSets, Resources, build_collections
from gippa.function.gsea import SIMPLE_ENGINE, EngineCall, EngineResult, GseaEngine, gsea_table
from gippa.function.ora import detection_enrichment
from gippa.function.running import gsea_running
from gippa.function.summary import category_summary
from gippa.genetable.eligibility import eligible
from gippa.io.fgsea_r import METHOD as R_METHOD
from gippa.io.fgsea_r import RFgseaEngine, r_fgsea_engine, resolve_launcher
from gippa.io.function_resources import read_resources
from gippa.io.function_tables import (
    write_category_summary_tsv,
    write_enrichment_detection_tsv,
    write_enrichment_gsea_tsv,
)
from gippa.io.genes_tsv import read_genes_parquet
from gippa.io.plot_support import write_support
from gippa.rng import for_step
from gippa.setupdb.installed import FUNCTION_RESOURCES
from gippa.steps.common import one_path
from gippa.steps.registry import StepContext, step

__all__ = ["read_folds", "run_function"]

NO_ELIGIBLE = "no_eligible_genes"
_FLAGS_EVENT = "function.run_flags"
_COLL_EVENT = "function.collections"
_ENGINE_EVENT = "function.gsea_engine"
#: M9-17: optional directory output; ``gsea_running`` goes to its ``support/`` (ADR-0043)
SUPPORT = "function_support"


def _empty(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series([], dtype=object) for c in columns})


def _union(*parts: Sequence[str]) -> list[str]:
    return list(dict.fromkeys(f for p in parts for f in p))


def _finite_or_none(x: float) -> float | None:
    v = float(x)
    return v if math.isfinite(v) else None


def run_function(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    folds: Mapping[str, int],
    res: Resources,
    cfg: GippaConfig,
    *,
    engine: GseaEngine = SIMPLE_ENGINE,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    """(category_summary, enrichment_detection, enrichment_gsea, run_flags). Pure.

    ``folds`` empty means no fitted folds (``folds=None`` to the ORA, ``tu_recal_skipped``).

    Raises:
        ValueError: an eligible gene without a finite ``p_detect`` or a ``tu_id``.
    """
    summary, detection, gsea, flags, _ = _run(genes, replicons, folds, res, cfg, engine)
    return summary, detection, gsea, flags


def _run(
    genes: pd.DataFrame,
    replicons: pd.DataFrame,
    folds: Mapping[str, int],
    res: Resources,
    cfg: GippaConfig,
    engine: GseaEngine,
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
                      for_step("function.gsea", cfg.seed), replicons, engine=engine)  # fmt: skip
    run_flags = _union(flags, detection.attrs.get("run_flags", []),
                       gsea.attrs.get("run_flags", []))  # fmt: skip
    return summary, detection, gsea, run_flags, colls


class _LazyRFgsea:
    """The R engine, built on its first call (§9); logs one ``function.gsea_engine`` event per
    call."""

    method = R_METHOD

    def __init__(self, log: logging.Logger) -> None:
        self._log = log
        self._engine: RFgseaEngine | None = None

    def __call__(self, call: EngineCall) -> EngineResult:
        if self._engine is None:
            self._engine = r_fgsea_engine(resolve_launcher(os.environ), tmp_root=Path.cwd())
        result = self._engine(call)
        prov = self._engine.provenance
        kind = "retest" if all(j.job_id.startswith("retest:") for j in call.jobs) else "main"
        data = {
            "call": kind,
            "jobs": len(call.jobs),
            "r": prov.r if prov else None,
            "fgsea": prov.fgsea if prov else None,
            "launcher": self._engine.launcher.kind,
            "warnings": [{"job_id": j, "message": m} for j, m in (prov.warnings if prov else ())],
        }
        self._log.info("function: R fgsea %s call, %d jobs", kind, len(call.jobs),
                       extra={"event": _ENGINE_EVENT, "data": data})  # fmt: skip
        return result


def read_folds(path: Path) -> dict[str, int]:
    """``genetable/tu_folds.parquet`` → ``{tu_id: fold}`` (``{}`` when header-only)."""
    df = pd.read_parquet(path)
    return {str(t): int(f) for t, f in zip(df["tu_id"], df["fold"], strict=True)}


@step(
    "function",
    module="function",
    inputs=["gene_table", "replicons", "tu_folds"],
    optional_inputs=list(FUNCTION_RESOURCES),
    outputs=["category_summary", "enrichment_detection", "enrichment_gsea", SUPPORT],
    optional_outputs=[SUPPORT],
)
def function_step(ctx: StepContext) -> None:
    """Category summary, detection ORA and GSEA tables (Module 7)."""
    genes = read_genes_parquet(one_path(ctx.inputs["gene_table"]))
    replicons = pd.read_parquet(one_path(ctx.inputs["replicons"]))
    folds = read_folds(one_path(ctx.inputs["tu_folds"]))
    paths = {k: one_path(v) for k in FUNCTION_RESOURCES if (v := ctx.inputs.get(k)) is not None}
    res = read_resources(paths)
    engine = _LazyRFgsea(ctx.log)
    summary, detection, gsea, flags, colls = _run(genes, replicons, folds, res, ctx.cfg, engine)
    for flag in flags:
        ctx.log.warning("function: run flag %s", flag)
    ctx.log.info("function: %d run flags", len(flags),
                 extra={"event": _FLAGS_EVENT, "data": {"flags": list(flags)}})  # fmt: skip
    # strict JSON: a coverage over no eligible gene is NaN, written as null (#181 review N1)
    data = {c.collection: {"n_terms": len(c.sets), "coverage": _finite_or_none(c.coverage),
                           "n_unknown_ids": int(c.n_unknown_ids)} for c in colls}  # fmt: skip
    ctx.log.info("function: %d collections built", len(data),
                 extra={"event": _COLL_EVENT, "data": data})  # fmt: skip
    write_category_summary_tsv(ctx.outputs["category_summary"], summary)
    write_enrichment_detection_tsv(ctx.outputs["enrichment_detection"], detection)
    write_enrichment_gsea_tsv(ctx.outputs["enrichment_gsea"], gsea)
    if SUPPORT in ctx.outputs:
        ctx.outputs[SUPPORT].mkdir(parents=True)
        running = gsea_running(genes, colls, gsea)
        if running is not None:
            write_support("gsea_running", running, ctx.outputs[SUPPORT])
