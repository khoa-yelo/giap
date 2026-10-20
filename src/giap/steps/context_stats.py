"""Step ``context_stats`` (GIAP_SPEC.md §6.2; M4-18): gene table → ``tables/regions_stats.tsv``
and ``tables/detection_summary.tsv`` (the M4 parts; M5 adds autocorrelation, islands, ori-ter).

Inputs: the final gene table (M4-15 ``work/genetable/genes.parquet``, key ``gene_table``), the
M4-06 contract tables C1-C14 (validated with ``load_m4_inputs``; an ``error`` issue stops the
step), optionally the ``--replicon-depth`` TSV (§5.3) and M4-15's ``detectability_summary``
parquet (detectability diagnostics rows; locally ``genetable/detection_summary_parts/
detectability.parquet``). Random streams, all from
``config.seed``: ``context_stats`` (M4-17 region tests), ``context_stats.plasmids`` (parent of the
per-plasmid streams) and ``context_stats.detection_summary`` (basic O - E CIs), so no part's draws
depend on another part's content.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from giap.config.model import GiapConfig
from giap.context.regions import test_all_regions
from giap.context.regions_table import depth_gaps, plasmid_rows, regions_stats, subset_rows
from giap.context.summary import basic_rows, detection_summary, qc_rows
from giap.io.depth import read_replicon_depth
from giap.io.regions_tsv import write_regions_stats
from giap.io.summary_tsv import write_detection_summary
from giap.io.work_tables import M4_INPUT_KEYS, OPTIONAL_KEYS, default_paths, load_m4_inputs
from giap.rng import for_step
from giap.steps.registry import InputPath, StepContext, step

__all__ = ["context_stats", "run_context_stats"]

_LOG = logging.getLogger("giap.steps.context_stats")
_REQUIRED = [k for k in M4_INPUT_KEYS if k not in OPTIONAL_KEYS]


def _path(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one path, got {len(p)}")
    return p


def _context_stats(
    inputs: Mapping[str, Path], outputs: Mapping[str, Path], cfg: GiapConfig,
    log: logging.Logger,
) -> None:  # fmt: skip
    m4, issues = load_m4_inputs({k: v for k, v in inputs.items() if k in M4_INPUT_KEYS})
    errors = [i for i in issues if i.severity == "error"]
    if errors:
        raise ValueError(f"context_stats: contract errors: {[i.message for i in errors[:5]]}")
    genes = pd.read_parquet(inputs["gene_table"])
    depth = None
    if "replicon_depth" in inputs:
        depth, dissues = read_replicon_depth(inputs["replicon_depth"])
        if depth is None:
            raise ValueError(f"replicon depth: {[i.message for i in dissues]}")
        gaps = depth_gaps(genes, m4.replicons, depth, cfg)
        if gaps:
            log.warning("replicon_depth: no depth for %s; dosage-adjusted plasmid rows omitted",
                        ", ".join(gaps))  # fmt: skip
    else:
        log.info("replicon_depth: not given; no dosage-adjusted plasmid rows")
    seed = cfg.seed
    regions = m4.regions[m4.regions["region_type"].isin(cfg.context.region_types)]
    tests = test_all_regions(genes, regions, m4.region_genes, m4.replicons, cfg,
                             for_step("context_stats", seed))  # fmt: skip
    plasmids = plasmid_rows(genes, m4.replicons, depth, cfg, for_step("context_stats.plasmids",
                                                                       seed))  # fmt: skip
    subsets = subset_rows(genes, regions, m4.region_genes)
    write_regions_stats(outputs["regions_stats"], regions_stats(tests, plasmids, subsets, regions))

    parts = [basic_rows(genes, cfg.stats.bootstrap, for_step("context_stats.detection_summary",
                                                               seed)),
             qc_rows(m4.peptides, m4.peptide_hits, m4.unmapped)]  # fmt: skip
    if "detectability_summary" in inputs:  # M4-15 detectability diagnostics rows
        parts.append(pd.read_parquet(inputs["detectability_summary"]))
    write_detection_summary(outputs["detection_summary"], detection_summary(parts))


def run_context_stats(work_dir: Path, cfg: GiapConfig, replicon_depth: Path | None) -> None:
    """Run the step on the local work-directory layout: contract tables directly under
    ``work_dir``, the gene table at ``genetable/genes.parquet``, outputs under ``tables/``."""
    inputs = {k: p for k, p in default_paths(work_dir).items()
              if k not in OPTIONAL_KEYS or p.exists()}  # fmt: skip
    inputs["gene_table"] = work_dir / "genetable" / "genes.parquet"
    part = work_dir / "genetable" / "detection_summary_parts" / "detectability.parquet"
    if part.exists():
        inputs["detectability_summary"] = part
    if replicon_depth is not None:
        inputs["replicon_depth"] = replicon_depth
    tables = work_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    outputs = {"regions_stats": tables / "regions_stats.tsv",
               "detection_summary": tables / "detection_summary.tsv"}  # fmt: skip
    _context_stats(inputs, outputs, cfg, _LOG)


@step(
    "context_stats",
    inputs=["gene_table", *_REQUIRED],
    optional_inputs=[*sorted(OPTIONAL_KEYS), "detectability_summary", "replicon_depth"],
    outputs=["regions_stats", "detection_summary"],
)
def context_stats(ctx: StepContext) -> None:
    """Region statistics and the detection summary from the gene table."""
    inputs = {k: _path(v) for k, v in ctx.inputs.items()}
    _context_stats(inputs, ctx.outputs, ctx.cfg, ctx.log)
