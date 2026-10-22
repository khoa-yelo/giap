"""Step ``context_stats`` (GIAP_SPEC.md §6.2; M4-18): gene table → ``tables/regions_stats.tsv``
and ``tables/detection_summary.tsv`` (the M4 parts; M5 adds autocorrelation, islands, ori-ter).

Inputs: the final gene table (M4-15 ``work/genetable/genes.parquet``, key ``gene_table``), the
M4-06 contract tables C1-C14 (validated with ``load_m4_inputs``; an ``error`` issue stops the
step), optionally the ``--replicon-depth`` TSV (§5.3), M4-15's ``detectability_summary``
parquet (detectability diagnostics rows; locally ``genetable/detection_summary_parts/
detectability.parquet``), ``tier_summary`` of ``gene_table`` (M4-20 decision 3; locally
``genetable/detection_summary_parts/tier.parquet``), the ADR-0019 ``stratum0_classes`` of
``placement`` (the ``stratum0_unmapped_share`` QC rows) and C15 ``features`` (M5-04), whose ``oric``
rows place the origin (without it: ``no_oric`` rows and ``layer_missing:features``).

Outputs: ``tables/regions_stats.tsv``, ``detection_summary.tsv``, ``autocorrelation.tsv``,
``islands.tsv``, ``genes.tsv`` (here until §6.2's ``finalize``, M5-09 lead decision 3) and
``work/context/{gene_context,island_diagnostics}.parquet``. Streams from ``config.seed``:
``context_stats``, ``.plasmids``, ``.detection_summary``, ``.autocorrelation``, ``.islands`` and
``.oriter``, so no part's draws depend on another's and M5 knobs never move the M4 tables.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.config.model import GiapConfig
from giap.context.autocorrelation import autocorrelation_rows
from giap.context.diagnostics import DIAGNOSTICS_COLUMNS, island_diagnostics
from giap.context.islands import islands_frame, scan_replicons, tu_series
from giap.context.oriter import oriter_rows
from giap.context.regions import replicon_index, test_all_regions
from giap.context.regions_table import depth_gaps, plasmid_rows, regions_stats, subset_rows
from giap.context.summary import basic_rows, detection_summary, qc_rows
from giap.context.tables import (
    GENE_CONTEXT_COLUMNS,
    autocorrelation_table,
    gene_context_table,
    islands_table,
)
from giap.genetable.assemble import GENES_COLUMNS
from giap.io.context_tsv import write_autocorrelation, write_islands
from giap.io.depth import read_replicon_depth
from giap.io.genes_tsv import write_genes_tsv
from giap.io.parquet import write_parquet_deterministic
from giap.io.regions_tsv import write_regions_stats
from giap.io.summary_tsv import write_detection_summary
from giap.io.work_tables import (
    M4_INPUT_KEYS,
    OPTIONAL_KEYS,
    M4Inputs,
    default_paths,
    load_m4_inputs,
)
from giap.rng import for_step
from giap.steps.gene_table import contract_paths
from giap.steps.registry import InputPath, StepContext, step

__all__ = ["context_stats", "run_context_stats"]

_LOG = logging.getLogger("giap.steps.context_stats")
_REQUIRED = [k for k in M4_INPUT_KEYS if k not in OPTIONAL_KEYS]
#: the M5 half's output keys, all optional at the registry level (wired together or not at all)
M5_OUTPUTS: tuple[str, ...] = (
    "autocorrelation", "islands", "gene_context", "island_diagnostics", "genes_tsv",
)  # fmt: skip


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

    s0 = inputs.get("stratum0_classes")  # ADR-0019 class table (M3-25 decision 3)
    parts = [basic_rows(genes, cfg.stats.bootstrap, for_step("context_stats.detection_summary",
                                                               seed)),
             qc_rows(m4.peptides, m4.peptide_hits, m4.unmapped,
                     stratum0=None if s0 is None else pd.read_parquet(s0))]  # fmt: skip
    if "tier_summary" in inputs:  # M4-12 tier counts and q-threshold curve (decision 3)
        parts.append(pd.read_parquet(inputs["tier_summary"]))
    if "detectability_summary" in inputs:  # M4-15 detectability diagnostics rows
        parts.append(pd.read_parquet(inputs["detectability_summary"]))
    write_detection_summary(outputs["detection_summary"], detection_summary(parts))

    _m5_context(inputs, outputs, cfg, log, m4, genes)


def _features(inputs: Mapping[str, Path], log: logging.Logger) -> pd.DataFrame | None:
    """The C15 ``features`` table, or None with the ``layer_missing:features`` run flag."""
    path = inputs.get("features")
    if path is None:
        log.warning(
            "layer_missing:features: no C15 features table; the ori-ter rows are no_oric and "
            "the genes.tsv replichore columns stay NA",
            extra={
                "event": "context_stats.layer_missing",
                "data": {"layer": "features", "flag": "layer_missing:features"},
            },
        )
        return None
    return pd.read_parquet(path)


def _m5_context(
    inputs: Mapping[str, Path], outputs: Mapping[str, Path], cfg: GiapConfig,
    log: logging.Logger, m4: M4Inputs, genes: pd.DataFrame,
) -> None:  # fmt: skip
    """The M5 half: autocorrelation, islands, ori-ter and genes.tsv (§7.3.5(c)-(e)).

    The M5 keys are ``optional_outputs``, so a caller wiring only the M4 half runs nothing
    here; ``run_context_stats`` and the CONTEXT subworkflow wire all of them.
    """
    missing = [k for k in M5_OUTPUTS if k not in outputs]
    if len(missing) == len(M5_OUTPUTS):
        log.info("M5 context outputs not wired; writing the M4 tables only")
        return
    if missing:
        raise ValueError(f"context_stats: M5 output(s) not given: {', '.join(missing)}")
    seed = cfg.seed
    features = _features(inputs, log)

    rows = autocorrelation_rows(
        genes, m4.replicons, cfg, for_step("context_stats.autocorrelation", seed)
    )
    oriter = oriter_rows(
        genes, m4.tus, m4.replicons, features, cfg, for_step("context_stats.oriter", seed)
    )
    write_autocorrelation(
        outputs["autocorrelation"], autocorrelation_table([rows, oriter], m4.replicons)
    )

    island_rng = for_step("context_stats.islands", seed)  # scan_replicons, then islands_frame
    span = float(cfg.context.islands.detrend_span_frac)
    series = [
        tu_series(
            replicon_index(genes, str(contig)),
            m4.tus,
            contig=str(contig),
            circular=str(topology) == "circular",
            length=int(length),
            span_frac=span,
        )
        for contig, topology, length in zip(
            m4.replicons["contig"], m4.replicons["topology"], m4.replicons["length"], strict=True
        )
    ]
    scan = scan_replicons(series, cfg, island_rng)
    islands = islands_table(islands_frame(scan, genes, cfg, island_rng), m4.replicons)
    write_islands(outputs["islands"], islands)
    _write_diagnostics(outputs["island_diagnostics"], island_diagnostics(scan, cfg))

    context = gene_context_table(genes, islands, m4.replicons, features)
    _write_gene_context(outputs["gene_context"], context)
    write_genes_tsv(outputs["genes_tsv"], _with_context(genes, context))


def _write_gene_context(path: Path, df: pd.DataFrame) -> None:
    """``work/context/gene_context.parquet`` (byte-reproducible, TE7)."""
    cols = [("locus_tag", pa.string()), ("island_id", pa.string()),
            ("replichore", pa.string()), ("replichore_position", pa.float64())]  # fmt: skip
    rows = [
        (str(tag), _text(island), _text(arm), float(pos))
        for tag, island, arm, pos in df[list(GENE_CONTEXT_COLUMNS)].itertuples(
            index=False, name=None
        )
    ]
    write_parquet_deterministic(path, cols, rows)


def _write_diagnostics(path: Path, df: pd.DataFrame) -> None:
    """``work/context/island_diagnostics.parquet`` (ADR-0021 Amendment 1)."""
    cols = [("contig", pa.string()), ("scanned", pa.bool_()), ("block_rule", pa.string()),
            ("l_hat", pa.int64()), ("block_len", pa.int64()),
            ("morans_i_lag1", pa.float64()), ("flags", pa.list_(pa.string()))]  # fmt: skip
    rows = [
        (str(contig), bool(scanned), str(rule), _count(l_hat), _count(block_len),
         float(moran), [str(f) for f in flags])
        for contig, scanned, rule, l_hat, block_len, moran, flags in df[
            list(DIAGNOSTICS_COLUMNS)
        ].itertuples(index=False, name=None)
    ]  # fmt: skip
    write_parquet_deterministic(path, cols, rows)


def _text(value: object) -> str | None:
    """A nullable string cell of a context column."""
    return None if value is None or pd.isna(value) else str(value)


def _count(value: object) -> int | None:
    """A nullable integer cell (``l_hat`` / ``block_len``: NA, never 0)."""
    return None if value is None or pd.isna(value) else int(value)  # type: ignore[call-overload]


def _with_context(genes: pd.DataFrame, context: pd.DataFrame) -> pd.DataFrame:
    """The gene table with the three §9.3 context columns joined in; absent
    ``GENES_COLUMNS`` are NA, so a partially assembled table still writes."""
    out = genes.copy()
    out.index = pd.RangeIndex(len(out))
    for column in ("island_id", "replichore", "replichore_position"):
        out[column] = context[column].to_numpy()
    for column in GENES_COLUMNS:
        if column not in out.columns:
            out[column] = pd.Series([None] * len(out), dtype=object)
    return out[list(GENES_COLUMNS)]


def run_context_stats(work_dir: Path, cfg: GiapConfig, replicon_depth: Path | None) -> None:
    """Run the step on a work directory: contract tables found by ``contract_paths`` (the flat
    fixture layout or the keyed ``giap run`` layout, M4-21 lead decision 3) plus any optional
    contract table present (C15 ``features.parquet``), the gene table at
    ``genetable/genes.parquet``, the tables under ``tables/`` and the two intermediate parquet
    files under ``context/``."""
    inputs = contract_paths(work_dir)
    inputs |= {k: p for k, p in default_paths(work_dir).items()
               if k in OPTIONAL_KEYS and k not in inputs and p.exists()}  # fmt: skip
    inputs["gene_table"] = work_dir / "genetable" / "genes.parquet"
    parts = work_dir / "genetable" / "detection_summary_parts"
    for key, name in (("detectability_summary", "detectability"), ("tier_summary", "tier")):
        if (parts / f"{name}.parquet").exists():
            inputs[key] = parts / f"{name}.parquet"
    classes = work_dir / "mapping" / "stratum0_classes.parquet"
    if classes.exists():
        inputs["stratum0_classes"] = classes
    features = work_dir / "features.parquet"  # C15 (M5-04), optional
    if features.exists():
        inputs["features"] = features
    if replicon_depth is not None:
        inputs["replicon_depth"] = replicon_depth
    tables = work_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    context = work_dir / "context"
    context.mkdir(parents=True, exist_ok=True)
    outputs = {"regions_stats": tables / "regions_stats.tsv",
               "detection_summary": tables / "detection_summary.tsv",
               "autocorrelation": tables / "autocorrelation.tsv",
               "islands": tables / "islands.tsv",
               "genes_tsv": tables / "genes.tsv",
               "gene_context": context / "gene_context.parquet",
               "island_diagnostics": context / "island_diagnostics.parquet"}  # fmt: skip
    _context_stats(inputs, outputs, cfg, _LOG)


@step(
    "context_stats",
    inputs=["gene_table", *_REQUIRED],
    optional_inputs=[
        *sorted(OPTIONAL_KEYS),
        "detectability_summary",
        "features",
        "replicon_depth",
        "stratum0_classes",
        "tier_summary",
    ],
    outputs=["regions_stats", "detection_summary", *M5_OUTPUTS],
    optional_outputs=list(M5_OUTPUTS),
)
def context_stats(ctx: StepContext) -> None:
    """Region statistics, autocorrelation, islands, ori-ter and genes.tsv from the gene table."""
    inputs = {k: _path(v) for k, v in ctx.inputs.items()}
    _context_stats(inputs, ctx.outputs, ctx.cfg, ctx.log)
