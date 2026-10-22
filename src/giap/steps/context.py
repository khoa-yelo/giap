"""Step ``context_layers`` (GIAP_SPEC.md §6.2, §7.1.8, §9.1; M2-20): parsed layers →
``layers.parquet`` and the ``annotation/`` directory.

Required: the ``parse_bakta`` gene model, ``pairs_raw`` / ``tus`` and the genome FASTA
(``gc_local_dev``). Every other layer is optional; a missing one gives NA columns and a
``layer_missing:<layer>`` warning. A precomputed layer (``annotation.precomputed.<key>`` set)
arrives as its M2-02 TSV under the producer's key (local executor) or ``promoters`` (Nextflow);
``--modules`` without ``mobile`` / ``bgc`` drops prophage / BGC rows of precomputed regions.
Under a precomputed regions layer ``region_layer`` still runs (M4-21 lead decision 2), so
``regions`` / ``region_genes`` may be its Parquet (already filtered) or the layer TSV itself.
``terminators_crosscheck.tsv`` and ``id_map.tsv`` are byte copies of their inputs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.layers import (
    INT_COLUMNS,
    LAYER_COLUMNS,
    LIST_COLUMNS,
    STR_COLUMNS,
    assign_tus,
    gene_layers,
    genes_of,
    layer_flags,
    merged_gff3_features,
    pairs_of,
    regions_of,
    replicons_of,
    terminator_records,
    terminators_of,
    tu_sources,
    tus_of,
    write_annotation,
)
from giap.annotate.precomputed import load_precomputed
from giap.annotate.promoters import promoter_records
from giap.annotate.regions import MODULE_REGION_TYPES
from giap.io.fasta import read_genome_fasta
from giap.io.parquet import write_parquet_deterministic
from giap.model.coords import Replicon
from giap.steps.common import one_path as _path
from giap.steps.registry import StepContext, step

_TYPES = {c: pa.int64() for c in INT_COLUMNS} | {c: pa.list_(pa.string()) for c in LIST_COLUMNS}
_TYPES |= {c: pa.string() for c in STR_COLUMNS}
_COLS = [(c, _TYPES.get(c, pa.float64())) for c in LAYER_COLUMNS]
_PARQUET_MAGIC = b"PAR1"


def _frame(ctx: StepContext, key: str, pre: str, reps: dict[str, Replicon]) -> pd.DataFrame | None:
    """The layer under ``key`` (parquet, or the M2-02 TSV when ``precomputed.<pre>`` is set)."""
    path = ctx.inputs.get(key)
    if path is None:
        return None
    if getattr(ctx.cfg.annotation.precomputed, pre) is not None:
        return load_precomputed(pre, _path(path), reps)
    return pd.read_parquet(_path(path))


def _is_parquet(path: Path) -> bool:
    with path.open("rb") as fh:
        return fh.read(4) == _PARQUET_MAGIC


def _regions(
    ctx: StepContext, reps: dict[str, Replicon]
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    given = ctx.inputs.get("regions")
    if given is not None and _is_parquet(_path(given)):  # region_layer's C7 (and C8)
        rg = ctx.inputs.get("region_genes")
        return pd.read_parquet(_path(given)), None if rg is None else pd.read_parquet(_path(rg))
    regions = _frame(ctx, "regions", "regions", reps)
    if regions is None or ctx.cfg.annotation.precomputed.regions is None:
        rg = ctx.inputs.get("region_genes")
        return regions, None if rg is None else pd.read_parquet(_path(rg))
    off = [t for m, t in MODULE_REGION_TYPES.items() if m not in ctx.cfg.modules]
    regions = regions[~regions["region_type"].isin(off)].reset_index(drop=True)
    members = [(r, g) for r, gs in zip(regions["region_id"], regions["genes"], strict=True)
               for g in gs]  # fmt: skip
    return regions, pd.DataFrame(members, columns=["region_id", "locus_tag"])


def _cell(c: str, v: Any) -> Any:
    if c in LIST_COLUMNS:
        return [str(x) for x in v]
    if v is None or v is pd.NA or (isinstance(v, float) and v != v):
        return None
    return int(v) if c in INT_COLUMNS else v


@step(
    "context_layers",
    inputs=["genes", "replicons", "pairs_raw", "tus", "genome"],
    optional_inputs=[
        "features", "promoter_windows", "promoters", "tx_rates", "terminators",
        "terminators_crosscheck", "rbs", "codon", "regions", "region_genes", "id_map", "tool_dirs",
    ],
    outputs=["annotation", "layers"],
)  # fmt: skip
def context_layers(ctx: StepContext) -> None:
    """Per-gene layers (``layers.parquet``) and every ``annotation/*`` file."""
    cfg, ins = ctx.cfg, ctx.inputs
    for flag in layer_flags(cfg, set(ins)):
        layer = flag.split(":", 1)[1]
        ctx.log.warning(
            "%s: layer not available; its columns are NA", flag,
            extra={"event": "context_layers.layer_missing", "data": {"layer": layer, "flag": flag}},
        )  # fmt: skip
    rep_df = pd.read_parquet(_path(ins["replicons"]))
    reps = replicons_of(rep_df)
    genes = pd.read_parquet(_path(ins["genes"]))
    tus, pairs = pd.read_parquet(_path(ins["tus"])), pd.read_parquet(_path(ins["pairs_raw"]))
    prom_key = "promoters" if "promoters" in ins else "promoter_windows"
    prom = _frame(ctx, prom_key, "promoters", reps)
    if prom is not None:
        prom = prom[["contig", "strand", "end3p", "score"]]
    tx = _frame(ctx, "tx_rates", "tx_rates", reps)
    tu_objs = tus_of(tus)
    if tx is not None and cfg.annotation.precomputed.tx_rates is not None:
        tx = assign_tus(tx, tu_objs, reps)
    terms = _frame(ctx, "terminators", "terminators", reps)
    rbs = _frame(ctx, "rbs", "rbs", reps)
    codon = None if "codon" not in ins else pd.read_parquet(_path(ins["codon"]))
    regions, region_genes = _regions(ctx, reps)
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ins["genome"]))}
    yes = cfg.annotation.promoters.promotech_yes
    k = cfg.annotation.promoters.window_3p_range
    df = gene_layers(genes, tus, pairs, prom, tx, terms, rbs, codon, region_genes, rep_df,
                     genome=genome, yes=yes, k_range=(k[0], k[1]))  # fmt: skip
    rows = [tuple(_cell(c, v) for c, v in zip(LAYER_COLUMNS, r, strict=True))
            for r in df.itertuples(index=False, name=None)]  # fmt: skip
    write_parquet_deterministic(ctx.outputs["layers"], _COLS, rows)

    pair_objs = pairs_of(pairs)
    regs = None if regions is None else regions_of(regions)
    prec = None if prom is None else promoter_records(prom, reps, yes)
    pbed = None if prec is None else [b for b, _ in prec]
    terms_rec = None if terms is None else terminator_records(terminators_of(terms), reps)
    feats = merged_gff3_features(
        genes_of(genes),
        tu_objs,
        prec or [],
        terms_rec or [],
        regs or [],
        tu_source=tu_sources(tu_objs, pair_objs),
    )
    copies: dict[str, Path] = {}
    if "terminators_crosscheck" in ins and cfg.annotation.precomputed.terminators is None:
        copies["terminators_crosscheck.tsv"] = _path(ins["terminators_crosscheck"])
    if "id_map" in ins:
        copies["id_map.tsv"] = _path(ins["id_map"])
    tool = ins.get("tool_dirs")
    tool_dirs = {p.name: p for p in (tool if isinstance(tool, list) else [tool] if tool else [])}
    write_annotation(
        ctx.outputs["annotation"], replicons=reps, features=feats, pairs=pair_objs, tus=tu_objs,
        regions=regs, promoters_bed=pbed,
        term_bed=None if terms_rec is None else [b for b, _ in terms_rec],
        copies=copies, tool_dirs=tool_dirs,
    )  # fmt: skip
