"""Promoter steps (GIAP_SPEC.md §6.2, §7.1.4; M2-13): ``promoter_windows``, ``promotech`` (ext)
and ``promoter_layer``.

``promoter_windows`` writes the deduplicated 40-nt windows of every CDS and sORF start as
``windows/windows_NNN.fa`` (``w<index>`` headers) plus ``windows_map.parquet``; the external
``promotech`` step scores each chunk (``promotech_NNN/sequences_predictions.csv``, never imported,
GPL-3.0); ``promoter_layer`` joins the scores back to unique positions. All three are replaced by
``annotation.precomputed.promoters``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.promotech_cmd import parse_promotech
from giap.annotate.promoters import dedupe_and_chunk, gene_windows, promoter_layer, write_chunks
from giap.annotate.uniop_adapter import genome_order
from giap.io.fasta import read_genome_fasta
from giap.io.parquet import write_parquet_deterministic
from giap.steps.operons import _genes, _path, _replicons
from giap.steps.registry import InputPath, StepContext, external, step

_STR, _INT, _F64 = pa.string(), pa.int64(), pa.float64()
_MAP_COLS = [
    ("header", _STR), ("contig", _STR), ("strand", _STR), ("end3p", _INT), ("owner", _STR),
    ("k", _INT),
]  # fmt: skip
_LAYER_COLS = [("contig", _STR), ("strand", _STR), ("end3p", _INT), ("score", _F64)]


def _rows(df: pd.DataFrame) -> list[tuple[object, ...]]:
    return list(df.itertuples(index=False, name=None))


@step(
    "promoter_windows",
    inputs=["genes", "replicons", "genome"],
    outputs=["windows", "windows_map"],
    precomputed="promoters",
)
def promoter_windows(ctx: StepContext) -> None:
    """Windows of every CDS/sORF start (§6.5 gene order, k ascending), deduplicated and chunked."""
    cfg = ctx.cfg.annotation.promoters
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ctx.inputs["genome"]))}
    genes = [g for g in _genes(_path(ctx.inputs["genes"])) if g.feature_type != "pseudogene"]
    k_range = (cfg.window_3p_range[0], cfg.window_3p_range[1])
    windows = [
        w
        for contig, gs in genome_order(genes, list(reps)).items()
        for g in gs
        for w in gene_windows(g, reps[contig], genome[contig], k_range, cfg.window_nt)
    ]
    chunks, wmap = dedupe_and_chunk(windows)
    write_chunks(ctx.outputs["windows"], chunks)
    write_parquet_deterministic(ctx.outputs["windows_map"], _MAP_COLS, _rows(wmap))
    ctx.log.info("promoter_windows: %d windows, %d unique", len(wmap), sum(map(len, chunks)))


external("promotech", inputs=["windows"], outputs=["promotech"], precomputed="promoters")


def _prediction_files(p: InputPath) -> list[Path]:
    paths = p if isinstance(p, list) else [p]
    out: list[Path] = []
    for path in sorted(paths):
        out += sorted(path.rglob("sequences_predictions.csv")) if path.is_dir() else [path]
    return out


@step(
    "promoter_layer",
    inputs=["windows_map", "promotech"],
    outputs=["promoter_windows"],
    precomputed="promoters",
)
def promoter_layer_step(ctx: StepContext) -> None:
    """Promotech scores joined to unique (contig, strand, end3p) window positions."""
    wmap = pd.read_parquet(_path(ctx.inputs["windows_map"]))
    scores: dict[str, float] = {}
    for path in _prediction_files(ctx.inputs["promotech"]):
        scores.update(parse_promotech(path))
    layer = promoter_layer(wmap, scores)
    write_parquet_deterministic(ctx.outputs["promoter_windows"], _LAYER_COLS, _rows(layer))
