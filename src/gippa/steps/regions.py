"""Step ``region_layer`` (GIPPA_SPEC.md §7.1.7; M2-19): genes + geNOMAD + antiSMASH → regions.

Required inputs are the ``parse_bakta`` tables; ``genomad`` (a geNOMAD output directory holding
``<prefix>_summary/``) and ``antismash`` (a directory with one top-level ``*.json``) are optional.
A missing optional layer is logged as ``layer_missing:<tool>``. Outputs: ``regions.parquet``
(``Region`` fields, 0-based half-open, nulls for AMR coordinates and NA scores) and
``region_genes.parquet`` (``RegionGene`` fields).

With ``annotation.precomputed.regions`` (read from the config, as ``operon_pairs`` reads its
layer) geNOMAD and antiSMASH are replaced, but the step still runs: it turns the M2-02 TSV into
the same two Parquet tables, the M4-06 C7/C8 that ``context_stats`` and ``sensitivity`` require
(M4-21 lead decision 2). Regions of a type whose module is off (``MODULE_REGION_TYPES``) are
dropped, and members come from the layer's ``genes`` lists (``precomputed_members``).
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from gippa.annotate.antismash import BgcRegion, parse_antismash_json
from gippa.annotate.genomad import GenomadResult, parse_genomad
from gippa.annotate.layers import regions_of
from gippa.annotate.precomputed import load_precomputed
from gippa.annotate.regions import (
    MODULE_REGION_TYPES,
    Region,
    RegionGene,
    build_regions,
    precomputed_members,
)
from gippa.config import GippaConfig
from gippa.io.parquet import write_parquet_deterministic
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, NcFeature
from gippa.steps.common import one_path as _path
from gippa.steps.common import read_replicons as _replicons
from gippa.steps.registry import StepContext, step

_STR, _INT, _F64 = pa.string(), pa.int64(), pa.float64()
_REGION_COLS = [
    ("region_id", _STR), ("region_type", _STR), ("contig", _STR), ("start", _INT), ("end", _INT),
    ("source_tool", _STR), ("score", _F64), ("genes", pa.list_(_STR)), ("notes", _STR),
]  # fmt: skip
_GENE_COLS = [
    ("region_id", _STR), ("locus_tag", _STR), ("overlap_frac", _F64), ("subset", _STR),
    ("core_biosynthetic", pa.bool_()),
]  # fmt: skip


def _genes(path: Path) -> list[Gene]:
    cols = ["locus_tag", "gene", "product", "feature_type", "contig", "start", "end", "strand"]
    return [
        Gene(
            locus_tag=r["locus_tag"],
            feature_type=r["feature_type"],
            interval=Interval(r["contig"], r["start"], r["end"], r["strand"]),
            aa="",
            nt=None,
            start_codon=None,
            gene=r["gene"],
            product=r["product"] or "",
            expert=("amrfinder",) if r["amr"] else (),
        )
        for r in pq.read_table(path, columns=[*cols, "amr"]).to_pylist()
    ]


def _features(path: Path) -> list[NcFeature]:
    cols = ["feature_id", "type", "contig", "start", "end", "strand"]
    return [
        NcFeature(
            r["feature_id"], r["type"], Interval(r["contig"], r["start"], r["end"], r["strand"])
        )
        for r in pq.read_table(path, columns=cols).to_pylist()
    ]


def _genomad(outdir: Path, reps: dict[str, Replicon]) -> GenomadResult:
    dirs = sorted(p for p in outdir.glob("*_summary") if p.is_dir())
    if len(dirs) != 1:
        raise ValueError(f"geNOMAD: expected one *_summary directory in {outdir}, got {len(dirs)}")
    return parse_genomad(outdir, dirs[0].name.removesuffix("_summary"), reps)


def _antismash(outdir: Path, reps: dict[str, Replicon]) -> list[BgcRegion]:
    jsons = sorted(outdir.glob("*.json"))
    if len(jsons) != 1:
        raise ValueError(f"antiSMASH: expected one top-level *.json in {outdir}, got {len(jsons)}")
    return parse_antismash_json(jsons[0], reps)


def _precomputed(
    layer: Path, cfg: GippaConfig, genes: list[Gene], reps: dict[str, Replicon]
) -> tuple[list[Region], list[RegionGene]]:
    off = {t for m, t in MODULE_REGION_TYPES.items() if m not in cfg.modules}
    df = load_precomputed("regions", layer, reps)
    regions = [r for r in regions_of(df) if r.region_type not in off]
    return regions, precomputed_members(regions, genes, reps)


def _layer_missing(layer: str) -> dict[str, object]:
    """The ``extra`` of a ``layer_missing:<layer>`` run-flag event (MVP-03 D-E)."""
    return {"event": "region_layer.layer_missing",
            "data": {"layer": layer, "flag": f"layer_missing:{layer}"}}  # fmt: skip


def _from_tools(
    ctx: StepContext, genes: list[Gene], reps: dict[str, Replicon]
) -> tuple[list[Region], list[RegionGene]]:
    gm = bgc = None
    if "genomad" in ctx.inputs:
        gm = _genomad(_path(ctx.inputs["genomad"]), reps)
    else:
        ctx.log.warning(
            "layer_missing:genomad: no geNOMAD output; no prophage regions",
            extra=_layer_missing("genomad"),
        )
    if "antismash" in ctx.inputs:
        bgc = _antismash(_path(ctx.inputs["antismash"]), reps)
    else:
        ctx.log.warning(
            "layer_missing:antismash: no antiSMASH output; no BGC regions",
            extra=_layer_missing("antismash"),
        )
    feats = _features(_path(ctx.inputs["features"]))
    return build_regions(genes, feats, reps, gm, bgc, ctx.cfg)


@step(
    "region_layer",
    inputs=["genes", "features", "replicons"],
    optional_inputs=["genomad", "antismash"],
    outputs=["regions", "region_genes"],
    precomputed="regions",
)
def region_layer(ctx: StepContext) -> None:
    """Prophage, plasmid, BGC, CRISPR and AMR regions with their member genes."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    layer = ctx.cfg.annotation.precomputed.regions
    if layer is not None:
        ctx.log.info("region_layer: regions from annotation.precomputed.regions (%s)", layer)
        regions, members = _precomputed(Path(layer), ctx.cfg, genes, reps)
    else:
        regions, members = _from_tools(ctx, genes, reps)
    rows = [
        (r.region_id, r.region_type, r.contig, r.interval.start if r.interval else None,
         r.interval.end if r.interval else None, r.source_tool, r.score, list(r.genes), r.notes)
        for r in regions
    ]  # fmt: skip
    write_parquet_deterministic(ctx.outputs["regions"], _REGION_COLS, rows)
    write_parquet_deterministic(
        ctx.outputs["region_genes"],
        _GENE_COLS,
        [
            (m.region_id, m.locus_tag, m.overlap_frac, m.subset, m.core_biosynthetic)
            for m in members
        ],
    )
