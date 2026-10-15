"""Step ``region_layer`` (GIAP_SPEC.md §7.1.7; M2-19): genes + geNOMAD + antiSMASH → regions.

Required inputs are the ``parse_bakta`` tables; ``genomad`` (a geNOMAD output directory holding
``<prefix>_summary/``) and ``antismash`` (a directory with one top-level ``*.json``) are optional.
A missing optional layer is logged as ``layer_missing:<tool>``. Outputs: ``regions.parquet``
(``Region`` fields, 0-based half-open, nulls for AMR coordinates and NA scores) and
``region_genes.parquet`` (``RegionGene`` fields).
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.annotate.antismash import BgcRegion, parse_antismash_json
from giap.annotate.genomad import GenomadResult, parse_genomad
from giap.annotate.regions import build_regions
from giap.io.parquet import write_parquet_deterministic
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, NcFeature
from giap.steps.operons import _path, _replicons
from giap.steps.registry import StepContext, step

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
    gm = bgc = None
    if "genomad" in ctx.inputs:
        gm = _genomad(_path(ctx.inputs["genomad"]), reps)
    else:
        ctx.log.warning("layer_missing:genomad: no geNOMAD output; no prophage regions")
    if "antismash" in ctx.inputs:
        bgc = _antismash(_path(ctx.inputs["antismash"]), reps)
    else:
        ctx.log.warning("layer_missing:antismash: no antiSMASH output; no BGC regions")
    genes = _genes(_path(ctx.inputs["genes"]))
    feats = _features(_path(ctx.inputs["features"]))
    regions, members = build_regions(genes, feats, reps, gm, bgc, ctx.cfg)
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
