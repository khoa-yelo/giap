"""Operon steps (GIAP_SPEC.md §6.2, §7.1.3; M2-12): ``uniop_input``, ``uniop`` (ext) and
``operon_pairs``.

``uniop_input`` writes the Prodigal-style GFF and its ID map; the external ``uniop`` step runs
UniOP (never imported, ADR-0008). ``operon_pairs`` merges ``uniop.pred`` (row numbers resolved
through ``uniop_cds``) with the distance model (parameters from ``calibration.yaml``), builds TUs
and their upstream FASTAs and records the pair source in ``operon_meta.json``. With
``annotation.precomputed.operon_pairs`` the layer is read directly and its ``source`` kept.
Work-directory tables are 0-based half-open; ``tus.parquet`` ``genes`` is a native list column.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.annotate.operons import DISTANCE_FLAG, build_pairs, build_tus, write_tu_regions
from giap.annotate.precomputed import load_precomputed
from giap.annotate.uniop_adapter import PairProb, parse_uniop_pred, write_prodigal_gff
from giap.data import calibration
from giap.io.fasta import read_genome_fasta
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import read_tsv, write_tsv
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.steps.registry import InputPath, StepContext, external, step

_STR, _INT, _F64 = pa.string(), pa.int64(), pa.float64()
_PAIR_COLS = [("gene_a", _STR), ("gene_b", _STR), ("uniop_p", _F64), ("source", _STR)]
_TU_COLS = [
    ("tu_id", _STR), ("contig", _STR), ("start", _INT), ("end", _INT), ("strand", _STR),
    ("genes", pa.list_(_STR)), ("upstream_start", _INT), ("upstream_end", _INT),
]  # fmt: skip


def _path(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one file, got {len(p)}")
    return p


def _replicons(path: Path) -> dict[str, Replicon]:
    rows = pq.read_table(path).to_pylist()
    return {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}


def _genes(path: Path) -> list[Gene]:
    cols = ["locus_tag", "feature_type", "contig", "start", "end", "strand", "aa", "nt"]
    cols.append("start_codon")
    return [
        Gene(
            locus_tag=r["locus_tag"],
            feature_type=r["feature_type"],
            interval=Interval(r["contig"], r["start"], r["end"], r["strand"]),
            aa=r["aa"] or "",
            nt=r["nt"],
            start_codon=r["start_codon"],
        )
        for r in pq.read_table(path, columns=cols).to_pylist()
    ]


@step(
    "uniop_input",
    inputs=["genes", "replicons"],
    outputs=["uniop_cds", "uniop_map"],
    precomputed="operon_pairs",
)
def uniop_input(ctx: StepContext) -> None:
    """genes → Prodigal-style GFF (``uniop_cds``) and ``prodigal_id → locus_tag`` map."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    id_map = write_prodigal_gff(ctx.outputs["uniop_cds"], genes, list(reps.values()))
    write_tsv(ctx.outputs["uniop_map"], ["prodigal_id", "locus_tag"], id_map)


external("uniop", inputs=["uniop_cds"], outputs=["uniop_pred"], precomputed="operon_pairs")


def _uniop_pairs(ctx: StepContext, reps: dict[str, Replicon]) -> list[PairProb] | None:
    layer = ctx.cfg.annotation.precomputed.operon_pairs
    if layer is not None:
        df = load_precomputed("operon_pairs", layer, reps)
        return [
            PairProb(str(a), str(b), float(p), src)
            for a, b, p, src in df[["gene_a", "gene_b", "uniop_p", "source"]].itertuples(
                index=False
            )
        ]
    pred = ctx.inputs.get("uniop_pred")
    if pred is None:
        return None
    _, rows = read_tsv(_path(ctx.inputs["uniop_map"]))
    id_map = {str(r[0]): str(r[1]) for r in rows}
    return parse_uniop_pred(_path(pred), _path(ctx.inputs["uniop_cds"]), id_map)


@step(
    "operon_pairs",
    inputs=["genes", "uniop_map", "uniop_cds", "replicons", "genome"],
    optional_inputs=["uniop_pred"],
    outputs=["pairs_raw", "tus", "tu_upstream", "operon_meta"],
)
def operon_pairs(ctx: StepContext) -> None:
    """UniOP / distance-model pairs, TUs (threshold ``annotation.uniop.pair_threshold``)."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    model = calibration.load()["operon_distance"]
    uniop = _uniop_pairs(ctx, reps)
    pairs, flags = build_pairs(genes, reps, uniop, float(model["d0"]), float(model["s"]))
    if ctx.cfg.annotation.precomputed.operon_pairs is not None:
        meta = {"operon_source": "precomputed", "reason": "annotation.precomputed.operon_pairs"}
    elif DISTANCE_FLAG in flags:
        why = "no uniop.pred" if "uniop_pred" not in ctx.inputs else "uniop.pred has no pair"
        meta = {"operon_source": "distance_model", "reason": f"UniOP failed ({why})"}
        ctx.log.warning("operon_pairs: %s; distance model used, run flag %s", why, DISTANCE_FLAG)
    else:
        meta = {"operon_source": "uniop", "reason": "UniOP probabilities (pin 3ed14dfc)"}
    ctx.outputs["operon_meta"].write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    ann = ctx.cfg.annotation
    up_bp = (ann.tu_upstream_bp[0], ann.tu_upstream_bp[1])
    tus = build_tus(genes, pairs, ann.uniop.pair_threshold, reps, up_bp)
    out = ctx.outputs
    pair_rows = [(p.gene_a, p.gene_b, p.uniop_p, p.source) for p in pairs]
    write_parquet_deterministic(out["pairs_raw"], _PAIR_COLS, pair_rows)
    tu_rows = [
        (t.tu_id, t.contig, t.interval.start, t.interval.end, t.strand, list(t.genes),
         t.upstream.start, t.upstream.end)
        for t in tus
    ]  # fmt: skip
    write_parquet_deterministic(out["tus"], _TU_COLS, tu_rows)
    genome = {r.id: r.sequence for r in read_genome_fasta(_path(ctx.inputs["genome"]))}
    write_tu_regions(out["tu_upstream"], tus, genome)
