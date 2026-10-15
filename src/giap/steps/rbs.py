"""RBS steps (GIAP_SPEC.md §6.2, §7.1.6; M2-16): ``ostir_input``, ``ostir`` (ext), ``rbs_layer``.

``ostir_input`` writes one OSTIR row per gene of ``genes.parquet`` (§6.5 order) with the anti-SD
from ``annotation_meta.json``; the external ``ostir`` step runs OSTIR (GPL, never imported);
``rbs_layer`` turns ``ostir.csv`` into the precomputed ``rbs`` layer: one row per gene with an
OSTIR row at its requested start, in input order, ``start_codon`` taken from the input window.
"""

from __future__ import annotations

import json

import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.ostir import ostir_row, parse_ostir, read_ostir_input, write_ostir_input
from giap.io.fasta import read_genome_fasta
from giap.io.parquet import write_parquet_deterministic
from giap.steps.operons import _genes, _path, _replicons
from giap.steps.registry import StepContext, external, step

_F64 = pa.float64()
_RBS_COLS = [
    ("locus_tag", pa.string()), ("tir_ostir", _F64), ("sd_dg", _F64), ("mrna_dg", _F64),
    ("spacing_dg", _F64), ("standby_dg", _F64), ("spacer_bp", pa.int64()),
    ("start_codon", pa.string()),
]  # fmt: skip


@step(
    "ostir_input",
    inputs=["genes", "replicons", "genome", "annotation_meta"],
    outputs=["ostir_input"],
    precomputed="rbs",
)
def ostir_input(ctx: StepContext) -> None:
    """genes → ``ostir_input.csv`` (windows from ``annotation.ostir``)."""
    reps = _replicons(_path(ctx.inputs["replicons"]))
    genes = _genes(_path(ctx.inputs["genes"]))
    seqs = {r.id: r.sequence for r in read_genome_fasta(_path(ctx.inputs["genome"]))}
    anti_sd = str(json.loads(_path(ctx.inputs["annotation_meta"]).read_text())["anti_sd"])
    cfg = ctx.cfg.annotation.ostir
    rows = [
        ostir_row(
            g, reps[g.interval.contig], seqs[g.interval.contig], anti_sd,
            up=cfg.upstream_nt, down=cfg.downstream_nt,
        )
        for g in genes
    ]  # fmt: skip
    write_ostir_input(ctx.outputs["ostir_input"], rows)


external("ostir", inputs=["ostir_input"], outputs=["ostir"], precomputed="rbs")


@step("rbs_layer", inputs=["ostir", "ostir_input"], outputs=["rbs"], precomputed="rbs")
def rbs_layer(ctx: StepContext) -> None:
    """``ostir.csv`` → ``rbs.parquet`` (the precomputed ``rbs`` columns)."""
    rows = read_ostir_input(_path(ctx.inputs["ostir_input"]))
    feats = parse_ostir(_path(ctx.inputs["ostir"]), {r.id: r.start for r in rows})
    out = [
        (r.id, f.tir_ostir, f.sd_dg, f.mrna_dg, f.spacing_dg, f.standby_dg, f.spacer_bp,
         r.seq[r.start - 1 : r.start + 2])
        for r in rows
        if (f := feats.get(r.id)) is not None
    ]  # fmt: skip
    write_parquet_deterministic(ctx.outputs["rbs"], _RBS_COLS, out)
