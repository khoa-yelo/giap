"""Annotation steps (GIAP_SPEC.md §6.2): ``parse_bakta`` (M2-04).

Bakta JSON plus the normalised genome FASTA become the internal work-directory tables
(0-based half-open coordinates, §6.5). Parquet is written without statistics and with fixed
metadata (``giap.io.parquet``) so that reruns are byte-identical (TE7).
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]

from giap.annotate.bakta import BaktaAnnotation, interval_nt, parse_bakta_json
from giap.annotate.bakta_layers import (
    annotation_source,
    anti_sd,
    is_ribosomal_protein,
    rrna16s,
    split_xrefs,
    trna_counts,
)
from giap.io.fasta import read_genome_fasta, write_fasta
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import write_tsv
from giap.model.coords import Replicon, interval_sort_key
from giap.model.genes import Gene, NcFeature
from giap.steps.registry import InputPath, StepContext, step

_STR, _INT, _BOOL = pa.string(), pa.int64(), pa.bool_()


def _path(p: InputPath) -> Path:
    if isinstance(p, list):
        raise ValueError(f"expected one file, got {len(p)}")
    return p


def _ordered(
    ann: BaktaAnnotation, order: Sequence[str]
) -> tuple[dict[str, Replicon], list[Gene], list[NcFeature]]:
    """Replicons, genes and features re-sorted by the genome FASTA order (§6.5)."""
    reps = ann.model.replicons
    if set(order) != set(reps):
        raise ValueError(
            f"genome FASTA contigs {sorted(order)} != Bakta JSON sequences {sorted(reps)}"
        )
    rank = {c: k for k, c in enumerate(order)}
    genes = sorted(ann.model.genes, key=lambda g: interval_sort_key(g.interval, rank, g.locus_tag))
    feats = sorted(ann.model.features, key=lambda f: interval_sort_key(f.interval, rank, f.id))
    return {c: reps[c] for c in order}, genes, feats


def _gene_row(g: Gene, rep: Replicon) -> tuple[Any, ...]:
    iv = g.interval
    amr = any("amrfinder" in e.lower() for e in g.expert)
    return (
        g.locus_tag, g.gene, g.product, g.feature_type, iv.contig, iv.start, iv.end, iv.strand,
        g.partial, iv.crosses_origin(rep), g.hypothetical, g.aa, g.nt, g.start_codon,
        is_ribosomal_protein(g), amr,
    )  # fmt: skip


def _feature_row(f: NcFeature) -> tuple[Any, ...]:
    iv = f.interval
    return (
        f.id, f.type, iv.contig, iv.start, iv.end, iv.strand, f.gene, f.product,
        f.amino_acid, f.anticodon,
    )  # fmt: skip


_GENE_COLS = [
    ("locus_tag", _STR), ("gene", _STR), ("product", _STR), ("feature_type", _STR),
    ("contig", _STR), ("start", _INT), ("end", _INT), ("strand", _STR), ("partial", _BOOL),
    ("crosses_origin", _BOOL), ("hypothetical", _BOOL), ("aa", _STR), ("nt", _STR),
    ("start_codon", _STR), ("ribosomal_protein", _BOOL), ("amr", _BOOL),
]  # fmt: skip
_FEATURE_COLS = [
    ("feature_id", _STR), ("type", _STR), ("contig", _STR), ("start", _INT), ("end", _INT),
    ("strand", _STR), ("gene", _STR), ("product", _STR), ("amino_acid", _STR),
    ("anticodon", _STR),
]  # fmt: skip
_REPLICON_COLS = [("contig", _STR), ("length", _INT), ("topology", _STR), ("kind", _STR)]
_XREF_COLS = [("locus_tag", _STR), ("db", _STR), ("id", _STR)]


def _trna_rows(feats: Sequence[NcFeature]) -> list[tuple[str, str | None, int]]:
    amino: dict[str, str | None] = {}
    for f in feats:
        if f.type == "trna" and f.anticodon:
            amino.setdefault(f.anticodon, f.amino_acid)
    return [(ac, amino[ac], n) for ac, n in trna_counts(feats).items()]


@step(
    "parse_bakta",
    inputs=["bakta_json", "genome"],
    outputs=[
        "genes",
        "features",
        "replicons",
        "rrna16s",
        "trna_counts",
        "xrefs",
        "annotation_meta",
    ],
    precomputed=None,
)
def parse_bakta(ctx: StepContext) -> None:
    """Bakta JSON → genes, features, replicons, 16S FASTA, tRNA counts, xrefs, meta (§6.2)."""
    records = read_genome_fasta(_path(ctx.inputs["genome"]))
    genome = {r.id: r.sequence for r in records}
    ann = parse_bakta_json(_path(ctx.inputs["bakta_json"]), genome=genome)
    reps, genes, feats = _ordered(ann, [r.id for r in records])
    for t, k in sorted(ann.skipped_types.items()):
        ctx.log.warning("parse_bakta: skipped %d Bakta feature(s) of unknown type %r", k, t)

    s16 = [
        (f.id, f.nt or interval_nt(ann.sequences[f.interval.contig], f.interval))
        for f in rrna16s(feats)
    ]
    anti, sd_flags = anti_sd([s for _, s in s16], ctx.cfg.annotation.ostir.anti_sd)
    flags = sorted({*sd_flags, *(["bakta_db_light"] if ann.db_type == "light" else [])})

    out = ctx.outputs
    write_parquet_deterministic(
        out["genes"], _GENE_COLS, [_gene_row(g, reps[g.interval.contig]) for g in genes]
    )
    write_parquet_deterministic(out["features"], _FEATURE_COLS, [_feature_row(f) for f in feats])
    rep_rows = [(r.id, r.length, r.topology, r.kind) for r in reps.values()]
    write_parquet_deterministic(out["replicons"], _REPLICON_COLS, rep_rows)
    write_fasta(out["rrna16s"], s16)
    write_tsv(out["trna_counts"], ["anticodon", "amino_acid", "n_copies"], _trna_rows(feats))
    write_parquet_deterministic(out["xrefs"], _XREF_COLS, split_xrefs(genes))
    meta: Mapping[str, Any] = {
        "bakta_version": ann.bakta_version,
        "db_version": ann.db_version,
        "db_type": ann.db_type,
        "annotation_source": annotation_source(ann.db_type),
        "anti_sd": anti,
        "flags": flags,
        "n_genes": len(genes),
        "n_features": len(feats),
        "skipped_types": dict(sorted(ann.skipped_types.items())),
    }
    out["annotation_meta"].write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
