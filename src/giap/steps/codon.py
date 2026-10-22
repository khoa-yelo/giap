"""Codon-usage step (GIAP_SPEC.md §6.2, §7.1.6; M2-17): ``codon_usage``.

``genes.parquet`` + ``trna_counts.tsv`` → ``codon.parquet``: one row per gene in the input order,
``locus_tag``, the eight ``CodonFeatures`` fields and the gene flag ``cai_reference``; missing
values are nulls. Run flags (``cai_reference_missing``) are logged and stored in the file metadata
under ``giap_flags`` (JSON list, ``[]`` when none); suppressed codon-bias warnings are logged once
as a count.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from giap.annotate.bakta_layers import is_ribosomal_protein
from giap.annotate.codon_usage import CodonFeatures, codon_usage_counted, scored_sequence
from giap.io.parquet import write_parquet_deterministic
from giap.io.tsv import read_tsv
from giap.model.coords import Interval
from giap.model.genes import Gene
from giap.steps.common import one_path as _path
from giap.steps.registry import StepContext, step

_F64 = pa.float64()
_FIELDS = [f.name for f in dataclasses.fields(CodonFeatures)]
_COLS = [
    ("locus_tag", pa.string()),
    *((n, pa.int64() if n == "rare_codons_2_10" else _F64) for n in _FIELDS),
    ("cai_reference", pa.bool_()),
]


def _read_genes(path: Path) -> list[Gene]:
    return [
        Gene(
            locus_tag=r["locus_tag"],
            feature_type=r["feature_type"],
            interval=Interval(r["contig"], r["start"], r["end"], r["strand"]),
            aa=r["aa"] or "",
            nt=r["nt"],
            start_codon=r["start_codon"],
            gene=r["gene"],
            product=r["product"] or "",
            partial=bool(r["partial"]),
        )
        for r in pq.read_table(path).to_pylist()
    ]


def _read_trna_counts(path: Path) -> dict[str, int]:
    header, rows = read_tsv(path)
    ac, n = header.index("anticodon"), header.index("n_copies")
    return {str(r[ac]): int(str(r[n])) for r in rows if r[ac] is not None}


@step("codon_usage", inputs=["genes", "trna_counts"], outputs=["codon"])
def codon_usage_step(ctx: StepContext) -> None:
    """genes + tRNA copy numbers → ``codon.parquet``."""
    genes = _read_genes(_path(ctx.inputs["genes"]))
    feats, flags, n_warned = codon_usage_counted(
        genes,
        _read_trna_counts(_path(ctx.inputs["trna_counts"])),
        ctx.cfg.annotation.translation_table,
    )
    for flag in flags:
        ctx.log.warning("codon_usage: run flag %s", flag)
    if n_warned:
        ctx.log.info("codon_usage: suppressed codon-bias RuntimeWarnings for %d genes", n_warned)
    rows = [
        (
            g.locus_tag,
            *(getattr(feats[g.locus_tag], n) for n in _FIELDS),
            is_ribosomal_protein(g) and scored_sequence(g) is not None,
        )
        for g in genes
    ]
    write_parquet_deterministic(
        ctx.outputs["codon"], _COLS, rows, metadata={"giap_flags": json.dumps(list(flags))}
    )
