"""Input helpers shared by the step modules (M2-16 review, M6-08 item 8): one declared path,
the C1 replicons and the C2 gene model read from their work-directory Parquet files, and the
operon source of ``operon_meta.json`` (M11-14d)."""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow.parquet as pq  # type: ignore[import-untyped]

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene
from gippa.steps.registry import InputPath

__all__ = ["one_path", "operon_source", "read_genes", "read_replicons"]


def one_path(p: InputPath) -> Path:
    """The single file of a step input; a list (a scattered key) raises ``ValueError``."""
    if isinstance(p, list):
        raise ValueError(f"expected one file, got {len(p)}")
    return p


def operon_source(path: InputPath | None) -> str | None:
    """``operon_source`` of ``operon_meta.json`` (``uniop`` | ``distance_model`` |
    ``precomputed``); ``None`` without the input. An unreadable file raises."""
    if path is None:
        return None
    data = json.loads(one_path(path).read_text(encoding="utf-8"))
    value = data.get("operon_source") if isinstance(data, dict) else None
    return value if isinstance(value, str) else None


def read_replicons(path: Path) -> dict[str, Replicon]:
    """C1 ``replicons.parquet`` → replicons by contig, in file order."""
    rows = pq.read_table(path).to_pylist()
    return {r["contig"]: Replicon(r["contig"], r["length"], r["topology"], r["kind"]) for r in rows}


def read_genes(path: Path) -> list[Gene]:
    """C2 ``genes.parquet`` → genes, in file order."""
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
