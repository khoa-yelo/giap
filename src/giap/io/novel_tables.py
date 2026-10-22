"""``tables/novel_orfs.tsv``, ``novel_orfs.gff3``, ``novel_orfs.faa`` and
``annotation_patch.gff3`` (docs/specs/M7-07-novel-orf-outputs.md; GIAP_SPEC.md §9.1, §9.3).

The TSV goes through M1-04 ``write_tsv``; the three positions become 1-based (``+ 1``) here and
``read_novel_orfs`` is the inverse. The GFF3 files use the M1-05 writer, the FASTA the M1-03 one.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pandas as pd

from giap.io.gff3 import write_gff3
from giap.io.protein_fasta import write_protein_fasta
from giap.io.tsv import read_tsv, write_tsv
from giap.model.coords import Replicon
from giap.model.genes import GeneModel
from giap.novel.outputs import (
    NOVEL_ORFS_COLUMNS,
    NOVEL_ORFS_TYPES,
    POSITION_COLUMNS,
    Kind,
    gff_features,
    is_na,
    novel_frame,
    patch_features,
    protein_records,
)

__all__ = ["read_novel_orfs", "write_annotation_patch", "write_novel_faa", "write_novel_gff3",
           "write_novel_orfs"]  # fmt: skip

_POS = frozenset(POSITION_COLUMNS)


def write_novel_orfs(path: Path, table: pd.DataFrame) -> None:
    """Write the ``orf_table`` frame as ``novel_orfs.tsv`` (positions 1-based)."""
    rows = []
    for r in table[list(NOVEL_ORFS_COLUMNS)].itertuples(index=False, name=None):
        rows.append([
            None if is_na(v) else int(v) + 1 if c in _POS else v
            for c, v in zip(NOVEL_ORFS_COLUMNS, r, strict=True)
        ])  # fmt: skip
    write_tsv(path, NOVEL_ORFS_COLUMNS, rows)


def _parse(v: str | None, kind: Kind, pos: bool) -> object:
    if v is None:
        return None
    if kind == "int":
        return int(v) - (1 if pos else 0)
    if kind == "float":
        return float(v)
    if kind == "bool":
        return v == "true"
    return v


def read_novel_orfs(path: Path) -> pd.DataFrame:
    """Read ``novel_orfs.tsv`` back into the ``orf_table`` representation."""
    header, rows = read_tsv(path)
    if tuple(header) != NOVEL_ORFS_COLUMNS:
        raise ValueError(f"{path}: header is not the novel_orfs.tsv column list")
    kinds = [(c, NOVEL_ORFS_TYPES[c], c in _POS) for c in header]
    recs = [{c: _parse(v, k, p) for (c, k, p), v in zip(kinds, r, strict=True)} for r in rows]
    return novel_frame(recs)


def write_novel_gff3(path: Path, table: pd.DataFrame, replicons: Mapping[str, Replicon]) -> None:
    write_gff3(path, gff_features(table), list(replicons.values()))


def write_annotation_patch(path: Path, table: pd.DataFrame, gm: GeneModel) -> None:
    write_gff3(path, patch_features(table, gm), list(gm.replicons.values()))


def write_novel_faa(
    path: Path, table: pd.DataFrame, genome: Mapping[str, str], replicons: Mapping[str, Replicon]
) -> None:
    recs = protein_records(table, genome, replicons)
    write_protein_fasta(path, [(r.header, r.sequence) for r in recs])
