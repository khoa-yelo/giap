"""``tables/pileup.tsv.gz`` and ``tables/variants.tsv`` (GIPPA_SPEC.md §7.5.5, §9.3, §6.5;
docs/specs/M7-03-variant-tables-and-vcf.md rule 10).

Frames hold 0-based ``codon_start`` wrapped into ``[0, length)``; the files hold it 1-based
closed. Only the §9.3 columns are written, by name, so internal pileup columns never leak.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from gippa.io.tsv import read_tsv, write_tsv
from gippa.variants.pileup import PILEUP_COLUMNS
from gippa.variants.vcf_records import UNESTIMATED, VARIANTS_COLUMNS

__all__ = ["read_variants", "write_pileup", "write_variants"]

_INT = ("codon_start", "n_ref", "n_alt")
_NULLABLE_INT = ("aa_pos", "nt_changes", "class_error_n")


def _write(path: Path, frame: pd.DataFrame, columns: Sequence[str]) -> None:
    sel = frame[list(columns)].copy()
    sel["codon_start"] = sel["codon_start"].astype("int64") + 1
    rows = sel.astype(object).where(sel.notna(), None).itertuples(index=False, name=None)
    write_tsv(path, columns, rows)


def write_pileup(path: Path, pile: pd.DataFrame) -> None:
    """Write ``pileup.tsv.gz`` (gzip with ``mtime = 0`` when ``path`` ends in ``.gz``)."""
    _write(path, pile, PILEUP_COLUMNS)


def write_variants(path: Path, sites: pd.DataFrame) -> None:
    """Write ``variants.tsv`` from the table returned by ``vcf_records``."""
    _write(path, sites, VARIANTS_COLUMNS)


# required residue columns: a two-residue block can read "NA" (Asn-Ala); the schema gives them
# no missing values (field-level ``missingValues: []``)
_RESIDUES = ("ref_aa", "alt_aa")


def _rate(v: str | None) -> object:
    if v is None:
        return pd.NA
    return v if v == UNESTIMATED else float(v)


def read_variants(path: Path) -> pd.DataFrame:
    """Read ``variants.tsv`` back into the frame :func:`write_variants` takes (its inverse)."""
    header, rows = read_tsv(path)
    if tuple(header) != VARIANTS_COLUMNS:
        raise ValueError(f"{path}: columns {header} are not the §9.3 variants columns")
    df = pd.DataFrame(rows, columns=list(VARIANTS_COLUMNS), dtype=object)
    for col in _INT:
        df[col] = df[col].astype("int64")
    df["codon_start"] -= 1
    for col in _NULLABLE_INT:
        df[col] = pd.array([pd.NA if v is None else int(v) for v in df[col]], dtype="Int64")
    df["class_error_rate"] = [_rate(v) for v in df["class_error_rate"]]
    df["flags"] = [() if v is None else tuple(v.split(";")) for v in df["flags"]]
    for col in ("gene", "alt_codon_vcf"):
        df[col] = [pd.NA if v is None else v for v in df[col]]
    for col in _RESIDUES:  # never missing: "NA" here is Asn-Ala, not the ADR-0006 token
        df[col] = ["NA" if v is None else v for v in df[col]]
    return df
