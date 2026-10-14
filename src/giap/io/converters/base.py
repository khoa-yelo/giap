"""Shared result type and helpers for converters (GIAP_SPEC.md §5.4)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import pandas as pd

from giap.io.issues import Issue
from giap.io.peptides import PEPTIDE_COLUMNS

__all__ = ["SORT_KEYS", "ConvertResult", "finish_peptides"]

SORT_KEYS: tuple[str, ...] = ("run", "peptide", "modified_peptide", "charge")


@dataclass(frozen=True, slots=True)
class ConvertResult:
    """Converter output: the generic peptide table, an optional abundance table, the issues.

    ``peptides`` holds generic columns only, in §5.3 order; ``abundance`` has the columns
    ``protein_id`` and ``maxlfq``. On error issues the tables may be empty.
    """

    peptides: pd.DataFrame
    abundance: pd.DataFrame | None
    issues: list[Issue]


def finish_peptides(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Keep ``columns`` (in §5.3 order) and sort rows stably by run, peptide, modified, charge."""
    cols = [c for c in PEPTIDE_COLUMNS if c in columns]
    keys = [c for c in SORT_KEYS if c in cols]
    out = df[cols]
    if keys:
        out = out.sort_values(keys, kind="mergesort", na_position="last")
    return out.reset_index(drop=True)
