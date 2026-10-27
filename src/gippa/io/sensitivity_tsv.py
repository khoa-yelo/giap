"""Writer for ``tables/sensitivity.tsv`` (GIPPA_SPEC.md §8.6, §2.4 G9; docs/specs/M4-19).

Columns ``gippa.context.sensitivity.COLUMNS``; counts as integers, ``jaccard`` / ``agreement``
through M1-04 ``write_tsv`` (6 significant digits), missing values ``NA``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from gippa.context.sensitivity import COLUMNS
from gippa.io.tsv import write_tsv

__all__ = ["write_sensitivity_tsv"]


def _float(v: object) -> float | None:
    return None if v is None or pd.isna(v) else float(v)  # type: ignore[arg-type]


def write_sensitivity_tsv(path: Path, df: pd.DataFrame) -> None:
    """Write ``sensitivity_rows`` output to ``path``."""
    rows = [
        [str(r.family), str(r.variant), int(r.n_primary), int(r.n_variant), int(r.n_shared),
         _float(r.jaccard), _float(r.agreement)]
        for r in df[list(COLUMNS)].itertuples(index=False)
    ]  # fmt: skip
    write_tsv(path, COLUMNS, rows)
