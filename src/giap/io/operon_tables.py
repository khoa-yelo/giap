"""Writers of ``tables/pairs.tsv``, ``tables/operons.tsv`` and ``tables/operon_flags.tsv``
(GIAP_SPEC.md §6.5, §9.2, §9.3; M6-08 item 3).

They take the internal 0-based half-open frames of ``merge_tests``, ``operons_table`` and
``flag_pairs``; cells follow ``giap.io.tsv`` (``NA``, ``%.6g``, ``true``/``false``, ``;``-joined
lists). Coordinates are converted here only: ``operons.tsv`` writes ``start + 1`` with ``end``
unchanged (an origin-crossing TU keeps ``end`` > replicon length). The pair tables carry no
coordinates. ``pairs.tsv`` and ``operon_flags.tsv`` are sorted by the number of ``PAIR_<n>``
(the contig/start order); ``operons.tsv`` keeps the order it is given. The inputs are not modified.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.io.tsv import write_tsv
from giap.operons.flags import FLAG_COLUMNS
from giap.operons.pairs import PAIR_COLUMNS
from giap.operons.tu_table import OPERON_COLUMNS

__all__ = ["write_operon_flags_tsv", "write_operons_tsv", "write_pairs_tsv"]


def _cell(v: object) -> object:
    if isinstance(v, np.ndarray | list | tuple):
        return [str(x) for x in v]
    return v


def _by_pair_number(df: pd.DataFrame) -> pd.DataFrame:
    num = df["pair_id"].astype(str).str.rsplit("_", n=1).str[-1].astype(int)
    return df.iloc[np.argsort(num.to_numpy(), kind="stable")]


def _write(path: Path, df: pd.DataFrame, cols: Sequence[str]) -> None:
    rows = [[_cell(v) for v in r] for r in df[list(cols)].itertuples(index=False, name=None)]
    write_tsv(path, cols, rows)


def write_pairs_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/pairs.tsv`` from the merged pair frame (an internal ``which`` column is dropped)."""
    _write(path, _by_pair_number(df) if len(df) else df, PAIR_COLUMNS)


def write_operons_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/operons.tsv`` from ``operons_table``: 1-based closed coordinates."""
    out = df.copy()
    if len(out):
        out["start"] = out["start"].astype("int64") + 1
    _write(path, out, OPERON_COLUMNS)


def write_operon_flags_tsv(path: Path, df: pd.DataFrame) -> None:
    """``tables/operon_flags.tsv`` from ``flag_pairs`` (header-only when no pair is flagged)."""
    _write(path, _by_pair_number(df) if len(df) else df, FLAG_COLUMNS)
