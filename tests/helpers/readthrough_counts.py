"""Read-through row counts for the M7-16 calibration checks (stressed arm, decoy-balance check).

A hit row *holds a cell* when its window ``[aa_start, aa_end)`` of frame ``frame_id`` contains a
consumable cell: a target cell for ``search == "target"``, a mirrored cell of the row's
``decoy_type`` for ``search == "decoy"`` (decoy rows carry decoy-frame aa coordinates). Frameshift
rows (``blocks`` set) never consume a cell and are never counted.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd


def holds_cell(hits: pd.DataFrame, cells: object) -> np.ndarray:
    """Boolean per row of ``hits`` (``gippa.mapping.hits`` columns) for a ``StopCells``."""
    out = np.zeros(len(hits), dtype=bool)
    if not len(hits):
        return out
    target: Sequence[np.ndarray] = cells.target  # type: ignore[attr-defined]
    decoys = cells.decoys  # type: ignore[attr-defined]
    blocks = hits["blocks"] if "blocks" in hits.columns else pd.Series([None] * len(hits))
    for i, (search, kind, fid, a, b, blk) in enumerate(
        zip(
            hits["search"],
            hits["decoy_type"],
            hits["frame_id"],
            hits["aa_start"],
            hits["aa_end"],
            blocks,
            strict=True,
        )
    ):
        if isinstance(blk, str) and blk:
            continue
        arr = target[int(fid)] if search == "target" else decoys[str(kind)][int(fid)]
        j = int(np.searchsorted(arr, int(a), side="left"))
        out[i] = j < len(arr) and int(arr[j]) < int(b)
    return out


def best_cost_rows(hits: pd.DataFrame) -> pd.DataFrame:
    """Rows at their peptide's minimal ``n_diff`` within each search."""
    if not len(hits):
        return hits
    m = hits.groupby(["peptide", "search"])["n_diff"].transform("min")
    return hits[hits["n_diff"] == m]
