"""Hit table shared by the exact and tolerant routes (GIPPA_SPEC.md §7.2; issue M3-04).

One row per peptide x locus x search. ``start``/``end`` are 0-based half-open and unwrapped
(§6.5); decoy rows carry the mirrored *target* interval, while ``frame_id``/``aa_start``/
``aa_end`` stay in the decoy frame.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import pandas as pd

HIT_COLUMNS: tuple[str, ...] = (
    "peptide",
    "query",
    "search",
    "decoy_type",
    "contig",
    "frame_id",
    "aa_start",
    "aa_end",
    "start",
    "end",
    "strand",
    "frame",
    "route",
    "n_diff",
    "diff_types",
    "blocks",
)

INT_COLUMNS: frozenset[str] = frozenset(
    {"frame_id", "aa_start", "aa_end", "start", "end", "frame", "n_diff"}
)


def hits_frame(rows: Iterable[Mapping[str, object]] | None = None) -> pd.DataFrame:
    """Build a hit table with ``HIT_COLUMNS`` in order; integer columns are ``int64``."""
    df = pd.DataFrame(list(rows or ()), columns=list(HIT_COLUMNS))
    for col in HIT_COLUMNS:
        if col in INT_COLUMNS:
            df[col] = df[col].astype("int64")
        else:
            df[col] = df[col].astype(object)
    return df
