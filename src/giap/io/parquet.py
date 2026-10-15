"""Deterministic Parquet writer (GIAP_SPEC.md §6.6, TE7).

Tables are written without column statistics and with the schema metadata given by the caller
only, so two writes of the same rows give identical bytes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]


def write_parquet_deterministic(
    path: Path,
    cols: Sequence[tuple[str, pa.DataType]],
    rows: Sequence[Sequence[Any]],
    metadata: Mapping[str, str] | None = None,
) -> None:
    """Write ``rows`` (tuples in ``cols`` order) to ``path`` as a byte-reproducible Parquet file.

    ``metadata`` (optional) is stored as key/value schema metadata in the file footer.
    """
    schema = pa.schema(
        [pa.field(name, typ) for name, typ in cols],
        metadata=dict(sorted(metadata.items())) if metadata else None,
    )
    arrays = [pa.array([r[i] for r in rows], type=t) for i, (_, t) in enumerate(cols)]
    pq.write_table(pa.Table.from_arrays(arrays, schema=schema), path, write_statistics=False)
