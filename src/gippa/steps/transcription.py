"""Transcription-rate steps (GIPPA_SPEC.md §6.2, §7.1.1, §7.1.4; M2-14): ``promoter_calculator``
(ext) and ``tx_layer``.

``promoter_calculator`` runs Promoter Calculator (GPL, never imported) once per
``tu_upstream/<tu_id>.fa`` and writes ``<tu_id>.csv``; ``tx_layer`` maps the kept hits of every CSV
(the file stem is the TU id) to genome TSSs and writes ``tx_rates.parquet`` (the precomputed
``tx_rates`` columns plus ``tu_id``), in TU order whatever the order of the CSVs.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from gippa.annotate.operons import Tu
from gippa.annotate.promoter_calculator import parse_pc_csv, tx_layer
from gippa.io.parquet import write_parquet_deterministic
from gippa.model.coords import Interval
from gippa.steps.common import one_path as _path
from gippa.steps.common import read_replicons as _replicons
from gippa.steps.registry import StepContext, external, step

_STR = pa.string()
_TX_COLS = [
    ("contig", _STR), ("strand", _STR), ("tss", pa.int64()), ("tx_rate", pa.float64()),
    ("tu_id", _STR),
]  # fmt: skip


def _tus(path: Path) -> list[Tu]:
    return [
        Tu(
            tu_id=r["tu_id"],
            contig=r["contig"],
            interval=Interval(r["contig"], r["start"], r["end"], r["strand"]),
            strand=r["strand"],
            genes=tuple(r["genes"]),
            upstream=Interval(r["contig"], r["upstream_start"], r["upstream_end"], r["strand"]),
        )
        for r in pq.read_table(path).to_pylist()
    ]


external("promoter_calculator", inputs=["tu_upstream"], outputs=["pc"], precomputed="tx_rates")


@step("tx_layer", inputs=["pc", "tus", "replicons"], outputs=["tx_rates"], precomputed="tx_rates")
def tx_layer_step(ctx: StepContext) -> None:
    """``pc/<tu_id>.csv`` → ``tx_rates.parquet``."""
    pc = ctx.inputs["pc"]
    hits = {p.stem: parse_pc_csv(p) for p in sorted(pc if isinstance(pc, list) else [pc])}
    tx = tx_layer(_tus(_path(ctx.inputs["tus"])), hits, _replicons(_path(ctx.inputs["replicons"])))
    write_parquet_deterministic(ctx.outputs["tx_rates"], _TX_COLS, list(tx.itertuples(index=False)))
