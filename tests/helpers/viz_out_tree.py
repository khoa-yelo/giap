"""Small synthetic ``OUT/`` trees for the M9 plot framework tests (M9-01).

Tables are written through ``gippa.io.tsv.write_tsv`` with the columns of their Frictionless
schema (``src/gippa/schemas/tables/<name>.schema.json``); a cell left out of a row dict is ``NA``.
Nothing here touches ``$GIPPA_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml

from gippa.io.tsv import write_tsv

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "src" / "gippa" / "schemas" / "tables"


def schema_fields(schema: str) -> list[dict[str, Any]]:
    """The ``fields`` of ``<schema>.schema.json``."""
    data = json.loads((SCHEMAS / f"{schema}.schema.json").read_text(encoding="utf-8"))
    fields: list[dict[str, Any]] = data["fields"]
    return fields


def schema_columns(schema: str) -> list[str]:
    return [f["name"] for f in schema_fields(schema)]


def write_schema_table(
    path: Path,
    schema: str,
    rows: Sequence[Mapping[str, object]] = (),
    drop: Iterable[str] = (),
) -> Path:
    """Write ``rows`` under the schema's header (minus ``drop``); absent cells are ``NA``."""
    dropped = set(drop)
    cols = [c for c in schema_columns(schema) if c not in dropped]
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, cols, [[row.get(c) for c in cols] for row in rows])
    return path


def write_events(path: Path, records: Iterable[Mapping[str, Any] | str]) -> Path:
    """One JSON object per line; a ``str`` record is written verbatim (for ``-Infinity``)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [r if isinstance(r, str) else json.dumps(dict(r), sort_keys=True) for r in records]
    path.write_text("".join(line + "\n" for line in lines), encoding="utf-8")
    return path


def write_resolved_config(path: Path, modules: Sequence[str]) -> Path:
    """A full resolved config (as ``gippa.runner.nextflow`` writes it) with ``modules``."""
    from gippa.config import resolve_config

    cfg = resolve_config(overrides=[f"modules=[{', '.join(modules)}]"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg.model_dump(mode="json"), sort_keys=True), encoding="utf-8")
    return path


def make_out(tmp: Path, name: str = "OUT") -> Path:
    """An ``OUT/`` that counts as a GIPPA run: it has ``tables/``."""
    out = tmp / name
    (out / "tables").mkdir(parents=True, exist_ok=True)
    return out


ISLAND_ROW: dict[str, object] = {
    "island_id": "ISL0001",
    "direction": "high",
    "contig": "chr1",
    "start": 1001,
    "end": 5400,
    "tu_first": "TU0001",
    "tu_last": "TU0003",
    "n_tus": 3,
    "n_genes": 5,
    "genes": "g1;g2;g3;g4;g5",
    "mean_abundance_log2": 21.5,
    "detected_frac": 0.8,
    "scan_s": 4.25,
    "block_len": 3,
    "stat": "scan_s",
    "effect": 1.5,
    "ci_low": None,
    "ci_high": None,
    "p_value": 0.001,
    "q_value": None,
    "method": "circular_shift_scan",
    "null_model": "circular_shift",
    "n_perm": 999,
    "flags": None,
}


def gene_row(**overrides: object) -> dict[str, object]:
    """One ``genes.tsv`` row: every boolean column filled, the rest ``NA`` unless overridden."""
    row: dict[str, object] = {}
    for f in schema_fields("genes"):
        if f["type"] == "boolean":
            row[f["name"]] = False
    row.update(
        locus_tag="g1",
        contig="chr1",
        start=101,
        end=400,
        strand="+",
        feature_type="CDS",
        detection_status="detected",
        detected=True,
        cog_category="KL",
    )
    row.update(overrides)
    return row
