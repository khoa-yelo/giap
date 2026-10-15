"""``DIR/db_manifest.json`` (GIAP_SPEC.md §7.9.1) and the prefetch log it adopts provenance from.

The manifest is ``{"schema_version": "1.0", "items": [DbRecord, ...]}`` sorted by name; writes
merge with the existing manifest by name and are atomic (temporary file, then ``os.replace``).
"""

from __future__ import annotations

import contextlib
import csv
import dataclasses
import json
import os
import tempfile
from collections.abc import Sequence
from pathlib import Path

from giap.setupdb.items import DbRecord

__all__ = ["MANIFEST", "PREFETCH_LOG", "SCHEMA_VERSION", "read_prefetch", "read_records",
           "write_manifest"]  # fmt: skip

MANIFEST = "db_manifest.json"
PREFETCH_LOG = "prefetch_manifest.tsv"
SCHEMA_VERSION = "1.0"
_FIELDS = {f.name for f in dataclasses.fields(DbRecord)}


def read_records(db_dir: Path) -> dict[str, DbRecord]:
    """Records of an existing manifest by name; ``{}`` when missing or not in this format."""
    try:
        data = json.loads((db_dir / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    items = data.get("items") if isinstance(data, dict) else None
    out: dict[str, DbRecord] = {}
    for raw in items if isinstance(items, list) else []:
        if isinstance(raw, dict) and set(raw) <= _FIELDS:
            with contextlib.suppress(TypeError):
                rec = DbRecord(**raw)
                out[rec.name] = rec
    return out


def write_manifest(db_dir: Path, records: Sequence[DbRecord]) -> None:
    """Merge ``records`` into ``db_dir/db_manifest.json`` (same name: the new record wins)."""
    merged = read_records(db_dir) | {r.name: r for r in records}
    doc = {
        "schema_version": SCHEMA_VERSION,
        "items": [dataclasses.asdict(merged[n]) for n in sorted(merged)],
    }
    text = json.dumps(doc, indent=2, ensure_ascii=False) + "\n"
    fd, tmp = tempfile.mkstemp(prefix=f".{MANIFEST}.", suffix=".tmp", dir=db_dir)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, db_dir / MANIFEST)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def read_prefetch(db_dir: Path) -> dict[str, tuple[str, str]]:
    """``{item: (checksum, ISO date)}`` from ``prefetch_manifest.tsv``; ``{}`` when missing."""
    try:
        text = (db_dir / PREFETCH_LOG).read_text(encoding="utf-8")
    except OSError:
        return {}
    out: dict[str, tuple[str, str]] = {}
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        name, checksum, date = row.get("item"), row.get("checksum"), row.get("date")
        if name and checksum and date:
            out[name] = (checksum, date[:10])
    return out
