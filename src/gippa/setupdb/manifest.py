"""``DIR/db_manifest.json`` (GIPPA_SPEC.md §7.9.1) and the prefetch log it adopts provenance from.

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
from collections.abc import Mapping, Sequence
from pathlib import Path
from urllib.parse import urlparse

from gippa.setupdb.items import DbRecord

__all__ = ["MANIFEST", "PREFETCH_LOG", "SCHEMA_VERSION", "read_prefetch", "read_records",
           "write_manifest"]  # fmt: skip

MANIFEST = "db_manifest.json"
PREFETCH_LOG = "prefetch_manifest.tsv"
SCHEMA_VERSION = "1.0"
_FIELDS = {f.name for f in dataclasses.fields(DbRecord)}
LOCAL_BUILD = "local_build"  # optional object of a locally built image (M11-14 item 1b)


def _read_items(db_dir: Path) -> dict[str, tuple[DbRecord, dict[str, str] | None]]:
    """Valid items of an existing manifest by name, with their ``local_build`` object."""
    try:
        data = json.loads((db_dir / MANIFEST).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    items = data.get("items") if isinstance(data, dict) else None
    out: dict[str, tuple[DbRecord, dict[str, str] | None]] = {}
    for raw in items if isinstance(items, list) else []:
        if not isinstance(raw, dict) or not set(raw) <= _FIELDS | {LOCAL_BUILD}:
            continue
        build = raw.get(LOCAL_BUILD)
        if build is not None and not (
            isinstance(build, dict) and all(isinstance(v, str) for v in build.values())
        ):
            continue
        with contextlib.suppress(TypeError):
            rec = DbRecord(**{k: v for k, v in raw.items() if k != LOCAL_BUILD})
            out[rec.name] = (rec, build)
    return out


def read_records(db_dir: Path) -> dict[str, DbRecord]:
    """Records of an existing manifest by name; ``{}`` when missing or not in this format."""
    return {name: rec for name, (rec, _build) in _read_items(db_dir).items()}


def _item(rec: DbRecord, build: Mapping[str, str] | None) -> dict[str, object]:
    out: dict[str, object] = dataclasses.asdict(rec)
    return out if build is None else out | {LOCAL_BUILD: dict(build)}


def write_manifest(
    db_dir: Path,
    records: Sequence[DbRecord],
    local_builds: Mapping[str, Mapping[str, str]] | None = None,
) -> None:
    """Merge ``records`` into ``db_dir/db_manifest.json`` (same name: the new record wins).
    ``local_builds`` gives a record its ``local_build`` object; a kept item keeps its own."""
    builds = local_builds or {}
    merged: dict[str, tuple[DbRecord, Mapping[str, str] | None]] = {
        **_read_items(db_dir), **{r.name: (r, builds.get(r.name)) for r in records}}  # fmt: skip
    doc = {
        "schema_version": SCHEMA_VERSION,
        "items": [_item(*merged[n]) for n in sorted(merged)],
    }
    text = json.dumps(doc, indent=2, ensure_ascii=False) + "\n"
    fd, tmp = tempfile.mkstemp(prefix=f".{MANIFEST}.", suffix=".tmp", dir=db_dir)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.chmod(tmp, 0o644)  # shared DB directory, whatever the umask (M2-09 review 6)
        os.replace(tmp, db_dir / MANIFEST)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


def read_prefetch(db_dir: Path) -> dict[str, tuple[str, str, str]]:
    """``{item: (checksum, ISO date, source host)}`` from ``prefetch_manifest.tsv``; ``{}`` when
    missing."""
    try:
        text = (db_dir / PREFETCH_LOG).read_text(encoding="utf-8")
    except OSError:
        return {}
    out: dict[str, tuple[str, str, str]] = {}
    for row in csv.DictReader(text.splitlines(), delimiter="\t"):
        name, checksum, date = row.get("item"), row.get("checksum"), row.get("date")
        if name and checksum and date:
            out[name] = (checksum, date[:10], urlparse(row.get("source") or "").hostname or "")
    return out
