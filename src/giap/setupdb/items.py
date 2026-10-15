"""Database items of ``giap setup``: official sources, the plan, and offline verification.

GIAP_SPEC.md §3.2 (sources), §3.9 (checksums), §7.9.1; ADR-0014 (pins). ``verify`` never touches
the network: it reads the version files the downloaders leave and hashes pinned archives.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path
from types import MappingProxyType
from typing import Literal
from urllib.parse import urlparse

import yaml

__all__ = [
    "SOURCES",
    "DbItem",
    "DbRecord",
    "SetupOptions",
    "Status",
    "check_dir",
    "file_digest",
    "plan",
    "tree_size",
    "verify",
]

Status = Literal["ok", "absent", "mismatch"]
BaktaChoice = Literal["light", "full", "skip"]
PROMOTECH_ARCHIVE = "RF-HOT.zip"
PROMOTECH_MODEL = "RF-HOT.model"


@dataclass(frozen=True, slots=True)
class DbItem:
    """One downloadable database or model (a row of ``giap/data/db_sources.yaml``)."""

    name: str
    version: str
    url: str
    host: str
    checksum: str  # "md5:<hex>" or "sha256:<hex>" of the archive as published; "" if none
    path: str  # relative to db_dir
    tool: str | None  # TOOLS key whose environment runs the downloader; None = HTTPS fetch
    check: str = ""  # verification rule: "bakta:light", "bakta:full", "genomad", ...
    size_gib: int = 0  # free disk needed before downloading
    prefetch: str | None = None  # row name in prefetch_manifest.tsv


@dataclass(frozen=True, slots=True)
class DbRecord:
    """One ``db_manifest.json`` item."""

    name: str
    version: str
    url: str
    checksum: str
    size_bytes: int
    retrieved: str  # ISO date, "" when absent
    path: str
    status: Status
    message: str = ""


@dataclass(frozen=True, slots=True)
class SetupOptions:
    """``giap setup`` options (GIAP_SPEC.md §5.1)."""

    db_dir: Path
    bakta: BaktaChoice = "light"
    genomad: bool = True
    antismash: bool = False
    promotech_models: bool = True
    check: bool = False
    go: bool = False
    kegg: bool = False
    cog: bool = False
    contaminants: bool = False


def _load_sources() -> Mapping[str, DbItem]:
    raw = yaml.safe_load(files("giap.data").joinpath("db_sources.yaml").read_text("utf-8"))
    out: dict[str, DbItem] = {}
    for name, row in raw.items():
        url = str(row["url"])
        out[name] = DbItem(
            name=name,
            version=str(row["version"]),
            url=url,
            host=urlparse(url).hostname or "",
            checksum=str(row.get("checksum") or ""),
            path=str(row["path"]),
            tool=row.get("tool"),
            check=str(row["check"]),
            size_gib=int(row.get("size_gib", 0)),
            prefetch=row.get("prefetch"),
        )
    return MappingProxyType(out)


SOURCES: Mapping[str, DbItem] = _load_sources()


def plan(opts: SetupOptions) -> list[DbItem]:
    """Items the options ask for, in a fixed order."""
    names = [f"bakta_{opts.bakta}"] if opts.bakta != "skip" else []
    names += [n for n, on in (("genomad", opts.genomad), ("antismash", opts.antismash),
                              ("promotech_models", opts.promotech_models)) if on]  # fmt: skip
    return [SOURCES[n] for n in names]


def today() -> str:
    return datetime.now(UTC).date().isoformat()


def file_digest(path: Path, algo: str) -> str:
    """``"<algo>:<hex>"`` of a file, read in 8 MiB blocks."""
    h = hashlib.new(algo)
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return f"{algo}:{h.hexdigest()}"


def tree_size(path: Path) -> int:
    """Bytes of the regular files under ``path`` (symlinks not followed)."""
    if path.is_file():
        return path.stat().st_size
    return sum(p.lstat().st_size for p in path.rglob("*") if p.is_file() and not p.is_symlink())


def _record(item: DbItem, status: Status, version: str, message: str, d: Path) -> DbRecord:
    ok = status == "ok"
    return DbRecord(item.name, version, item.url, item.checksum, tree_size(d) if ok else 0,
                    today() if ok else "", item.path, status, message)  # fmt: skip


def _check_bakta(item: DbItem, d: Path) -> tuple[Status, str, str]:
    try:
        info = json.loads((d / "version.json").read_text(encoding="utf-8"))
        major, minor, kind = int(info["major"]), int(info["minor"]), str(info["type"])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return "absent", item.version, f"no readable {item.path}/version.json ({exc})"
    want = item.check.split(":", 1)[1]
    if major != int(item.version.split(".")[0]) or kind != want:
        return (
            "absent",
            f"{major}.{minor}",
            f"found v{major}.{minor} {kind}, need v{item.version} {want}",
        )
    return "ok", f"{major}.{minor}", ""


def _check_genomad(item: DbItem, d: Path) -> tuple[Status, str, str]:
    try:
        found = (d / "version.txt").read_text(encoding="utf-8").strip()
    except OSError:
        return "absent", item.version, f"no {item.path}/version.txt"
    if found != item.version and not found.startswith(item.version + "."):
        return "absent", found, f"found v{found}, need v{item.version}"
    return "ok", found, ""


def _check_promotech(item: DbItem, d: Path) -> tuple[Status, str, str]:
    archive = d / PROMOTECH_ARCHIVE
    if not (d / PROMOTECH_MODEL).is_file() or not archive.is_file():
        return "absent", item.version, f"no {PROMOTECH_MODEL} with its {PROMOTECH_ARCHIVE}"
    got = file_digest(archive, item.checksum.split(":", 1)[0])
    if got != item.checksum:
        return "mismatch", item.version, f"{PROMOTECH_ARCHIVE} is {got}, pinned {item.checksum}"
    return "ok", item.version, ""


def _check_antismash(item: DbItem, d: Path) -> tuple[Status, str, str]:
    # Presence only; `antismash --check-prereqs` in its environment is M2-09.
    if d.is_dir() and any(d.iterdir()):
        return "ok", item.version, ""
    return "absent", item.version, f"no {item.path}"


_CHECKS = {
    "bakta": _check_bakta,
    "genomad": _check_genomad,
    "promotech": _check_promotech,
    "antismash": _check_antismash,
}


def check_dir(item: DbItem, d: Path) -> DbRecord:
    """Verify an installed copy of ``item`` in directory ``d`` (no network)."""
    status, version, message = _CHECKS[item.check.split(":", 1)[0]](item, d)
    return _record(item, status, version, message, d)


def verify(item: DbItem, db_dir: Path) -> DbRecord:
    """Verify ``item`` where it lives under ``db_dir`` (no network)."""
    return check_dir(item, db_dir / item.path)
