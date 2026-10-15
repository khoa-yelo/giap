"""Database items of ``giap setup``: official sources, the plan, and offline verification.

GIAP_SPEC.md §3.2 (sources), §3.9 (checksums), §7.9.1; ADR-0014 (pins). ``verify`` never touches
the network: it reads the version files the downloaders leave and hashes pinned archives.
"""

from __future__ import annotations

import hashlib
import json
import re
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
    "check_file",
    "count_fasta_records",
    "file_digest",
    "plan",
    "tree_size",
    "verify",
    "verify_go",
    "verify_kegg_links",
]

Status = Literal["ok", "absent", "mismatch"]
BaktaChoice = Literal["light", "full", "skip"]
PROMOTECH_ARCHIVE = "RF-HOT.zip"
PROMOTECH_MODEL = "RF-HOT.model"
FILE_CHECKS = frozenset({"go", "kegg", "sha256", "fasta"})  # single-file items (M2-09)
RESOURCES = (
    ("go", ("go_basic", "goslim_generic")),
    ("kegg", ("kegg_ko_pathway", "kegg_pathway_names")),
    ("cog", ("cog2020_def", "cog2020_fun")),
    ("contaminants", ("contaminants",)),
)


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
            version=str(row.get("version") or ""),
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
    names += [n for flag, group in RESOURCES if getattr(opts, flag) for n in group]
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


_GO_VERSION = re.compile(r"data-version: (?:go/)?releases/(\d{4}-\d{2}-\d{2})(?:/\S+)?")
_KO_LINK = re.compile(r"ko:K\d{5}\tpath:(?:map|ko)\d{5}")
_PATHWAY = re.compile(r"map\d{5}\t\S.*")


def verify_go(path: Path, release: str) -> bool:
    """True when the OBO header announces GO release ``release`` (either live header form)."""
    try:
        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip() or line.startswith("["):
                    return False  # end of the header
                if line.startswith("data-version:"):
                    m = _GO_VERSION.fullmatch(line.strip())
                    return m is not None and m.group(1) == release
    except OSError:
        return False
    return False


def _rows(path: Path) -> list[str]:
    try:
        return [ln for ln in path.read_text(encoding="utf-8", errors="replace").splitlines() if ln]
    except OSError:
        return []


def verify_kegg_links(path: Path) -> bool:
    """True for a non-empty KO -> pathway table with at least one ``path:map`` link."""
    rows = _rows(path)
    return all(_KO_LINK.fullmatch(r) for r in rows) and any("\tpath:map" in r for r in rows)


def count_fasta_records(path: Path) -> int:
    """Number of ``>`` header lines."""
    with path.open("rb") as fh:
        return sum(1 for line in fh if line.startswith(b">"))


def _content_problem(item: DbItem, f: Path) -> str:
    kind, _, arg = item.check.partition(":")
    if kind == "go" and not verify_go(f, item.version):
        return f"version mismatch (expected releases/{item.version})"
    if kind == "kegg" and arg == "links" and not verify_kegg_links(f):
        return "not a KEGG KO -> pathway link table"
    if (
        kind == "kegg"
        and arg == "names"
        and not ((rows := _rows(f)) and all(map(_PATHWAY.fullmatch, rows)))
    ):
        return "not a KEGG pathway list"
    if kind == "fasta" and (n := count_fasta_records(f)) != int(arg):
        return f"{n} records, expected {arg} records"
    return ""


def check_file(item: DbItem, f: Path, expected: str, version: str) -> DbRecord:
    """Verify one file item at ``f``: present, sha256 equal to ``expected`` (when given), and its
    content rule (GO release header, KEGG table shape, FASTA record count)."""
    if not f.is_file():
        return DbRecord(item.name, version, item.url, item.checksum, 0, "", item.path, "absent",
                        f"no {item.path}")  # fmt: skip
    digest = file_digest(f, "sha256")
    problem = (f"{item.path} is {digest}, expected {expected}" if expected not in ("", digest)
               else _content_problem(item, f))  # fmt: skip
    if problem:
        return DbRecord(item.name, version, item.url, item.checksum, 0, "", item.path, "mismatch",
                        problem)  # fmt: skip
    return DbRecord(item.name, version, item.url, digest, f.stat().st_size, today(), item.path,
                    "ok")  # fmt: skip


def verify(item: DbItem, db_dir: Path) -> DbRecord:
    """Verify ``item`` where it lives under ``db_dir`` (no network). Unpinned files (KEGG) are
    checked against the checksum recorded in ``db_manifest.json`` when they were retrieved."""
    if item.check.partition(":")[0] not in FILE_CHECKS:
        return check_dir(item, db_dir / item.path)
    if item.checksum:
        return check_file(item, db_dir / item.path, item.checksum, item.version)
    from giap.setupdb.manifest import read_records  # manifest imports this module

    prev = read_records(db_dir).get(item.name)
    if prev is None or prev.status != "ok" or (prev.url, prev.path) != (item.url, item.path):
        return DbRecord(item.name, item.version, item.url, "", 0, "", item.path, "absent",
                        f"no recorded checksum for {item.path}")  # fmt: skip
    return check_file(item, db_dir / item.path, prev.checksum, prev.version)
