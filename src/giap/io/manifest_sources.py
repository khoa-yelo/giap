"""File readers and the strict writer behind ``run_manifest.json`` (MVP-03 D-B, D-C; ADR-0040).

Everything here touches the file system; :mod:`giap.manifest` stays pure. MVP-03c adds the
``versions.yml`` (``tools``), ``db_manifest.json`` (``databases``), ``sensitivity.tsv`` and
``island_diagnostics.parquet`` readers.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq  # type: ignore[import-untyped]
import yaml

from giap.manifest import canonical_json
from giap.runner.nextflow import _task_dir

LOGGER = logging.getLogger(__name__)

_DIGEST = re.compile(r"sha256[:-]([0-9a-f]{64})")
_COUNTS = ("n_primary", "n_variant", "n_shared")


def rotate_events(logs_dir: Path) -> None:
    """``events.jsonl`` → ``events.previous.jsonl`` and ``trace.tsv`` → ``trace.previous.tsv``
    (one generation kept; lead decision 5, review B1)."""
    for cur in (logs_dir / "events.jsonl", logs_dir / "trace.tsv"):
        if cur.is_file():
            cur.replace(cur.with_stem(f"{cur.stem}.previous"))


def read_events(path: Path) -> list[dict[str, Any]]:
    """The events of ``path`` in file order (``[]`` if absent); a non-object line raises
    ``ValueError`` (``json.JSONDecodeError`` is one): events.jsonl is GIAP's own output."""
    if not path.is_file():
        return []
    rows = [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not all(isinstance(r, dict) for r in rows):
        raise ValueError(f"{path}: a line is not a JSON object")
    return rows


def read_json(path: Path | None) -> Any:
    """Parsed JSON of ``path``; None when absent or empty (a stub-run output)."""
    if path is None or not path.is_file() or not path.read_text(encoding="utf-8").strip():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def read_trace(path: Path) -> list[dict[str, str]] | None:
    """``trace.tsv`` rows (human-readable units), or None when absent."""
    if not path.is_file():
        return None
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def input_row(role: str, path: Path) -> dict[str, Any]:
    """``{role, file, path, sha256, size}`` of one input file."""
    return {"role": role, "file": path.name, "path": str(path), "sha256": _sha256(path),
            "size": path.stat().st_size}  # fmt: skip


def git_info(root: Path) -> tuple[str | None, bool | None]:
    """``(HEAD, dirty)`` when ``root`` is the top of a git checkout, else ``(None, None)``."""

    def git(*args: str) -> str | None:
        try:
            res = subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                                 text=True, check=False)  # fmt: skip
        except OSError:
            return None
        return res.stdout.strip() if res.returncode == 0 else None

    top = git("rev-parse", "--show-toplevel")
    if top is None or Path(top).resolve() != root.resolve():
        return None, None
    head, status = git("rev-parse", "HEAD"), git("status", "--porcelain")
    if head is None or status is None:
        return None, None
    return head, bool(status)


def write_manifest(path: Path, manifest: Mapping[str, Any]) -> None:
    """Write ``canonical_json(manifest)``; NaN/±∞ raise ``ValueError`` before any write."""
    text = canonical_json(manifest)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _nonempty(path: Path) -> bool:
    """False for a missing or 0-byte file (a ``-stub-run`` ``touch`` output)."""
    return path.is_file() and path.stat().st_size > 0


def read_sensitivity(path: Path | None) -> list[dict[str, Any]] | None:
    """Rows of ``tables/sensitivity.tsv`` (counts ``int``, other numbers ``float``, ``NA`` →
    None); ``[]`` when header-only, None when absent or empty."""
    if path is None or not _nonempty(path):
        return None
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in rows:
        num = {k: None if r[k] in ("", "NA") else (int if k in _COUNTS else float)(r[k])
               for k in (*_COUNTS, "jaccard", "agreement")}  # fmt: skip
        out.append({"family": r["family"], "variant": r["variant"], **num})
    return out


def read_island_diagnostics(path: Path | None) -> list[dict[str, Any]] | None:
    """Rows of ``context/island_diagnostics.parquet`` in file (C1) order; None if absent/empty."""
    if path is None or not _nonempty(path):
        return None
    rows: list[dict[str, Any]] = pq.read_table(path).to_pylist()
    return rows


def read_db_manifest(db_dir: Path | None) -> dict[str, Any] | None:
    """``<db_dir>/db_manifest.json`` (either shape), None without a directory or file."""
    doc = read_json(db_dir / "db_manifest.json") if db_dir is not None else None
    return doc if isinstance(doc, dict) else None


def _yaml(path: Path) -> Any:
    """``versions.yml`` with every scalar a string (``1.10`` stays ``1.10``); ``None`` when it is
    unreadable (a tool's stderr in its ``--version`` line), so the run tail never fails on it."""
    try:
        return yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    except (yaml.YAMLError, UnicodeDecodeError) as exc:
        LOGGER.warning("unreadable versions.yml skipped: %s (%s)", path, type(exc).__name__)
        return None


def tool_reports(
    trace: Sequence[Mapping[str, str]], work: Path, locks: Path
) -> list[dict[str, Any]]:
    """One ``{process, name, version, container_digest, conda_lock_hash}`` per tool line of the
    ``versions.yml`` of each COMPLETED/CACHED task (lead decision 9); a missing file is skipped."""
    out: list[dict[str, Any]] = []
    for row in trace:
        versions = Path(_task_dir(row, work)) / "versions.yml"
        done = (row.get("status") or "").strip() in {"COMPLETED", "CACHED"}
        if not done or not versions.is_file():
            continue
        m = _DIGEST.search(row.get("container") or "")
        doc = _yaml(versions)
        for process, tools in (doc if isinstance(doc, dict) else {}).items():
            for name, version in (tools if isinstance(tools, dict) else {}).items():
                lock = locks / f"{name}.explicit.txt"
                out.append({"process": str(process), "name": str(name), "version": str(version),
                            "container_digest": f"sha256:{m.group(1)}" if m else None,
                            "conda_lock_hash": _sha256(lock) if not m and lock.is_file()
                            else None})  # fmt: skip
    return out
