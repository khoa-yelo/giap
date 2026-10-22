"""File readers and the strict writer behind ``run_manifest.json`` (MVP-03 D-B, D-C; ADR-0038).

Everything here touches the file system; :mod:`giap.manifest` stays pure.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import yaml

from giap.manifest import canonical_json
from giap.runner.nextflow import _task_dir

_DIGEST = re.compile(r"@sha256[:-]([0-9a-f]{64})")
_EMPTY = {"", "-"}


def rotate_events(logs_dir: Path) -> None:
    """``events.jsonl`` → ``events.previous.jsonl`` (one generation kept; lead decision 5)."""
    cur = logs_dir / "events.jsonl"
    if cur.is_file():
        cur.replace(logs_dir / "events.previous.jsonl")


def read_events(path: Path) -> list[dict[str, Any]]:
    """The events of ``path`` in file order (``[]`` if absent); a non-object line raises."""
    if not path.is_file():
        return []
    out: list[dict[str, Any]] = []
    for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            obj = json.loads(line)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{path}:{n}: not JSON: {exc}") from exc
        if not isinstance(obj, dict):
            raise ValueError(f"{path}:{n}: not a JSON object")
        out.append(obj)
    return out


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


def tool_reports(
    trace_rows: Sequence[Mapping[str, str]], work_dir: Path, locks_dir: Path
) -> list[dict[str, Any]]:
    """One ``{process, name, version, container_digest, conda_lock_hash}`` per tool line of the
    ``versions.yml`` of each COMPLETED/CACHED task (lead decision 9)."""
    out: list[dict[str, Any]] = []
    for row in trace_rows:
        if (row.get("status") or "").strip() not in {"COMPLETED", "CACHED"}:
            continue
        versions = Path(_task_dir(row, work_dir)) / "versions.yml"
        if not versions.is_file():
            continue
        doc = yaml.load(versions.read_text(encoding="utf-8"), Loader=yaml.BaseLoader) or {}
        container = (row.get("container") or "").strip()
        m = _DIGEST.search(container)
        digest = f"sha256:{m.group(1)}" if m else None
        for process, tools in doc.items():
            for name, version in (tools or {}).items():
                lock = locks_dir / f"{name}.explicit.txt"
                lock_hash = _sha256(lock) if container in _EMPTY and lock.is_file() else None
                out.append({"process": process, "name": name, "version": str(version),
                            "container_digest": digest, "conda_lock_hash": lock_hash})  # fmt: skip
    return out


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
