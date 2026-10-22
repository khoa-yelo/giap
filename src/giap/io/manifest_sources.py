"""File readers and the strict writer behind ``run_manifest.json`` (MVP-03 D-B, D-C; ADR-0040).

Everything here touches the file system; :mod:`giap.manifest` stays pure. The ``versions.yml``
(``tools``) and ``db_manifest.json`` (``databases``) readers are MVP-03c.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from giap.manifest import canonical_json


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
