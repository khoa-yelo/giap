"""``RunView``: every file read of the plot framework (M9-01 behaviour 3-5; M9-index decision 2).

One table, ``PATHS``, maps a key to its ``OUT/``-relative path: the ``giap.runner.local``
``PUBLISHED_FILES`` entries plus the files inside the ``PUBLISHED`` directories, preflight, logs
and the resolved config. ``from_outdir`` resolves keys under ``OUT/``; ``from_inputs`` resolves
them from a step's input map (file keys, or a directory key such as ``mapping`` whose
``OUT/`` directory holds the file), so both see the same keys and read the same frames.

TSVs that have a Frictionless schema (``giap/schemas/tables/<stem>.schema.json``) are read with
its types and its ``missingValues``; a schema column missing from the file raises ``ValueError``.
List columns stay ``;``-joined strings (M9-01 lead decision 2). Coordinates are left as written
(1-based closed, §9.2).
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from functools import cache
from importlib import resources
from pathlib import Path
from typing import Any, get_args

import pandas as pd
import yaml

from giap.config.model import Module
from giap.io._open import open_text
from giap.runner.local import PUBLISHED, PUBLISHED_FILES

_EXTRA: dict[str, str] = {
    "preflight": "preflight/preflight.json",
    "peptide_hits": "mapping/peptide_hits.tsv",
    "protein_hits": "mapping/protein_hits.tsv",
    "unmapped": "mapping/unmapped.tsv",
    "mapping_fdr": "mapping/mapping_fdr.tsv",
    "tus": "annotation/tus.tsv",
    "regions": "annotation/regions.tsv",
    "promoters": "annotation/promoters.bed",
    "terminators": "annotation/terminators.bed",
    "events": "logs/events.jsonl",
    "resolved_config": "resolved_config.yaml",
}
#: key -> ``OUT/``-relative path (behaviour 3)
PATHS: Mapping[str, str] = {**PUBLISHED_FILES, **_EXTRA}
#: further ``OUT/``-relative candidates tried after ``PATHS[key]``
_FALLBACK: dict[str, tuple[str, ...]] = {"resolved_config": ("work/resolved_config.yaml",)}
STEP_EVENTS = "logs/steps"
SUPPORT = "plots/support"
_BED_COLUMNS = ("chrom", "start", "end", "name", "score", "strand")


class NotAGiapRun(ValueError):
    """``OUT/`` has neither ``tables/`` nor ``preflight/``."""


@cache
def _schema(stem: str) -> dict[str, Any] | None:
    res = resources.files("giap.schemas").joinpath("tables", f"{stem}.schema.json")
    if not res.is_file():
        return None
    data: dict[str, Any] = json.loads(res.read_text(encoding="utf-8"))
    return data


def _typed(values: list[str], kind: str, missing: frozenset[str], where: str) -> pd.Series:
    cells = [None if v in missing else v for v in values]
    try:
        if kind == "integer":
            return pd.Series([None if v is None else int(v) for v in cells], dtype="Int64")
        if kind == "number":
            return pd.Series([float("nan") if v is None else float(v) for v in cells], dtype=float)
        if kind == "boolean":
            bools = {"true": True, "false": False}
            return pd.Series([None if v is None else bools[v] for v in cells], dtype="boolean")
    except (KeyError, ValueError) as err:
        raise ValueError(f"{where}: a value is not of type {kind}: {err}") from err
    return pd.Series(cells, dtype=object)


def _read_lines(path: Path) -> tuple[list[str], list[list[str]]]:
    with open_text(path, newline="") as fh:
        lines = fh.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise ValueError(f"{path}: empty TSV (no header)")
    header = lines[0].rstrip("\r").split("\t")
    rows = [line.rstrip("\r").split("\t") for line in lines[1:]]
    for i, row in enumerate(rows, start=2):
        if len(row) != len(header):
            raise ValueError(f"{path}: line {i} has {len(row)} fields, header has {len(header)}")
    return header, rows


def read_schema_tsv(path: Path, table: str, schema: Mapping[str, Any]) -> pd.DataFrame:
    """A TSV typed by its Frictionless ``schema``; schema columns first, in schema order."""
    header, rows = _read_lines(path)
    fields = {f["name"]: str(f.get("type", "string")) for f in schema["fields"]}
    for name in fields:
        if name not in header:
            raise ValueError(f"table {table} ({path}): schema column {name!r} is missing")
    missing = frozenset(schema.get("missingValues", [""]))
    cols = {
        name: _typed([r[i] for r in rows], fields.get(name, "string"), missing, f"{table}.{name}")
        for i, name in enumerate(header)
    }
    order = [*fields, *(c for c in header if c not in fields)]
    return pd.DataFrame(cols, columns=order)


def _read_plain_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", na_values=["NA"], keep_default_na=False, compression="infer")


class RunView:
    """Read-only access to one run's published files, by ``PATHS`` key."""

    def __init__(self, files: Mapping[str, Path], dirs: Mapping[str, Path]) -> None:
        self._files = dict(files)
        self._dirs = dict(dirs)  # STEP_EVENTS / SUPPORT -> directory

    @classmethod
    def from_outdir(cls, out: Path) -> RunView:
        if not ((out / "tables").is_dir() or (out / "preflight").is_dir()):
            raise NotAGiapRun(f"{out} is not a GIAP output directory (no tables/ or preflight/)")
        files = {}
        for key, rel in PATHS.items():
            for cand in (rel, *_FALLBACK.get(key, ())):
                if (out / cand).is_file():
                    files[key] = out / cand
                    break
        return cls(files, {d: out / d for d in (STEP_EVENTS, SUPPORT)})

    @classmethod
    def from_inputs(cls, inputs: Mapping[str, Path]) -> RunView:
        """From a step's input map: a file key, or a ``PUBLISHED`` directory key holding it;
        the directory keys ``step_events`` and ``plots_support`` stand for ``logs/steps`` and
        ``plots/support``."""
        dirs = {sub: inputs[k] for k, sub in PUBLISHED.items() if k in inputs}
        files: dict[str, Path] = {}
        for key, rel in PATHS.items():
            head, _, rest = rel.partition("/")
            if key in inputs:
                cand = inputs[key]
            elif head in dirs:
                cand = dirs[head] / rest
            else:
                continue
            if cand.is_file():
                files[key] = cand
        extra = {STEP_EVENTS: "step_events", SUPPORT: "plots_support"}
        return cls(files, {d: inputs[k] for d, k in extra.items() if k in inputs})

    def path(self, key: str) -> Path | None:
        """The file behind ``key``, or ``None`` when it is absent; unknown keys raise."""
        if key not in PATHS:
            raise KeyError(f"unknown run file key {key!r}")
        return self._files.get(key)

    def has(self, key: str) -> bool:
        return self.path(key) is not None

    def table(self, key: str) -> pd.DataFrame | None:
        """The table behind ``key``; ``None`` when absent, zero rows when header-only."""
        path = self.path(key)
        if path is None:
            return None
        rel = PATHS[key]
        if rel.endswith(".parquet"):
            return pd.read_parquet(path)
        if rel.endswith(".bed"):
            frame = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)
            frame.columns = list(_BED_COLUMNS[: frame.shape[1]]) + [
                f"col{i}" for i in range(len(_BED_COLUMNS), frame.shape[1])
            ]
            return frame.astype({"start": "int64", "end": "int64"})
        if not rel.endswith((".tsv", ".tsv.gz")):
            raise ValueError(f"{key} ({rel}) is not a table")
        stem = Path(rel).name.split(".")[0]
        schema = _schema(stem)
        return _read_plain_tsv(path) if schema is None else read_schema_tsv(path, stem, schema)

    def preflight(self) -> Mapping[str, Any] | None:
        path = self.path("preflight")
        if path is None:
            return None
        with open_text(path) as fh:
            data: dict[str, Any] = json.load(fh)
        return data

    def support(self, name: str) -> pd.DataFrame | None:
        """``plots/support/<name>.tsv[.gz]`` (M9-03), or ``None``."""
        base = self._dirs.get(SUPPORT)
        if base is None:
            return None
        for suffix in (".tsv", ".tsv.gz"):
            path = base / f"{name}{suffix}"
            if path.is_file():
                return _read_plain_tsv(path)
        return None

    def _event_files(self) -> list[Path]:
        found = [p for p in (self._files.get("events"),) if p is not None]
        steps = self._dirs.get(STEP_EVENTS)
        if steps is not None and steps.is_dir():
            found += sorted(steps.glob("*.events.jsonl"))
        return found

    def run_flags(self) -> frozenset[str]:
        """Every ``data.flag`` and ``data.flags`` of every events file (behaviour 5)."""
        flags: set[str] = set()
        for path in self._event_files():
            with open_text(path) as fh:
                for line in fh:
                    if not line.strip():
                        continue
                    data = json.loads(line).get("data")
                    if not isinstance(data, dict):
                        continue
                    if isinstance(data.get("flag"), str):
                        flags.add(data["flag"])
                    if isinstance(data.get("flags"), list):
                        flags.update(str(f) for f in data["flags"])
        return frozenset(flags)

    def modules_off(self) -> frozenset[str]:
        """The optional modules a readable ``resolved_config.yaml`` leaves out; ∅ otherwise."""
        path = self._files.get("resolved_config")
        if path is None:
            return frozenset()
        with open_text(path) as fh:
            cfg = yaml.safe_load(fh)
        modules = cfg.get("modules") if isinstance(cfg, dict) else None
        if not isinstance(modules, list):
            return frozenset()
        return frozenset(set(get_args(Module)) - set(modules) - {"core"})
