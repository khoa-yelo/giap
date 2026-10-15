"""Loader for the M4 input contract (docs/specs/M4-06-input-contract.md; checks in
``giap.genetable.contract``).

M2's Parquet intermediates are read as written (0-based half-open, Arrow list columns become Python
lists; a string where a list is expected is left as is and reported). M3's final ``mapping/*.tsv``
tables (1-based closed, ``;`` lists, ``NA``/``INF`` per ADR-0006) are typed from ``TABLES`` and
converted to 0-based (``start - 1``) on read. Content problems are returned as issues, never raised.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from giap.genetable.contract import TABLES, Column, validate_inputs, validate_table
from giap.io.issues import Issue, Severity
from giap.io.tsv import read_tsv

__all__ = ["FILE_NAMES", "M4_INPUT_KEYS", "M4Inputs", "default_paths", "load_m4_inputs"]


@dataclass(frozen=True, slots=True)
class M4Inputs:
    """Every table M4 reads (internal representation); optional tables are ``None`` when absent."""

    replicons: pd.DataFrame
    genes: pd.DataFrame
    xrefs: pd.DataFrame
    annotation_source: str
    tus: pd.DataFrame
    layers: pd.DataFrame
    regions: pd.DataFrame
    region_genes: pd.DataFrame
    peptides: pd.DataFrame | None
    abundance: pd.DataFrame | None
    peptide_hits: pd.DataFrame
    protein_hits: pd.DataFrame
    unmapped: pd.DataFrame
    id_map: pd.DataFrame | None


FILE_NAMES: Mapping[str, str] = {
    "replicons": "replicons.parquet",
    "genes": "genes.parquet",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
    "tus": "tus.parquet",
    "layers": "layers.parquet",
    "regions": "regions.parquet",
    "region_genes": "region_genes.parquet",
    "peptides": "peptides.norm.parquet",
    "abundance": "abundance.norm.parquet",
    "peptide_hits": "peptide_hits.tsv",
    "protein_hits": "protein_hits.tsv",
    "unmapped": "unmapped.tsv",
    "id_map": "id_map.tsv",
}
M4_INPUT_KEYS: tuple[str, ...] = tuple(FILE_NAMES)
OPTIONAL_KEYS = frozenset({"peptides", "abundance", "id_map"})
_ONE_BASED = frozenset({"peptide_hits", "protein_hits"})  # M3 TSVs with coordinates


def default_paths(work_dir: Path) -> dict[str, Path]:
    """Key → path in the fixture / local-executor layout (every file directly in ``work_dir``)."""
    return {key: work_dir / name for key, name in FILE_NAMES.items()}


def _issue(key: str, rule: str, message: str, severity: Severity = "error") -> Issue:
    return Issue(f"contract.{key}.{rule}", severity, message)


def _convert(values: list[str | None], col: Column) -> pd.Series:
    """Type one TSV column per ``col``; unparseable cells leave the raw strings (reported later)."""
    if col.dtype == "list":
        return pd.Series([[] if v is None else v.split(";") for v in values], dtype=object)
    try:
        if col.dtype == "int":
            ints = [None if v is None else int(v) for v in values]
            return pd.Series(ints, dtype="Int64" if None in ints else "int64")
        if col.dtype == "float":
            return pd.Series([np.nan if v is None else float(v) for v in values], dtype="float64")
        if col.dtype == "bool":
            table = {"true": True, "false": False}
            bools = [None if v is None else table[v] for v in values]
            return pd.Series(bools, dtype="boolean" if None in bools else bool)
    except (ValueError, KeyError):
        pass
    return pd.Series(values, dtype=object)


def _read_tsv_table(key: str, path: Path) -> pd.DataFrame:
    header, rows = read_tsv(path)
    cols = {c.name: c for c in TABLES[key]}
    data: dict[str, pd.Series] = {}
    for j, name in enumerate(header):
        values = [r[j] for r in rows]
        s = _convert(values, cols[name]) if name in cols else pd.Series(values, dtype=object)
        if key in _ONE_BASED and name == "start" and pd.api.types.is_integer_dtype(s.dtype):
            s = s - 1
        data[name] = s
    return pd.DataFrame(data, columns=header)


def _read_parquet_table(key: str, path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path, engine="pyarrow")
    for c in TABLES[key]:
        if c.dtype == "list" and c.name in df.columns:
            df[c.name] = pd.Series(
                [v.tolist() if isinstance(v, np.ndarray) else v for v in df[c.name].tolist()],
                index=df.index,
                dtype=object,
            )
    return df


def _empty(key: str) -> pd.DataFrame:
    return pd.DataFrame({c.name: pd.Series([], dtype=object) for c in TABLES[key]})


def _read_meta(path: Path | None, issues: list[Issue]) -> str:
    if path is None or not path.is_file():
        where = FILE_NAMES["annotation_meta"] if path is None else f"{path.name} ({path})"
        issues.append(_issue("annotation_meta", "missing_file", f"missing {where}"))
        return ""
    try:
        meta = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        issues.append(_issue("annotation_meta", "unreadable", f"{path.name}: {exc}"))
        return ""
    source = meta.get("annotation_source") if isinstance(meta, dict) else None
    if not isinstance(source, str) or not source:
        msg = f"{path.name} has no annotation_source string"
        issues.append(_issue("annotation_meta", "annotation_source", msg))
        return ""
    return source


def _load_table(key: str, path: Path | None, issues: list[Issue]) -> pd.DataFrame | None:
    optional = key in OPTIONAL_KEYS
    if path is None or not path.is_file():
        where = FILE_NAMES[key] if path is None else f"{path.name} ({path})"
        sev: Severity = "warn" if optional else "error"
        issues.append(_issue(key, "missing_file", f"missing {'optional ' * optional}{where}", sev))
        return None if optional else _empty(key)
    try:
        if path.name.endswith(".parquet"):
            df = _read_parquet_table(key, path)
        else:
            df = _read_tsv_table(key, path)
    except Exception as exc:  # any unreadable file is reported, never raised (spec)
        issues.append(_issue(key, "unreadable", f"{path.name}: {type(exc).__name__}: {exc}"))
        return None if optional else _empty(key)
    issues.extend(validate_table(key, df))
    return df


def load_m4_inputs(paths: Mapping[str, Path]) -> tuple[M4Inputs, list[Issue]]:
    """Read every contract file in ``paths`` (keys of ``M4_INPUT_KEYS``) and validate it.

    Missing required files are errors (an empty frame stands in), missing optional files warnings
    (``None``); per-table and cross-table issues follow. Never raises on content.
    """
    issues: list[Issue] = []
    frames: dict[str, pd.DataFrame | None] = {
        key: _load_table(key, paths.get(key), issues) for key in M4_INPUT_KEYS if key in TABLES
    }
    source = _read_meta(paths.get("annotation_meta"), issues)

    def req(key: str) -> pd.DataFrame:
        df = frames[key]
        return df if df is not None else _empty(key)

    inputs = M4Inputs(
        replicons=req("replicons"),
        genes=req("genes"),
        xrefs=req("xrefs"),
        annotation_source=source,
        tus=req("tus"),
        layers=req("layers"),
        regions=req("regions"),
        region_genes=req("region_genes"),
        peptides=frames["peptides"],
        abundance=frames["abundance"],
        peptide_hits=req("peptide_hits"),
        protein_hits=req("protein_hits"),
        unmapped=req("unmapped"),
        id_map=frames["id_map"],
    )
    issues.extend(validate_inputs(inputs))
    return inputs, issues
