"""Generic converter: validate, normalise and report (GIAP_SPEC.md §5.3, §5.4).

Several inputs are concatenated in argument order. ``peptide`` values with modification syntax are
stripped (``stripped``, info), booleans are canonicalised, unknown columns are dropped
(``unknown_column``, warn) and ``--source``/``--search-db`` fill missing cells.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.config.model import InputConfig
from giap.io.converters.base import ConvertResult, finish_peptides
from giap.io.converters.strip import strip_modifications
from giap.io.issues import Issue, row_issue
from giap.io.peptides import (
    _BAD_TYPE,
    PEPTIDE_COLUMNS,
    _empty_column,
    _rows,
    read_peptide_table,
    validate_peptides,
)

__all__ = ["convert_generic"]


def _read_all(paths: Sequence[Path]) -> tuple[list[pd.DataFrame], list[str], list[Issue]]:
    frames: list[pd.DataFrame] = []
    unknown: list[str] = []
    issues: list[Issue] = []
    for path in paths:
        try:
            df, found = read_peptide_table(path)
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            issues.append(Issue("unreadable", "error", f"{path}: {exc}"))
            continue
        for issue in found:
            if issue.code != "unknown_column":
                issues.append(dataclasses.replace(issue, message=f"{path.name}: {issue.message}"))
        unknown += [c for c in df.columns if c not in PEPTIDE_COLUMNS and c != _BAD_TYPE]
        frames.append(df)
    return frames, list(dict.fromkeys(unknown)), issues


def convert_generic(
    paths: Sequence[Path], *, source: str | None, search_db: str | None, cfg: InputConfig
) -> ConvertResult:
    """Convert generic peptide TSVs (§5.3) into one normalised generic table."""
    frames, unknown, issues = _read_all(paths)
    if unknown:
        issues.append(Issue("unknown_column", "warn", f"unknown columns dropped: {unknown}"))
    present = [c for c in PEPTIDE_COLUMNS if any(c in f.columns for f in frames)]
    for flag, col in ((source, "source"), (search_db, "search_db")):
        if flag is not None and col not in present:
            present.append(col)
    keep = [*present, _BAD_TYPE]
    aligned = []
    for f in frames:
        f = f.copy()
        for col in present:
            if col not in f.columns:
                f[col] = _empty_column(col, len(f))
        if _BAD_TYPE not in f.columns:
            f[_BAD_TYPE] = np.zeros(len(f), dtype=bool)
        aligned.append(f[keep])
    if not aligned:
        return ConvertResult(pd.DataFrame(columns=present), None, issues)
    df = pd.concat(aligned, ignore_index=True)
    for flag, col in ((source, "source"), (search_db, "search_db")):
        if flag is not None:
            df[col] = df[col].fillna(flag)

    if "peptide" in df.columns:
        pep = df["peptide"]
        changed = pep.notna() & ~pep.str.fullmatch(r"[A-Z]+").fillna(False).astype(bool)
        mask = changed.to_numpy(dtype=bool)
        if mask.any():
            df.loc[mask, "peptide"] = pep[mask].map(strip_modifications)
        stripped = row_issue("stripped", "info", "modifications stripped from peptide", _rows(mask))
        if stripped is not None:
            issues.append(stripped)

    _, found = validate_peptides(df, cfg)
    issues += found
    return ConvertResult(finish_peptides(df, present), None, issues)
