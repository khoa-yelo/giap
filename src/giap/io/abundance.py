"""Abundance table: read, detect the scale, count missing values (GIAP_SPEC.md §5.3, §7.0).

The first line is the header; column 1 holds protein IDs (UniProt ``sp|ACC|NAME`` resolves to
``ACC`` as in ``proteins.faa``, ADR-0005 item 2), column 2 the values, later columns are ignored
with a warning. ``NA`` and empty cells are missing; ``inf``/``nan`` text is an error.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from giap.io.issues import Issue, row_issue
from giap.io.protein_fasta import resolve_protein_id
from giap.io.tsv import read_tsv

__all__ = [
    "AbundanceStats",
    "AbundanceTable",
    "Scale",
    "abundance_stats",
    "detect_scale",
    "read_abundance",
    "resolve_scale",
    "validate_abundance",
]

Scale = Literal["log2", "linear"]

# Values at or below this (with no negatives) are taken as log2 (§5.3).
LOG2_MAX = 50.0
_NAN_TEXT = frozenset({"nan", "+nan", "-nan"})


@dataclass(frozen=True, slots=True)
class AbundanceTable:
    """Parsed abundance table; ``values[i]`` is ``None`` when row ``i`` is ``NA``/empty."""

    id_column: str
    value_column: str
    ids: tuple[str, ...]
    values: tuple[float | None, ...]
    ignored_columns: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class AbundanceStats:
    n_rows: int
    n_na: int
    n_zero: int
    n_negative: int
    missing_frac: float


def _rows(mask: np.ndarray) -> list[int]:
    """1-based data-row numbers where ``mask`` is true."""
    return [int(i) + 1 for i in np.flatnonzero(mask)]


def read_abundance(path: Path) -> tuple[AbundanceTable | None, list[Issue]]:
    """Read an abundance TSV; the table is ``None`` when its shape or values are unusable."""
    try:
        header, rows = read_tsv(path)
    except ValueError as exc:
        return None, [Issue("malformed_table", "error", str(exc))]
    if len(header) < 2:
        msg = "abundance table needs an ID column and a value column"
        return None, [Issue("no_value_column", "error", msg)]
    if not rows:
        return None, [Issue("empty_table", "error", "abundance table has a header but no rows")]

    issues: list[Issue] = []
    ignored = tuple(header[2:])
    if ignored:
        msg = f"only one value column is used; ignored columns: {', '.join(ignored)}"
        issues.append(Issue("extra_columns", "warn", msg))

    raw_ids = pd.Series([r[0] for r in rows], dtype=object)
    ids = raw_ids.map(lambda s: "" if s is None else resolve_protein_id(s)[0])
    empty = (ids == "").to_numpy()
    dup = (ids.duplicated(keep="first") & ~ids.eq("")).to_numpy()
    for code, mask, msg in (
        ("empty_id", empty, "empty protein ID"),
        ("duplicate_ids", dup, "protein ID repeats an earlier row"),
    ):
        issue = row_issue(code, "error", msg, _rows(mask))
        if issue is not None:
            issues.append(issue)

    text = pd.Series([r[1] for r in rows], dtype=object)
    missing = text.isna().to_numpy()
    parsed = pd.to_numeric(text, errors="coerce").to_numpy(dtype=np.float64)
    nan_text = text.str.strip().str.lower().isin(_NAN_TEXT).to_numpy()
    unparsed = ~missing & np.isnan(parsed) & ~nan_text
    non_finite = ~missing & ~unparsed & ~np.isfinite(parsed)
    bad_values = False
    for code, mask, msg in (
        ("non_numeric_values", unparsed, f"non-numeric value in {header[1]!r}"),
        ("non_finite", non_finite, f"inf/nan value in {header[1]!r} (write NA for missing)"),
    ):
        issue = row_issue(code, "error", msg, _rows(mask))
        if issue is not None:
            issues.append(issue)
            bad_values = True
    if bad_values:
        return None, issues

    # ``float`` on the original text keeps exact round-trip parsing.
    values = tuple(None if s is None else float(s) for s in text)
    table = AbundanceTable(header[0], header[1], tuple(ids), values, ignored)
    return table, issues


def _present(values: Sequence[float | None]) -> np.ndarray:
    return np.array([v for v in values if v is not None], dtype=np.float64)


def detect_scale(values: Sequence[float | None]) -> Scale:
    """``log2`` if any value is negative or all are ≤ 50; NA ignored; all-NA → ``linear``."""
    arr = _present(values)
    if arr.size and ((arr < 0).any() or (arr <= LOG2_MAX).all()):
        return "log2"
    return "linear"


def resolve_scale(
    values: Sequence[float | None], requested: Literal["auto", "log2", "linear"]
) -> Scale:
    """The requested scale, or the detected one for ``auto``."""
    return detect_scale(values) if requested == "auto" else requested


def abundance_stats(table: AbundanceTable) -> AbundanceStats:
    n_rows = len(table.values)
    arr = _present(table.values)
    n_na = n_rows - arr.size
    return AbundanceStats(
        n_rows=n_rows,
        n_na=n_na,
        n_zero=int((arr == 0).sum()),
        n_negative=int((arr < 0).sum()),
        missing_frac=n_na / n_rows if n_rows else 0.0,
    )


def validate_abundance(
    table: AbundanceTable, scale: Scale, missing_warn_frac: float
) -> list[Issue]:
    """Scale-dependent checks: undetermined scale, negatives on linear, too many missing."""
    issues: list[Issue] = []
    stats = abundance_stats(table)
    if stats.n_na == stats.n_rows:
        msg = "every abundance value is missing; scale cannot be detected (linear assumed)"
        issues.append(Issue("scale_undetermined", "warn", msg, stats.n_rows))
    if scale == "linear":
        vals = np.array([np.nan if v is None else v for v in table.values], dtype=np.float64)
        issue = row_issue(
            "negative_on_linear", "error", "negative abundance on a linear scale", _rows(vals < 0)
        )
        if issue is not None:
            issues.append(issue)
    if stats.missing_frac > missing_warn_frac:
        msg = (
            f"{stats.missing_frac:.1%} of abundance values are missing (> {missing_warn_frac:.0%})"
        )
        issues.append(Issue("too_many_missing", "warn", msg, stats.n_na))
    return issues
