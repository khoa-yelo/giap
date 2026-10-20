"""Reader for ``--replicon-depth`` (GIAP_SPEC.md §5.3: TSV ``contig``, ``mean_depth``; M4-18)."""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from giap.io.issues import Issue, row_issue
from giap.io.tsv import read_tsv

__all__ = ["read_replicon_depth"]

_COLUMNS = ("contig", "mean_depth")


def _depth(v: str | None) -> float:
    try:
        return float(v) if v is not None else math.nan
    except ValueError:
        return math.nan


def read_replicon_depth(path: Path) -> tuple[pd.DataFrame | None, list[Issue]]:
    """(``contig``, ``mean_depth``) frame and issues. Missing columns, non-numeric or non-finite
    depths, depths <= 0 (log2 undefined) and duplicate contigs are errors; the frame is then
    ``None``."""
    header, rows = read_tsv(path)
    missing = [c for c in _COLUMNS if c not in header]
    if missing:
        return None, [Issue("replicon_depth.columns", "error",
                            f"{path}: missing column(s) {', '.join(missing)}")]  # fmt: skip
    ic, idep = header.index("contig"), header.index("mean_depth")
    contigs = [r[ic] for r in rows]
    raw = [r[idep] for r in rows]
    depth = [_depth(v) for v in raw]
    seen: set[str | None] = set()
    dup: list[int] = []
    for i, c in enumerate(contigs, 1):
        if c in seen:
            dup.append(i)
        seen.add(c)
    checks = [
        ("replicon_depth.contig_missing", "contig is missing",
         [i for i, c in enumerate(contigs, 1) if c is None]),
        ("replicon_depth.not_numeric", "mean_depth is not a finite number",
         [i for i, d in enumerate(depth, 1) if not math.isfinite(d)]),
        ("replicon_depth.nonpositive", "mean_depth <= 0 (log2 undefined)",
         [i for i, d in enumerate(depth, 1) if math.isfinite(d) and d <= 0]),
        ("replicon_depth.duplicate_contig", "contig listed twice", dup),
    ]  # fmt: skip
    issues = [x for code, msg, idx in checks if (x := row_issue(code, "error", msg, idx))]
    if issues:
        return None, issues
    return pd.DataFrame({"contig": [str(c) for c in contigs], "mean_depth": depth}), []
