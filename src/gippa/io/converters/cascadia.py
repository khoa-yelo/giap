"""Cascadia 0.0.7 converter (BiblioSpec-style ``.ssl``; GIPPA_SPEC.md §5.4, M10-09).

``peptide`` ← stripped ``sequence``; ``modified_peptide`` ← ``sequence``; ``score`` ← ``score``;
``run`` ← the stem of ``file``; ``charge`` ← ``charge``; ``spectrum_id`` ← ``cascadia:<scan>``
(``scan`` is a prediction index, not a scan number); ``source`` = denovo. No q-value or PEP is
derived from ``score`` (it is not calibrated). The RT columns are never read. Several inputs are
concatenated in argument order; issue rows are 1-based over that concatenation.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from gippa.config.model import InputConfig
from gippa.io.converters.base import ConvertResult, finish_peptides
from gippa.io.converters.strip import strip_modifications
from gippa.io.issues import Issue, row_issue
from gippa.io.peptides import _rows, _typed, validate_peptides
from gippa.io.tsv import read_tsv

__all__ = ["CASCADIA_DEFAULT_CUTOFF", "convert_cascadia"]

CASCADIA_DEFAULT_CUTOFF = 0.8
"""Cascadia 0.0.7's ``-t/--score_threshold`` default; it writes a row only if ``score > 0.8``
(``cascadia/utils.py::write_results``). Counted (info), never applied by GIPPA."""

_STR = pd.StringDtype()
_REQUIRED = ("file", "scan", "charge", "sequence", "score-type", "score")
_KNOWN = frozenset((*_REQUIRED, "retention-time", "start-time", "end-time"))
_SCAN = re.compile(r"[0-9]+")
_PLAIN = re.compile(r"[A-Z]+")
_COLUMNS = ("peptide", "modified_peptide", "score", "run", "charge", "spectrum_id", "source")


def _stem(file: str | None) -> str | None:
    """Text after the last ``/``, minus one ``.gz`` then one ``.mzML`` (case-insensitive)."""
    name = (file or "").rsplit("/", 1)[-1]
    for suffix in (".gz", ".mzml"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
    return name or None


def _score(cell: str | None) -> float | None:
    try:
        value = float(cell) if cell is not None else math.nan
    except ValueError:
        return None
    return None if math.isnan(value) else value


def _fail(issues: list[Issue]) -> ConvertResult:
    return ConvertResult(pd.DataFrame(), None, issues)


def _read(
    paths: Sequence[Path],
) -> tuple[dict[str, list[str | None]], list[int], list[Issue]]:
    """The required columns of all inputs, concatenated; the input index of each row; issues."""
    cols: dict[str, list[str | None]] = {c: [] for c in _REQUIRED}
    origin: list[int] = []
    issues: list[Issue] = []
    for k, path in enumerate(paths):
        try:
            header, rows = read_tsv(path)
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            msg = str(exc) if str(path) in str(exc) else f"{path}: {exc}"
            issues.append(Issue("unreadable", "error", msg))
            continue
        missing = [c for c in _REQUIRED if c not in header]
        for c in missing:
            issues.append(Issue(f"missing_column:{c}", "error", f"{path}: no {c!r} column"))
        unknown = [c for c in header if c not in _KNOWN]
        if unknown:
            issues.append(Issue("unknown_column", "warn", f"{path}: columns ignored: {unknown}"))
        if missing:
            continue
        if not rows:
            issues.append(Issue("empty_input", "info", f"{path}: header only, no data rows"))
        for c in _REQUIRED:
            j = header.index(c)
            cols[c] += [r[j] for r in rows]
        origin += [k] * len(rows)
    return cols, origin, issues


def convert_cascadia(
    paths: Sequence[Path], *, search_db: str | None, cfg: InputConfig
) -> ConvertResult:
    """Convert Cascadia 0.0.7 ``.ssl`` files into the generic peptide table (de novo rows)."""
    if search_db is not None:
        msg = f"--search-db {search_db} is invalid with --from cascadia (de novo rows)"
        return _fail([Issue("search_db_invalid", "error", msg)])
    cols, origin, issues = _read(paths)
    if any(i.severity == "error" for i in issues):
        return _fail(issues)

    def check(code: str, message: str, bad: Sequence[bool] | np.ndarray) -> None:
        issue = row_issue(code, "error", message, _rows(np.asarray(bad, dtype=bool)))
        if issue is not None:
            issues.append(issue)

    runs = [_stem(f) for f in cols["file"]]
    check("bad_value:file", "file: missing, or its stem is empty", [r is None for r in runs])
    sources: dict[str, list[int]] = {}
    for run, k in zip(runs, origin, strict=True):
        if run is not None and k not in sources.setdefault(run, []):
            sources[run].append(k)
    for run, ks in sources.items():
        if len(ks) > 1:
            named = ", ".join(str(paths[k]) for k in ks)
            issues.append(
                Issue("run_in_several_inputs", "error", f"run {run!r} comes from: {named}")
            )
    scans = cols["scan"]
    check(
        "bad_type:scan",
        "scan: not a non-negative integer",
        [s is None or _SCAN.fullmatch(s) is None for s in scans],
    )
    types = cols["score-type"]
    wrong = [t != "UNKNOWN" for t in types]
    if any(wrong):
        first = types[wrong.index(True)]
        check("unsupported_score_type", f"score-type {first!r} is not 'UNKNOWN'", wrong)
    parsed = [_score(s) for s in cols["score"]]
    check("bad_type:score", "score: missing, not a number or NaN", [s is None for s in parsed])
    score = np.array([math.nan if s is None else s for s in parsed], dtype="float64")
    check("out_of_range:score", "score outside [0, 1]", (score < 0) | (score > 1))
    charge, bad_charge = _typed(pd.Series(cols["charge"], dtype=object), "charge")
    check("bad_type:charge", "charge: value is not a valid Int64", bad_charge)
    if any(i.severity == "error" for i in issues):
        return _fail(issues)

    seqs = cols["sequence"]
    strip = [s is not None and _PLAIN.fullmatch(s) is None for s in seqs]
    peptide = [
        (strip_modifications(s) or None) if s is not None and t else s
        for s, t in zip(seqs, strip, strict=True)
    ]
    n = len(seqs)
    df = pd.DataFrame(
        {
            "peptide": pd.Series(peptide, dtype=_STR),
            "modified_peptide": pd.Series(seqs, dtype=_STR),
            "score": score,
            "run": pd.Series(runs, dtype=_STR),
            "charge": charge,
            "spectrum_id": pd.Series([f"cascadia:{s}" for s in scans], dtype=_STR),
            "source": pd.Series(["denovo"] * n, dtype=_STR),
        },
        index=pd.RangeIndex(n),
    )
    for found in (
        row_issue("stripped", "info", "modifications stripped from peptide", _rows(strip)),
        row_issue(
            "score_at_or_below_default_cutoff",
            "info",
            f"score <= {CASCADIA_DEFAULT_CUTOFF} (Cascadia's default cut-off; kept)",
            _rows(score <= CASCADIA_DEFAULT_CUTOFF),
        ),
    ):
        if found is not None:
            issues.append(found)
    _, found_v = validate_peptides(df, cfg)
    return ConvertResult(finish_peptides(df, _COLUMNS), None, issues + found_v)
