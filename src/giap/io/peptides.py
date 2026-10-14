"""Generic peptide table: read, validate, aggregate (GIAP_SPEC.md §5.3; ADR-0005 item 2, ADR-0006).

``read_peptide_table`` types the columns, ``validate_peptides`` (pure) keeps evidence-eligible rows
and derives ``is_decoy``/``is_contaminant``/``source``/``search_db``, ``aggregate_peptides`` groups
by stripped ``peptide`` (I/L kept as written) into one row per peptide plus per-run detail.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.config.model import InputConfig
from giap.io.issues import Issue, Severity, row_issue
from giap.io.tsv import read_tsv

__all__ = [
    "AGG_COLUMNS",
    "PEPTIDE_COLUMNS",
    "PEPTIDE_RUN_COLUMNS",
    "Issue",
    "Severity",
    "aggregate_peptides",
    "read_peptide_table",
    "validate_peptides",
]

PEPTIDE_COLUMNS: tuple[str, ...] = (
    "peptide",
    "modified_peptide",
    "protein_ids",
    "q_value",
    "pep",
    "score",
    "is_decoy",
    "intensity",
    "run",
    "charge",
    "spectrum_id",
    "source",
    "search_db",
)
AGG_COLUMNS: tuple[str, ...] = (
    "peptide",
    "q_value",
    "pep",
    "score",
    "is_decoy",
    "is_contaminant",
    "protein_ids",
    "modified_peptides",
    "n_rows",
    "n_runs",
    "intensity",
    "source",
    "search_db",
)
PEPTIDE_RUN_COLUMNS: tuple[str, ...] = ("peptide", "run", "intensity", "n_rows")

_FLOAT_COLS = ("q_value", "pep", "score", "intensity")
_INT_COLS = ("charge",)
_BOOL_COLS = ("is_decoy",)
_STR = pd.StringDtype()
_TRUE = frozenset({"true", "TRUE", "True", "1"})
_FALSE = frozenset({"false", "FALSE", "False", "0"})
_ENUMS = {"source": ("database", "denovo"), "search_db": ("same", "other")}
# Aggregated source/search_db: the first value present in this order wins.
_PRIORITY = {"source": ("database", "denovo", "auto"), "search_db": ("same", "other", "auto")}


def _parse_float(s: str) -> float:
    return float(s)


def _parse_int(s: str) -> int:
    f = float(s)
    if not f.is_integer():
        raise ValueError(s)
    return int(f)


def _parse_bool(s: str) -> bool:
    if s in _TRUE:
        return True
    if s in _FALSE:
        return False
    raise ValueError(s)


def _typed(
    values: list[str | None], parse: Callable[[str], object], dtype: str, col: str
) -> tuple[pd.Series, Issue | None]:
    out: list[object] = []
    bad: list[int] = []
    for i, v in enumerate(values, start=1):
        if v is None:
            out.append(None)
            continue
        try:
            out.append(parse(v))
        except ValueError:
            out.append(None)
            bad.append(i)
    issue = row_issue(f"bad_type:{col}", "error", f"{col}: value is not a valid {dtype}", bad)
    if dtype == "float64":
        return pd.Series([np.nan if x is None else x for x in out], dtype="float64"), issue
    return pd.Series(out, dtype=dtype), issue


def read_peptide_table(path: Path) -> tuple[pd.DataFrame, list[Issue]]:
    """Read a (possibly gzipped) peptide TSV; ``NA`` and empty cells are missing.

    Known columns are typed (float64, Int64, boolean, string); unparseable values become missing
    with a ``bad_type:<col>`` error. Unknown columns are kept as strings (``unknown_column`` info).
    """
    header, rows = read_tsv(path)
    issues: list[Issue] = []
    data: dict[str, pd.Series] = {}
    for j, col in enumerate(header):
        values = [r[j] for r in rows]
        issue: Issue | None = None
        if col in _FLOAT_COLS:
            data[col], issue = _typed(values, _parse_float, "float64", col)
        elif col in _INT_COLS:
            data[col], issue = _typed(values, _parse_int, "Int64", col)
        elif col in _BOOL_COLS:
            data[col], issue = _typed(values, _parse_bool, "boolean", col)
        else:
            data[col] = pd.Series(values, dtype=_STR)
        if issue is not None:
            issues.append(issue)
    unknown = [c for c in header if c not in PEPTIDE_COLUMNS]
    if unknown:
        issues.append(Issue("unknown_column", "info", f"unknown columns kept: {unknown}"))
    return pd.DataFrame(data, index=pd.RangeIndex(len(rows))), issues


def _empty_column(col: str, n: int) -> pd.Series:
    if col in _FLOAT_COLS:
        return pd.Series(np.full(n, np.nan), dtype="float64")
    if col in _INT_COLS:
        return pd.Series([None] * n, dtype="Int64")
    if col in _BOOL_COLS:
        return pd.Series([None] * n, dtype="boolean")
    return pd.Series([None] * n, dtype=_STR)


def _all_prefixed(ids: object, prefixes: Sequence[str]) -> bool:
    """True when every ``;``-separated entry of ``ids`` starts with one of ``prefixes``."""
    if not isinstance(ids, str):
        return False
    entries = [e for e in ids.split(";") if e]
    return bool(entries) and all(e.startswith(tuple(prefixes)) for e in entries)


def _rows(mask: Iterable[bool]) -> list[int]:
    return [i for i, m in enumerate(mask, start=1) if m]


def validate_peptides(df_raw: pd.DataFrame, cfg: InputConfig) -> tuple[pd.DataFrame, list[Issue]]:
    """Return the evidence-eligible rows (fresh frame, ``RangeIndex``) and the issues found.

    Row numbers in issue examples are 1-based positions in ``df_raw``. Rows above
    ``cfg.q_value_max`` are kept (needed for FDR estimation).
    """
    issues: list[Issue] = []

    def add(issue: Issue | None) -> None:
        if issue is not None:
            issues.append(issue)

    df = df_raw.reset_index(drop=True).copy()
    n = len(df)
    if "peptide" not in df.columns:
        issues.append(Issue("missing_peptide_column", "error", "no 'peptide' column"))
    elif n == 0:
        issues.append(Issue("empty_table", "error", "peptide table has no data rows"))
    for col in PEPTIDE_COLUMNS:
        if col not in df.columns:
            df[col] = _empty_column(col, n)
    if issues:
        return _finish(df.iloc[0:0], cfg), issues

    for col, lo, hi in (("q_value", 0, 1), ("pep", 0, 1), ("intensity", 0, math.inf)):
        v = df[col]
        add(
            row_issue(
                f"out_of_range:{col}",
                "error",
                f"{col} outside [{lo}, {hi}]",
                _rows(v.notna() & ~v.between(lo, hi)),
            )
        )
    charge = df["charge"]
    add(
        row_issue(
            "out_of_range:charge",
            "error",
            "charge < 1",
            _rows((charge.notna() & (charge < 1)).fillna(False)),
        )
    )
    for col, allowed in _ENUMS.items():
        v = df[col]
        add(
            row_issue(
                f"bad_enum:{col}",
                "error",
                f"{col} not in {list(allowed)}",
                _rows(v.notna() & ~v.isin(allowed)),
            )
        )

    sid = df[["peptide", "run", "spectrum_id"]].copy()
    sid["run"] = sid["run"].fillna("")
    sid = sid[sid["spectrum_id"].notna()]
    n_pep = sid.groupby(["run", "spectrum_id"])["peptide"].transform("nunique")
    dup = pd.Series(False, index=df.index)
    dup.loc[n_pep.index[n_pep > 1]] = True
    add(
        row_issue(
            "duplicate_spectrum_id",
            "error",
            "(run, spectrum_id) repeated with different peptides",
            _rows(dup),
        )
    )

    pep = df["peptide"]
    letters = pep.str.fullmatch(r"[A-Z]+").fillna(False).astype(bool)
    length = pep.str.len()
    lo_len, hi_len = cfg.peptide_len[0], cfg.peptide_len[1]
    short = letters & (length < lo_len).fillna(False).astype(bool)
    long = letters & (length > hi_len).fillna(False).astype(bool)
    add(
        row_issue(
            "non_letters", "warn", "peptide is not uppercase A-Z only (dropped)", _rows(~letters)
        )
    )
    add(
        row_issue("short_peptide", "warn", f"peptide shorter than {lo_len} (dropped)", _rows(short))
    )
    add(row_issue("long_peptide", "warn", f"peptide longer than {hi_len} (dropped)", _rows(long)))
    keep = letters & ~short & ~long

    decoy_no_ids = keep & df["is_decoy"].fillna(False).astype(bool) & df["protein_ids"].isna()
    add(
        row_issue(
            "decoy_without_protein_ids",
            "warn",
            "is_decoy is true but protein_ids is missing",
            _rows(decoy_no_ids),
        )
    )

    sources = df.loc[keep, ["peptide", "source"]].copy()
    sources["source"] = sources["source"].fillna(cfg.peptide_source)
    n_src = sources.groupby("peptide")["source"].transform(
        lambda s: len(set(s) & {"database", "denovo"})
    )
    mixed = pd.Series(False, index=df.index)
    mixed.loc[n_src.index[n_src > 1]] = True
    add(
        row_issue(
            "mixed_source",
            "warn",
            "peptide has both database and denovo rows (database wins)",
            _rows(mixed),
        )
    )

    return _finish(df[keep], cfg), issues


def _finish(df: pd.DataFrame, cfg: InputConfig) -> pd.DataFrame:
    """Derive ``is_decoy``/``is_contaminant`` and fill ``source``/``search_db`` from config."""
    df = df.reset_index(drop=True).copy()
    ids = df["protein_ids"].astype(object)
    by_prefix = ids.map(lambda x: _all_prefixed(x, cfg.decoy_prefixes)).astype(bool)
    df["is_decoy"] = df["is_decoy"].fillna(by_prefix).astype(bool)
    df["is_contaminant"] = ids.map(lambda x: _all_prefixed(x, cfg.contaminant_prefixes)).astype(
        bool
    )
    df["source"] = df["source"].fillna(cfg.peptide_source).astype(_STR)
    df["search_db"] = df["search_db"].fillna(cfg.search_db).astype(_STR)
    extra = [c for c in df.columns if c not in PEPTIDE_COLUMNS and c != "is_contaminant"]
    return df[[*PEPTIDE_COLUMNS, "is_contaminant", *extra]]


def _union(values: pd.Series, sep: str | None) -> list[str]:
    out: set[str] = set()
    for v in values:
        if isinstance(v, str):
            out.update(e for e in (v.split(sep) if sep else [v]) if e)
    return sorted(out)


def _first_by_priority(values: pd.Series, order: Sequence[str]) -> str:
    present = set(values.dropna())
    for v in order:
        if v in present:
            return v
    return min(present) if present else ""


def aggregate_peptides(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate validated rows to ``(peptides, peptide_runs)``, both sorted by peptide (then run).

    Per peptide: ``q_value``/``pep`` min, ``score`` max, ``is_decoy``/``is_contaminant`` all rows,
    ``protein_ids``/``modified_peptides`` sorted unions, ``n_rows``, ``n_runs``, and ``intensity``
    = median over runs of the per-run sums (runs with no intensity excluded; none → NA). A missing
    ``run`` is the run ``""``. ``peptide_runs`` holds the per-run sums and row counts.
    """
    work = df.copy()
    work["run"] = work["run"].fillna("").astype(_STR)
    work["peptide"] = work["peptide"].astype(_STR)

    by_run = work.groupby(["peptide", "run"], sort=True, observed=True)
    runs = pd.DataFrame(
        {
            "intensity": by_run["intensity"].sum(min_count=1).astype("float64"),
            "n_rows": by_run.size().astype("int64"),
        }
    ).reset_index()
    runs = runs[list(PEPTIDE_RUN_COLUMNS)].reset_index(drop=True)

    g = work.groupby("peptide", sort=True, observed=True)
    agg = pd.DataFrame(
        {
            "q_value": g["q_value"].min().astype("float64"),
            "pep": g["pep"].min().astype("float64"),
            "score": g["score"].max().astype("float64"),
            "is_decoy": g["is_decoy"].all().astype(bool),
            "is_contaminant": g["is_contaminant"].all().astype(bool),
            "protein_ids": g["protein_ids"].agg(lambda s: _union(s, ";")).astype(object),
            "modified_peptides": g["modified_peptide"]
            .agg(lambda s: _union(s, None))
            .astype(object),
            "n_rows": g.size().astype("int64"),
            "n_runs": g["run"].nunique().astype("int64"),
            "intensity": runs.groupby("peptide", sort=True, observed=True)["intensity"]
            .median()
            .astype("float64"),
            "source": g["source"]
            .agg(lambda s: _first_by_priority(s, _PRIORITY["source"]))
            .astype(_STR),
            "search_db": g["search_db"]
            .agg(lambda s: _first_by_priority(s, _PRIORITY["search_db"]))
            .astype(_STR),
        }
    )
    agg.index.name = "peptide"
    agg = agg.reset_index()
    agg["peptide"] = agg["peptide"].astype(_STR)
    return agg[list(AGG_COLUMNS)].reset_index(drop=True), runs
