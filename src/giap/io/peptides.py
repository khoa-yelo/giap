"""Generic peptide table: read, validate, aggregate (GIAP_SPEC.md §5.3; ADR-0005 item 2, ADR-0006).

``read_peptide_table`` types the columns, ``validate_peptides`` (pure) keeps evidence-eligible rows
and derives ``is_decoy``/``is_contaminant``/``source``/``search_db``, ``aggregate_peptides`` groups
by stripped ``peptide`` (I/L kept as written) into one row per peptide plus per-run detail.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
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
# ``nan`` text is an error here (it would otherwise pass as a missing value, ADR-0005).
_NAN_IS_ERROR = ("q_value", "pep", "intensity")
_NAN_TEXT = frozenset({"nan", "+nan", "-nan"})
_INT_COLS = ("charge",)
_BOOL_COLS = ("is_decoy",)
_STR = pd.StringDtype()
_BOOL_TEXT = {
    **dict.fromkeys(("true", "TRUE", "True", "1"), True),
    **dict.fromkeys(("false", "FALSE", "False", "0"), False),
}
_ENUMS = {"source": ("database", "denovo"), "search_db": ("same", "other")}
# Aggregated source/search_db: the first value present in this order wins.
_PRIORITY = {"source": ("database", "denovo", "auto"), "search_db": ("same", "other", "auto")}
# Private marker column set by ``read_peptide_table``: the row had an unparseable value.
_BAD_TYPE = "_bad_type"


def _typed(raw: pd.Series, col: str) -> tuple[pd.Series, np.ndarray]:
    """Parse one known column; return the typed series and the mask of unparseable cells."""
    present = raw.notna().to_numpy()
    if col in _BOOL_COLS:
        mapped = raw.map(_BOOL_TEXT)
        bad = present & mapped.isna().to_numpy()
        return mapped.astype("boolean"), bad
    num = pd.to_numeric(raw, errors="coerce").astype("float64")
    parsed = num.notna().to_numpy()
    if col in _INT_COLS:
        values = num.to_numpy()
        integral = parsed & np.isfinite(values)
        integral[integral] = values[integral] == np.round(values[integral])
        bad = present & ~integral
        out = pd.array(np.where(integral, values, 0).astype("int64"), dtype="Int64")
        out[~integral] = pd.NA
        return pd.Series(out), bad
    bad = present & ~parsed
    if col not in _NAN_IS_ERROR:
        bad &= ~raw.str.lower().isin(_NAN_TEXT).to_numpy(dtype=bool, na_value=False)
    return num, bad


def read_peptide_table(path: Path) -> tuple[pd.DataFrame, list[Issue]]:
    """Read a (possibly gzipped) peptide TSV; ``NA`` and empty cells are missing.

    Known columns are typed (float64, Int64, boolean, string); unparseable values (and ``nan``
    text in q_value/pep/intensity) become missing with a ``bad_type:<col>`` error, and the row is
    marked so that ``validate_peptides`` drops it. Unknown columns are kept as strings
    (``unknown_column`` info).
    """
    header, rows = read_tsv(path)
    n = len(rows)
    columns = list(zip(*rows, strict=True)) if rows else [() for _ in header]
    issues: list[Issue] = []
    data: dict[str, pd.Series] = {}
    any_bad = np.zeros(n, dtype=bool)
    for col, values in zip(header, columns, strict=True):
        if col in _FLOAT_COLS or col in _INT_COLS or col in _BOOL_COLS:
            raw = pd.Series(np.asarray(values, dtype=object), dtype=object)
            data[col], bad = _typed(raw, col)
            kind = "float64" if col in _FLOAT_COLS else "Int64" if col in _INT_COLS else "boolean"
            issue = row_issue(
                f"bad_type:{col}", "error", f"{col}: value is not a valid {kind}", _rows(bad)
            )
            if issue is not None:
                issues.append(issue)
            any_bad |= bad
        else:
            data[col] = pd.Series(np.asarray(values, dtype=object), dtype=_STR)
    unknown = [c for c in header if c not in PEPTIDE_COLUMNS]
    if unknown:
        issues.append(Issue("unknown_column", "info", f"unknown columns kept: {unknown}"))
    df = pd.DataFrame(data, index=pd.RangeIndex(n))
    if any_bad.any():
        df[_BAD_TYPE] = any_bad
    return df, issues


def _empty_column(col: str, n: int) -> pd.Series:
    if col in _FLOAT_COLS:
        return pd.Series(np.full(n, np.nan), dtype="float64")
    if col in _INT_COLS:
        return pd.Series([None] * n, dtype="Int64")
    if col in _BOOL_COLS:
        return pd.Series([None] * n, dtype="boolean")
    return pd.Series([None] * n, dtype=_STR)


def _all_prefixed(ids: pd.Series, prefixes: Sequence[str]) -> np.ndarray:
    """True where every ``;``-separated entry of normalised ``ids`` starts with a prefix."""
    if not prefixes:
        return np.zeros(len(ids), dtype=bool)
    alt = "|".join(re.escape(p) for p in prefixes)
    pattern = rf"(?:{alt})[^;]*(?:;(?:{alt})[^;]*)*"
    return np.asarray(ids.str.fullmatch(pattern).to_numpy(dtype=bool, na_value=False), dtype=bool)


def _normalise_ids(ids: pd.Series) -> pd.Series:
    """Strip whitespace around ``;``-separated entries and drop empty entries (empty → NA)."""
    out = (
        ids.astype(_STR)
        .str.replace(r"\s*;\s*", ";", regex=True)
        .str.replace(r";{2,}", ";", regex=True)
        .str.strip()
        .str.strip(";")
    )
    return out.mask(out == "")


def _rows(mask: np.ndarray | pd.Series) -> np.ndarray:
    """1-based row numbers where ``mask`` is true."""
    return np.flatnonzero(np.asarray(mask, dtype=bool)) + 1


def _mask(values: pd.Series) -> np.ndarray:
    return np.asarray(values.to_numpy(dtype=bool, na_value=False), dtype=bool)


def validate_peptides(df_raw: pd.DataFrame, cfg: InputConfig) -> tuple[pd.DataFrame, list[Issue]]:
    """Return the evidence-eligible rows (fresh frame, ``RangeIndex``) and the issues found.

    Row numbers in issue examples are 1-based positions in ``df_raw``. Rows above
    ``cfg.q_value_max`` are kept (needed for FDR estimation). Rows with any error-level issue
    are dropped (the issue is still reported). Rows sharing (run, spectrum_id, peptide, charge)
    are one spectrum's evidence: only the first is kept (``duplicate_spectrum_row`` warning).
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
    error = np.zeros(n, dtype=bool)
    if _BAD_TYPE in df.columns:
        error |= _mask(df.pop(_BAD_TYPE))
    if issues:
        return _finish(df.iloc[0:0], cfg), issues

    def add_error(code: str, message: str, mask: np.ndarray) -> None:
        nonlocal error
        add(row_issue(code, "error", message, _rows(mask)))
        error |= mask

    for col, hi in (("q_value", "1]"), ("pep", "1]"), ("intensity", "inf)")):
        v = df[col].to_numpy(dtype="float64", na_value=np.nan)
        ok = np.isnan(v) | ((v >= 0) & (v <= 1 if hi == "1]" else np.isfinite(v)))
        add_error(f"out_of_range:{col}", f"{col} outside [0, {hi}", ~ok)
    add_error("out_of_range:charge", "charge < 1", _mask(df["charge"] < 1))
    for col, allowed in _ENUMS.items():
        v = df[col]
        add_error(
            f"bad_enum:{col}", f"{col} not in {list(allowed)}", _mask(v.notna() & ~v.isin(allowed))
        )

    run = df["run"].fillna("")
    has_sid = _mask(df["spectrum_id"].notna())
    n_pep = (
        pd.DataFrame({"run": run, "sid": df["spectrum_id"], "peptide": df["peptide"]})[has_sid]
        .groupby(["run", "sid"], sort=False)["peptide"]
        .transform("nunique")
    )
    conflict = np.zeros(n, dtype=bool)
    conflict[n_pep.index[n_pep.to_numpy() > 1]] = True
    add_error(
        "duplicate_spectrum_id",
        "(run, spectrum_id) repeated with different peptides (rows dropped)",
        conflict,
    )

    pep = df["peptide"]
    letters = _mask(pep.str.fullmatch(r"[A-Z]+"))
    length = pep.str.len()
    lo_len, hi_len = cfg.peptide_len[0], cfg.peptide_len[1]
    short = letters & _mask(length < lo_len)
    long = letters & _mask(length > hi_len)
    add(
        row_issue(
            "non_letters", "warn", "peptide is not uppercase A-Z only (dropped)", _rows(~letters)
        )
    )
    add(
        row_issue("short_peptide", "warn", f"peptide shorter than {lo_len} (dropped)", _rows(short))
    )
    add(row_issue("long_peptide", "warn", f"peptide longer than {hi_len} (dropped)", _rows(long)))
    keep = letters & ~short & ~long & ~error

    decoy_no_ids = keep & _mask(df["is_decoy"]) & _mask(df["protein_ids"].isna())
    add(
        row_issue(
            "decoy_without_protein_ids",
            "warn",
            "is_decoy is true but protein_ids is missing",
            _rows(decoy_no_ids),
        )
    )

    key = pd.DataFrame(
        {"run": run, "sid": df["spectrum_id"], "peptide": pep, "charge": df["charge"]}
    )
    cand = keep & has_sid
    dup = np.zeros(n, dtype=bool)
    dup[cand] = key[cand].duplicated(keep="first").to_numpy()
    add(
        row_issue(
            "duplicate_spectrum_row",
            "warn",
            "(run, spectrum_id, peptide, charge) repeated; first row kept",
            _rows(dup),
        )
    )
    keep &= ~dup

    src = df["source"].fillna(cfg.peptide_source)
    named = keep & _mask(src.isin(("database", "denovo")))
    n_src = (
        pd.DataFrame({"peptide": pep, "source": src})[named]
        .groupby("peptide", sort=False)["source"]
        .transform("nunique")
    )
    mixed = np.zeros(n, dtype=bool)
    mixed[n_src.index[n_src.to_numpy() > 1]] = True
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
    df["protein_ids"] = _normalise_ids(df["protein_ids"])
    ids = df["protein_ids"]
    by_prefix = pd.Series(_all_prefixed(ids, cfg.decoy_prefixes), index=df.index)
    df["is_decoy"] = df["is_decoy"].astype("boolean").fillna(by_prefix).astype(bool)
    df["is_contaminant"] = _all_prefixed(ids, cfg.contaminant_prefixes)
    df["source"] = df["source"].fillna(cfg.peptide_source).astype(_STR)
    df["search_db"] = df["search_db"].fillna(cfg.search_db).astype(_STR)
    extra = [
        c for c in df.columns if c not in PEPTIDE_COLUMNS and c not in ("is_contaminant", _BAD_TYPE)
    ]
    return df[[*PEPTIDE_COLUMNS, "is_contaminant", *extra]]


def _sorted_union(
    peptide: pd.Series, values: pd.Series, sep: str | None, peptides: pd.Index
) -> pd.Series:
    """Per peptide (in ``peptides`` order), the sorted union of the non-missing entries."""
    pairs = pd.DataFrame({"peptide": peptide.to_numpy(), "v": values.astype(_STR).to_numpy()})
    pairs = pairs[pairs["v"].notna()]
    if sep is not None:
        pairs = pairs.assign(v=pairs["v"].str.split(sep)).explode("v")
        pairs = pairs[pairs["v"].notna() & (pairs["v"] != "")]
    pairs = pairs.drop_duplicates().sort_values(["peptide", "v"])
    keys = pairs["peptide"].to_numpy(dtype=object)
    vals = pairs["v"].to_numpy(dtype=object)
    if not len(keys):
        return pd.Series([[] for _ in peptides], index=peptides, dtype=object)
    starts = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1]])
    chunks = np.split(vals, starts[1:])
    found = dict(zip(keys[starts], (c.tolist() for c in chunks), strict=True))
    return pd.Series([found.get(p, []) for p in peptides], index=peptides, dtype=object)


def _first_by_priority(
    peptide: pd.Series, values: pd.Series, order: Sequence[str], peptides: pd.Index
) -> pd.Series:
    """Per peptide, the present value earliest in ``order`` (unknown values after, sorted)."""
    present = values.dropna().unique()
    categories = [*order, *sorted(set(present) - set(order))]
    codes = pd.Categorical(values, categories=categories).codes.astype("int64")
    codes[codes < 0] = len(categories)
    best = pd.Series(codes).groupby(peptide.to_numpy(), sort=True).min().reindex(peptides)
    lookup = np.array([*categories, ""], dtype=object)
    return pd.Series(lookup[best.to_numpy()], index=peptides).astype(_STR)


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
    n_rows = g.size().astype("int64")
    peptides = n_rows.index
    pep = work["peptide"]
    agg = pd.DataFrame(
        {
            "q_value": g["q_value"].min().astype("float64"),
            "pep": g["pep"].min().astype("float64"),
            "score": g["score"].max().astype("float64"),
            "is_decoy": g["is_decoy"].all().astype(bool),
            "is_contaminant": g["is_contaminant"].all().astype(bool),
            "protein_ids": _sorted_union(pep, work["protein_ids"], ";", peptides),
            "modified_peptides": _sorted_union(pep, work["modified_peptide"], None, peptides),
            "n_rows": n_rows,
            "n_runs": g["run"].nunique().astype("int64"),
            "intensity": runs.groupby("peptide", sort=True, observed=True)["intensity"]
            .median()
            .astype("float64"),
            "source": _first_by_priority(pep, work["source"], _PRIORITY["source"], peptides),
            "search_db": _first_by_priority(
                pep, work["search_db"], _PRIORITY["search_db"], peptides
            ),
        }
    )
    agg.index.name = "peptide"
    agg = agg.reset_index()
    agg["peptide"] = agg["peptide"].astype(_STR)
    return agg[list(AGG_COLUMNS)].reset_index(drop=True), runs
