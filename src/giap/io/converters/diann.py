"""DIA-NN 2.x converter (``report.parquet`` or ``report.tsv``; GIAP_SPEC.md §5.4).

``peptide`` ← ``Stripped.Sequence``; ``modified_peptide`` ← ``Modified.Sequence``; ``protein_ids``
← ``Protein.Ids``; ``q_value`` ← ``Global.Q.Value``, else ``Q.Value``; ``pep`` ← ``PEP``;
``is_decoy`` ← ``Decoy`` (0/1); ``intensity`` ← ``Precursor.Normalised``, else
``Precursor.Quantity`` (0 → missing); ``run`` ← ``Run``; ``charge`` ← ``Precursor.Charge``;
``source`` = database. ``peptide`` values with modification syntax are stripped as in the generic
converter. Abundance: per ``Protein.Group``, the median over runs of ``PG.MaxLFQ`` (missing and 0
excluded, decoy rows excluded, groups with no value dropped), column ``maxlfq``.

Only the mapped columns are loaded; parquet numbers are converted directly.
"""

from __future__ import annotations

import csv
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from giap.config.model import InputConfig
from giap.io._open import open_text
from giap.io.converters.base import ConvertResult, finish_peptides
from giap.io.converters.strip import strip_modifications
from giap.io.issues import Issue, row_issue
from giap.io.peptides import _rows, _typed, validate_peptides
from giap.io.tsv import NA

__all__ = ["convert_diann"]

_STR = pd.StringDtype()
# generic column ← DIA-NN column(s), first present wins; required when ``True``.
_MAP: tuple[tuple[str, tuple[str, ...], bool], ...] = (
    ("peptide", ("Stripped.Sequence",), True),
    ("modified_peptide", ("Modified.Sequence",), True),
    ("protein_ids", ("Protein.Ids",), True),
    ("q_value", ("Global.Q.Value", "Q.Value"), True),
    ("pep", ("PEP",), False),
    ("is_decoy", ("Decoy",), False),
    ("intensity", ("Precursor.Normalised", "Precursor.Quantity"), False),
    ("run", ("Run",), True),
    ("charge", ("Precursor.Charge",), True),
)
_TEXT = ("peptide", "modified_peptide", "protein_ids", "run")
_WANTED = frozenset([n for _, names, _ in _MAP for n in names] + ["Protein.Group", "PG.MaxLFQ"])


def _read(path: Path) -> pd.DataFrame:
    """The mapped columns of one report: parquet columns keep their types, TSV cells are ``str``
    (``None`` for ``NA`` or empty)."""
    if path.suffix == ".parquet":
        names = [c for c in pq.read_schema(path).names if c in _WANTED]
        return pd.read_parquet(path, engine="pyarrow", columns=names)
    with open_text(path, encoding="utf-8-sig", newline="") as fh:
        raw = pd.read_csv(
            fh,
            sep="\t",
            dtype=str,
            usecols=lambda c: c.rstrip("\r") in _WANTED,
            keep_default_na=False,
            na_filter=False,
            quoting=csv.QUOTE_NONE,
            lineterminator="\n",
            engine="c",
        )
    raw.columns = [c.rstrip("\r") for c in raw.columns]
    raw = raw.astype(object)
    for col in raw.columns:
        cells = raw[col].str.rstrip("\r")
        raw[col] = cells.where(~cells.isin(("", NA)), None)
    return raw.reset_index(drop=True)


def _numeric(cells: pd.Series, col: str) -> tuple[pd.Series, np.ndarray]:
    """Type a numeric column: native numbers directly, text through the peptide-table parser."""
    if not is_numeric_dtype(cells.dtype) or is_bool_dtype(cells.dtype):
        return _typed(cells.astype(object), col)
    values = cells.to_numpy(dtype="float64", na_value=np.nan)
    if col != "charge":
        return pd.Series(values, dtype="float64"), np.zeros(len(values), dtype=bool)
    present = ~np.isnan(values)
    integral = present & np.isfinite(values)
    integral[integral] = values[integral] == np.round(values[integral])
    out = pd.array(np.where(integral, values, 0).astype("int64"), dtype="Int64")
    out[~integral] = pd.NA
    return pd.Series(out), present & ~integral


def _decoy(raw: pd.Series) -> tuple[pd.Series, np.ndarray]:
    num = pd.to_numeric(raw, errors="coerce").to_numpy(dtype="float64", na_value=np.nan)
    ok = (num == 0) | (num == 1)
    bad = raw.notna().to_numpy(dtype=bool) & ~ok
    out = pd.array(num == 1, dtype="boolean")
    out[~ok] = pd.NA
    return pd.Series(out), bad


def _strip(pep: pd.Series) -> pd.Series:
    plain = pep.str.fullmatch(r"[A-Z]+").fillna(True).to_numpy(dtype=bool)
    if plain.all():
        return pep
    pep = pep.copy()
    pep[~plain] = pep[~plain].map(strip_modifications)
    return pep


def _abundance(raw: pd.DataFrame) -> pd.DataFrame | None:
    if "Protein.Group" not in raw.columns or "PG.MaxLFQ" not in raw.columns:
        return None
    keep = raw["Protein.Group"].notna().to_numpy(dtype=bool, copy=True)
    if "Decoy" in raw.columns:
        decoy = pd.to_numeric(raw["Decoy"], errors="coerce").to_numpy(dtype="float64", na_value=0)
        keep &= decoy != 1
    value = pd.to_numeric(raw["PG.MaxLFQ"], errors="coerce").to_numpy(
        dtype="float64", na_value=np.nan
    )
    run = raw["Run"] if "Run" in raw.columns else pd.Series([""] * len(raw))
    frame = pd.DataFrame(
        {
            "protein_id": raw["Protein.Group"].astype(_STR).to_numpy(),
            "run": run.astype(_STR).fillna("").to_numpy(),
            "v": np.where(value == 0, np.nan, value),
        }
    )[keep]
    per_run = frame.groupby(["protein_id", "run"], sort=True)["v"].first()
    median = per_run.groupby(level="protein_id", sort=True).median().astype("float64").dropna()
    out = pd.DataFrame({"protein_id": median.index.astype(_STR), "maxlfq": median.to_numpy()})
    return out.reset_index(drop=True)


def convert_diann(
    paths: Sequence[Path], *, search_db: str | None, cfg: InputConfig
) -> ConvertResult:
    """Convert DIA-NN 2.x reports into the generic peptide table (and MaxLFQ abundance)."""
    issues: list[Issue] = []
    frames = []
    for path in paths:
        try:
            frames.append(_read(path))
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            issues.append(Issue("unreadable", "error", f"{path}: {exc}"))
    if issues:
        return ConvertResult(pd.DataFrame(), None, issues)
    raw = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    out: dict[str, pd.Series] = {}
    for col, names, required in _MAP:
        name = next((n for n in names if n in raw.columns), None)
        if name is None:
            if required:
                issues.append(
                    Issue(f"missing_column:{names[-1]}", "error", f"no {names[-1]!r} column")
                )
            continue
        cells = raw[name].reset_index(drop=True)
        if col in _TEXT:
            out[col] = _strip(cells.astype(_STR)) if col == "peptide" else cells.astype(_STR)
            continue
        typed, bad = _decoy(cells) if col == "is_decoy" else _numeric(cells, col)
        if col == "intensity":
            typed = typed.where(typed != 0)
        out[col] = typed
        issue = row_issue(f"bad_type:{name}", "error", f"{name}: invalid value", _rows(bad))
        if issue is not None:
            issues.append(issue)
    if any(i.severity == "error" for i in issues):
        return ConvertResult(pd.DataFrame(), None, issues)

    df = pd.DataFrame(out, index=pd.RangeIndex(len(raw)))
    df["source"] = pd.Series(["database"] * len(df), dtype=_STR)
    if search_db is not None:
        df["search_db"] = pd.Series([search_db] * len(df), dtype=_STR)
    _, found = validate_peptides(df, cfg)
    issues += found
    return ConvertResult(finish_peptides(df, list(df.columns)), _abundance(raw), issues)
