"""DIA-NN 2.x converter (``report.parquet`` or ``report.tsv``; GIAP_SPEC.md §5.4).

``peptide`` ← ``Stripped.Sequence``; ``modified_peptide`` ← ``Modified.Sequence``; ``protein_ids``
← ``Protein.Ids``; ``q_value`` ← ``Global.Q.Value``, else ``Q.Value``; ``pep`` ← ``PEP``;
``is_decoy`` ← ``Decoy`` (0/1); ``intensity`` ← ``Precursor.Normalised``, else
``Precursor.Quantity`` (0 → missing); ``run`` ← ``Run``; ``charge`` ← ``Precursor.Charge``;
``source`` = database. Abundance: per ``Protein.Group``, the median over runs of ``PG.MaxLFQ``
(missing and 0 excluded), column ``maxlfq``.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from giap.config.model import InputConfig
from giap.io.converters.base import ConvertResult, finish_peptides
from giap.io.issues import Issue, row_issue
from giap.io.peptides import _rows, _typed, validate_peptides
from giap.io.tsv import read_tsv

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


def _read(path: Path) -> pd.DataFrame:
    """One report as an object frame of ``str`` cells (``None`` for missing), any format."""
    if path.suffix == ".parquet":
        raw = pd.read_parquet(path, engine="pyarrow")
        return pd.DataFrame(
            {c: [None if pd.isna(v) else str(v) for v in raw[c]] for c in raw.columns},
            index=pd.RangeIndex(len(raw)),
            dtype=object,
        )
    header, rows = read_tsv(path)
    return pd.DataFrame(rows, columns=header, dtype=object)


def _decoy(raw: pd.Series) -> tuple[pd.Series, np.ndarray]:
    num = pd.to_numeric(raw, errors="coerce")
    ok = num.isin((0, 1)).to_numpy(dtype=bool)
    bad = raw.notna().to_numpy(dtype=bool) & ~ok
    out = pd.Series(pd.array([None] * len(raw), dtype="boolean"))
    out[ok] = num[ok].to_numpy() == 1
    return out, bad


def _abundance(raw: pd.DataFrame) -> pd.DataFrame | None:
    if "Protein.Group" not in raw.columns or "PG.MaxLFQ" not in raw.columns:
        return None
    pg = raw["Protein.Group"].astype(_STR)
    value = pd.to_numeric(raw["PG.MaxLFQ"], errors="coerce").astype("float64")
    run = raw["Run"].astype(_STR).fillna("") if "Run" in raw.columns else pd.Series([""] * len(raw))
    frame = pd.DataFrame({"protein_id": pg, "run": run, "v": value.where(value != 0)})
    frame = frame[frame["protein_id"].notna()]
    per_run = frame.groupby(["protein_id", "run"], sort=True)["v"].first()
    median = per_run.groupby(level="protein_id", sort=True).median().astype("float64")
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
    raw = raw.astype(object).where(raw.notna(), None)

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
            out[col] = cells.astype(_STR)
            continue
        typed, bad = _decoy(cells) if col == "is_decoy" else _typed(cells, col)
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
