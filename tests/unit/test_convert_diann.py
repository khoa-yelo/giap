"""DIA-NN 2.x converter (issue M1-07; GIAP_SPEC.md §5.4; ADR-0005 item 2).

Fixtures are literal report rows using the DIA-NN 2.x column names, written to ``tmp_path`` as
``report.tsv`` and (from the same frame, with pyarrow) ``report.parquet``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
import pytest

from giap.config.model import InputConfig
from giap.io.converters import CONVERTERS
from giap.io.converters.base import ConvertResult
from giap.io.converters.diann import convert_diann
from giap.io.issues import Issue
from giap.io.peptides import PEPTIDE_COLUMNS

# Example 5 of the spec.
ROW5: dict[str, object] = {
    "Run": "r1",
    "Protein.Group": "P1;P2",
    "Protein.Ids": "P1;P2",
    "Stripped.Sequence": "PEPTIDEK",
    "Modified.Sequence": "PEPTIDEK",
    "Precursor.Charge": 2,
    "Q.Value": 0.002,
    "Global.Q.Value": 0.004,
    "PEP": 0.01,
    "Precursor.Quantity": 1400.0,
    "Precursor.Normalised": 1500.0,
    "PG.MaxLFQ": 10.0,
}


def _frame(rows: Sequence[dict[str, object]], drop: Sequence[str] = ()) -> pd.DataFrame:
    df = pd.DataFrame(list(rows))
    return df.drop(columns=list(drop))


def _tsv(tmp_path: Path, df: pd.DataFrame, name: str = "report.tsv") -> Path:
    p = tmp_path / name
    df.to_csv(p, sep="\t", index=False, lineterminator="\n")
    return p


def _parquet(tmp_path: Path, df: pd.DataFrame, name: str = "report.parquet") -> Path:
    p = tmp_path / name
    df.to_parquet(p, engine="pyarrow", index=False)
    return p


def _run(paths: Sequence[Path], search_db: str | None = None) -> ConvertResult:
    return convert_diann(list(paths), search_db=search_db, cfg=InputConfig())


def _codes(issues: list[Issue]) -> set[str]:
    return {i.code for i in issues}


def _errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def test_registered_in_converters() -> None:
    assert CONVERTERS["diann"] is convert_diann
    assert "generic" in CONVERTERS


def test_example5_mapping(tmp_path: Path) -> None:
    res = _run([_tsv(tmp_path, _frame([ROW5]))])
    assert isinstance(res, ConvertResult)
    assert not _errors(res.issues)
    df = res.peptides
    assert list(df.columns) == [
        "peptide",
        "modified_peptide",
        "protein_ids",
        "q_value",
        "pep",
        "intensity",
        "run",
        "charge",
        "source",
    ]
    assert len(df) == 1
    row = df.iloc[0]
    assert row["peptide"] == "PEPTIDEK"
    assert row["modified_peptide"] == "PEPTIDEK"
    assert row["protein_ids"] == "P1;P2"
    assert row["q_value"] == pytest.approx(0.004)
    assert row["pep"] == pytest.approx(0.01)
    assert row["intensity"] == pytest.approx(1500)
    assert row["run"] == "r1"
    assert int(row["charge"]) == 2
    assert row["source"] == "database"


def test_example6_fallback_columns(tmp_path: Path) -> None:
    df_in = _frame([ROW5], drop=["Global.Q.Value", "Precursor.Normalised"])
    res = _run([_tsv(tmp_path, df_in)])
    assert not _errors(res.issues)
    row = res.peptides.iloc[0]
    assert row["q_value"] == pytest.approx(0.002)
    assert row["intensity"] == pytest.approx(1400)


def test_columns_are_generic_and_in_spec_order(tmp_path: Path) -> None:
    res = _run([_tsv(tmp_path, _frame([{**ROW5, "Decoy": 0}]))], search_db="same")
    cols = list(res.peptides.columns)
    assert set(cols) <= set(PEPTIDE_COLUMNS)
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]
    assert "is_decoy" in cols
    assert "search_db" in cols
    assert set(res.peptides["search_db"]) == {"same"}


def test_decoy_column_maps_to_bool(tmp_path: Path) -> None:
    rows = [
        {**ROW5, "Decoy": 0},
        {**ROW5, "Run": "r2", "Stripped.Sequence": "KEDITPEP", "Decoy": 1},
    ]
    res = _run([_tsv(tmp_path, _frame(rows))])
    df = res.peptides.set_index("peptide")
    assert bool(df.loc["PEPTIDEK", "is_decoy"]) is False
    assert bool(df.loc["KEDITPEP", "is_decoy"]) is True


def test_decoy_absent_means_no_is_decoy_column(tmp_path: Path) -> None:
    res = _run([_tsv(tmp_path, _frame([ROW5]))])
    assert "is_decoy" not in res.peptides.columns


def test_zero_intensity_becomes_na(tmp_path: Path) -> None:
    rows = [ROW5, {**ROW5, "Run": "r2", "Precursor.Normalised": 0.0}]
    res = _run([_tsv(tmp_path, _frame(rows))])
    df = res.peptides.set_index("run")
    assert df.loc["r1", "intensity"] == pytest.approx(1500)
    assert _na(df.loc["r2", "intensity"])


def test_modified_sequence_kept_and_peptide_from_stripped(tmp_path: Path) -> None:
    row = {
        **ROW5,
        "Stripped.Sequence": "PEPMIDEK",
        "Modified.Sequence": "PEPM(UniMod:35)IDEK",
    }
    res = _run([_tsv(tmp_path, _frame([row]))])
    r = res.peptides.iloc[0]
    assert r["peptide"] == "PEPMIDEK"
    assert r["modified_peptide"] == "PEPM(UniMod:35)IDEK"


def test_rows_sorted_by_run_peptide_modified_charge(tmp_path: Path) -> None:
    rows = [
        {**ROW5, "Run": "r2", "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK"},
        {**ROW5, "Run": "r1", "Stripped.Sequence": "PEPK", "Modified.Sequence": "PEPK",
         "Precursor.Charge": 3},
        {**ROW5, "Run": "r1", "Stripped.Sequence": "PEPK", "Modified.Sequence": "PEPK",
         "Precursor.Charge": 2},
        {**ROW5, "Run": "r1", "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK"},
    ]  # fmt: skip
    res = _run([_tsv(tmp_path, _frame(rows))])
    got = [
        (r.run, r.peptide, int(r.charge))
        for r in res.peptides[["run", "peptide", "charge"]].itertuples()
    ]
    assert got == [("r1", "AAAK", 2), ("r1", "PEPK", 2), ("r1", "PEPK", 3), ("r2", "AAAK", 2)]


def test_parquet_matches_tsv(tmp_path: Path) -> None:
    rows = [
        {**ROW5, "Decoy": 0},
        {**ROW5, "Run": "r2", "Precursor.Normalised": 0.0, "Decoy": 0, "PG.MaxLFQ": 30.0},
    ]
    df_in = _frame(rows)
    a = _run([_tsv(tmp_path, df_in)])
    b = _run([_parquet(tmp_path, df_in)])
    pd.testing.assert_frame_equal(
        a.peptides.reset_index(drop=True), b.peptides.reset_index(drop=True)
    )
    assert a.abundance is not None and b.abundance is not None
    pd.testing.assert_frame_equal(
        a.abundance.reset_index(drop=True), b.abundance.reset_index(drop=True)
    )


def test_several_inputs_concatenated_in_order_and_not_deduplicated(tmp_path: Path) -> None:
    p = _tsv(tmp_path, _frame([ROW5]))
    res = _run([p, p])
    assert len(res.peptides) == 2


def test_example7_maxlfq_median_over_runs_zero_excluded(tmp_path: Path) -> None:
    rows = [
        {**ROW5, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW5, "Run": "r2", "PG.MaxLFQ": 30.0},
        {**ROW5, "Run": "r3", "PG.MaxLFQ": 0.0},
    ]
    res = _run([_tsv(tmp_path, _frame(rows))])
    ab = res.abundance
    assert ab is not None
    assert list(ab.columns) == ["protein_id", "maxlfq"]
    assert ab["protein_id"].tolist() == ["P1;P2"]
    assert ab["maxlfq"].iloc[0] == pytest.approx(20)


def test_maxlfq_one_value_per_run_not_per_precursor(tmp_path: Path) -> None:
    # PG.MaxLFQ is repeated on every precursor row of a run; the median is over runs.
    rows = [
        {**ROW5, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW5, "Run": "r1", "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK",
         "PG.MaxLFQ": 10.0},
        {**ROW5, "Run": "r2", "PG.MaxLFQ": 30.0},
    ]  # fmt: skip
    res = _run([_tsv(tmp_path, _frame(rows))])
    assert res.abundance is not None
    assert res.abundance["maxlfq"].iloc[0] == pytest.approx(20)


def test_maxlfq_one_row_per_protein_group_na_excluded(tmp_path: Path) -> None:
    rows = [
        {**ROW5, "Run": "r1", "PG.MaxLFQ": 10.0},
        {**ROW5, "Run": "r2", "PG.MaxLFQ": float("nan")},
        {**ROW5, "Run": "r1", "Protein.Group": "Q9", "Protein.Ids": "Q9",
         "Stripped.Sequence": "AAAK", "Modified.Sequence": "AAAK", "PG.MaxLFQ": 7.0},
    ]  # fmt: skip
    res = _run([_tsv(tmp_path, _frame(rows))])
    ab = res.abundance
    assert ab is not None
    assert sorted(ab["protein_id"]) == ["P1;P2", "Q9"]
    vals = dict(zip(ab["protein_id"], ab["maxlfq"], strict=True))
    assert vals["P1;P2"] == pytest.approx(10)
    assert vals["Q9"] == pytest.approx(7)


def test_example8_missing_stripped_sequence(tmp_path: Path) -> None:
    res = _run([_tsv(tmp_path, _frame([ROW5], drop=["Stripped.Sequence"]))])
    errs = _errors(res.issues)
    assert "missing_column:Stripped.Sequence" in _codes(errs)


def test_empty_report_is_error(tmp_path: Path) -> None:
    res = _run([_tsv(tmp_path, _frame([ROW5]).iloc[0:0])])
    assert "empty_table" in _codes(_errors(res.issues))
