"""C16 ``peptide_runs`` (work/peptide_runs.parquet, M1-06/M1-12) joins the M4-06 contract as an
optional table (docs/specs/M6-02-detection-limit-model.md lead decision 1; M6-index decision 9).

Pinned: ``TABLES["peptide_runs"]`` has exactly the M1-06 columns ``PEPTIDE_RUN_COLUMNS``
(``peptide`` str, ``run`` str, ``intensity`` float nullable, ``n_rows`` int); ``validate_table``
accepts M1-06's own output and reports bad frames as ``contract.peptide_runs.*`` issues; the
loader reads it through ``load_tables`` like C15 (a missing file is a warning, never an error:
the table is optional). ``M4Inputs`` and ``M4_INPUT_KEYS`` are unchanged (as for C15).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from gippa.io.issues import Issue
from gippa.io.peptides import PEPTIDE_RUN_COLUMNS, aggregate_peptides


def _errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def _runs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": pd.Series(["AAAAAAK", "AAAAAAK", "CCCCCCK"], dtype=object),
            "run": pd.Series(["r1", "r2", "r1"], dtype=object),
            "intensity": pd.Series([10.0, None, 3.5], dtype="float64"),
            "n_rows": pd.Series([1, 2, 1], dtype="int64"),
        }
    )


def test_tables_has_c16_with_the_m1_06_columns() -> None:
    from gippa.genetable.contract import CONTRACT_COLUMNS, TABLES

    assert "peptide_runs" in TABLES
    assert CONTRACT_COLUMNS["peptide_runs"] == tuple(PEPTIDE_RUN_COLUMNS)
    cols = {c.name: c for c in TABLES["peptide_runs"]}
    assert cols["peptide"].dtype == "str" and not cols["peptide"].nullable
    assert cols["run"].dtype == "str" and not cols["run"].nullable
    assert cols["intensity"].dtype == "float" and cols["intensity"].nullable
    assert cols["n_rows"].dtype == "int" and not cols["n_rows"].nullable


def test_valid_frame_has_no_issues() -> None:
    from gippa.genetable.contract import validate_table

    assert validate_table("peptide_runs", _runs()) == []


def test_m1_06_output_validates(tmp_path: Path) -> None:
    """The real producer's frame (M1-06 read -> validate -> ``aggregate_peptides``) passes C16,
    including the missing-run convention (run ``""``)."""
    from gippa.config.model import InputConfig
    from gippa.genetable.contract import validate_table
    from gippa.io.peptides import read_peptide_table, validate_peptides

    path = tmp_path / "peptides.tsv"
    lines = [
        "peptide|q_value|run|charge|intensity",
        "PEPTIDEK|0.004|r1|2|100",
        "PEPTIDEK|0.001|r1|3|50",
        "PEPTIDEK|0.02|r2|2|90",
        "SAMPLERK|0.01|r2|2|",
    ]
    path.write_text("\n".join(x.replace("|", "\t") for x in lines) + "\n", encoding="utf-8")
    raw, _ = read_peptide_table(path)
    df, issues = validate_peptides(raw, InputConfig())
    assert _errors(issues) == []
    _, runs = aggregate_peptides(df)
    assert _errors(validate_table("peptide_runs", runs)) == []
    no_run = tmp_path / "norun.tsv"
    no_run.write_text("peptide\nPEPTIDEK\n", encoding="utf-8")
    raw2, _ = read_peptide_table(no_run)
    df2, _ = validate_peptides(raw2, InputConfig())
    _, runs2 = aggregate_peptides(df2)
    assert list(runs2["run"]) == [""]
    assert _errors(validate_table("peptide_runs", runs2)) == []


@pytest.mark.parametrize("drop", PEPTIDE_RUN_COLUMNS)
def test_missing_column_is_error(drop: str) -> None:
    from gippa.genetable.contract import validate_table

    issues = validate_table("peptide_runs", _runs().drop(columns=[drop]))
    assert any(
        i.code == "contract.peptide_runs.missing_column" and i.severity == "error" for i in issues
    )


def test_bad_dtype_and_na_are_errors() -> None:
    from gippa.genetable.contract import validate_table

    bad = _runs().assign(n_rows=["one", "two", "three"])
    assert any(
        i.code == "contract.peptide_runs.dtype"
        for i in _errors(validate_table("peptide_runs", bad))
    )
    na = _runs()
    na.loc[0, "peptide"] = None
    assert any(
        i.code.startswith("contract.peptide_runs.")
        for i in _errors(validate_table("peptide_runs", na))
    )


def test_duplicate_peptide_run_is_error() -> None:
    """M1-06 aggregates per (peptide, run): a repeated pair is a malformed C16."""
    from gippa.genetable.contract import validate_table

    dup = pd.concat([_runs(), _runs().iloc[[0]]], ignore_index=True)
    errs = _errors(validate_table("peptide_runs", dup))
    assert any(i.code.startswith("contract.peptide_runs.") for i in errs)
    assert _errors(validate_table("peptide_runs", _runs())) == []


def test_load_tables_reads_c16(tmp_path: Path) -> None:
    from gippa.io.work_tables import load_tables

    path = tmp_path / "peptide_runs.parquet"
    _runs().to_parquet(path, engine="pyarrow", index=False)
    frames, issues = load_tables({"peptide_runs": path}, ["peptide_runs"])
    assert _errors(issues) == []
    pd.testing.assert_frame_equal(
        frames["peptide_runs"].reset_index(drop=True), _runs(), check_dtype=False
    )


def test_missing_c16_is_a_warning_not_an_error(tmp_path: Path) -> None:
    from gippa.io.work_tables import load_tables

    frames, issues = load_tables({"peptide_runs": tmp_path / "absent.parquet"}, ["peptide_runs"])
    assert _errors(issues) == []
    assert any(
        i.severity == "warn" and i.code == "contract.peptide_runs.missing_file" for i in issues
    )
    assert len(frames["peptide_runs"]) == 0
    assert set(PEPTIDE_RUN_COLUMNS) <= set(frames["peptide_runs"].columns)
