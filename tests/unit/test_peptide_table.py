"""Generic peptide table: read, validate, aggregate (issue M1-06; GIPPA_SPEC.md §5.3; ADR-0005/6).

Tables are literal TSV text written to ``tmp_path``. ``_load`` runs ``read_peptide_table`` then
``validate_peptides`` and returns the issues of both steps, so a check may live in either one.
"""

from __future__ import annotations

import dataclasses
import gzip
import math
from pathlib import Path

import pandas as pd
import pytest

from gippa.config.model import InputConfig
from gippa.io.issues import Issue as IssuesIssue
from gippa.io.peptides import (
    PEPTIDE_COLUMNS,
    Issue,
    aggregate_peptides,
    read_peptide_table,
    validate_peptides,
)


def _write(tmp_path: Path, text: str, name: str = "peptides.tsv") -> Path:
    p = tmp_path / name
    data = text.encode("utf-8")
    if name.endswith(".gz"):
        data = gzip.compress(data)
    p.write_bytes(data)
    return p


def _load(
    tmp_path: Path, text: str, cfg: InputConfig | None = None, name: str = "peptides.tsv"
) -> tuple[pd.DataFrame, list[Issue]]:
    df_raw, issues_read = read_peptide_table(_write(tmp_path, text, name))
    df, issues_val = validate_peptides(df_raw, cfg or InputConfig())
    return df, [*issues_read, *issues_val]


def _codes(issues: list[Issue]) -> set[str]:
    return {i.code for i in issues}


def _issue(issues: list[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert found, f"no issue {code!r} in {[i.code for i in issues]}"
    return found[0]


def _errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def _row(agg: pd.DataFrame, peptide: str) -> pd.Series:
    sel = agg[agg["peptide"] == peptide]
    assert len(sel) == 1, f"expected one row for {peptide}, got {len(sel)}"
    return sel.iloc[0]


def _tsv(*lines: str) -> str:
    return "\n".join(line.replace("|", "\t") for line in lines) + "\n"


# --- interface ----------------------------------------------------------------------------------


def test_peptide_columns_in_spec_order() -> None:
    assert PEPTIDE_COLUMNS == (
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


def test_issue_is_shared_frozen_dataclass() -> None:
    assert Issue is IssuesIssue
    i = Issue(code="x", severity="warn", message="m")
    assert i.n_rows == 0 and i.examples == ()
    with pytest.raises(dataclasses.FrozenInstanceError):
        i.code = "y"  # type: ignore[misc]


# --- Examples 1-10 ------------------------------------------------------------------------------


def test_example_01_two_peptides_il_kept(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide", "PEPTIDEK", "LLLLIIIK"))
    assert [i for i in issues if i.severity in ("error", "warn")] == []
    agg, _ = aggregate_peptides(df)
    assert list(agg["peptide"]) == ["LLLLIIIK", "PEPTIDEK"]


def test_example_02_missing_peptide_column(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("sequence|q_value", "PEPTIDEK|0.01"))
    assert _issue(issues, "missing_peptide_column").severity == "error"


def test_example_03_q_value_out_of_range(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("peptide|q_value", "PEPTIDEK|1.2"))
    iss = _issue(issues, "out_of_range:q_value")
    assert iss.severity == "error"
    assert iss.examples == ("1",)
    assert iss.n_rows == 1


def test_example_04_intensity_sum_within_run_median_across(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|q_value|run|charge|intensity",
        "PEPTIDEK|0.004|r1|2|100",
        "PEPTIDEK|0.001|r1|3|50",
        "PEPTIDEK|0.02|r2|2|90",
    )
    df, issues = _load(tmp_path, text)
    assert _errors(issues) == []
    agg, runs = aggregate_peptides(df)
    assert len(agg) == 1
    r = _row(agg, "PEPTIDEK")
    assert r["q_value"] == pytest.approx(0.001)
    assert r["intensity"] == pytest.approx(120.0)
    assert r["n_runs"] == 2
    assert r["n_rows"] == 3
    assert list(runs.columns[:4]) == ["peptide", "run", "intensity", "n_rows"]
    assert list(runs["run"]) == ["r1", "r2"]
    assert list(runs["intensity"]) == pytest.approx([150.0, 90.0])
    assert list(runs["n_rows"]) == [2, 1]


def test_example_05_is_decoy_true_only_if_all_rows(tmp_path: Path) -> None:
    text = _tsv("peptide|protein_ids|is_decoy", "PEPTIDEK|rev_P1|true", "PEPTIDEK|P1|false")
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert not bool(_row(agg, "PEPTIDEK")["is_decoy"])


def test_example_06_protein_ids_union(tmp_path: Path) -> None:
    text = _tsv("peptide|protein_ids", "PEPTIDEK|P1;P2", "PEPTIDEK|P2;P3")
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert list(_row(agg, "PEPTIDEK")["protein_ids"]) == ["P1", "P2", "P3"]


def test_example_07_decoy_from_prefix(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|protein_ids", "AAAAKPEP|rev_P1"))
    assert bool(df["is_decoy"].iloc[0])
    agg, _ = aggregate_peptides(df)
    assert bool(_row(agg, "AAAAKPEP")["is_decoy"])


def test_example_08_short_peptide_dropped(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide", "PEPK"), InputConfig(peptide_len=[7, 40]))
    iss = _issue(issues, "short_peptide")
    assert iss.severity == "warn"
    assert iss.n_rows == 1
    assert len(df) == 0


def test_example_09_modification_syntax_dropped(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide", "PEPT(UniMod:21)IDEK", "PEPTIDEK"))
    iss = _issue(issues, "non_letters")
    assert iss.severity == "warn"
    assert iss.n_rows == 1
    assert iss.examples == ("1",)
    assert list(df["peptide"]) == ["PEPTIDEK"]


def test_example_10_header_only_is_empty_table(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("peptide|q_value"))
    assert _issue(issues, "empty_table").severity == "error"


# --- validation rules ---------------------------------------------------------------------------


def test_long_peptide_dropped_boundary_kept(tmp_path: Path) -> None:
    p40, p41 = "A" * 39 + "K", "A" * 40 + "K"
    df, issues = _load(tmp_path, _tsv("peptide", p40, p41))
    iss = _issue(issues, "long_peptide")
    assert iss.severity == "warn"
    assert iss.n_rows == 1
    assert iss.examples == ("2",)
    assert list(df["peptide"]) == [p40]


def test_lowercase_peptide_is_non_letters(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide", "peptidek", "PEPTIDEK"))
    assert _issue(issues, "non_letters").severity == "warn"
    assert list(df["peptide"]) == ["PEPTIDEK"]


def test_peptide_len_from_config(tmp_path: Path) -> None:
    df, issues = _load(
        tmp_path, _tsv("peptide", "PEPK", "PEPTIDEK"), InputConfig(peptide_len=[4, 6])
    )
    assert "short_peptide" not in _codes(issues)
    assert _issue(issues, "long_peptide").n_rows == 1
    assert list(df["peptide"]) == ["PEPK"]


@pytest.mark.parametrize("col", ["q_value", "pep"])
@pytest.mark.parametrize("value", ["-0.1", "1.5", "INF", "-INF"])
def test_probability_out_of_range(tmp_path: Path, col: str, value: str) -> None:
    _, issues = _load(tmp_path, _tsv(f"peptide|{col}", "PEPTIDEK|0.5", f"PEPTIDEK|{value}"))
    iss = _issue(issues, f"out_of_range:{col}")
    assert iss.severity == "error"
    assert iss.examples == ("2",)


@pytest.mark.parametrize("value", ["0", "1"])
def test_probability_bounds_inclusive(tmp_path: Path, value: str) -> None:
    _, issues = _load(tmp_path, _tsv("peptide|q_value|pep", f"PEPTIDEK|{value}|{value}"))
    assert _errors(issues) == []


def test_negative_intensity_and_zero_charge_out_of_range(tmp_path: Path) -> None:
    text = _tsv("peptide|intensity|charge", "PEPTIDEK|-1|2", "PEPTIDEK|5|0")
    _, issues = _load(tmp_path, text)
    assert _issue(issues, "out_of_range:intensity").examples == ("1",)
    assert _issue(issues, "out_of_range:charge").examples == ("2",)
    assert all(i.severity == "error" for i in issues if i.code.startswith("out_of_range"))


@pytest.mark.parametrize(
    ("col", "value"),
    [
        ("q_value", "abc"),
        ("intensity", "lots"),
        ("score", "x1"),
        ("charge", "2.5"),
        ("charge", "two"),
        ("is_decoy", "yes"),
        ("is_decoy", "2"),
    ],
)
def test_bad_type(tmp_path: Path, col: str, value: str) -> None:
    _, issues = _load(tmp_path, _tsv(f"peptide|{col}", f"PEPTIDEK|{value}"))
    iss = _issue(issues, f"bad_type:{col}")
    assert iss.severity == "error"
    assert iss.examples == ("1",)


def test_integral_float_charge_accepted(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide|charge", "PEPTIDEK|2.0"))
    assert "bad_type:charge" not in _codes(issues)
    assert int(df["charge"].iloc[0]) == 2


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        *[(v, True) for v in ("true", "TRUE", "True", "1")],
        *[(v, False) for v in ("false", "FALSE", "False", "0")],
    ],
)
def test_bool_spellings(tmp_path: Path, value: str, expected: bool) -> None:
    df, issues = _load(tmp_path, _tsv("peptide|protein_ids|is_decoy", f"PEPTIDEK|P1|{value}"))
    assert _errors(issues) == []
    assert bool(df["is_decoy"].iloc[0]) is expected


def test_adr0006_na_and_inf_encoding(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|q_value|score|intensity|run",
        "PEPTIDEK|NA|INF|NA|r1",
        "PEPTIDEK||-INF||r2",
        "AAAAKPEP|0.01|-INF|10|r1",
    )
    df, issues = _load(tmp_path, text)
    assert _errors(issues) == []
    assert df["q_value"].isna().sum() == 2
    agg, _ = aggregate_peptides(df)
    r = _row(agg, "PEPTIDEK")
    assert r["score"] == math.inf
    assert pd.isna(r["q_value"])
    assert pd.isna(r["intensity"])
    assert _row(agg, "AAAAKPEP")["score"] == -math.inf


@pytest.mark.parametrize(("col", "value"), [("source", "manual"), ("search_db", "maybe")])
def test_bad_enum(tmp_path: Path, col: str, value: str) -> None:
    _, issues = _load(tmp_path, _tsv(f"peptide|{col}", f"PEPTIDEK|{value}"))
    iss = _issue(issues, f"bad_enum:{col}")
    assert iss.severity == "error"
    assert iss.examples == ("1",)


def test_duplicate_spectrum_id_within_run(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|spectrum_id",
        "PEPTIDEK|r1|s1",
        "AAAAKPEP|r1|s1",
        "AAAAKPEP|r2|s1",
    )
    _, issues = _load(tmp_path, text)
    iss = _issue(issues, "duplicate_spectrum_id")
    assert iss.severity == "error"
    assert iss.n_rows >= 1


def test_same_spectrum_id_same_peptide_is_fine(tmp_path: Path) -> None:
    text = _tsv("peptide|run|spectrum_id|charge", "PEPTIDEK|r1|s1|2", "PEPTIDEK|r1|s1|3")
    _, issues = _load(tmp_path, text)
    assert "duplicate_spectrum_id" not in _codes(issues)


def test_spectrum_id_without_run_checked_over_whole_table(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("peptide|spectrum_id", "PEPTIDEK|s1", "AAAAKPEP|s1"))
    assert _issue(issues, "duplicate_spectrum_id").severity == "error"


def test_decoy_without_protein_ids(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("peptide|protein_ids|is_decoy", "PEPTIDEK|NA|true"))
    iss = _issue(issues, "decoy_without_protein_ids")
    assert iss.severity == "warn"
    assert iss.examples == ("1",)


def test_unknown_column_kept_with_info(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide|my_note", "PEPTIDEK|hello"))
    assert _issue(issues, "unknown_column").severity == "info"
    assert "my_note" in df.columns
    assert df["my_note"].iloc[0] == "hello"


def test_examples_capped_at_five(tmp_path: Path) -> None:
    rows = ["PEPTIDEK|2"] * 8
    _, issues = _load(tmp_path, _tsv("peptide|q_value", *rows))
    iss = _issue(issues, "out_of_range:q_value")
    assert iss.n_rows == 8
    assert iss.examples == ("1", "2", "3", "4", "5")


def test_rows_above_q_value_max_are_kept(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|q_value", "PEPTIDEK|0.5"), InputConfig(q_value_max=0.01))
    assert len(df) == 1


def test_gzipped_input(tmp_path: Path) -> None:
    df, issues = _load(tmp_path, _tsv("peptide", "PEPTIDEK"), name="peptides.tsv.gz")
    assert _errors(issues) == []
    assert list(df["peptide"]) == ["PEPTIDEK"]


def test_validate_is_pure(tmp_path: Path) -> None:
    df_raw, _ = read_peptide_table(_write(tmp_path, _tsv("peptide", "PEPK", "PEPTIDEK")))
    before = df_raw.copy(deep=True)
    validate_peptides(df_raw, InputConfig())
    pd.testing.assert_frame_equal(df_raw, before)


# --- derived columns from config ----------------------------------------------------------------


def test_decoy_and_contaminant_prefixes_from_config(tmp_path: Path) -> None:
    cfg = InputConfig(decoy_prefixes=["DEC_"], contaminant_prefixes=["KONT_"])
    text = _tsv(
        "peptide|protein_ids",
        "AAAAAAAK|DEC_P1;DEC_P2",
        "CCCCCCCK|DEC_P1;P2",
        "DDDDDDDK|rev_P1",
        "EEEEEEEK|KONT_X;KONT_Y",
        "FFFFFFFK|KONT_X;P9",
    )
    df, _ = _load(tmp_path, text, cfg)
    by = df.set_index("peptide")
    assert bool(by.loc["AAAAAAAK", "is_decoy"])
    assert not bool(by.loc["CCCCCCCK", "is_decoy"])
    assert not bool(by.loc["DDDDDDDK", "is_decoy"])  # rev_ is not in this config
    assert bool(by.loc["EEEEEEEK", "is_contaminant"])
    assert not bool(by.loc["FFFFFFFK", "is_contaminant"])
    assert not bool(by.loc["AAAAAAAK", "is_contaminant"])


def test_default_contaminant_prefix(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|protein_ids", "PEPTIDEK|CON__P1"))
    assert bool(df["is_contaminant"].iloc[0])
    assert not bool(df["is_decoy"].iloc[0])


def test_explicit_is_decoy_overrides_prefix(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide|protein_ids|is_decoy", "PEPTIDEK|rev_P1|false"))
    assert not bool(df["is_decoy"].iloc[0])


def test_no_protein_ids_means_not_decoy(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide", "PEPTIDEK"))
    assert not bool(df["is_decoy"].iloc[0])
    assert not bool(df["is_contaminant"].iloc[0])


@pytest.mark.parametrize("source", ["database", "denovo", "auto"])
def test_source_default_from_config(tmp_path: Path, source: str) -> None:
    text = _tsv("peptide|source", "PEPTIDEK|NA", "AAAAKPEP|denovo")
    df, _ = _load(tmp_path, text, InputConfig(peptide_source=source))
    by = df.set_index("peptide")
    assert by.loc["PEPTIDEK", "source"] == source
    assert by.loc["AAAAKPEP", "source"] == "denovo"


@pytest.mark.parametrize("search_db", ["same", "other", "auto"])
def test_search_db_default_from_config(tmp_path: Path, search_db: str) -> None:
    df, _ = _load(tmp_path, _tsv("peptide", "PEPTIDEK"), InputConfig(search_db=search_db))
    assert df["search_db"].iloc[0] == search_db


def test_mixed_source_database_wins(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|source",
        "PEPTIDEK|denovo",
        "PEPTIDEK|database",
        "AAAAKPEP|denovo",
        "AAAAKPEP|denovo",
    )
    df, issues = _load(tmp_path, text)
    iss = _issue(issues, "mixed_source")
    assert iss.severity == "warn"
    assert iss.n_rows >= 1
    agg, _ = aggregate_peptides(df)
    assert _row(agg, "PEPTIDEK")["source"] == "database"
    assert _row(agg, "AAAAKPEP")["source"] == "denovo"


# --- aggregation --------------------------------------------------------------------------------


def test_aggregation_min_max_rules(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|q_value|pep|score",
        "PEPTIDEK|0.01|0.2|10",
        "PEPTIDEK|0.003|0.5|42",
        "PEPTIDEK|NA|0.05|7",
    )
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    r = _row(agg, "PEPTIDEK")
    assert r["q_value"] == pytest.approx(0.003)
    assert r["pep"] == pytest.approx(0.05)
    assert r["score"] == pytest.approx(42.0)
    assert r["n_rows"] == 3


def test_is_decoy_all_rows_true(tmp_path: Path) -> None:
    text = _tsv("peptide|protein_ids|is_decoy", "PEPTIDEK|rev_P1|true", "PEPTIDEK|rev_P2|true")
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert bool(_row(agg, "PEPTIDEK")["is_decoy"])


def test_is_contaminant_all_rows(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|protein_ids",
        "PEPTIDEK|CON__A",
        "PEPTIDEK|P1",
        "AAAAKPEP|CON__A",
        "AAAAKPEP|cRAP01",
    )
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert not bool(_row(agg, "PEPTIDEK")["is_contaminant"])
    assert bool(_row(agg, "AAAAKPEP")["is_contaminant"])


def test_modified_peptides_sorted_union(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|modified_peptide",
        "PEPTIDEK|PEPT(UniMod:21)IDEK",
        "PEPTIDEK|PEPTIDEK",
        "PEPTIDEK|PEPT(UniMod:21)IDEK",
    )
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert list(_row(agg, "PEPTIDEK")["modified_peptides"]) == ["PEPT(UniMod:21)IDEK", "PEPTIDEK"]


def test_intensity_summed_over_modified_forms_within_run(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|modified_peptide|run|intensity",
        "PEPTIDEK|PEPTIDEK|r1|10",
        "PEPTIDEK|PEPT(UniMod:21)IDEK|r1|5",
        "PEPTIDEK|PEPTIDEK|r2|30",
        "PEPTIDEK|PEPTIDEK|r3|40",
    )
    df, _ = _load(tmp_path, text)
    agg, runs = aggregate_peptides(df)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(30.0)  # median(15, 30, 40)
    assert list(runs["intensity"]) == pytest.approx([15.0, 30.0, 40.0])


def test_median_uses_only_quantified_runs(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|intensity",
        "PEPTIDEK|r1|100",
        "PEPTIDEK|r2|NA",
        "PEPTIDEK|r3|50",
        "AAAAKPEP|r1|7",
    )
    df, _ = _load(tmp_path, text)
    agg, runs = aggregate_peptides(df)
    r = _row(agg, "PEPTIDEK")
    assert r["intensity"] == pytest.approx(75.0)  # not median(100, 0, 50) = 50
    assert r["n_runs"] == 3
    pr = runs[runs["peptide"] == "PEPTIDEK"].set_index("run")
    assert list(pr.index) == ["r1", "r2", "r3"]
    assert pd.isna(pr.loc["r2", "intensity"])
    assert pr.loc["r2", "n_rows"] == 1
    assert _row(agg, "AAAAKPEP")["intensity"] == pytest.approx(7.0)


def test_zero_intensity_is_a_value(tmp_path: Path) -> None:
    text = _tsv("peptide|run|intensity", "PEPTIDEK|r1|0", "PEPTIDEK|r2|10")
    df, _ = _load(tmp_path, text)
    agg, runs = aggregate_peptides(df)
    assert _row(agg, "PEPTIDEK")["intensity"] == pytest.approx(5.0)
    assert list(runs["intensity"]) == pytest.approx([0.0, 10.0])


def test_all_intensities_missing_gives_na(tmp_path: Path) -> None:
    text = _tsv("peptide|run|intensity", "PEPTIDEK|r1|NA", "PEPTIDEK|r2|NA")
    df, _ = _load(tmp_path, text)
    agg, _ = aggregate_peptides(df)
    assert pd.isna(_row(agg, "PEPTIDEK")["intensity"])


def test_missing_run_is_single_empty_run(tmp_path: Path) -> None:
    text = _tsv("peptide|intensity", "PEPTIDEK|10", "PEPTIDEK|20")
    df, _ = _load(tmp_path, text)
    agg, runs = aggregate_peptides(df)
    r = _row(agg, "PEPTIDEK")
    assert r["n_runs"] == 1
    assert r["intensity"] == pytest.approx(30.0)
    assert list(runs["run"]) == [""]
    assert list(runs["n_rows"]) == [2]


def test_outputs_sorted(tmp_path: Path) -> None:
    text = _tsv(
        "peptide|run|intensity",
        "YYYYYYYK|r2|1",
        "AAAAAAAK|r2|1",
        "YYYYYYYK|r1|1",
        "AAAAAAAK|r1|1",
    )
    df, _ = _load(tmp_path, text)
    agg, runs = aggregate_peptides(df)
    assert list(agg["peptide"]) == ["AAAAAAAK", "YYYYYYYK"]
    assert list(zip(runs["peptide"], runs["run"], strict=True)) == [
        ("AAAAAAAK", "r1"),
        ("AAAAAAAK", "r2"),
        ("YYYYYYYK", "r1"),
        ("YYYYYYYK", "r2"),
    ]


def test_il_not_merged(tmp_path: Path) -> None:
    df, _ = _load(tmp_path, _tsv("peptide", "PEPTLDEK", "PEPTIDEK"))
    agg, _ = aggregate_peptides(df)
    assert list(agg["peptide"]) == ["PEPTIDEK", "PEPTLDEK"]
