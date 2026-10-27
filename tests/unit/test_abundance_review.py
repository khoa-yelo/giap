"""Lead decisions from the M1-08 code review (abundance table; GIPPA_SPEC.md §5.3, §7.0; ADR-0005/6).

1. a table with valid shape is always returned, even with ``non_numeric_values``/``non_finite``
   errors: bad cells are ``None`` and flagged in ``bad_rows`` (one bool per row, row order);
   counts and the detected scale stay available. Only shape errors return ``table=None``;
2. ``AbundanceTable.raw_ids`` holds the IDs as written, ``ids`` the resolved IDs;
3. ``scale_undetermined`` only when the requested scale is ``auto``
   (``validate_abundance`` takes the requested scale);
4. protein groups keep written order (``P2;P3`` and ``P3;P2`` are not duplicates);
5. duplicate-ID examples name the later row(s); two blank IDs give only ``empty_id``.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import pytest

from gippa.io.abundance import (
    AbundanceTable,
    abundance_stats,
    detect_scale,
    read_abundance,
    validate_abundance,
)
from gippa.io.issues import Issue

WARN_FRAC = 0.20


def _tsv(*lines: str) -> str:
    return "\n".join(line.replace("|", "\t") for line in lines) + "\n"


def _read(tmp_path: Path, text: str) -> tuple[AbundanceTable | None, list[Issue]]:
    p = tmp_path / "abundance.tsv"
    p.write_text(text, encoding="utf-8")
    return read_abundance(p)


def _codes(issues: list[Issue]) -> list[str]:
    return [i.code for i in issues]


def _issue(issues: list[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert found, f"no issue {code!r} in {_codes(issues)}"
    return found[0]


# --- 1. table returned despite value errors ------------------------------------------------------


def test_table_fields_include_raw_ids_and_bad_rows() -> None:
    names = {f.name for f in dataclasses.fields(AbundanceTable)}
    assert {"id_column", "value_column", "raw_ids", "ids", "values", "bad_rows"} <= names
    assert "ignored_columns" in names


def test_non_numeric_value_still_returns_table_with_bad_cell_none(tmp_path: Path) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1e6", "P2|abc", "P3|2e6"))
    assert _issue(issues, "non_numeric_values").severity == "error"
    assert table is not None
    assert table.ids == ("P1", "P2", "P3")
    assert table.values == (1e6, None, 2e6)
    assert tuple(table.bad_rows) == (False, True, False)


@pytest.mark.parametrize("text", ["inf", "-INF", "nan", "NaN", "1e400"])
def test_non_finite_value_still_returns_table_with_bad_cell_none(tmp_path: Path, text: str) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", f"P1|{text}", "P2|1e6"))
    nf = _issue(issues, "non_finite")
    assert nf.severity == "error"
    assert nf.examples == ("1",)
    assert table is not None
    assert table.values == (None, 1e6)
    assert tuple(table.bad_rows) == (True, False)


def test_bad_rows_all_false_for_clean_table(tmp_path: Path) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1e6", "P2|NA", "P3|0"))
    assert table is not None
    assert [i for i in issues if i.severity == "error"] == []
    assert tuple(table.bad_rows) == (False, False, False)
    assert all(isinstance(b, bool) for b in table.bad_rows)


def test_counts_and_scale_available_with_bad_cells(tmp_path: Path) -> None:
    text = _tsv("protein_id|v", "P1|1e6", "P2|abc", "P3|NA", "P4|0", "P5|inf", "P6|2e6")
    table, issues = _read(tmp_path, text)
    assert {"non_numeric_values", "non_finite"} <= set(_codes(issues))
    assert table is not None
    assert tuple(table.bad_rows) == (False, True, False, False, True, False)
    stats = abundance_stats(table)
    assert stats.n_rows == 6
    # Bad cells are errors, not missing: only the NA row counts.
    assert stats.n_na == 1
    assert stats.n_zero == 1
    assert detect_scale(table.values) == "linear"


@pytest.mark.parametrize(
    "text",
    [
        "",  # empty file
        _tsv("protein_id|v"),  # header only
        _tsv("protein_id", "P1", "P2"),  # one column
        _tsv("protein_id|v", "P1|1|extra", "P2|2"),  # ragged
    ],
    ids=["empty", "header_only", "one_column", "ragged"],
)
def test_shape_errors_return_no_table(tmp_path: Path, text: str) -> None:
    table, issues = _read(tmp_path, text)
    assert table is None
    assert [i for i in issues if i.severity == "error"]


# --- 2. raw_ids vs resolved ids ------------------------------------------------------------------


def test_raw_ids_keep_written_form(tmp_path: Path) -> None:
    text = "protein_id\tv\nsp|P0A7V8|RS4_ECOLI\t1e6\nQ1\t2e6\n"
    table, _ = _read(tmp_path, text)
    assert table is not None
    assert table.raw_ids == ("sp|P0A7V8|RS4_ECOLI", "Q1")
    assert table.ids == ("P0A7V8", "Q1")
    assert len(table.raw_ids) == len(table.ids) == len(table.values) == len(table.bad_rows)


# --- 3. scale_undetermined only for requested auto ----------------------------------------------


def _all_na(tmp_path: Path) -> AbundanceTable:
    table, _ = _read(tmp_path, _tsv("protein_id|v", "P1|NA", "P2|"))
    assert table is not None
    return table


def test_scale_undetermined_when_auto_requested(tmp_path: Path) -> None:
    issues = validate_abundance(_all_na(tmp_path), "auto", WARN_FRAC)
    assert _issue(issues, "scale_undetermined").severity == "warn"


@pytest.mark.parametrize("requested", ["log2", "linear"])
def test_no_scale_undetermined_when_scale_given(tmp_path: Path, requested: str) -> None:
    issues = validate_abundance(_all_na(tmp_path), requested, WARN_FRAC)  # type: ignore[arg-type]
    assert "scale_undetermined" not in _codes(issues)


def test_auto_requested_resolves_scale_for_negative_check(tmp_path: Path) -> None:
    table, _ = _read(tmp_path, _tsv("protein_id|v", "P1|-1", "P2|100"))
    assert table is not None
    # auto → detected log2 (a negative is present): no negative_on_linear.
    assert "negative_on_linear" not in _codes(validate_abundance(table, "auto", WARN_FRAC))
    assert "negative_on_linear" in _codes(validate_abundance(table, "linear", WARN_FRAC))


# --- 4. protein groups keep written order --------------------------------------------------------


def test_protein_group_order_is_kept_and_not_duplicate(tmp_path: Path) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", "P2;P3|1e6", "P3;P2|2e6"))
    assert table is not None
    assert table.ids == ("P2;P3", "P3;P2")
    assert "duplicate_ids" not in _codes(issues)


# --- 5. duplicate / empty ID examples ------------------------------------------------------------


def test_duplicate_id_examples_name_later_row(tmp_path: Path) -> None:
    _, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1e6", "P1|2e6"))
    dup = _issue(issues, "duplicate_ids")
    assert dup.examples == ("2",)
    assert dup.n_rows == 1


def test_duplicate_id_examples_name_every_later_row(tmp_path: Path) -> None:
    _, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1", "P2|2", "P1|3", "P1|4", "P2|5"))
    dup = _issue(issues, "duplicate_ids")
    assert dup.examples == ("3", "4", "5")


def test_two_blank_ids_give_only_empty_id(tmp_path: Path) -> None:
    _, issues = _read(tmp_path, _tsv("protein_id|v", "|1e6", "P1|2e6", "|3e6"))
    empty = _issue(issues, "empty_id")
    assert empty.examples == ("1", "3")
    assert "duplicate_ids" not in _codes(issues)
