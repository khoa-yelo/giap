"""Abundance table reader and scale detection (issue M1-08; GIAP_SPEC.md §5.3, §7.0; ADR-0005/6).

Tables are literal TSV text written to ``tmp_path`` (``|`` in the helpers stands for a tab).
``_load`` runs ``read_abundance`` and, when a table comes back, ``validate_abundance`` on the
auto-resolved scale, and returns the issues of both steps, so a check may live in either one.
"""

from __future__ import annotations

import dataclasses
import gzip
import math
from pathlib import Path

import pytest

from giap.config import resolve_config
from giap.io.abundance import (
    AbundanceStats,
    AbundanceTable,
    abundance_stats,
    detect_scale,
    read_abundance,
    resolve_scale,
    validate_abundance,
)
from giap.io.issues import Issue

WARN_FRAC = 0.20


def _tsv(*lines: str) -> str:
    return "\n".join(line.replace("|", "\t") for line in lines) + "\n"


def _write(tmp_path: Path, text: str, name: str = "abundance.tsv") -> Path:
    p = tmp_path / name
    data = text.encode("utf-8")
    if name.endswith(".gz"):
        data = gzip.compress(data)
    p.write_bytes(data)
    return p


def _read(
    tmp_path: Path, text: str, name: str = "abundance.tsv"
) -> tuple[AbundanceTable | None, list[Issue]]:
    return read_abundance(_write(tmp_path, text, name))


def _load(
    tmp_path: Path, text: str, requested: str = "auto", name: str = "abundance.tsv"
) -> tuple[AbundanceTable | None, list[Issue]]:
    table, issues = _read(tmp_path, text, name)
    if table is None:
        return None, issues
    scale = resolve_scale(table.values, requested)  # type: ignore[arg-type]
    return table, [*issues, *validate_abundance(table, scale, WARN_FRAC)]


def _table(values: list[float | None], ids: list[str] | None = None) -> AbundanceTable:
    ids = ids or [f"P{i}" for i in range(1, len(values) + 1)]
    return AbundanceTable(
        id_column="protein_id",
        value_column="v",
        ids=tuple(ids),
        values=tuple(values),
        ignored_columns=(),
    )


def _codes(issues: list[Issue]) -> set[str]:
    return {i.code for i in issues}


def _issue(issues: list[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert found, f"no issue {code!r} in {[i.code for i in issues]}"
    return found[0]


def _errors(issues: list[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


# --- interface ----------------------------------------------------------------------------------


def test_dataclasses_are_frozen_with_spec_fields() -> None:
    assert [f.name for f in dataclasses.fields(AbundanceTable)] == [
        "id_column",
        "value_column",
        "ids",
        "values",
        "ignored_columns",
    ]
    assert [f.name for f in dataclasses.fields(AbundanceStats)] == [
        "n_rows",
        "n_na",
        "n_zero",
        "n_negative",
        "missing_frac",
    ]
    t = _table([1.0])
    with pytest.raises(dataclasses.FrozenInstanceError):
        t.ids = ()  # type: ignore[misc]


def test_issues_use_shared_issue_type(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("protein_id|v|w", "P1|3|x"))
    assert issues
    assert all(isinstance(i, Issue) for i in issues)


def test_config_defaults_for_abundance() -> None:
    cfg = resolve_config()
    assert cfg.abundance.scale == "auto"
    assert cfg.preflight.abundance_missing_warn_frac == pytest.approx(0.20)


# --- spec examples 1-10 -------------------------------------------------------------------------


def test_example_1_linear_ibaq_table_ok(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|ibaq", "P1|1e6", "P2|2.5e5"))
    assert table is not None
    assert table.id_column == "protein_id"
    assert table.value_column == "ibaq"
    assert table.ids == ("P1", "P2")
    assert table.values == (1e6, 2.5e5)
    assert table.ignored_columns == ()
    assert detect_scale(table.values) == "linear"
    assert issues == []


def test_example_2_negative_value_means_log2(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|log2_abund", "P1|21.3", "P2|-0.5"))
    assert table is not None
    assert table.values == (21.3, -0.5)
    assert detect_scale(table.values) == "log2"
    assert _errors(issues) == []


def test_example_3_all_values_at_most_50_is_log2() -> None:
    assert detect_scale([12.1, 30.0, 50]) == "log2"


def test_example_4_value_above_50_is_linear() -> None:
    assert detect_scale([12.1, 50.01]) == "linear"


def test_example_5_extra_columns_warn_and_are_ignored(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|v|w", "P1|3|x"))
    assert table is not None
    assert table.value_column == "v"
    assert table.values == (3.0,)
    assert table.ignored_columns == ("w",)
    extra = _issue(issues, "extra_columns")
    assert extra.severity == "warn"
    assert "w" in extra.message
    assert _errors(issues) == []


def test_example_6_non_numeric_value_is_error_with_example_row(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("protein_id|v", "P1|abc"))
    bad = _issue(issues, "non_numeric_values")
    assert bad.severity == "error"
    assert bad.examples == ("1",)


def test_example_7_negative_on_requested_linear_is_error() -> None:
    t = _table([-1.0, 100.0])
    scale = resolve_scale(t.values, "linear")
    assert scale == "linear"
    bad = _issue(validate_abundance(t, scale, WARN_FRAC), "negative_on_linear")
    assert bad.severity == "error"


def test_example_8_negative_with_auto_is_log2_without_error() -> None:
    t = _table([-1.0, 100.0])
    scale = resolve_scale(t.values, "auto")
    assert scale == "log2"
    issues = validate_abundance(t, scale, WARN_FRAC)
    assert "negative_on_linear" not in _codes(issues)
    assert _errors(issues) == []


def test_example_9_missing_and_zero_counts_warn_above_20_percent(tmp_path: Path) -> None:
    rows = ["P1|NA", "P2|NA", "P3|NA", "P4|0"] + [f"P{i}|{1000 * i}" for i in range(5, 11)]
    table, issues = _load(tmp_path, _tsv("protein_id|ibaq", *rows))
    assert table is not None
    stats = abundance_stats(table)
    assert stats.n_rows == 10
    assert stats.n_na == 3
    assert stats.n_zero == 1
    assert stats.n_negative == 0
    assert stats.missing_frac == pytest.approx(0.3)
    warn = _issue(issues, "too_many_missing")
    assert warn.severity == "warn"


def test_example_10_exactly_20_percent_missing_does_not_warn(tmp_path: Path) -> None:
    rows = ["P1|NA", "P2|NA"] + [f"P{i}|{1000 * i}" for i in range(3, 11)]
    table, issues = _load(tmp_path, _tsv("protein_id|ibaq", *rows))
    assert table is not None
    assert abundance_stats(table).missing_frac == pytest.approx(0.2)
    assert "too_many_missing" not in _codes(issues)


# --- detect_scale / resolve_scale ---------------------------------------------------------------


def test_detect_scale_boundary_exactly_50_is_log2() -> None:
    assert detect_scale([50.0]) == "log2"
    assert detect_scale([1.0, 50.0]) == "log2"


def test_detect_scale_boundary_just_above_50_is_linear() -> None:
    assert detect_scale([1.0, 50.000001]) == "linear"


def test_detect_scale_negative_wins_over_large_values() -> None:
    assert detect_scale([1e9, -0.001]) == "log2"


def test_detect_scale_ignores_na() -> None:
    assert detect_scale([None, 10.0, None]) == "log2"
    assert detect_scale([None, 1e6]) == "linear"


def test_detect_scale_all_na_is_linear() -> None:
    assert detect_scale([None, None]) == "linear"


def test_all_na_table_reports_scale_undetermined_warn(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|v", "P1|NA", "P2|"))
    assert table is not None
    assert table.values == (None, None)
    und = _issue(issues, "scale_undetermined")
    assert und.severity == "warn"


@pytest.mark.parametrize("requested", ["log2", "linear"])
def test_resolve_scale_explicit_overrides_detection(requested: str) -> None:
    assert resolve_scale([1e6, 2e6], requested) == requested  # type: ignore[arg-type]
    assert resolve_scale([1.0, 2.0], requested) == requested  # type: ignore[arg-type]


def test_resolve_scale_auto_delegates_to_detection() -> None:
    assert resolve_scale([1e6, 2e6], "auto") == "linear"
    assert resolve_scale([1.0, 2.0], "auto") == "log2"


# --- validate_abundance / stats -----------------------------------------------------------------


def test_zero_on_log2_is_valid() -> None:
    t = _table([0.0, 12.0, -1.0])
    issues = validate_abundance(t, "log2", WARN_FRAC)
    assert _errors(issues) == []
    assert "too_many_missing" not in _codes(issues)


def test_linear_zeros_counted_but_not_missing() -> None:
    t = _table([0.0, 0.0, 0.0, 1e6])
    stats = abundance_stats(t)
    assert stats.n_zero == 3
    assert stats.n_na == 0
    assert stats.missing_frac == 0.0
    issues = validate_abundance(t, "linear", WARN_FRAC)
    assert "too_many_missing" not in _codes(issues)
    assert _errors(issues) == []


def test_stats_count_negatives() -> None:
    stats = abundance_stats(_table([-1.0, -2.0, 3.0, None]))
    assert stats.n_negative == 2
    assert stats.n_na == 1
    assert stats.missing_frac == pytest.approx(0.25)


def test_negative_on_linear_reports_rows() -> None:
    bad = _issue(
        validate_abundance(_table([1e6, -5.0, 2e6, -1.0]), "linear", WARN_FRAC),
        "negative_on_linear",
    )
    assert bad.n_rows == 2
    assert bad.examples == ("2", "4")


def test_missing_warn_frac_threshold_is_respected() -> None:
    t = _table([None, 1.0, 2.0, 3.0])  # 25 % missing
    assert "too_many_missing" in _codes(validate_abundance(t, "log2", 0.20))
    assert "too_many_missing" not in _codes(validate_abundance(t, "log2", 0.30))


# --- reader edge cases --------------------------------------------------------------------------


def test_header_only_is_empty_table_error(tmp_path: Path) -> None:
    _, issues = _read(tmp_path, _tsv("protein_id|v"))
    err = _issue(issues, "empty_table")
    assert err.severity == "error"


def test_single_column_is_no_value_column_error(tmp_path: Path) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id", "P1", "P2"))
    assert table is None
    err = _issue(issues, "no_value_column")
    assert err.severity == "error"


def test_empty_file_is_error_not_exception(tmp_path: Path) -> None:
    table, issues = _read(tmp_path, "")
    assert table is None
    assert _errors(issues)


def test_duplicate_ids_are_error(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("protein_id|v", "P1|1e6", "P2|2e6", "P1|3e6"))
    err = _issue(issues, "duplicate_ids")
    assert err.severity == "error"
    assert "3" in err.examples


def test_empty_id_is_error(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("protein_id|v", "P1|1e6", "|2e6"))
    err = _issue(issues, "empty_id")
    assert err.severity == "error"
    assert err.examples == ("2",)


def test_scientific_notation_parses(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|v", "P1|1.5E+06", "P2|3e-2"))
    assert table is not None
    assert table.values == (1.5e6, 3e-2)
    assert _errors(issues) == []


@pytest.mark.parametrize("text", ["inf", "INF", "-inf", "-INF", "nan", "NaN"])
def test_non_finite_text_is_error_not_float_or_missing(tmp_path: Path, text: str) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1e6", f"P2|{text}", "P3|2e6"))
    bad = [i for i in issues if i.code in {"bad_value", "non_finite"}]
    assert bad, f"no bad_value/non_finite issue in {[i.code for i in issues]}"
    assert bad[0].severity == "error"
    assert bad[0].examples == ("2",)
    if table is not None:
        # Neither kept as a non-finite float nor counted as missing (NA).
        assert all(v is None or math.isfinite(v) for v in table.values)
        assert abundance_stats(table).n_na == 0


@pytest.mark.parametrize("missing", ["NA", ""])
def test_na_and_empty_are_missing(tmp_path: Path, missing: str) -> None:
    table, issues = _read(tmp_path, _tsv("protein_id|v", f"P1|{missing}", "P2|1e6"))
    assert table is not None
    assert table.values == (None, 1e6)
    assert _errors(issues) == []


def test_non_numeric_examples_list_every_bad_row(tmp_path: Path) -> None:
    _, issues = _read(tmp_path, _tsv("protein_id|v", "P1|1", "P2|x", "P3|2", "P4|y"))
    bad = _issue(issues, "non_numeric_values")
    assert bad.n_rows == 2
    assert bad.examples == ("2", "4")


def test_uniprot_ids_resolve_like_proteins_faa(tmp_path: Path) -> None:
    table, issues = _load(tmp_path, _tsv("protein_id|v", "sp|P0A7V8|RS4_ECOLI|1e6", "tr|Q1|X|2e6"))
    assert table is not None
    assert table.ids == ("P0A7V8", "Q1")
    assert _errors(issues) == []


def test_ids_after_resolution_collide_as_duplicates(tmp_path: Path) -> None:
    _, issues = _load(tmp_path, _tsv("protein_id|v", "sp|P1|A_ECOLI|1e6", "P1|2e6"))
    assert "duplicate_ids" in _codes(issues)


def test_header_row_is_required_first_line_never_data(tmp_path: Path) -> None:
    # A file whose first line looks like data is still read as a header.
    table, _ = _read(tmp_path, _tsv("P1|1e6", "P2|2e6"))
    assert table is not None
    assert table.id_column == "P1"
    assert table.value_column == "1e6"
    assert table.ids == ("P2",)


def test_bom_and_crlf_are_tolerated(tmp_path: Path) -> None:
    p = tmp_path / "abundance.tsv"
    p.write_bytes(b"\xef\xbb\xbfprotein_id\tibaq\r\nP1\t1e6\r\nP2\t2.5e5\r\n")
    table, issues = read_abundance(p)
    assert table is not None
    assert table.id_column == "protein_id"
    assert table.value_column == "ibaq"
    assert table.ids == ("P1", "P2")
    assert table.values == (1e6, 2.5e5)
    assert issues == []


def test_gzipped_input_reads_identically(tmp_path: Path) -> None:
    text = _tsv("protein_id|ibaq|extra", "P1|1e6|a", "P2|NA|b", "P3|0|c")
    plain = _read(tmp_path, text, "abundance.tsv")
    gz = _read(tmp_path, text, "abundance.tsv.gz")
    assert plain[0] is not None
    assert gz == plain
