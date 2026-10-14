"""Generic converter: validate, normalise, report (issue M1-07; GIAP_SPEC.md §5.3, §5.4)."""

from __future__ import annotations

from pathlib import Path

import pytest

from giap.config.model import InputConfig
from giap.io.converters import CONVERTERS
from giap.io.converters.base import ConvertResult
from giap.io.converters.generic import convert_generic
from giap.io.issues import Issue
from giap.io.peptides import PEPTIDE_COLUMNS


def _write(tmp_path: Path, text: str, name: str = "in.tsv") -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _run(
    paths: list[Path], source: str | None = None, search_db: str | None = None
) -> ConvertResult:
    return convert_generic(paths, source=source, search_db=search_db, cfg=InputConfig())


def _issue(issues: list[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert found, f"no issue {code!r} in {[i.code for i in issues]}"
    return found[0]


def test_registered() -> None:
    assert CONVERTERS["generic"] is convert_generic


def test_example9_strip_and_canonical_bool(tmp_path: Path) -> None:
    p = _write(tmp_path, "peptide\tis_decoy\nPEPT[+79.966]IDEK\t1\nAAAAAAAK\t0\n")
    res = _run([p])
    assert not [i for i in res.issues if i.severity == "error"]
    df = res.peptides.set_index("peptide")
    assert set(df.index) == {"PEPTIDEK", "AAAAAAAK"}
    assert bool(df.loc["PEPTIDEK", "is_decoy"]) is True
    assert bool(df.loc["AAAAAAAK", "is_decoy"]) is False
    issue = _issue(res.issues, "stripped")
    assert issue.severity == "info"
    assert issue.n_rows == 1 or "1 row" in issue.message


def test_no_stripped_issue_for_clean_input(tmp_path: Path) -> None:
    res = _run([_write(tmp_path, "peptide\nPEPTIDEK\n")])
    assert "stripped" not in {i.code for i in res.issues}
    assert res.abundance is None


def test_unknown_column_dropped_with_warning(tmp_path: Path) -> None:
    p = _write(tmp_path, "peptide\tfoo\tq_value\nPEPTIDEK\tx\t0.001\n")
    res = _run([p])
    assert "foo" not in res.peptides.columns
    assert list(res.peptides.columns) == ["peptide", "q_value"]
    assert _issue(res.issues, "unknown_column").severity == "warn"


def test_columns_in_spec_order(tmp_path: Path) -> None:
    p = _write(tmp_path, "run\tq_value\tpeptide\nr1\t0.001\tPEPTIDEK\n")
    cols = list(_run([p]).peptides.columns)
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]


def test_source_and_search_db_flags_fill_columns(tmp_path: Path) -> None:
    p = _write(tmp_path, "peptide\nPEPTIDEK\n")
    res = _run([p], source="denovo", search_db="other")
    assert set(res.peptides["source"]) == {"denovo"}
    assert set(res.peptides["search_db"]) == {"other"}


def test_error_issue_is_reported(tmp_path: Path) -> None:
    p = _write(tmp_path, "peptide\tq_value\nPEPTIDEK\tabc\n")
    res = _run([p])
    assert [i for i in res.issues if i.severity == "error"]


def test_several_inputs_concatenated_and_sorted(tmp_path: Path) -> None:
    a = _write(tmp_path, "peptide\trun\nPEPTIDEK\tr2\n", "a.tsv")
    b = _write(tmp_path, "peptide\trun\nAAAAAAAK\tr1\n", "b.tsv")
    res = _run([a, b, a])
    assert res.peptides["run"].tolist() == ["r1", "r2", "r2"]
    assert res.peptides["peptide"].tolist() == ["AAAAAAAK", "PEPTIDEK", "PEPTIDEK"]


@pytest.mark.parametrize("text", ["TRUE", "True", "true", "1"])
def test_boolean_spellings_canonicalised(tmp_path: Path, text: str) -> None:
    p = _write(tmp_path, f"peptide\tis_decoy\nPEPTIDEK\t{text}\n")
    assert bool(_run([p]).peptides["is_decoy"].iloc[0]) is True
