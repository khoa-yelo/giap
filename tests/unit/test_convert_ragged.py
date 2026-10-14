"""Ragged rows in a DIA-NN ``report.tsv`` are errors (GIAP_SPEC.md §5.4; lead decision).

A data row with fewer fields than the header (truncated file, corrupt row) is rejected exactly like
a row with too many fields: it is never padded with missing values, because a padded row would
carry a missing q-value, which means "no threshold applies". Empty fields with their tabs present
are well formed and stay accepted.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from pathlib import Path

import pytest
from click.testing import CliRunner

from giap.cli import cli
from giap.config.model import InputConfig
from giap.exitcodes import ExitCode
from giap.io.converters.diann import convert_diann

HEADER = [
    "Run",
    "Protein.Group",
    "Protein.Ids",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Q.Value",
    "Global.Q.Value",
    "PEP",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "PG.MaxLFQ",
]


def _row(run: str) -> list[str]:
    lead = [run, "P1;P2", "P1;P2", "PEPTIDEK", "PEPTIDEK", "2", "0.002", "0.004", "0.01"]
    return [*lead, "1400", "1500", "10"]


def _short_cut_q(cells: list[str]) -> list[str]:
    return cells[:6]  # truncated before Q.Value / Global.Q.Value


def _short_one(cells: list[str]) -> list[str]:
    return cells[:-1]  # only the last field missing


def _long(cells: list[str]) -> list[str]:
    return [*cells, "extra"]


RAGGED: dict[str, Callable[[list[str]], list[str]]] = {
    "short_cut_q": _short_cut_q,
    "short_one": _short_one,
    "long": _long,
}


def _write(tmp_path: Path, rows: list[list[str]]) -> Path:
    p = tmp_path / "report.tsv"
    lines = ["\t".join(HEADER)] + ["\t".join(r) for r in rows]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def _ragged_report(tmp_path: Path, kind: str) -> Path:
    """Three data rows (file lines 2-4); the one on line 3 is ragged."""
    return _write(tmp_path, [_row("r1"), RAGGED[kind](_row("r2")), _row("r3")])


def _errors(path: Path) -> list[str]:
    res = convert_diann([path], search_db=None, cfg=InputConfig())
    return [f"{i.code}: {i.message}" for i in res.issues if i.severity == "error"]


@pytest.mark.parametrize("kind", sorted(RAGGED))
def test_ragged_row_is_error(tmp_path: Path, kind: str) -> None:
    res = convert_diann([_ragged_report(tmp_path, kind)], search_db=None, cfg=InputConfig())
    codes = {i.code for i in res.issues if i.severity == "error"}
    assert codes & {"unreadable", "malformed_table"}, res.issues
    assert res.peptides.empty


@pytest.mark.parametrize("kind", sorted(RAGGED))
def test_ragged_error_names_line_number(tmp_path: Path, kind: str) -> None:
    errors = _errors(_ragged_report(tmp_path, kind))
    assert any(re.search(r"\bline 3\b", e) for e in errors), errors


@pytest.mark.parametrize("kind", sorted(RAGGED))
def test_cli_ragged_exits_1_writes_nothing(tmp_path: Path, kind: str) -> None:
    src = _ragged_report(tmp_path, kind)
    out = tmp_path / "peptides.tsv"
    ab = tmp_path / "abundance.tsv"
    argv = ["convert", "--from", "diann", str(src), "--out", str(out), "--abundance-out", str(ab)]
    result = CliRunner().invoke(cli, argv)
    assert result.exit_code == ExitCode.INVALID_INPUT, result.output
    assert not out.exists()
    assert not ab.exists()


def test_empty_trailing_fields_are_well_formed(tmp_path: Path) -> None:
    blank = _row("r2")
    blank[-2:] = ["", ""]  # Precursor.Normalised and PG.MaxLFQ empty, tabs present
    src = _write(tmp_path, [_row("r1"), blank, _row("r3")])
    assert src.read_text().splitlines()[2].endswith("\t\t")
    res = convert_diann([src], search_db=None, cfg=InputConfig())
    assert not [i for i in res.issues if i.severity == "error"], res.issues
    assert len(res.peptides) == 3
    assert res.peptides["q_value"].notna().all()
    assert res.peptides["intensity"].isna().tolist() == [False, True, False]  # column-level map

    out = tmp_path / "peptides.tsv"
    result = CliRunner().invoke(cli, ["convert", "--from", "diann", str(src), "--out", str(out)])
    assert result.exit_code == ExitCode.SUCCESS, result.output
    assert out.exists()
