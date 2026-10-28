"""Cascadia 0.0.7 ``.ssl`` converter (M10-09; GIPPA_SPEC.md §5.3, §5.4).

Spec: docs/specs/M10-09-cascadia-convert.md (facts W1-W9, decisions D1-D7, examples 1-11 and 13,
acceptance criteria 2-7; lead decisions 2026-10-24: Q1-Q13 take their defaults).

Every ``.ssl`` fixture is hand-made text in the 0.0.7 shape (W1-W7: rounded mass shifts,
``UNKNOWN`` score-type, absolute ``file`` paths, gapped ``scan`` values) written into
``tmp_path``; no ``.ssl`` file is committed (criterion 9).
"""

from __future__ import annotations

import ast
import gzip
import itertools
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.config.model import GippaConfig, InputConfig
from gippa.io.converters.base import ConvertResult
from gippa.io.converters.strip import strip_modifications
from gippa.io.issues import Issue
from gippa.io.peptides import PEPTIDE_COLUMNS, read_peptide_table
from gippa.io.tsv import write_tsv

REPO = Path(__file__).resolve().parents[2]
MODULE = REPO / "src" / "gippa" / "io" / "converters" / "cascadia.py"

# W1: the header 0.0.7 writes, literally.
W1 = (
    "file",
    "scan",
    "charge",
    "sequence",
    "score-type",
    "score",
    "retention-time",
    "start-time",
    "end-time",
)
FILE_WT1 = "/x/mzml/WT_1.mzML"
OUT_HEADER = "peptide\tmodified_peptide\tscore\trun\tcharge\tspectrum_id\tsource\n"
OUT_COLUMNS = ["peptide", "modified_peptide", "score", "run", "charge", "spectrum_id", "source"]
DEFAULTS: dict[str, str] = {
    "file": FILE_WT1,
    "score-type": "UNKNOWN",
    "retention-time": "4.5",
    "start-time": "4.4",
    "end-time": "4.6",
}
# The W5 mass-shift forms of 0.0.7 (and one unresolved name, which 0.0.7 leaves as the name).
EXAMPLE2 = {
    "[+42.01]-SAMPLLER": "SAMPLLER",
    "[-17.03]-QLLDPEKR": "QLLDPEKR",
    "AC[+57.02]DN[+0.98]LLK": "ACDNLLK",
    "[+25.980265]-GLLDPEKR": "GLLDPEKR",
    "Q[+0.98]VLLDPEK": "QVLLDPEK",
    "M[Oxidation]VLLDPEK": "MVLLDPEK",
}
ABSENT = ("q_value", "pep", "protein_ids", "is_decoy", "intensity", "search_db")


# ------------------------------------------------------------------ helpers


def ssl_text(
    rows: Sequence[Mapping[str, str | None]],
    header: Sequence[str] = W1,
    *,
    newline: str = "\n",
) -> str:
    """The 0.0.7 text: the header, then one tab-joined line per row (defaults per column)."""
    lines = ["\t".join(header)]
    for row in rows:
        cells = []
        for col in header:
            value = row[col] if col in row else DEFAULTS.get(col, "")
            cells.append("" if value is None else value)
        lines.append("\t".join(cells))
    return newline.join(lines) + newline


def write_ssl(
    tmp_path: Path,
    rows: Sequence[Mapping[str, str | None]],
    header: Sequence[str] = W1,
    *,
    name: str = "WT_1.ssl",
) -> Path:
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(ssl_text(rows, header), encoding="utf-8")
    return p


def row(
    scan: str, charge: str, sequence: str, score: str, **extra: str | None
) -> dict[str, str | None]:
    """One data row; ``score_type=...`` sets the ``score-type`` cell, ``None`` an empty cell."""
    out: dict[str, str | None] = {
        "scan": scan,
        "charge": charge,
        "sequence": sequence,
        "score": score,
    }
    out.update({k.replace("_", "-"): v for k, v in extra.items()})
    return out


def convert(paths: Sequence[Path], search_db: str | None = None) -> ConvertResult:
    from gippa.io.converters.cascadia import convert_cascadia

    return convert_cascadia(list(paths), search_db=search_db, cfg=InputConfig())


def to_bytes(res: ConvertResult, out: Path) -> bytes:
    """The bytes the CLI would write for ``res`` (``write_tsv``, §6.6)."""
    df = res.peptides
    write_tsv(out, list(df.columns), df.itertuples(index=False, name=None))
    return out.read_bytes()


def errors(issues: Sequence[Issue]) -> list[Issue]:
    return [i for i in issues if i.severity == "error"]


def by_code(issues: Sequence[Issue], code: str) -> Issue:
    found = [i for i in issues if i.code == code]
    assert len(found) == 1, (code, issues)
    return found[0]


def codes(issues: Sequence[Issue]) -> set[str]:
    return {i.code for i in issues}


EX1_ROWS = [
    row("0", "2", "PEPTLDEK", "0.95123"),
    row("3", "3", "M[+15.99]PEPTLDEKR", "0.88120001"),
]
EX1_OUT = (
    OUT_HEADER + "MPEPTLDEKR\tM[+15.99]PEPTLDEKR\t0.8812\tWT_1\t3\tcascadia:3\tdenovo\n"
    "PEPTLDEK\tPEPTLDEK\t0.95123\tWT_1\t2\tcascadia:0\tdenovo\n"
)


# ------------------------------------------------------------------ criterion 1 (registry)


def test_registered_in_converters() -> None:
    import gippa.io.converters as conv
    from gippa.io.converters.cascadia import convert_cascadia

    assert conv.CONVERTERS["cascadia"] is convert_cascadia
    assert "convert_cascadia" in conv.__all__
    assert conv.convert_cascadia is convert_cascadia


def test_default_cutoff_constant_is_named() -> None:
    from gippa.io.converters import cascadia

    assert cascadia.CASCADIA_DEFAULT_CUTOFF == 0.8


# ------------------------------------------------------------------ examples 1-3


def test_example1_mapping_and_bytes(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, EX1_ROWS)])
    assert not errors(res.issues), res.issues
    assert res.abundance is None
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT
    stripped = by_code(res.issues, "stripped")
    assert stripped.severity == "info"
    assert stripped.n_rows == 1
    assert stripped.examples == ("2",)
    # Example 1: the issue list is the single info ``stripped``. A standard W1 header (RT
    # columns included) gives no ``unknown_column`` warning.
    assert [(i.code, i.severity) for i in res.issues] == [("stripped", "info")]


@pytest.mark.parametrize(("sequence", "peptide"), list(EXAMPLE2.items()))
def test_example2_modifications_strip(tmp_path: Path, sequence: str, peptide: str) -> None:
    res = convert([write_ssl(tmp_path, [row("1", "2", sequence, "0.9")])])
    assert not errors(res.issues), res.issues
    assert res.peptides["peptide"].tolist() == [peptide]
    assert res.peptides["modified_peptide"].tolist() == [sequence]


def test_example2_all_rows_in_one_file(tmp_path: Path) -> None:
    rows = [row(str(2 * i + 1), "2", s, "0.9") for i, s in enumerate(EXAMPLE2)]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    got = dict(zip(res.peptides["modified_peptide"], res.peptides["peptide"], strict=True))
    assert got == EXAMPLE2
    assert by_code(res.issues, "stripped").n_rows == len(EXAMPLE2)


def test_sequence_that_strips_to_empty_becomes_missing(tmp_path: Path) -> None:
    """Step 8: an empty stripped result is missing, not ``""``; the validator warns, row kept."""
    rows = [row("0", "2", "PEPTLDEK", "0.9"), row("1", "2", "[+42.01]-", "0.9")]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert len(res.peptides) == 2
    df = res.peptides.set_index("spectrum_id")
    assert pd.isna(df.loc["cascadia:1", "peptide"])
    assert df.loc["cascadia:1", "modified_peptide"] == "[+42.01]-"
    assert "" not in res.peptides["peptide"].dropna().tolist()
    assert by_code(res.issues, "stripped").examples == ("2",)
    warn = by_code(res.issues, "non_letters")
    assert warn.severity == "warn"
    assert warn.examples == ("2",)
    line = [x for x in to_bytes(res, tmp_path / "p.tsv").decode().splitlines() if "cascadia:1" in x]
    assert line[0].split("\t")[:2] == ["NA", "[+42.01]-"]


def test_example3_i_and_l_kept_as_written(tmp_path: Path) -> None:
    rows = [row("4", "2", "PEPTIDEK", "0.9"), row("5", "2", "PEPTLDEK", "0.85")]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert res.peptides["peptide"].tolist() == ["PEPTIDEK", "PEPTLDEK"]
    assert res.peptides["modified_peptide"].tolist() == ["PEPTIDEK", "PEPTLDEK"]
    assert res.peptides["spectrum_id"].tolist() == ["cascadia:4", "cascadia:5"]
    assert "stripped" not in codes(res.issues)


# ------------------------------------------------------------------ examples 4-6 (score, type)


def test_example4_low_scores_kept_and_counted(tmp_path: Path) -> None:
    rows = [
        row("0", "2", "PEPTLDEKA", "0.8"),
        row("1", "2", "PEPTLDEKC", "0.5"),
        row("2", "2", "PEPTLDEKD", "0.80000001"),
        row("3", "2", "PEPTLDEKE", "1.0"),
    ]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert len(res.peptides) == 4
    low = by_code(res.issues, "score_at_or_below_default_cutoff")
    assert low.severity == "info"
    assert low.n_rows == 2
    assert low.examples == ("1", "2")
    text = to_bytes(res, tmp_path / "p.tsv").decode().splitlines()[1:]
    assert [line.split("\t")[2] for line in text] == ["0.8", "0.5", "0.8", "1"]


def test_no_low_score_issue_when_all_above_cutoff(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, EX1_ROWS)])
    assert "score_at_or_below_default_cutoff" not in codes(res.issues)


def test_example5_bad_and_out_of_range_scores(tmp_path: Path) -> None:
    cells = ["nan", "NA", "abc", "", "inf", "1.2", "-0.1", "0.9"]
    rows = [row(str(i), "2", "PEPTLDEK" + "ACDEFGHK"[i], s) for i, s in enumerate(cells)]
    res = convert([write_ssl(tmp_path, rows)])
    assert res.peptides.empty
    assert res.abundance is None
    bad = by_code(res.issues, "bad_type:score")
    assert bad.severity == "error"
    assert (bad.n_rows, bad.examples) == (4, ("1", "2", "3", "4"))
    out = by_code(res.issues, "out_of_range:score")
    assert out.severity == "error"
    assert (out.n_rows, out.examples) == (3, ("5", "6", "7"))


@pytest.mark.parametrize("cell", ["-inf", "+nan", "NaN", "1.0000001", "-0"])
def test_more_score_cells(tmp_path: Path, cell: str) -> None:
    """``-0`` is 0, inside [0, 1]: valid; the others are errors (D5)."""
    res = convert([write_ssl(tmp_path, [row("0", "2", "PEPTLDEK", cell)])])
    if cell == "-0":
        assert not errors(res.issues), res.issues
        assert "score_at_or_below_default_cutoff" in codes(res.issues)
    else:
        assert errors(res.issues)
        assert codes(errors(res.issues)) <= {"bad_type:score", "out_of_range:score"}


def test_example6_unsupported_score_type(tmp_path: Path) -> None:
    rows = [
        row("0", "2", "PEPTLDEK", "0.9", score_type="PERCOLATOR QVALUE"),
        row("1", "2", "SAMPLLER", "0.9"),
    ]
    res = convert([write_ssl(tmp_path, rows)])
    assert res.peptides.empty
    issue = by_code(res.issues, "unsupported_score_type")
    assert issue.severity == "error"
    assert issue.examples == ("1",)
    assert "PERCOLATOR QVALUE" in issue.message


def test_unsupported_score_type_message_names_first_offending_value(tmp_path: Path) -> None:
    """Step 5: the message quotes the first offending value, not the first row's value."""
    rows = [
        row("0", "2", "PEPTLDEK", "0.9"),
        row("1", "2", "SAMPLLER", "0.9", score_type="XCORR"),
        row("2", "2", "SAMPLLERK", "0.9", score_type="PERCOLATOR QVALUE"),
    ]
    res = convert([write_ssl(tmp_path, rows)])
    assert res.peptides.empty
    issue = by_code(res.issues, "unsupported_score_type")
    assert issue.severity == "error"
    assert (issue.n_rows, issue.examples) == (2, ("2", "3"))
    assert "'XCORR'" in issue.message
    assert "PERCOLATOR" not in issue.message


def test_missing_score_type_is_the_same_error(tmp_path: Path) -> None:
    rows = [row("0", "2", "PEPTLDEK", "0.9"), row("1", "2", "SAMPLLER", "0.9", score_type="NA")]
    res = convert([write_ssl(tmp_path, rows)])
    issue = by_code(res.issues, "unsupported_score_type")
    assert issue.severity == "error"
    assert issue.examples == ("2",)


# ------------------------------------------------------------------ example 7 (duplicates)


def test_example7a_duplicate_spectrum_row_is_a_warning(tmp_path: Path) -> None:
    rows = [row("0", "2", "PEPTLDEK", "0.9"), row("0", "2", "PEPTLDEK", "0.85")]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert len(res.peptides) == 2
    assert by_code(res.issues, "duplicate_spectrum_row").severity == "warn"


def test_example7b_duplicate_spectrum_id_is_an_error(tmp_path: Path) -> None:
    rows = [row("0", "2", "PEPTLDEK", "0.9"), row("0", "2", "SAMPLLER", "0.85")]
    res = convert([write_ssl(tmp_path, rows)])
    assert by_code(res.issues, "duplicate_spectrum_id").severity == "error"


# ------------------------------------------------------------------ examples 8-9 (runs)

EX8_OUT = (
    OUT_HEADER + "PEPTLDEK\tPEPTLDEK\t0.95\tWT_1\t2\tcascadia:0\tdenovo\n"
    "PEPTLDEK\tPEPTLDEK\t0.9\tWT_2\t2\tcascadia:0\tdenovo\n"
)


def _ex8(tmp_path: Path) -> tuple[Path, Path]:
    wt2 = write_ssl(
        tmp_path, [row("0", "2", "PEPTLDEK", "0.9", file="/x/mzml/WT_2.mzML")], name="WT_2.ssl"
    )
    wt1 = write_ssl(tmp_path, [row("0", "2", "PEPTLDEK", "0.95")], name="WT_1.ssl")
    return wt2, wt1


def test_example8_two_runs_sorted_by_run(tmp_path: Path) -> None:
    wt2, wt1 = _ex8(tmp_path)
    res = convert([wt2, wt1])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX8_OUT


def test_argument_order_does_not_change_bytes(tmp_path: Path) -> None:
    """Criterion 3: distinct runs in a different argument order give the same bytes."""
    wt2, wt1 = _ex8(tmp_path)
    a = to_bytes(convert([wt2, wt1]), tmp_path / "a.tsv")
    b = to_bytes(convert([wt1, wt2]), tmp_path / "b.tsv")
    assert a == b


def test_example9_same_stem_in_two_directories(tmp_path: Path) -> None:
    a = write_ssl(tmp_path / "a", [row("0", "2", "PEPTLDEK", "0.9", file="/a/WT_1.mzML")])
    b = write_ssl(tmp_path / "b", [row("1", "2", "SAMPLLER", "0.9", file="/b/WT_1.mzML")])
    res = convert([a, b])
    assert res.peptides.empty
    issue = by_code(res.issues, "run_in_several_inputs")
    assert issue.severity == "error"
    assert str(a) in issue.message
    assert str(b) in issue.message


def test_example9_same_file_twice(tmp_path: Path) -> None:
    p = write_ssl(tmp_path, EX1_ROWS)
    res = convert([p, p])
    assert res.peptides.empty
    issue = by_code(res.issues, "run_in_several_inputs")
    assert issue.severity == "error"
    assert str(p) in issue.message


# ------------------------------------------------------------------ example 10 (empty, unreadable)


def test_example10a_header_only_alone(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, [])])
    assert by_code(res.issues, "empty_table").severity == "error"
    empty = by_code(res.issues, "empty_input")
    assert empty.severity == "info"
    assert "WT_1.ssl" in empty.message


def test_example10b_header_only_among_others(tmp_path: Path) -> None:
    wt1 = write_ssl(tmp_path, [], name="WT_1.ssl")
    wt2 = write_ssl(
        tmp_path, [row("0", "2", "PEPTLDEK", "0.9", file="/x/mzml/WT_2.mzML")], name="WT_2.ssl"
    )
    res = convert([wt1, wt2])
    assert not errors(res.issues), res.issues
    assert res.peptides["run"].tolist() == ["WT_2"]
    empty = by_code(res.issues, "empty_input")
    assert empty.severity == "info"
    assert "WT_1.ssl" in empty.message
    assert "WT_2.ssl" not in empty.message


def test_example10c_zero_byte_file(tmp_path: Path) -> None:
    p = tmp_path / "WT_1.ssl"
    p.write_bytes(b"")
    res = convert([p])
    issue = by_code(res.issues, "unreadable")
    assert issue.severity == "error"
    assert "empty TSV (no header)" in issue.message
    assert str(p) in issue.message


@pytest.mark.parametrize("missing", ["file", "scan", "charge", "sequence", "score-type", "score"])
def test_example10d_missing_required_column(tmp_path: Path, missing: str) -> None:
    header = [c for c in W1 if c != missing]
    res = convert([write_ssl(tmp_path, EX1_ROWS, header)])
    assert res.peptides.empty
    assert by_code(res.issues, f"missing_column:{missing}").severity == "error"


def test_ragged_row_is_unreadable_with_line_number(tmp_path: Path) -> None:
    p = tmp_path / "WT_1.ssl"
    good = ssl_text(EX1_ROWS)
    p.write_text(good + "x\t1\t2\n", encoding="utf-8")
    res = convert([p])
    issue = by_code(res.issues, "unreadable")
    assert issue.severity == "error"
    assert "line 4" in issue.message


def test_bad_encoding_is_unreadable(tmp_path: Path) -> None:
    p = tmp_path / "WT_1.ssl"
    p.write_bytes(ssl_text(EX1_ROWS).encode() + b"\xff\xfe\t\n")
    res = convert([p])
    issue = by_code(res.issues, "unreadable")
    assert issue.severity == "error"
    assert str(p) in issue.message


def test_missing_input_file_is_unreadable(tmp_path: Path) -> None:
    res = convert([tmp_path / "absent.ssl"])
    issue = by_code(res.issues, "unreadable")
    assert issue.severity == "error"
    assert str(tmp_path / "absent.ssl") in issue.message


# ------------------------------------------------------------------ example 11 (charge, length)


def test_example11_charge_zero(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, [row("0", "0", "PEPTLDEK", "0.9")])])
    assert by_code(res.issues, "out_of_range:charge").severity == "error"


def test_example11_charge_not_integral(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, [row("0", "2.5", "PEPTLDEK", "0.9")])])
    assert by_code(res.issues, "bad_type:charge").severity == "error"
    assert res.peptides.empty


def test_example11_charge_two_point_zero(tmp_path: Path) -> None:
    res = convert([write_ssl(tmp_path, [row("0", "2.0", "PEPTLDEK", "0.9")])])
    assert not errors(res.issues), res.issues
    line = to_bytes(res, tmp_path / "p.tsv").decode().splitlines()[1]
    assert line.split("\t")[4] == "2"


@pytest.mark.parametrize("cell", [None, "NA"])
def test_f7_missing_charge_is_left_to_the_validator(tmp_path: Path, cell: str | None) -> None:
    """Lead decision 2026-10-24 (D5): a missing charge converts to NA with no converter issue."""
    rows = [row("0", "2", "PEPTLDEK", "0.9"), {**row("1", "2", "SAMPLLER", "0.9"), "charge": cell}]
    res = convert([write_ssl(tmp_path, rows)])
    assert list(res.issues) == [], res.issues
    df = res.peptides.set_index("spectrum_id")
    assert len(df) == 2
    assert pd.isna(df.loc["cascadia:1", "charge"])
    assert df.loc["cascadia:0", "charge"] == 2
    assert str(res.peptides["charge"].dtype) == "Int64"


@pytest.mark.parametrize("cell", [None, "NA"])
def test_f7_missing_sequence_warns_and_row_is_written(tmp_path: Path, cell: str | None) -> None:
    """Lead decision 2026-10-24 (D5): a missing sequence gives only the warn ``non_letters``."""
    rows = [row("0", "2", "PEPTLDEK", "0.9"), {**row("1", "2", "X", "0.9"), "sequence": cell}]
    res = convert([write_ssl(tmp_path, rows)])
    assert [(i.code, i.severity, i.examples) for i in res.issues] == [
        ("non_letters", "warn", ("2",))
    ]
    df = res.peptides.set_index("spectrum_id")
    assert len(df) == 2
    assert pd.isna(df.loc["cascadia:1", "peptide"])
    assert pd.isna(df.loc["cascadia:1", "modified_peptide"])


def test_example11_short_and_long_peptides_written_with_warnings(tmp_path: Path) -> None:
    rows = [row("0", "2", "PEPK", "0.9"), row("1", "2", "A" * 41, "0.9")]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert sorted(res.peptides["peptide"]) == ["A" * 41, "PEPK"]
    assert by_code(res.issues, "short_peptide").severity == "warn"
    assert by_code(res.issues, "long_peptide").severity == "warn"


# ------------------------------------------------------------------ example 13 (RT, extra columns)


def test_example13_no_rt_columns(tmp_path: Path) -> None:
    header = W1[:6]
    res = convert([write_ssl(tmp_path, EX1_ROWS, header)])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


def test_example13_extra_column_warns(tmp_path: Path) -> None:
    header = [*W1, "ion-mobility"]
    rows = [{**r, "ion-mobility": "0.9"} for r in EX1_ROWS]
    res = convert([write_ssl(tmp_path, rows, header)])
    assert not errors(res.issues), res.issues
    issue = by_code(res.issues, "unknown_column")
    assert issue.severity == "warn"
    assert "ion-mobility" in issue.message
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


def test_unknown_column_warns_once_per_file(tmp_path: Path) -> None:
    header = [*W1, "ion-mobility"]
    a = write_ssl(tmp_path, [{**EX1_ROWS[0], "ion-mobility": "1"}], header, name="WT_1.ssl")
    b = write_ssl(
        tmp_path,
        [{**EX1_ROWS[0], "ion-mobility": "1", "file": "/x/mzml/WT_2.mzML"}],
        header,
        name="WT_2.ssl",
    )
    res = convert([a, b])
    assert not errors(res.issues), res.issues
    assert [i.code for i in res.issues].count("unknown_column") == 2


@pytest.mark.parametrize("col", ["retention-time", "start-time", "end-time"])
def test_example13_rt_values_never_read(tmp_path: Path, col: str) -> None:
    rows = [{**r, col: "abc"} for r in EX1_ROWS]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


def test_column_order_is_free(tmp_path: Path) -> None:
    header = list(reversed(W1))
    res = convert([write_ssl(tmp_path, EX1_ROWS, header)])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


# ------------------------------------------------------------------ D1 run, D2 scan


@pytest.mark.parametrize(
    ("file", "run"),
    [
        ("/data/mzml/WT_1.mzML", "WT_1"),
        ("WT_2.mzml.gz", "WT_2"),
        ("run3", "run3"),
        ("/a/b/X.MZML.GZ", "X"),
        ("/a/b/sample.mzML.mzML", "sample.mzML"),
        ("/a/b/run.raw", "run.raw"),
        ("/a/b/run.gz", "run"),
        (r"C:\data\WT_1.mzML", r"C:\data\WT_1"),
    ],
)
def test_d1_run_is_the_stem(tmp_path: Path, file: str, run: str) -> None:
    res = convert([write_ssl(tmp_path, [row("0", "2", "PEPTLDEK", "0.9", file=file)])])
    assert not errors(res.issues), res.issues
    assert res.peptides["run"].tolist() == [run]


@pytest.mark.parametrize("file", [None, "/x/mzml/", "/x/mzml/.mzML", ".mzML.gz"])
def test_d1_empty_file_or_stem_is_an_error(tmp_path: Path, file: str | None) -> None:
    """D1, step 3: a missing ``file`` or an empty stem is the error ``bad_value:file``."""
    res = convert([write_ssl(tmp_path, [row("0", "2", "PEPTLDEK", "0.9", file=file)])])
    issue = by_code(res.issues, "bad_value:file")
    assert issue.severity == "error"
    assert (issue.n_rows, issue.examples) == (1, ("1",))
    assert codes(errors(res.issues)) == {"bad_value:file"}
    assert res.peptides.empty


@pytest.mark.parametrize("scan", ["-1", "3.0", "abc", "1e3", None, "+3", " 3"])
def test_d2_scan_not_a_nonnegative_integer(tmp_path: Path, scan: str | None) -> None:
    rows = [row("0", "2", "SAMPLLER", "0.9"), {**row("x", "2", "PEPTLDEK", "0.9"), "scan": scan}]
    res = convert([write_ssl(tmp_path, rows)])
    issue = by_code(res.issues, "bad_type:scan")
    assert issue.severity == "error"
    assert issue.examples == ("2",)
    assert res.peptides.empty


def test_d2_scan_leading_zeros_are_distinct_ids(tmp_path: Path) -> None:
    """``[0-9]+``: ``0`` and ``00`` both pass and stay as written (0.0.7 writes ``str(idx)``)."""
    rows = [row("0", "2", "PEPTLDEK", "0.9"), row("00", "2", "SAMPLLER", "0.9")]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert sorted(res.peptides["spectrum_id"]) == ["cascadia:0", "cascadia:00"]


@pytest.mark.parametrize("scan", ["\u0663", "\uff13", "1\u0660"])
def test_d2_scan_non_ascii_digit_is_an_error(tmp_path: Path, scan: str) -> None:
    """``[0-9]+`` is ASCII only: Arabic-Indic or full-width digits are not a scan (D2)."""
    rows = [row("0", "2", "SAMPLLER", "0.9"), row(scan, "2", "PEPTLDEK", "0.9")]
    res = convert([write_ssl(tmp_path, rows)])
    issue = by_code(res.issues, "bad_type:scan")
    assert issue.severity == "error"
    assert issue.examples == ("2",)
    assert res.peptides.empty


def test_criterion6_spectrum_id_is_prefixed(tmp_path: Path) -> None:
    rows = [row(s, "2", f"PEPTLDEK{'ACDEF'[i]}", "0.9") for i, s in enumerate(["0", "7", "12"])]
    res = convert([write_ssl(tmp_path, rows)])
    sids = res.peptides["spectrum_id"].tolist()
    assert sids == ["cascadia:0", "cascadia:7", "cascadia:12"]
    assert not set(sids) & {"0", "7", "12"}


# ------------------------------------------------------------------ D7 search_db


@pytest.mark.parametrize("search_db", ["same", "other"])
def test_d7_search_db_is_invalid(tmp_path: Path, search_db: str) -> None:
    res = convert([write_ssl(tmp_path, EX1_ROWS)], search_db=search_db)
    assert by_code(res.issues, "search_db_invalid").severity == "error"
    assert res.peptides.empty


# ------------------------------------------------------------------ ties, criteria 3-5


def test_ties_keep_file_order(tmp_path: Path) -> None:
    """Same run, peptide, modified form and charge: file order is kept (stable sort)."""
    rows = [
        row("10", "2", "PEPTLDEK", "0.81"),
        row("9", "2", "PEPTLDEK", "0.99"),
        row("2", "2", "AAAPEPTLK", "0.9"),
    ]
    res = convert([write_ssl(tmp_path, rows)])
    assert not errors(res.issues), res.issues
    assert res.peptides["spectrum_id"].tolist() == ["cascadia:2", "cascadia:10", "cascadia:9"]


def test_criterion3_repeat_is_byte_identical(tmp_path: Path) -> None:
    p = write_ssl(tmp_path, EX1_ROWS)
    assert to_bytes(convert([p]), tmp_path / "a.tsv") == to_bytes(convert([p]), tmp_path / "b.tsv")
    ga = to_bytes(convert([p]), tmp_path / "a.tsv.gz")
    gb = to_bytes(convert([p]), tmp_path / "b.tsv.gz")
    assert ga == gb
    assert gzip.decompress(ga).decode() == EX1_OUT


def test_gzip_input_by_magic_bytes(tmp_path: Path) -> None:
    p = tmp_path / "WT_1.ssl"  # no .gz suffix: magic bytes decide (W1)
    p.write_bytes(gzip.compress(ssl_text(EX1_ROWS).encode(), mtime=0))
    res = convert([p])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


def test_crlf_and_bom_input(tmp_path: Path) -> None:
    p = tmp_path / "WT_1.ssl"
    p.write_bytes(b"\xef\xbb\xbf" + ssl_text(EX1_ROWS, newline="\r\n").encode())
    res = convert([p])
    assert not errors(res.issues), res.issues
    assert to_bytes(res, tmp_path / "p.tsv").decode() == EX1_OUT


def test_criterion4_matches_m10_08_stripping_rule() -> None:
    """M10-08's evaluator rule (docs/specs/M10-08-a3-feasibility.md "evaluate" step 2; not yet
    merged as validation/pilot/a3_feasibility.py, so encoded inline): delete [..] and (..)
    groups; delete signed masses and ``-``; delete characters outside A-Z; uppercase; fold I->L.
    """

    def m10_08(seq: str) -> str:
        seq = re.sub(r"\[[^\]]*\]|\([^)]*\)", "", seq)
        seq = re.sub(r"[+-]\d+(\.\d+)?", "", seq).replace("-", "")
        seq = re.sub(r"[^A-Z]", "", seq).upper()
        return seq.replace("I", "L")

    for seq, peptide in EXAMPLE2.items():
        assert strip_modifications(seq).replace("I", "L") == m10_08(seq) == peptide


def test_criterion5_columns_and_sources(tmp_path: Path) -> None:
    from gippa.steps.mapping import resolve_sources

    res = convert([write_ssl(tmp_path, EX1_ROWS)])
    cols = list(res.peptides.columns)
    assert cols == OUT_COLUMNS
    assert cols == [c for c in PEPTIDE_COLUMNS if c in cols]
    assert not set(ABSENT) & set(cols)
    assert set(res.peptides["source"]) == {"denovo"}
    out = tmp_path / "p.tsv"
    to_bytes(res, out)
    df, found = read_peptide_table(out)
    assert not errors(found), found
    cfg = GippaConfig()
    assert resolve_sources(df, cfg).tolist() == ["denovo", "denovo"]
    assert resolve_sources(df.drop(columns=["source"]), cfg).tolist() == ["denovo", "denovo"]


def test_criterion5_dtypes(tmp_path: Path) -> None:
    df = convert([write_ssl(tmp_path, EX1_ROWS)]).peptides
    assert str(df["charge"].dtype) == "Int64"
    assert df["score"].dtype == "float64"


def test_criterion7_no_tool_imports() -> None:
    tree = ast.parse(MODULE.read_text(encoding="utf-8"), filename=str(MODULE))
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.append(node.module)
    roots = {n.split(".")[0].lower() for n in names}
    assert not roots & {"cascadia", "pyteomics", "torch", "depthcharge"}, roots


def test_criterion7_check_imports_passes() -> None:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "check_imports", REPO / "scripts/ci/check_imports.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert MODULE.is_file()
    assert mod.violations(MODULE.parent) == []


# ------------------------------------------------------------------ properties

_TOKENS = ["[+15.99]", "[+0.98]", "[+57.02]"]
_NTERM = ["", "[+42.01]-", "[+43.01]-", "[-17.03]-", "[+25.980265]-"]
_RES = "ACDEFGHKLMNPQRSTVWY"  # no I: 0.0.7 folds I to L (W5)


@st.composite
def _sequence(draw: st.DrawFn) -> str:
    residues = draw(st.lists(st.sampled_from(_RES), min_size=7, max_size=30))
    out = [draw(st.sampled_from(_NTERM))]
    for r in residues:
        out.append(r + (draw(st.sampled_from(["", *_TOKENS])) if draw(st.booleans()) else ""))
    return "".join(out)


@st.composite
def _rows(draw: st.DrawFn) -> list[dict[str, str | None]]:
    n = draw(st.integers(min_value=1, max_value=12))
    gaps = draw(st.lists(st.integers(min_value=1, max_value=50), min_size=n, max_size=n))
    scans = list(itertools.accumulate(gaps, initial=draw(st.integers(0, 5))))[:n]
    out = []
    for s in scans:
        score = draw(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))
        out.append(row(str(s), str(draw(st.integers(1, 6))), draw(_sequence()), repr(score)))
    return out


_SETTINGS = settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)


@_SETTINGS
@given(rows=_rows(), perm=st.permutations(W1))
def test_property_mapping_invariants(
    tmp_path_factory: pytest.TempPathFactory, rows: list[dict[str, str | None]], perm: list[str]
) -> None:
    """Valid 0.0.7-shaped input: one row per input row, the multiset of (run, spectrum_id,
    modified_peptide, charge, score) equals the raw one, peptide = strip(sequence), all denovo;
    the column order of the input does not change the output bytes."""
    tmp = tmp_path_factory.mktemp("prop")
    p = write_ssl(tmp, rows)
    res = convert([p])
    assert not errors(res.issues), res.issues
    df = res.peptides
    assert len(df) == len(rows)
    got = sorted(
        zip(df["run"], df["spectrum_id"], df["modified_peptide"], df["charge"], df["score"],
            strict=True)
    )  # fmt: skip
    want = sorted(
        (
            "WT_1",
            f"cascadia:{r['scan']}",
            r["sequence"],
            int(str(r["charge"])),
            float(str(r["score"])),
        )
        for r in rows
    )
    assert got == want
    for seq, pep in zip(df["modified_peptide"], df["peptide"], strict=True):
        assert pep == strip_modifications(seq)
        assert re.fullmatch(r"[A-Z]+", pep)
    assert set(df["source"]) == {"denovo"}
    low = sum(float(str(r["score"])) <= 0.8 for r in rows)
    issue = [i for i in res.issues if i.code == "score_at_or_below_default_cutoff"]
    assert (issue[0].n_rows if issue else 0) == low
    q = write_ssl(tmp, rows, perm, name="perm/WT_1.ssl")
    assert to_bytes(convert([q]), tmp / "b.tsv") == to_bytes(res, tmp / "a.tsv")


@_SETTINGS
@given(rows=_rows())
def test_property_sorted_by_keys(
    tmp_path_factory: pytest.TempPathFactory, rows: list[dict[str, str | None]]
) -> None:
    tmp = tmp_path_factory.mktemp("prop")
    df = convert([write_ssl(tmp, rows)]).peptides
    keys = list(zip(df["run"], df["peptide"], df["modified_peptide"], df["charge"], strict=True))
    assert keys == sorted(keys)


@_SETTINGS
@given(score=st.floats(allow_nan=True, allow_infinity=True))
def test_property_score_range(tmp_path_factory: pytest.TempPathFactory, score: float) -> None:
    """Exactly the finite scores in [0, 1] convert (D5)."""
    tmp = tmp_path_factory.mktemp("prop")
    res = convert([write_ssl(tmp, [row("0", "2", "PEPTLDEK", repr(score))])])
    ok = math.isfinite(score) and 0.0 <= score <= 1.0
    assert (not errors(res.issues)) == ok, (score, res.issues)
