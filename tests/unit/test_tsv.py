"""TSV reader/writer conventions (issue M1-04; GIAP_SPEC.md §6.6, §9.2, §9.3).

Floats ``%.6g``, integers as integers, booleans ``true``/``false``, ``NA`` for missing, lists
``;``-joined; ``.gz`` output is byte-reproducible (``mtime=0``, no filename in the header).
"""

from __future__ import annotations

import gzip
import math

import numpy as np
import pandas as pd
import pytest

from giap.io.tsv import COLUMN_NAME, NA, format_value, read_tsv, write_tsv

# --- Examples 1-4, 9 ----------------------------------------------------------------------------


def test_example_01_floats_ints_bools() -> None:
    assert format_value(0.000123456789) == "0.000123457"
    assert format_value(1234567.0) == "1.23457e+06"
    assert format_value(3) == "3"
    assert format_value(True) == "true"


def test_example_02_missing_negative_zero_list() -> None:
    assert format_value(None) == "NA"
    assert format_value(float("nan")) == "NA"
    assert format_value(-0.0) == "0"
    assert format_value(["a", "b"]) == "a;b"


def test_example_03_write_tsv_bytes(tmp_path) -> None:
    p = tmp_path / "genes.tsv"
    write_tsv(p, ["gene_id", "length_aa"], [["g1", 300]])
    assert p.read_bytes() == b"gene_id\tlength_aa\ng1\t300\n"


def test_example_04_column_not_snake_case(tmp_path) -> None:
    with pytest.raises(ValueError):
        write_tsv(tmp_path / "x.tsv", ["GeneID"], [["g1"]])


def test_example_09_gzip_reproducible(tmp_path) -> None:
    rows = [["g1", 300, 0.5, True, None], ["g2", 12, 1e-9, False, ["x", "y"]]]
    cols = ["gene_id", "length_aa", "score", "is_novel", "flags"]
    a, b = tmp_path / "a" / "x.tsv.gz", tmp_path / "b" / "x.tsv.gz"
    a.parent.mkdir()
    b.parent.mkdir()
    write_tsv(a, cols, rows)
    write_tsv(b, cols, rows)
    assert a.read_bytes() == b.read_bytes()
    # Rewriting the same path also gives the same bytes.
    first = a.read_bytes()
    write_tsv(a, cols, rows)
    assert a.read_bytes() == first


def test_gzip_header_has_zero_mtime_and_no_filename(tmp_path) -> None:
    p = tmp_path / "x.tsv.gz"
    write_tsv(p, ["gene_id"], [["g1"]])
    data = p.read_bytes()
    assert data[:2] == b"\x1f\x8b"
    assert data[3] & 0x08 == 0, "FNAME flag set: filename stored in gzip header"
    assert data[4:8] == b"\x00\x00\x00\x00", "gzip mtime must be 0"
    assert gzip.decompress(data) == b"gene_id\ng1\n"


def test_gzip_bytes_independent_of_file_name(tmp_path) -> None:
    a, b = tmp_path / "one.tsv.gz", tmp_path / "two.tsv.gz"
    write_tsv(a, ["gene_id"], [["g1"]])
    write_tsv(b, ["gene_id"], [["g1"]])
    assert a.read_bytes() == b.read_bytes()


# --- format_value -------------------------------------------------------------------------------


def test_na_constant() -> None:
    assert NA == "NA"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (False, "false"),
        (True, "true"),
        (0, "0"),
        (-7, "-7"),
        (np.int64(42), "42"),
        (np.int32(-5), "-5"),
        (2**53 + 1, "9007199254740993"),
        (10**20, "100000000000000000000"),
        (0.5, "0.5"),
        (1.0, "1"),
        (100.0, "100"),
        (1e-9, "1e-09"),
        (123456.7, "123457"),
        (np.float32(0.1), "0.1"),
        (np.float64(2.5), "2.5"),
        (float("inf"), "INF"),  # ADR-0006
        (float("-inf"), "-INF"),  # ADR-0006
        (math.nan, "NA"),
        (np.nan, "NA"),
        (pd.NA, "NA"),
        (-0.0, "0"),
        ("abc", "abc"),
        ("", ""),
        ("NA", "NA"),
        (("a", "b", "c"), "a;b;c"),
        (["x"], "x"),
        ([], "NA"),
        ((), "NA"),
    ],
)
def test_format_value(value, expected: str) -> None:
    assert format_value(value) == expected


def test_bool_checked_before_int() -> None:
    # bool is a subclass of int; it must not be written as 1/0.
    assert format_value(True) != "1"
    assert format_value(False) != "0"


def test_list_items_formatted_recursively() -> None:
    assert format_value(["a", 1, 0.25, True]) == "a;1;0.25;true"


@pytest.mark.parametrize("bad", ["a\tb", "a\nb", "a\rb", "\t", "end\n"])
def test_format_value_rejects_tab_and_newlines(bad: str) -> None:
    with pytest.raises(ValueError):
        format_value(bad)


def test_format_value_rejects_tab_inside_list_item() -> None:
    with pytest.raises(ValueError):
        format_value(["ok", "bad\titem"])


# --- write_tsv ----------------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["gene_id", "a", "length_aa", "x1", "p_value_bh", "a_1_b"])
def test_column_name_pattern_accepts(name: str) -> None:
    assert COLUMN_NAME.match(name)


@pytest.mark.parametrize("name", ["GeneID", "1abc", "_x", "gene-id", "gene id", "", "Gene_id"])
def test_column_name_pattern_rejects(name: str) -> None:
    assert not COLUMN_NAME.match(name)


@pytest.mark.parametrize("name", ["gene-id", "_x", "1abc", "gene id", ""])
def test_write_tsv_rejects_bad_column(tmp_path, name: str) -> None:
    with pytest.raises(ValueError):
        write_tsv(tmp_path / "x.tsv", ["ok", name], [["a", "b"]])


def test_write_tsv_rejects_row_width_mismatch(tmp_path) -> None:
    with pytest.raises(ValueError):
        write_tsv(tmp_path / "x.tsv", ["a", "b"], [["1", "2"], ["3"]])
    with pytest.raises(ValueError):
        write_tsv(tmp_path / "y.tsv", ["a", "b"], [["1", "2", "3"]])


def test_write_tsv_zero_rows_header_only(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["gene_id", "length_aa"], [])
    assert p.read_bytes() == b"gene_id\tlength_aa\n"


def test_write_tsv_empty_string_is_empty_field(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["a", "b", "c"], [["", None, "z"]])
    assert p.read_bytes() == b"a\tb\tc\n\tNA\tz\n"


def test_write_tsv_formats_all_types(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(
        p,
        ["gene_id", "length_aa", "score", "is_novel", "flags", "missing"],
        [["g1", 300, 1234567.0, True, ["low_coverage", "edge"], float("nan")]],
    )
    assert p.read_text(encoding="utf-8") == (
        "gene_id\tlength_aa\tscore\tis_novel\tflags\tmissing\n"
        "g1\t300\t1.23457e+06\ttrue\tlow_coverage;edge\tNA\n"
    )


def test_write_tsv_lf_endings_no_trailing_tab(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["a", "b"], [["1", "2"], ["3", None]])
    data = p.read_bytes()
    assert b"\r" not in data
    assert data.endswith(b"\n")
    for line in data.split(b"\n")[:-1]:
        assert not line.endswith(b"\t")


def test_write_tsv_utf8_no_bom(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["name"], [["β-lactamase"]])
    data = p.read_bytes()
    assert not data.startswith(b"\xef\xbb\xbf")
    assert data == "name\nβ-lactamase\n".encode()


def test_write_tsv_rejects_tab_in_cell(tmp_path) -> None:
    with pytest.raises(ValueError):
        write_tsv(tmp_path / "x.tsv", ["a"], [["x\ty"]])


def test_write_tsv_accepts_generator_rows(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["n"], ([i] for i in range(3)))
    assert p.read_bytes() == b"n\n0\n1\n2\n"


# --- read_tsv -----------------------------------------------------------------------------------


def test_read_tsv_round_trip_na_and_empty_to_none(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["a", "b", "c", "d"], [["g1", 300, None, ""], ["g2", 0.5, True, ["x", "y"]]])
    header, rows = read_tsv(p)
    assert header == ["a", "b", "c", "d"]
    assert rows == [["g1", "300", None, None], ["g2", "0.5", "true", "x;y"]]


def test_read_tsv_gzip(tmp_path) -> None:
    p = tmp_path / "x.tsv.gz"
    write_tsv(p, ["a", "b"], [["1", None]])
    assert read_tsv(p) == (["a", "b"], [["1", None]])


def test_read_tsv_gzip_detected_by_magic(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    p.write_bytes(gzip.compress(b"a\tb\n1\tNA\n"))
    assert read_tsv(p) == (["a", "b"], [["1", None]])


def test_read_tsv_strips_bom(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    p.write_bytes(b"\xef\xbb\xbfgene_id\tn\ng1\t3\n")
    header, rows = read_tsv(p)
    assert header == ["gene_id", "n"]
    assert rows == [["g1", "3"]]


def test_read_tsv_zero_rows(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    write_tsv(p, ["a", "b"], [])
    assert read_tsv(p) == (["a", "b"], [])


def test_read_tsv_ragged_row_reports_line_number(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    p.write_bytes(b"a\tb\n1\t2\n3\n4\t5\n")
    with pytest.raises(ValueError, match=r"\b3\b"):
        read_tsv(p)


def test_read_tsv_ragged_row_too_wide(tmp_path) -> None:
    p = tmp_path / "x.tsv"
    p.write_bytes(b"a\tb\n1\t2\t3\n")
    with pytest.raises(ValueError, match=r"\b2\b"):
        read_tsv(p)
