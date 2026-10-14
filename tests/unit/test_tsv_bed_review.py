"""Review pins for TSV/BED/proBed writers (ADR-0006; GIAP_SPEC.md §6.5, §6.6, §9.2, §9.5).

Covers: CRLF/BOM tolerance of ``read_tsv`` (must survive consolidation onto ``giap.io._open``),
numpy scalar formatting, ±inf/NaN encoding, strict list cells, and strict BED12/proBed writing.
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pytest

from giap.io.bed import BedRecord, bed_lines
from giap.io.probed import PROBED_EXTRA, ProBedRecord, write_probed
from giap.io.tsv import format_value, read_tsv, write_tsv
from giap.model import Interval, Replicon

C = Replicon("c1", 1000, "circular")
L = Replicon("c1", 1000, "linear")
IV = Interval("c1", 100, 200, "+")  # length 100

# --- 1. read_tsv: CRLF line endings and a UTF-8 BOM ---------------------------------------------


@pytest.mark.parametrize("gz", [False, True])
def test_read_tsv_crlf_and_bom_match_lf(tmp_path: Path, gz: bool) -> None:
    lf = b"a\tb\n1\tx\n2\tNA\n"
    crlf = b"\xef\xbb\xbf" + lf.replace(b"\n", b"\r\n")
    p_lf, p_crlf = tmp_path / "lf.tsv", tmp_path / "crlf.tsv"
    p_lf.write_bytes(lf)
    p_crlf.write_bytes(gzip.compress(crlf, mtime=0) if gz else crlf)
    expected = (["a", "b"], [["1", "x"], ["2", None]])
    assert read_tsv(p_lf) == expected
    assert read_tsv(p_crlf) == expected  # no BOM in "a", no trailing "\r" in the last column


# --- 2. numpy scalars --------------------------------------------------------------------------


def test_numpy_bool_formats_as_true_false() -> None:
    assert format_value(np.bool_(True)) == "true"
    assert format_value(np.bool_(False)) == "false"


def test_numpy_int_and_float_format_like_python() -> None:
    assert format_value(np.int64(42)) == format_value(42) == "42"
    assert format_value(np.float64(0.1234567)) == format_value(0.1234567) == "0.123457"
    assert format_value(np.float64(3.0)) == format_value(3.0)


# --- 3. ADR-0006 item 1: infinities and NaN ------------------------------------------------------


def test_inf_and_nan_encoding() -> None:
    assert format_value(float("inf")) == "INF"
    assert format_value(float("-inf")) == "-INF"
    assert format_value(float("nan")) == "NA"
    assert format_value(np.float64("inf")) == "INF"
    assert format_value(np.float64("-inf")) == "-INF"


def test_inf_round_trip_through_tsv(tmp_path: Path) -> None:
    p = tmp_path / "t.tsv"
    write_tsv(p, ["v"], [[float("inf")], [float("-inf")], [float("nan")], [None]])
    assert p.read_text(encoding="utf-8") == "v\nINF\n-INF\nNA\nNA\n"
    _, rows = read_tsv(p)
    assert float(rows[0][0]) == math.inf
    assert float(rows[1][0]) == -math.inf
    assert rows[2][0] is None and rows[3][0] is None


# --- 4. ADR-0006 item 2: list cells ------------------------------------------------------------


def test_list_item_containing_separator_is_error() -> None:
    assert format_value(["a", "b"]) == "a;b"
    with pytest.raises(ValueError):
        format_value(["a;b", "c"])


def test_list_with_missing_item_is_error() -> None:
    with pytest.raises(ValueError):
        format_value(["a", None])


# --- 5. ADR-0006 item 3: strict BED ------------------------------------------------------------


@pytest.mark.parametrize(
    "blocks",
    [
        ((50, 10), (40, 60)),  # first block does not start at 0
        ((40, 60), (50, 0)),  # unsorted
        ((60, 0), (50, 50)),  # overlapping ([0,60) and [50,100))
        ((50, 0), (30, 60)),  # last block ends at 90, record length is 100
    ],
    ids=["not_at_0", "unsorted", "overlapping", "short_end"],
)
def test_bed12_bad_blocks_raise_value_error(blocks: tuple[tuple[int, int], ...]) -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, "g", blocks=blocks), L)


@pytest.mark.parametrize(
    "thick",
    [Interval("c1", 90, 150, "+"), Interval("c1", 150, 210, "+"), Interval("c1", 90, 210, "+")],
    ids=["left", "right", "both"],
)
def test_bed12_thick_outside_interval_raises(thick: Interval) -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, "g", thick=thick, blocks=((100, 0),)), L)


def test_bed12_wrapped_thick_split_across_origin() -> None:
    # Interval [950, 1050) on a 1000-bp circle; thick [980, 1020) in the same unwrapped frame.
    # part1 = [950, 1000): thick ∩ part1 = [980, 1000); one block of 1000-950 = 50 at offset 0.
    # part2 = [0, 1050-1000) = [0, 50): thick ∩ part2 = [1000, 1020) - 1000 = [0, 20);
    #         one block of 50 at offset 0.
    rec = BedRecord(
        Interval("c1", 950, 1050, "+"),
        "g",
        thick=Interval("c1", 980, 1020, "+"),
        blocks=((100, 0),),
    )
    assert sorted(bed_lines(rec, C)) == sorted(
        [
            "c1\t950\t1000\tg_part1\t0\t+\t980\t1000\t0\t1\t50\t0",
            "c1\t0\t50\tg_part2\t0\t+\t0\t20\t0\t1\t50\t0",
        ]
    )


def test_bed12_thick_only_in_part2_gives_empty_thick_on_part1() -> None:
    # Thick [1010, 1040) lies wholly after the origin. part1 [950, 1000) has no thick region,
    # written as thickStart = thickEnd = chromStart = 950 (BED convention for "no thick").
    # part2: [1010, 1040) - 1000 = [10, 40).
    rec = BedRecord(
        Interval("c1", 950, 1050, "+"),
        "g",
        thick=Interval("c1", 1010, 1040, "+"),
        blocks=((100, 0),),
    )
    assert sorted(bed_lines(rec, C)) == sorted(
        [
            "c1\t950\t1000\tg_part1\t0\t+\t950\t950\t0\t1\t50\t0",
            "c1\t0\t50\tg_part2\t0\t+\t10\t40\t0\t1\t50\t0",
        ]
    )


@pytest.mark.parametrize("name", ["a b", "a\tb", "a\nb", " a"])
def test_bed_name_with_whitespace_raises(name: str) -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, name), L)


@pytest.mark.parametrize("score", [True, False, 1.5, 500.0, "5", -1, 1001])
def test_bed_bad_score_raises(score: object) -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, "g", score=score), L)  # type: ignore[arg-type]


@pytest.mark.parametrize("score", [0, 1000])
def test_bed_score_bounds_accepted(score: int) -> None:
    assert bed_lines(BedRecord(IV, "g", score=score), L) == [f"c1\t100\t200\tg\t{score}\t+"]


def test_bed6_with_thick_raises() -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, "g", thick=Interval("c1", 120, 180, "+")), L)


def test_bed6_with_rgb_raises() -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(IV, "g", rgb="255,0,0"), L)


# --- 6. proBed extras --------------------------------------------------------------------------


def _probed_fields(tmp_path: Path, extra: dict[str, object]) -> dict[str, str]:
    p = tmp_path / "x.probed"
    write_probed(p, [ProBedRecord(BedRecord(IV, "pep"), extra)], {"c1": L}, "v1")  # type: ignore[arg-type]
    data = [ln for ln in p.read_text(encoding="utf-8").splitlines() if not ln.startswith("#")]
    assert len(data) == 1
    f = data[0].split("\t")
    assert len(f) == 25
    return dict(zip(PROBED_EXTRA, f[12:], strict=True))


def test_probed_modifications_list_joined_with_comma(tmp_path: Path) -> None:
    got = _probed_fields(tmp_path, {"modifications": ["3-Oxidation", "7-Phospho"]})
    assert got["modifications"] == "3-Oxidation,7-Phospho"


def test_probed_empty_string_extra_is_dot(tmp_path: Path) -> None:
    got = _probed_fields(tmp_path, {"uri": "", "datasetID": ""})
    assert got["uri"] == "." and got["datasetID"] == "."
