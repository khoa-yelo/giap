"""BED writer (issue M1-04; GIAP_SPEC.md §6.5, decision 16; ADR-0005 item 1).

Internal intervals are 0-based half-open and BED is too, so coordinates are written unchanged.
Origin-crossing features on circular replicons become two records ``_part1`` / ``_part2``.
"""

from __future__ import annotations

import pytest
from giap.io.bed import BedRecord, bed_lines, write_bed

from giap.model import Interval, Replicon

C = Replicon("c1", 1000, "circular")
L = Replicon("c1", 1000, "linear")


def _fields(line: str) -> list[str]:
    assert not line.endswith("\n")
    return line.split("\t")


def _ints(field: str) -> list[int]:
    """Parse a BED12 comma list, tolerating a trailing comma."""
    return [int(x) for x in field.split(",") if x != ""]


# --- Examples 5-7 -------------------------------------------------------------------------------


def test_example_05_bed6_minus_strand() -> None:
    rec = BedRecord(Interval("c1", 99, 200, "-"), "g1")
    assert bed_lines(rec, C) == ["c1\t99\t200\tg1\t0\t-"]


def test_example_06_origin_crossing_split() -> None:
    rec = BedRecord(Interval("c1", 950, 1050, "+"), "g2")
    lines = bed_lines(rec, C)
    assert sorted(lines) == sorted(["c1\t0\t50\tg2_part2\t0\t+", "c1\t950\t1000\tg2_part1\t0\t+"])


def test_example_06_write_bed_order(tmp_path) -> None:
    p = tmp_path / "x.bed"
    write_bed(p, [BedRecord(Interval("c1", 950, 1050, "+"), "g2")], {"c1": C})
    assert p.read_text() == "c1\t0\t50\tg2_part2\t0\t+\nc1\t950\t1000\tg2_part1\t0\t+\n"


def test_example_07_linear_overrun_raises() -> None:
    rec = BedRecord(Interval("c1", 950, 1050, "+"), "g2")
    with pytest.raises(ValueError):
        bed_lines(rec, L)


def test_example_07_write_bed_linear_overrun_raises(tmp_path) -> None:
    rec = BedRecord(Interval("c1", 950, 1050, "+"), "g2")
    with pytest.raises(ValueError):
        write_bed(tmp_path / "x.bed", [rec], {"c1": L})


# --- BED6 details -------------------------------------------------------------------------------


def test_bed6_unstranded_dot() -> None:
    assert bed_lines(BedRecord(Interval("c1", 0, 10, "."), "g"), C) == ["c1\t0\t10\tg\t0\t."]


def test_bed6_score_written() -> None:
    rec = BedRecord(Interval("c1", 0, 10, "+"), "g", score=1000)
    assert bed_lines(rec, C) == ["c1\t0\t10\tg\t1000\t+"]


def test_bed_full_length_linear_ok() -> None:
    assert bed_lines(BedRecord(Interval("c1", 0, 1000, "+"), "g"), L) == ["c1\t0\t1000\tg\t0\t+"]


def test_end_exactly_at_length_is_not_split() -> None:
    assert bed_lines(BedRecord(Interval("c1", 900, 1000, "+"), "g"), C) == [
        "c1\t900\t1000\tg\t0\t+"
    ]


@pytest.mark.parametrize("score", [-1, 1001, 5000])
def test_score_out_of_range_raises(score: int) -> None:
    rec = BedRecord(Interval("c1", 0, 10, "+"), "g", score=score)
    with pytest.raises(ValueError):
        bed_lines(rec, C)


def test_wrong_contig_raises() -> None:
    with pytest.raises(ValueError):
        bed_lines(BedRecord(Interval("c2", 0, 10, "+"), "g"), C)


# --- BED12 --------------------------------------------------------------------------------------


def test_bed12_linear_blocks_and_default_thick() -> None:
    rec = BedRecord(Interval("c1", 99, 200, "+"), "g1", blocks=((10, 0), (20, 81)))
    (line,) = bed_lines(rec, C)
    f = _fields(line)
    assert len(f) == 12
    assert f[:6] == ["c1", "99", "200", "g1", "0", "+"]
    assert f[6:8] == ["99", "200"]  # thick defaults to the full interval
    assert f[8] == "0"
    assert f[9] == "2"
    assert _ints(f[10]) == [10, 20]
    assert _ints(f[11]) == [0, 81]


def test_bed12_explicit_thick_and_rgb() -> None:
    rec = BedRecord(
        Interval("c1", 99, 200, "-"),
        "g1",
        score=5,
        thick=Interval("c1", 120, 180, "-"),
        rgb="255,0,0",
        blocks=((101, 0),),
    )
    (line,) = bed_lines(rec, C)
    f = _fields(line)
    assert f[:9] == ["c1", "99", "200", "g1", "5", "-", "120", "180", "255,0,0"]
    assert f[9] == "1"
    assert _ints(f[10]) == [101]
    assert _ints(f[11]) == [0]


def test_bed12_origin_crossing_blocks_split_and_rebased() -> None:
    # Interval [950, 1050) on C; block 1 = [950, 980), block 2 = [990, 1050) crosses the origin.
    rec = BedRecord(Interval("c1", 950, 1050, "+"), "g3", blocks=((30, 0), (60, 40)))
    lines = bed_lines(rec, C)
    assert len(lines) == 2
    by_name = {_fields(x)[3]: _fields(x) for x in lines}
    p1, p2 = by_name["g3_part1"], by_name["g3_part2"]
    assert len(p1) == len(p2) == 12

    assert p1[1:3] == ["950", "1000"]
    assert p1[6:8] == ["950", "1000"]
    assert p1[9] == "2"
    assert _ints(p1[10]) == [30, 10]
    assert _ints(p1[11]) == [0, 40]

    assert p2[1:3] == ["0", "50"]
    assert p2[6:8] == ["0", "50"]
    assert p2[9] == "1"
    assert _ints(p2[10]) == [50]
    assert _ints(p2[11]) == [0]


def test_bed12_origin_crossing_single_block() -> None:
    rec = BedRecord(Interval("c1", 950, 1050, "-"), "g4", blocks=((100, 0),))
    by_name = {_fields(x)[3]: _fields(x) for x in bed_lines(rec, C)}
    assert by_name["g4_part1"][9:12] and _ints(by_name["g4_part1"][10]) == [50]
    assert _ints(by_name["g4_part1"][11]) == [0]
    assert by_name["g4_part1"][9] == "1"
    assert by_name["g4_part2"][9] == "1"
    assert _ints(by_name["g4_part2"][10]) == [50]
    assert _ints(by_name["g4_part2"][11]) == [0]
    assert by_name["g4_part2"][5] == "-"


# --- write_bed ----------------------------------------------------------------------------------


def test_write_bed_sorted_by_contig_order_start_end_strand_name(tmp_path) -> None:
    reps = {"chr": Replicon("chr", 1000, "circular"), "p1": Replicon("p1", 500, "circular")}
    recs = [
        BedRecord(Interval("p1", 5, 10, "+"), "p"),
        BedRecord(Interval("chr", 100, 200, "."), "c"),
        BedRecord(Interval("chr", 100, 200, "-"), "b"),
        BedRecord(Interval("chr", 100, 200, "+"), "z"),
        BedRecord(Interval("chr", 100, 200, "+"), "a"),
        BedRecord(Interval("chr", 100, 150, "+"), "y"),
        BedRecord(Interval("chr", 50, 300, "+"), "x"),
    ]
    p = tmp_path / "x.bed"
    write_bed(p, recs, reps)
    names = [line.split("\t")[3] for line in p.read_text().splitlines()]
    assert names == ["x", "y", "a", "z", "b", "c", "p"]


def test_write_bed_contig_order_follows_mapping_not_alphabet(tmp_path) -> None:
    reps = {"zeta": Replicon("zeta", 100, "linear"), "alpha": Replicon("alpha", 100, "linear")}
    recs = [BedRecord(Interval("alpha", 0, 5, "+"), "a"), BedRecord(Interval("zeta", 0, 5), "z")]
    p = tmp_path / "x.bed"
    write_bed(p, recs, reps)
    assert [line.split("\t")[0] for line in p.read_text().splitlines()] == ["zeta", "alpha"]


def test_write_bed_empty_and_lf_endings(tmp_path) -> None:
    p = tmp_path / "x.bed"
    write_bed(p, [], {"c1": C})
    assert p.read_bytes() == b""
    write_bed(p, [BedRecord(Interval("c1", 1, 2, "+"), "g")], {"c1": C})
    assert p.read_bytes() == b"c1\t1\t2\tg\t0\t+\n"


def test_write_bed_deterministic(tmp_path) -> None:
    recs = [
        BedRecord(Interval("c1", 950, 1050, "+"), "g2"),
        BedRecord(Interval("c1", 99, 200, "-"), "g1"),
    ]
    a, b = tmp_path / "a.bed", tmp_path / "b.bed"
    write_bed(a, recs, {"c1": C})
    write_bed(b, list(reversed(recs)), {"c1": C})
    assert a.read_bytes() == b.read_bytes()


def test_bed_record_is_frozen() -> None:
    rec = BedRecord(Interval("c1", 0, 10, "+"), "g")
    assert rec.score == 0 and rec.thick is None and rec.rgb == "0" and rec.blocks == ()
    with pytest.raises(AttributeError):
        rec.name = "h"  # type: ignore[misc]
