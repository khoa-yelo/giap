"""proBed 1.0 writer (issue M1-04; GIPPA_SPEC.md §7.2.9, §6.5; ADR-0005 item 1).

Always BED12 + 13 columns (25 fields); missing proBed values are ``.`` (not ``NA``).
"""

from __future__ import annotations

import pytest

from gippa.io.bed import BedRecord
from gippa.io.probed import PROBED_EXTRA, ProBedRecord, write_probed
from gippa.model import Interval, Replicon

C = Replicon("c1", 1000, "circular")
L = Replicon("c1", 1000, "linear")
HEADER = "# proBed-version 1.0"


def _data_lines(path) -> list[list[str]]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    return [ln.split("\t") for ln in text.splitlines() if not ln.startswith("#")]


def test_probed_extra_columns() -> None:
    assert PROBED_EXTRA == (
        "proteinAccession",
        "peptideSequence",
        "uniqueness",
        "genomeReferenceVersion",
        "psmScore",
        "fdr",
        "modifications",
        "charge",
        "expMassToCharge",
        "calcMassToCharge",
        "psmRank",
        "datasetID",
        "uri",
    )


# --- Example 8 ----------------------------------------------------------------------------------


def test_example_08_probed_fields(tmp_path) -> None:
    p = tmp_path / "mappings.probed"
    rec = ProBedRecord(
        BedRecord(Interval("c1", 9, 30, "+"), "pep1"),
        {"peptideSequence": "PEPTIDEK", "fdr": 0.001},
    )
    write_probed(p, [rec], {"c1": C}, "sha256:abc")
    (f,) = _data_lines(p)
    assert len(f) == 25
    assert f[0:6] == ["c1", "9", "30", "pep1", "0", "+"]
    assert f[6:8] == ["9", "30"]  # thick defaults to the full interval
    assert f[9] == "1"  # field 10 blockCount
    assert f[10].rstrip(",") == "21"  # field 11 blockSizes
    assert f[11].rstrip(",") == "0"  # field 12 chromStarts
    assert f[13] == "PEPTIDEK"  # field 14 peptideSequence
    assert f[17] == "0.001"  # field 18 fdr
    unspecified = [i for i, k in enumerate(PROBED_EXTRA) if k not in {"peptideSequence", "fdr"}]
    for i in unspecified:
        if PROBED_EXTRA[i] == "genomeReferenceVersion":
            continue
        assert f[12 + i] == ".", PROBED_EXTRA[i]


def test_genome_version_written(tmp_path) -> None:
    p = tmp_path / "x.probed"
    rec = ProBedRecord(BedRecord(Interval("c1", 9, 30, "+"), "pep1"), {})
    write_probed(p, [rec], {"c1": C}, "sha256:abc")
    (f,) = _data_lines(p)
    assert f[12 + PROBED_EXTRA.index("genomeReferenceVersion")] == "sha256:abc"


def test_header_comment_first_line(tmp_path) -> None:
    p = tmp_path / "x.probed"
    write_probed(p, [ProBedRecord(BedRecord(Interval("c1", 0, 3), "p"), {})], {"c1": C}, "v1")
    assert p.read_text().splitlines()[0] == HEADER


def test_empty_file_has_header_only(tmp_path) -> None:
    p = tmp_path / "x.probed"
    write_probed(p, [], {"c1": C}, "v1")
    assert p.read_text() == HEADER + "\n"


# --- Missing values, formats, validation --------------------------------------------------------


def test_missing_values_are_dot_not_na(tmp_path) -> None:
    p = tmp_path / "x.probed"
    rec = ProBedRecord(
        BedRecord(Interval("c1", 0, 30, "-"), "p"),
        {"proteinAccession": None, "charge": 2, "psmScore": 12.3456789, "uniqueness": "unique"},
    )
    write_probed(p, [rec], {"c1": C}, "v1")
    (f,) = _data_lines(p)
    extra = dict(zip(PROBED_EXTRA, f[12:], strict=True))
    assert extra["proteinAccession"] == "."
    assert extra["charge"] == "2"
    assert extra["psmScore"] == "12.3457"
    assert extra["uniqueness"] == "unique"
    assert "NA" not in f


def test_unknown_extra_key_raises(tmp_path) -> None:
    rec = ProBedRecord(BedRecord(Interval("c1", 0, 30, "+"), "p"), {"notAColumn": 1})
    with pytest.raises((ValueError, KeyError)):
        write_probed(tmp_path / "x.probed", [rec], {"c1": C}, "v1")


def test_score_out_of_range_raises(tmp_path) -> None:
    rec = ProBedRecord(BedRecord(Interval("c1", 0, 30, "+"), "p", score=1001), {})
    with pytest.raises(ValueError):
        write_probed(tmp_path / "x.probed", [rec], {"c1": C}, "v1")


def test_explicit_blocks_kept(tmp_path) -> None:
    p = tmp_path / "x.probed"
    bed = BedRecord(Interval("c1", 100, 160, "+"), "p", blocks=((20, 0), (10, 50)))
    write_probed(p, [ProBedRecord(bed, {})], {"c1": C}, "v1")
    (f,) = _data_lines(p)
    assert len(f) == 25
    assert f[9] == "2"
    assert [int(x) for x in f[10].split(",") if x] == [20, 10]
    assert [int(x) for x in f[11].split(",") if x] == [0, 50]


# --- Origin crossing ----------------------------------------------------------------------------


def test_origin_crossing_two_records_part2_first(tmp_path) -> None:
    p = tmp_path / "x.probed"
    rec = ProBedRecord(
        BedRecord(Interval("c1", 990, 1020, "+"), "pep"), {"peptideSequence": "PEPTIDEKR"}
    )
    write_probed(p, [rec], {"c1": C}, "v1")
    rows = _data_lines(p)
    assert len(rows) == 2
    assert all(len(r) == 25 for r in rows)
    assert [r[:4] for r in rows] == [
        ["c1", "0", "20", "pep_part2"],
        ["c1", "990", "1000", "pep_part1"],
    ]
    assert [r[9] for r in rows] == ["1", "1"]
    assert [r[10].rstrip(",") for r in rows] == ["20", "10"]
    assert all(r[13] == "PEPTIDEKR" for r in rows)


def test_linear_overrun_raises(tmp_path) -> None:
    rec = ProBedRecord(BedRecord(Interval("c1", 990, 1020, "+"), "pep"), {})
    with pytest.raises(ValueError):
        write_probed(tmp_path / "x.probed", [rec], {"c1": L}, "v1")


def test_sorted_and_deterministic(tmp_path) -> None:
    recs = [
        ProBedRecord(BedRecord(Interval("c1", 500, 530, "+"), "b"), {"fdr": 0.01}),
        ProBedRecord(BedRecord(Interval("c1", 10, 40, "-"), "a"), {"fdr": 0.02}),
    ]
    a, b = tmp_path / "a.probed", tmp_path / "b.probed"
    write_probed(a, recs, {"c1": C}, "v1")
    write_probed(b, list(reversed(recs)), {"c1": C}, "v1")
    assert a.read_bytes() == b.read_bytes()
    assert [r[3] for r in _data_lines(a)] == ["a", "b"]
