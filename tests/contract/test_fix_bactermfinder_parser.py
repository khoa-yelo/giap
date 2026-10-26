"""Contract test: BacTermFinder parser on its real output table.

Spec: docs/specs/fix-bactermfinder-parser.md (GIAP_SPEC.md §5.4, §6.5, §7.1.1; backlog M2-15).

Recordings in ``tests/data/recorded/bactermfinder/024870447f62/`` (PROVENANCE.md in each):

* ``ref_genome_test/``: 120-row excerpt of the output that ships with the pinned commit for its
  bundled genome (CP054306.1, 630 bp, step 3);
* ``smoke_1895829/``: the 4 rows that FarmShare Slurm job 1895829 printed (contig ``tiny``,
  step 10).

``SampleName`` = ``<u_id>_<seq_id>_<start>_<end>_<strand>``; ``start``/``end`` are the Python
slice bounds ``seq[i:i + 101]`` (genome_scan.py lines 44-53), i.e. 0-based half-open, so the
internal interval is ``[start, end)`` with no shift. ``seq_id`` may contain ``_``: split it
from the right.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.terminators import Terminator, parse_bactermfinder
from giap.io.fasta import read_genome_fasta
from giap.model.coords import Interval, Replicon, check_interval

REC = Path(__file__).resolve().parents[1] / "data" / "recorded" / "bactermfinder" / "024870447f62"
REF = REC / "ref_genome_test"
SMOKE = REC / "smoke_1895829"
REF_CSV = REF / "outref-genome-test.fasta_mean.excerpt.csv"
SMOKE_CSV = SMOKE / "outtiny.fa_mean.head.csv"
HEADER = ",".join(
    ["SampleName"] + [f"probability_{k}" for k in ("binary", "ENAC", "PS2", "NCP", "mean")]
)
WINDOW = 101


def _reps(fasta: Path) -> dict[str, Replicon]:
    return {r.id: Replicon(r.id, len(r.sequence), "circular") for r in read_genome_fasta(fasta)}


def _expected(path: Path) -> list[tuple[str, int, int, str, float]]:
    """Independent oracle: split SampleName (u_id from the left, the rest from the right)."""
    out = []
    for row in csv.DictReader(path.read_text().splitlines()):
        _u, rest = row["SampleName"].split("_", 1)
        contig, start, end, strand = rest.rsplit("_", 3)
        out.append((contig, int(start), int(end), strand, float(row["probability_mean"])))
    return out


def _key(t: Terminator) -> tuple[str, int, int, str, float]:
    iv = t.interval
    return (iv.contig, iv.start, iv.end, iv.strand, t.confidence)


def _write(tmp: Path, names_probs: list[tuple[str, str]], header: str = HEADER) -> Path:
    p = tmp / "outgenome.fa_mean.csv"
    lines = [header] + [f"{n},0.1,0.2,0.3,0.4,{pr}" for n, pr in names_probs]
    p.write_text("\n".join(lines) + "\n")
    return p


# ------------------------------------------------------------------ recorded samples


def test_recorded_header_is_the_real_one() -> None:
    assert REF_CSV.read_text().splitlines()[0] == HEADER
    assert SMOKE_CSV.read_text().splitlines()[0] == HEADER


def test_recorded_excerpt_gives_the_expected_rows_in_file_order() -> None:
    reps = _reps(REF / "ref-genome-test.fasta")
    assert reps == {"CP054306.1": Replicon("CP054306.1", 630, "circular")}
    got = parse_bactermfinder(REF_CSV, reps)
    assert len(got) == 120
    assert [_key(t) for t in got] == _expected(REF_CSV)
    for t in got:
        check_interval(t.interval, reps["CP054306.1"])
        assert t.source == "BacTermFinder"
        assert 0.0 <= t.confidence <= 1.0
        assert t.interval.end - t.interval.start == WINDOW


def test_recorded_excerpt_spot_rows() -> None:
    got = parse_bactermfinder(REF_CSV, _reps(REF / "ref-genome-test.fasta"))
    by = {(t.interval.start, t.interval.strand): t for t in got}
    # first window: 1_CP054306.1_0_101_+  -> [0, 101)
    assert got[0].interval == Interval("CP054306.1", 0, 101, "+")
    assert got[0].confidence == pytest.approx(0.0300241185, abs=1e-12)
    assert got[1].interval == Interval("CP054306.1", 0, 101, "-")
    assert got[1].confidence == pytest.approx(0.031558088500000005, abs=1e-12)
    # 79_CP054306.1_234_335_+  (highest probability_mean of the table)
    assert by[(234, "+")].interval == Interval("CP054306.1", 234, 335, "+")
    assert by[(234, "+")].confidence == pytest.approx(0.190896124, abs=1e-12)
    # 33_CP054306.1_96_197_-
    assert by[(96, "-")].confidence == pytest.approx(0.1732537305, abs=1e-12)
    # last window 177_CP054306.1_528_629_- : ends inside the 630-bp replicon
    assert got[-1].interval == Interval("CP054306.1", 528, 629, "-")
    assert got[-1].confidence == pytest.approx(0.017605707249999998, abs=1e-12)


def test_coordinates_are_0_based_half_open_without_shift() -> None:
    """genome_scan.py: window u_id covers seq[i:i+101] with i = step * (u_id - 1), written as
    start=i, end=i+101. A 1-based reading (start - 1) would give start = -1 for the first window."""
    got = parse_bactermfinder(REF_CSV, _reps(REF / "ref-genome-test.fasta"))
    names = [r["SampleName"] for r in csv.DictReader(REF_CSV.read_text().splitlines())]
    for t, name in zip(got, names, strict=True):
        u_id = int(name.split("_", 1)[0])
        assert t.interval.start == 3 * (u_id - 1)
        assert t.interval.end == t.interval.start + WINDOW
    assert min(t.interval.start for t in got) == 0


def test_each_window_has_a_plus_then_minus_row_with_the_same_interval() -> None:
    got = parse_bactermfinder(REF_CSV, _reps(REF / "ref-genome-test.fasta"))
    for plus, minus in zip(got[::2], got[1::2], strict=True):
        assert (plus.interval.strand, minus.interval.strand) == ("+", "-")
        assert (plus.interval.start, plus.interval.end) == (
            minus.interval.start,
            minus.interval.end,
        )


def test_smoke_job_head() -> None:
    reps = _reps(SMOKE / "tiny.fa")
    assert reps == {"tiny": Replicon("tiny", 630, "circular")}
    got = parse_bactermfinder(SMOKE_CSV, reps)
    assert [_key(t) for t in got] == [
        ("tiny", 0, 101, "+", 0.030024133499999998),
        ("tiny", 0, 101, "-", 0.031558096),
        ("tiny", 10, 111, "+", 0.0380584085),
        ("tiny", 10, 111, "-", 0.0196476205),
    ]


# ------------------------------------------------------------------ SampleName splitting


@pytest.mark.parametrize(
    ("name", "contig", "start", "end", "strand"),
    [
        ("7_NZ_CP012345.1_18_119_-", "NZ_CP012345.1", 18, 119, "-"),
        ("1_contig_2_5_0_101_+", "contig_2_5", 0, 101, "+"),
        ("12_a_b_c_33_134_+", "a_b_c", 33, 134, "+"),
        ("3_chr1_6_107_-", "chr1", 6, 107, "-"),
    ],
)
def test_contig_ids_with_underscores_split_from_the_right(
    tmp_path: Path, name: str, contig: str, start: int, end: int, strand: str
) -> None:
    reps = {contig: Replicon(contig, 1_000, "circular")}
    (t,) = parse_bactermfinder(_write(tmp_path, [(name, "0.42")]), reps)
    assert t.interval == Interval(contig, start, end, strand)  # type: ignore[arg-type]
    assert t.confidence == 0.42 and t.source == "BacTermFinder"


def test_several_contigs_u_id_restarts(tmp_path: Path) -> None:
    rows = [
        ("1_NZ_A_1_0_101_+", "0.5"),
        ("1_NZ_A_1_0_101_-", "0.6"),
        ("1_plasmid_2_0_101_+", "0.7"),
        ("2_plasmid_2_3_104_-", "0.8"),
    ]
    reps = {"NZ_A_1": Replicon("NZ_A_1", 500, "circular"),
            "plasmid_2": Replicon("plasmid_2", 300, "circular")}  # fmt: skip
    got = parse_bactermfinder(_write(tmp_path, rows), reps)
    assert [_key(t) for t in got] == [
        ("NZ_A_1", 0, 101, "+", 0.5),
        ("NZ_A_1", 0, 101, "-", 0.6),
        ("plasmid_2", 0, 101, "+", 0.7),
        ("plasmid_2", 3, 104, "-", 0.8),
    ]


CONTIG = st.text(
    alphabet=st.sampled_from("ABCNZPabcxyz0123456789._|-"), min_size=1, max_size=20
).filter(lambda s: s.strip() == s)


@settings(max_examples=200, deadline=None)
@given(
    rows=st.lists(
        st.tuples(
            CONTIG,
            st.integers(0, 10_000),
            st.integers(1, 500),
            st.sampled_from(["+", "-"]),
            st.floats(0.0, 1.0, allow_nan=False),
        ),
        min_size=1,
        max_size=20,
    )
)
def test_property_sample_name_round_trip(
    tmp_path_factory: pytest.TempPathFactory,
    rows: list[tuple[str, int, int, str, float]],
) -> None:
    """Any SampleName genome_scan.py can write (seq_id with ``_``, ``-``, ``.``) parses back to
    exactly the (contig, start, end, strand) it encodes, and confidence = probability_mean."""
    reps = {c: Replicon(c, 20_000, "circular") for c, *_ in rows}
    names = [
        (f"{i + 1}_{c}_{s}_{s + n}_{strand}", repr(p))
        for i, (c, s, n, strand, p) in enumerate(rows)
    ]
    got = parse_bactermfinder(_write(tmp_path_factory.mktemp("btf"), names), reps)
    assert [_key(t) for t in got] == [(c, s, s + n, strand, p) for c, s, n, strand, p in rows]


# ------------------------------------------------------------------ edge cases and errors


def test_header_only_and_empty_files_give_no_calls(tmp_path: Path) -> None:
    reps = {"tiny": Replicon("tiny", 630, "circular")}
    assert parse_bactermfinder(_write(tmp_path, []), reps) == []
    empty = tmp_path / "empty.csv"
    empty.write_text("")
    assert parse_bactermfinder(empty, reps) == []


@pytest.mark.parametrize("missing", ["SampleName", "probability_mean"])
def test_missing_required_column_names_the_file(tmp_path: Path, missing: str) -> None:
    header = HEADER.replace(missing, "something_else")
    p = _write(tmp_path, [("1_tiny_0_101_+", "0.1")], header=header)
    with pytest.raises(ValueError, match=r"outgenome\.fa_mean\.csv"):
        parse_bactermfinder(p, {"tiny": Replicon("tiny", 630, "circular")})


@pytest.mark.parametrize(
    "name",
    [
        "1_tiny_0_101",  # no strand
        "tiny_0_101_+",  # no u_id
        "1_tiny_zero_101_+",  # non-integer start
        "1_tiny_0_101_.",  # strand not +/-
        "1_tiny_0_101_x",
    ],
)
def test_malformed_sample_name_names_file_and_line(tmp_path: Path, name: str) -> None:
    p = _write(tmp_path, [("1_tiny_0_101_+", "0.1"), (name, "0.2")])
    with pytest.raises(ValueError, match=r"outgenome\.fa_mean\.csv:3"):
        parse_bactermfinder(p, {"tiny": Replicon("tiny", 630, "circular")})


def test_unknown_contig_or_out_of_range_interval_raise(tmp_path: Path) -> None:
    reps = {"tiny": Replicon("tiny", 630, "circular")}
    with pytest.raises(ValueError, match=r"outgenome\.fa_mean\.csv:2.*other"):
        parse_bactermfinder(_write(tmp_path, [("1_other_0_101_+", "0.1")]), reps)
    with pytest.raises(ValueError):
        parse_bactermfinder(_write(tmp_path, [("1_tiny_700_801_+", "0.1")]), reps)
