"""Hypothesis properties for OSTIR rows (M2-16 acceptance criterion; GIPPA_SPEC.md §7.1.6).

On circular replicons ``len(seq) == up + down`` (100 by default) and
``seq[start-1:start+2] == gene.start_codon``; on linear contigs the same identity holds after
clipping. Rows equal a per-base brute-force oracle, ``-`` rows equal the ``+`` row of the mirrored
gene on the reverse-complement genome, the writer round-trips through ``csv`` and the parser keeps
exactly the rows at the requested start.
"""

from __future__ import annotations

import csv
from pathlib import Path
from tempfile import TemporaryDirectory

from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

ASD = "ACCTCCTTA"
_COMP = str.maketrans("ACGT", "TGCA")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


@st.composite
def cases(draw: st.DrawFn, topology: str) -> tuple[Gene, Replicon, str]:
    n = draw(st.integers(min_value=200, max_value=1500))
    genome = draw(st.text(alphabet="ACGT", min_size=n, max_size=n))
    length = 3 * draw(st.integers(min_value=1, max_value=min(300, (n - 1) // 3)))
    strand = draw(st.sampled_from("+-"))
    if topology == "circular":
        start = draw(st.integers(min_value=0, max_value=n - 1))  # may cross the origin
    else:
        start = draw(st.integers(min_value=0, max_value=n - length))
    end = start + length
    if strand == "+":
        codon = "".join(genome[(start + k) % n] for k in range(3))
    else:
        codon = rc("".join(genome[(end - 3 + k) % n] for k in range(3)))
    g = Gene(
        locus_tag="G1",
        feature_type="cds",
        interval=Interval("c", start, end, strand),  # type: ignore[arg-type]
        aa="",
        nt=None,
        start_codon=codon,
    )
    return g, Replicon("c", n, topology, "unknown"), genome  # type: ignore[arg-type]


def oracle(g: Gene, rep: Replicon, genome: str, up: int, down: int) -> tuple[str, int]:
    """Per-base window and the 1-based start, wrapping or clipping one position at a time."""
    n, iv = rep.length, g.interval
    out: list[str] = []
    kept_up = 0
    for i in range(-up, down):
        if iv.strand == "+":
            pos = iv.start + i
            base = lambda p: genome[p % n]  # noqa: E731
        else:
            pos = iv.end - 1 - i
            base = lambda p: genome[p % n].translate(_COMP)  # noqa: E731
        if rep.topology == "linear" and not 0 <= pos < n:
            continue
        out.append(base(pos))
        kept_up += i < 0
    return "".join(out), kept_up + 1


@settings(max_examples=200, deadline=None)
@given(cases("circular"))
def test_circular_rows_have_100_nt_and_the_start_codon_at_51(c: tuple[Gene, Replicon, str]) -> None:
    from gippa.annotate.ostir import ostir_row

    g, rep, genome = c
    row = ostir_row(g, rep, genome, ASD)
    assert len(row.seq) == 100
    assert row.start == row.end == 51
    assert row.seq[row.start - 1 : row.start + 2] == g.start_codon
    assert (row.id, row.anti_sd) == ("G1", ASD)
    assert row.seq == oracle(g, rep, genome, 50, 50)[0]


@settings(max_examples=200, deadline=None)
@given(cases("linear"))
def test_linear_rows_clip_and_keep_the_identity(c: tuple[Gene, Replicon, str]) -> None:
    from gippa.annotate.ostir import ostir_row

    g, rep, genome = c
    row = ostir_row(g, rep, genome, ASD)
    seq, start = oracle(g, rep, genome, 50, 50)
    assert (row.seq, row.start, row.end) == (seq, start, start)
    assert 1 <= row.start <= 51 and len(row.seq) <= 100
    assert row.seq[row.start - 1 : row.start + 2] == g.start_codon


@settings(max_examples=100, deadline=None)
@given(
    cases("circular") | cases("linear"),
    st.integers(min_value=0, max_value=80),
    st.integers(min_value=3, max_value=80),
)
def test_window_sizes_follow_up_and_down(c: tuple[Gene, Replicon, str], up: int, down: int) -> None:
    from gippa.annotate.ostir import ostir_row

    g, rep, genome = c
    row = ostir_row(g, rep, genome, ASD, up=up, down=down)
    seq, start = oracle(g, rep, genome, up, down)
    assert (row.seq, row.start, row.end) == (seq, start, start)
    if rep.topology == "circular":
        assert len(row.seq) == up + down and row.start == up + 1


@settings(max_examples=100, deadline=None)
@given(cases("circular") | cases("linear"))
def test_minus_row_equals_plus_row_of_the_mirrored_gene(c: tuple[Gene, Replicon, str]) -> None:
    from gippa.annotate.ostir import ostir_row

    g, rep, genome = c
    n, iv = rep.length, g.interval
    flip = "-" if iv.strand == "+" else "+"
    ms, me = n - iv.end, n - iv.start
    if ms < 0:  # origin-crossing: keep the unwrapped convention (0 <= start < length)
        ms, me = ms + n, me + n
    mirrored = Gene(
        locus_tag="G1",
        feature_type="cds",
        interval=Interval("c", ms, me, flip),  # type: ignore[arg-type]
        aa="",
        nt=None,
        start_codon=g.start_codon,
    )
    assert ostir_row(g, rep, genome, ASD) == ostir_row(mirrored, rep, rc(genome), ASD)


ids = st.text(alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789", min_size=1, max_size=12)


@settings(max_examples=100, deadline=None)
@given(
    st.lists(
        st.tuples(ids, st.text(alphabet="ACGT", min_size=3, max_size=120),
                  st.integers(min_value=1, max_value=51)),
        max_size=20,
        unique_by=lambda t: t[0],
    )
)  # fmt: skip
def test_writer_round_trips_through_csv(rows: list[tuple[str, str, int]]) -> None:
    from gippa.annotate.ostir import OstirRow, write_ostir_input

    objs = [OstirRow(i, s, p, p, ASD) for i, s, p in rows]
    with TemporaryDirectory() as d:
        path = Path(d) / "ostir_input.csv"
        write_ostir_input(path, objs)
        raw = path.read_bytes()
        with path.open(newline="") as fh:
            got = list(csv.reader(fh))
    assert b"\r" not in raw and raw.endswith(b"\n")
    assert got[0] == ["id", "seq", "start", "end", "anti-Shine-Dalgarno"]
    assert got[1:] == [[i, s, str(p), str(p), ASD] for i, s, p in rows]


HEADER = (
    "name,start_codon,start_position,expression,RBS_distance_bp,dG_total,dG_rRNA:mRNA,dG_mRNA,"
    "dG_spacing,dG_standby,dG_start_codon"
)


@settings(max_examples=100, deadline=None)
@given(
    st.dictionaries(ids, st.integers(min_value=1, max_value=51), max_size=12),
    st.lists(
        st.tuples(ids, st.integers(min_value=1, max_value=99),
                  st.floats(min_value=0, max_value=1e5, allow_nan=False),
                  st.integers(min_value=0, max_value=20)),
        max_size=30,
    ),
)  # fmt: skip
def test_parser_keeps_exactly_the_requested_positions(
    requested: dict[str, int], rows: list[tuple[str, int, float, int]]
) -> None:
    from gippa.annotate.ostir import parse_ostir

    # at most one row per (id, position): OSTIR scores each start once
    uniq: dict[tuple[str, int], tuple[float, int]] = {}
    for name, pos, tir, sp in rows:
        uniq.setdefault((name, pos), (tir, sp))
    lines = [
        f"{name},ATG,{pos},{tir!r},{sp},0.0,-1.0,-2.0,0.0,0.0,-1.194"
        for (name, pos), (tir, sp) in uniq.items()
    ]
    with TemporaryDirectory() as d:
        path = Path(d) / "ostir.csv"
        path.write_text("\n".join([HEADER, *lines]) + "\n")
        got = parse_ostir(path, requested)
    want = {i for i, p in requested.items() if (i, p) in uniq}
    assert set(got) == want
    for i in want:
        tir, sp = uniq[(i, requested[i])]
        assert got[i].tir_ostir == tir and got[i].spacer_bp == sp
