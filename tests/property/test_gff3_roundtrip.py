"""Hypothesis round trip for GFF3 write -> read (issue M1-05; GIAP_SPEC.md §6.5, §9.5)."""

from __future__ import annotations

import tempfile
from pathlib import Path

from hypothesis import given, settings
from hypothesis import strategies as st

from giap.io.gff3 import GffFeature, read_gff3, write_gff3
from giap.model import Interval, Replicon, interval_sort_key

FAST = settings(max_examples=150, deadline=None)

SEQID = st.from_regex(r"[A-Za-z0-9][A-Za-z0-9._-]{0,8}", fullmatch=True)
KEY = st.from_regex(r"[a-z][A-Za-z0-9_]{0,6}", fullmatch=True)
# printable ASCII and control characters, plus the GFF3 reserved characters ;=,&% (and %xx text)
VALUE = st.one_of(
    st.text(alphabet=st.characters(min_codepoint=1, max_codepoint=126), min_size=1, max_size=12),
    st.text(alphabet=";=,&% %3B%2Cab", min_size=1, max_size=12),
)
SCORE = st.one_of(
    st.none(),
    st.floats(allow_nan=False, allow_infinity=False).map(lambda x: float(f"{x:.6g}")),
)


@st.composite
def replicons(draw: st.DrawFn) -> list[Replicon]:
    ids = draw(st.lists(SEQID, min_size=1, max_size=3, unique=True))
    return [
        Replicon(i, draw(st.integers(1, 5000)), draw(st.sampled_from(["circular", "linear"])))
        for i in ids
    ]


@st.composite
def feature(draw: st.DrawFn, reps: list[Replicon], n: int) -> GffFeature:
    rep = draw(st.sampled_from(reps))
    start = draw(st.integers(0, rep.length - 1))
    max_len = rep.length if rep.topology == "circular" else rep.length - start
    end = start + draw(st.integers(1, max_len))
    extra = draw(st.dictionaries(KEY, VALUE, max_size=3))
    ftype = draw(st.sampled_from(["CDS", "gene", "tRNA", "ncRNA", "repeat_region"]))
    # The writer requires a phase on CDS (ADR-0005 item 1); other types may omit it.
    phase_st = st.integers(0, 2) if ftype == "CDS" else st.one_of(st.none(), st.integers(0, 2))
    return GffFeature(
        Interval(rep.id, start, end, draw(st.sampled_from(["+", "-", "."]))),
        draw(st.sampled_from(["giap", "Prodigal", "Bakta"])),
        ftype,
        score=draw(SCORE),
        phase=draw(phase_st),
        attributes={"ID": f"f{n}", **extra},
    )


@st.composite
def gff_input(draw: st.DrawFn) -> tuple[list[Replicon], list[GffFeature]]:
    reps = draw(replicons())
    n = draw(st.integers(0, 6))
    return reps, [draw(feature(reps, i)) for i in range(n)]


@FAST
@given(gff_input())
def test_gff3_write_read_round_trip(data: tuple[list[Replicon], list[GffFeature]]) -> None:
    reps, feats = data
    order = {r.id: i for i, r in enumerate(reps)}
    expected = sorted(
        feats, key=lambda f: interval_sort_key(f.interval, order, str(f.attributes["ID"]))
    )
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / "a.gff3", Path(d) / "b.gff3"
        write_gff3(a, feats, reps)
        write_gff3(b, list(reversed(feats)), reps)
        assert a.read_bytes() == b.read_bytes()
        got = read_gff3(a, {r.id: r for r in reps})
    assert got == expected
    for g, e in zip(got, expected, strict=True):
        assert list(g.attributes) == list(e.attributes)
