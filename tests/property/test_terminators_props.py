"""Hypothesis properties for the terminator layer (M2-15; GIPPA_SPEC.md §7.1.1, §7.1.5).

* parsers round-trip the real TerminatorNet (TAB, 11 columns, 1-based closed forward coordinates)
  and TransTermHP ``.tt`` formats into 0-based half-open intervals (confidence ``pred`` / CONF/100);
* the ``.coords`` writer is invertible (minus strand written start > end);
* ``crosscheck`` equals a brute-force oracle, is monotone in ``tol_bp`` and its ``all`` row sums;
* ``terminator_between`` equals a brute-force oracle on linear contigs, ignores BacTermFinder, does
  not depend on the order of calls, is monotone in ``yes`` and is mirror-symmetric (reflecting a
  ``+`` configuration onto the ``-`` strand gives the same answer);
* ``terminator_conf_downstream`` (no threshold) ≥ ``terminator_between`` whenever the latter is set.
"""

from __future__ import annotations

import random
from pathlib import Path
from tempfile import TemporaryDirectory

from hypothesis import assume, given, settings
from hypothesis import strategies as st

from gippa.annotate.terminators import (
    Terminator,
    crosscheck,
    parse_terminatornet,
    parse_transterm_tt,
    terminator_between,
    terminator_conf_downstream,
    write_transterm_coords,
)
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene

L = 5_000
LIN = Replicon("c", L, "linear")
REPS = {"c": LIN, "d": Replicon("d", L, "circular")}
SOURCES = ["TerminatorNet", "TransTermHP", "BacTermFinder"]


def _gene(tag: str, contig: str, start: int, end: int, strand: str, ft: str = "cds") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type=ft,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


@st.composite
def ivs(draw: st.DrawFn, contigs: tuple[str, ...] = ("c",), strands: str = "+-") -> Interval:
    s = draw(st.integers(0, L - 2))
    e = draw(st.integers(s + 1, min(L, s + 80)))
    return Interval(draw(st.sampled_from(contigs)), s, e, draw(st.sampled_from(strands)))  # type: ignore[arg-type]


confs = st.integers(0, 100).map(lambda k: k / 100)


@st.composite
def terms(draw: st.DrawFn, sources: list[str] | None = None, **kw: object) -> Terminator:
    src = draw(st.sampled_from(sources or ["TerminatorNet"]))
    return Terminator(draw(ivs(**kw)), draw(confs), src)  # type: ignore[arg-type]


# ------------------------------------------------------------------- parser round trips


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.tuples(ivs(contigs=("c", "d")), st.floats(0, 1, allow_nan=False), st.floats(-60, 0)),
        max_size=25,
    ),
    st.booleans(),
)
def test_terminatornet_round_trip(
    calls: list[tuple[Interval, float, float]], none_line: bool
) -> None:
    lines = []
    for i, (iv, pred, gh) in enumerate(calls, start=1):
        lines.append("\t".join(
            [iv.contig, f"TERM_{i}", "AAAAAAAA", "G" * iv.length_bp, "TTTTTTTT", "." * iv.length_bp,
             repr(gh), repr(pred), str(iv.start + 1), str(iv.end), iv.strand]
        ))  # fmt: skip
    if none_line:
        lines.append("d\tNO TERMINATORS FOUND")
    with TemporaryDirectory() as d:
        p = Path(d) / "t.tsv"
        p.write_text("\n".join(lines) + "\n")
        got = parse_terminatornet(p, REPS)
    want = [Terminator(iv, pred, "TerminatorNet", gh) for iv, pred, gh in calls]
    key = lambda t: (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand,  # noqa: E731
                     t.confidence, t.hairpin_dg)  # fmt: skip
    assert sorted(got, key=key) == sorted(want, key=key)


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.tuples(
            ivs(contigs=("c", "d")), st.integers(0, 100), st.sampled_from(["F", "T", "GF", "Tr"])
        ),
        max_size=25,
    )
)
def test_transterm_tt_round_trip(calls: list[tuple[Interval, int, str]]) -> None:
    out = [
        "TransTermHP v2.09 (built on Oct  8 2025)",
        "",
        "Genes are interspersed, and start the first column.",
        "",
    ]
    n = 0
    for contig in ("c", "d"):
        out += ["", f"SEQUENCE {contig}  (length {L})", ""]
        out.append(f"G_{contig}_1        10 - 300     + | ")
        for iv, conf, loc in (c for c in calls if c[0].contig == contig):
            n += 1
            a, b = (iv.start + 1, iv.end) if iv.strand == "+" else (iv.end, iv.start + 1)
            out.append(
                f"  TERM {n:<6}{a:>7} - {b:<7}  {iv.strand} {loc:<2} {conf:>4} -11.5 -3.8 | "
            )
            out.append("  AAAAAAAAAAAAAAA      GCCGCC TTTT GGCGGC     TTTTTTTTTTTTTTT")
    with TemporaryDirectory() as d:
        p = Path(d) / "t.tt"
        p.write_text("\n".join(out) + "\n")
        got = parse_transterm_tt(p, REPS)
    seen = sorted(
        (t.interval.contig, t.interval.start, t.interval.end, t.interval.strand,
         round(t.confidence * 100), t.source)
        for t in got
    )  # fmt: skip
    want = sorted(
        (iv.contig, iv.start, iv.end, iv.strand, conf, "TransTermHP") for iv, conf, _ in calls
    )
    assert seen == want


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.tuples(ivs(contigs=("c", "d")), st.sampled_from(["cds", "sorf", "pseudogene"])),
        max_size=30,
    )
)
def test_coords_writer_is_invertible(genes: list[tuple[Interval, str]]) -> None:
    gs = [
        _gene(f"T_{i:05d}", iv.contig, iv.start, iv.end, iv.strand, ft)
        for i, (iv, ft) in enumerate(genes)
    ]
    with TemporaryDirectory() as d:
        p = Path(d) / "g.coords"
        write_transterm_coords(p, gs)
        lines = p.read_text().splitlines()
    assert len(lines) == len(gs)
    back = {}
    for ln in lines:
        tag, a, b, contig = ln.split(" ")
        a1, b1 = int(a), int(b)
        back[tag] = Interval(contig, min(a1, b1) - 1, max(a1, b1), "+" if a1 < b1 else "-")
    for g in gs:
        iv = g.interval
        got = back[g.locus_tag]
        assert (got.contig, got.start, got.end) == (iv.contig, iv.start, iv.end)
        if iv.length_bp > 1:  # a 1-bp gene has start == end; strand is not encoded then
            assert got.strand == iv.strand


# ------------------------------------------------------------------- crosscheck


def _gap(a: Interval, b: Interval) -> int:
    return max(0, max(a.start, b.start) - min(a.end, b.end))


@settings(max_examples=80, deadline=None)
@given(
    st.lists(terms(contigs=("c", "d")), max_size=20),
    st.lists(terms(["TransTermHP"], contigs=("c", "d")), max_size=20),
    st.integers(0, 60),
)
def test_crosscheck_matches_oracle(prim: list[Terminator], oth: list[Terminator], tol: int) -> None:
    rows = crosscheck(prim, oth, tol_bp=tol)
    assert rows[-1].contig == "all"
    by = {r.contig: r for r in rows}
    for contig in {t.interval.contig for t in prim + oth}:
        p = [t for t in prim if t.interval.contig == contig]
        o = [t for t in oth if t.interval.contig == contig]
        m = sum(
            any(
                x.interval.strand == t.interval.strand and _gap(x.interval, t.interval) <= tol
                for x in o
            )
            for t in p
        )
        r = by[contig]
        assert (r.n_terminatornet, r.n_transtermhp, r.n_matched) == (len(p), len(o), m)
        if p:
            assert abs(r.share_matched - m / len(p)) < 1e-9
    per = [r for r in rows if r.contig != "all"]
    a = by["all"]
    assert a.n_terminatornet == sum(r.n_terminatornet for r in per) == len(prim)
    assert a.n_transtermhp == sum(r.n_transtermhp for r in per) == len(oth)
    assert a.n_matched == sum(r.n_matched for r in per)
    wider = {r.contig: r.n_matched for r in crosscheck(prim, oth, tol_bp=tol + 15)}
    assert all(wider[c] >= by[c].n_matched for c in by)


# ------------------------------------------------------------------- terminator_between


def _oracle(up: Gene, down: Gene, ts: list[Terminator], slack: int, yes: float) -> float | None:
    strand = up.interval.strand
    if strand == "+":
        lo, hi = up.interval.end - slack, down.interval.start
    else:
        lo, hi = down.interval.end, up.interval.start + slack
    vals = [
        t.confidence
        for t in ts
        if t.source != "BacTermFinder"
        and t.interval.contig == up.interval.contig
        and t.interval.strand == strand
        and lo <= t.interval.start
        and t.interval.end <= hi
        and t.confidence >= yes
    ]
    return max(vals) if vals else None


@st.composite
def pair(draw: st.DrawFn) -> tuple[Gene, Gene]:
    a = draw(st.integers(0, 2_000))
    e = draw(st.integers(a + 30, a + 1_500))
    s = draw(st.integers(max(a + 1, e - 200), e + 400))
    f = draw(st.integers(max(s + 30, e + 1), min(L, max(s + 30, e + 1) + 1_500)))
    assume(f <= L)
    return _gene("U", "c", a, e, "+"), _gene("D", "c", s, f, "+")


def _mirror_iv(iv: Interval) -> Interval:
    return Interval(iv.contig, L - iv.end, L - iv.start, "-" if iv.strand == "+" else "+")


def _mirror_gene(g: Gene) -> Gene:
    iv = _mirror_iv(g.interval)
    return _gene(g.locus_tag, iv.contig, iv.start, iv.end, iv.strand)


@settings(max_examples=150, deadline=None)
@given(
    pair(),
    st.lists(terms(SOURCES, contigs=("c",), strands="+-"), max_size=15),
    st.integers(0, 120),
    confs,
    st.randoms(use_true_random=False),
)
def test_between_oracle_order_mirror_and_monotone(
    ud: tuple[Gene, Gene], ts: list[Terminator], slack: int, yes: float, rnd: random.Random
) -> None:
    up, down = ud
    got = terminator_between(up, down, ts, LIN, slack_bp=slack, yes=yes)
    assert got == _oracle(up, down, ts, slack, yes)
    shuffled = ts[:]
    rnd.shuffle(shuffled)
    assert terminator_between(up, down, shuffled, LIN, slack_bp=slack, yes=yes) == got
    no_btf = [t for t in ts if t.source != "BacTermFinder"]
    assert terminator_between(up, down, no_btf, LIN, slack_bp=slack, yes=yes) == got
    mts = [Terminator(_mirror_iv(t.interval), t.confidence, t.source) for t in ts]
    mirrored = terminator_between(
        _mirror_gene(up), _mirror_gene(down), mts, LIN, slack_bp=slack, yes=yes
    )
    assert mirrored == got
    stricter = terminator_between(up, down, ts, LIN, slack_bp=slack, yes=min(1.0, yes + 0.2))
    assert stricter is None or (got is not None and stricter == got)
    if got is not None:
        dn = terminator_conf_downstream(up, down, ts, LIN, slack_bp=slack)
        assert dn is not None and dn >= got
