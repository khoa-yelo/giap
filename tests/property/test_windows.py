"""Hypothesis properties for promoter windows (issue M2-13 acceptance; GIAP_SPEC.md §7.1.4).

Module 1 tests: "window generation property tests: count, strand, 3'-end range, origin wrap on
circular replicons, no wrap on linear contigs". Plus: every window is 40 nt and equals a
brute-force oracle; ``+``/``-`` revcomp symmetry on the mirrored genome; dedupe/chunk invariants;
the layer reproduces every map row's score; ``promoters_bed`` has one record per window >= yes.
"""

from __future__ import annotations

import zlib

import pandas as pd
from helpers.promoter_oracle import api, oracle_windows, revcomp
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

SETTINGS = settings(max_examples=60, deadline=None)


def _gene(tag: str, contig: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


@st.composite
def cases(draw: st.DrawFn, topology: str) -> tuple[Gene, Replicon, str]:
    length = draw(st.integers(300, 1500))
    seq = draw(st.text("ACGT", min_size=length, max_size=length))
    glen = draw(st.integers(3, min(600, length - 1)))
    if topology == "circular":
        start = draw(st.integers(0, length - 1))  # may cross the origin (end > length)
    else:
        start = draw(st.integers(0, length - glen))
    strand = draw(st.sampled_from("+-"))
    rep = Replicon("c1", length, topology)  # type: ignore[arg-type]
    return _gene("G1", "c1", start, start + glen, strand), rep, seq


def _start_codon(g: Gene) -> int:
    """0-based first base of the start codon (``-``: E - 1), wrapped."""
    return g.interval.start if g.interval.strand == "+" else g.interval.end - 1


@SETTINGS
@given(cases("circular"))
def test_circular_count_strand_length_and_3p_range(case: tuple[Gene, Replicon, str]) -> None:
    g, rep, seq = case
    ws = api("gene_windows")(g, rep, seq)
    assert len(ws) == 221
    assert [w.k for w in ws] == list(range(-200, 21))
    n, s0 = rep.length, _start_codon(g)
    for w in ws:
        assert (w.contig, w.strand, w.owner) == ("c1", g.interval.strand, "G1")
        assert len(w.seq) == 40
        assert 1 <= w.end3p <= n
        # 3'-most base sits k bases downstream of the base just before the start codon
        sign = 1 if g.interval.strand == "+" else -1
        assert (w.end3p - 1 - (s0 - sign) - sign * w.k) % n == 0


@SETTINGS
@given(cases("circular"))
def test_circular_windows_equal_oracle(case: tuple[Gene, Replicon, str]) -> None:
    g, rep, seq = case
    ws = api("gene_windows")(g, rep, seq)
    iv = g.interval
    assert [(w.k, w.end3p, w.seq) for w in ws] == oracle_windows(
        iv.start, iv.end, iv.strand, seq, True
    )


@SETTINGS
@given(cases("linear"))
def test_linear_never_wraps(case: tuple[Gene, Replicon, str]) -> None:
    g, rep, seq = case
    ws = api("gene_windows")(g, rep, seq)
    iv = g.interval
    assert [(w.k, w.end3p, w.seq) for w in ws] == oracle_windows(
        iv.start, iv.end, iv.strand, seq, False
    )
    for w in ws:
        lo = w.end3p - 40 if w.strand == "+" else w.end3p - 1
        assert lo >= 0 and lo + 40 <= rep.length
        assert len(w.seq) == 40


@SETTINGS
@given(cases("circular"))
def test_revcomp_symmetry_with_mirrored_gene(case: tuple[Gene, Replicon, str]) -> None:
    g, rep, seq = case
    n, iv = rep.length, g.interval
    mirror_seq = revcomp(seq)
    other = "-" if iv.strand == "+" else "+"
    ms, me = n - iv.end, n - iv.start
    if ms < 0:  # origin-crossing gene: shift the mirrored interval into [0, n)
        ms, me = ms + n, me + n
    mg = _gene("G1", "c1", ms, me, other)
    a = api("gene_windows")(g, rep, seq)
    b = api("gene_windows")(mg, rep, mirror_seq)
    assert [w.seq for w in a] == [w.seq for w in b]
    assert [w.k for w in a] == [w.k for w in b]
    assert [n - w.end3p + 1 for w in a] == [w.end3p for w in b]


@SETTINGS
@given(cases("circular"))
def test_gene_windows_deterministic(case: tuple[Gene, Replicon, str]) -> None:
    g, rep, seq = case
    assert api("gene_windows")(g, rep, seq) == api("gene_windows")(g, rep, seq)


@st.composite
def window_lists(draw: st.DrawFn) -> list:
    window = api("Window")
    pool = draw(st.lists(st.text("ACGTN", min_size=40, max_size=40), min_size=1, max_size=12))
    n = draw(st.integers(1, 40))
    out = []
    for _ in range(n):
        out.append(
            window(
                draw(st.sampled_from(["c1", "c2"])),
                draw(st.sampled_from("+-")),
                draw(st.integers(1, 500)),
                draw(st.sampled_from(pool)),
                f"G{draw(st.integers(0, 5))}",
                draw(st.integers(-200, 20)),
            )
        )
    return out


@SETTINGS
@given(window_lists(), st.integers(1, 7))
def test_dedupe_and_chunk_invariants(ws: list, size: int) -> None:
    chunks, wmap = api("dedupe_and_chunk")(ws, chunk_size=size)
    records = [r for c in chunks for r in c]
    uniq = list(dict.fromkeys(w.seq for w in ws))  # first-seen order
    assert records == [(f"w{i}", s) for i, s in enumerate(uniq)]
    assert all(len(c) == size for c in chunks[:-1]) and 1 <= len(chunks[-1]) <= size
    assert list(wmap.columns) == ["header", "contig", "strand", "end3p", "owner", "k"]
    seq_of = dict(records)
    assert len(wmap) == len(ws)
    for w, r in zip(ws, wmap.itertuples(index=False), strict=True):
        assert seq_of[r.header] == w.seq
        assert (r.contig, r.strand, r.end3p, r.owner, r.k) == (
            w.contig, w.strand, w.end3p, w.owner, w.k,
        )  # fmt: skip


@SETTINGS
@given(cases("circular"), st.integers(0, 2**32 - 1))
def test_layer_carries_each_header_score(case: tuple[Gene, Replicon, str], salt: int) -> None:
    g, rep, seq = case
    chunks, wmap = api("dedupe_and_chunk")(api("gene_windows")(g, rep, seq))
    scores = {h: (zlib.crc32(f"{h}:{salt}".encode()) % 1000) / 1000 for c in chunks for h, _ in c}
    layer = api("promoter_layer")(wmap, scores)
    assert not layer.duplicated(["contig", "strand", "end3p"]).any()
    got = {(r.contig, r.strand, r.end3p): r.score for r in layer.itertuples(index=False)}
    assert set(got) == {(r.contig, r.strand, r.end3p) for r in wmap.itertuples(index=False)}
    for r in wmap.itertuples(index=False):
        assert got[(r.contig, r.strand, r.end3p)] == scores[r.header]


@st.composite
def score_frames(draw: st.DrawFn) -> pd.DataFrame:
    rows = draw(
        st.lists(
            st.tuples(st.sampled_from("+-"), st.integers(41, 1900), st.floats(0, 1)),
            min_size=0,
            max_size=25,
            unique_by=lambda t: (t[0], t[1]),
        )
    )
    return pd.DataFrame(
        [("c1", s, p, x) for s, p, x in rows], columns=["contig", "strand", "end3p", "score"]
    )


def _lo(end3p: int, strand: str) -> int:
    """Genome start of a window: ``+`` ``[p - 40, p)``, ``-`` ``[p - 1, p + 39)``."""
    return end3p - 40 if strand == "+" else end3p - 1


@SETTINGS
@given(score_frames(), st.sampled_from([0.3, 0.5, 0.8]))
def test_promoters_bed_one_record_per_window_at_or_above_yes(df: pd.DataFrame, yes: float) -> None:
    """Lead decision 3: windows are never merged; numbering follows §6.5 order."""
    rep = Replicon("c1", 2000, "linear")
    recs = api("promoters_bed")(df, {"c1": rep}, yes=yes)
    keep = df[df["score"] >= yes]
    want = sorted(
        (_lo(p, s), _lo(p, s) + 40, s, round(1000 * x))
        for s, p, x in zip(keep["strand"], keep["end3p"], keep["score"], strict=True)
    )
    got = [(r.interval.start, r.interval.end, r.interval.strand, r.score) for r in recs]
    assert sorted(got) == want
    by_name = sorted(recs, key=lambda r: int(r.name.removeprefix("promoter_")))
    assert [r.name for r in by_name] == [f"promoter_{i}" for i in range(1, len(recs) + 1)]
    keys = [(r.interval.start, r.interval.end, r.interval.strand != "+") for r in by_name]
    assert keys == sorted(keys)
