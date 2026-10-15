"""Hypothesis properties of tblastn HSP chaining (issue M3-13; GIAP_SPEC.md §7.2.3, §6.5).

- Merging is invariant to the input order of HSPs (acceptance criterion).
- Merged coverage <= 1, identity in [0, 1], ``q`` range inside the protein (acceptance criterion).
- Every merged interval is canonical (``0 <= start < L``, unwrapped ``end <= start + L``; linear
  contigs never wrap), so ``assign_alignments`` accepts the output (M3-12 review decision 4).
- Ranks per protein are exactly ``1..k``, ordered by non-increasing summed bitscore (the summed
  bitscore of a chain is not observable, so this is checked on generated chain layouts).
- A protein cut into collinear pieces with gaps <= 300 nt and no query overlap gives one chain
  covering the whole span, on either strand and across the origin of a circular replicon;
  pieces > 300 nt apart on a linear contig stay separate chains.
- Chains never mix proteins, contigs or strands.
"""

from __future__ import annotations

import pytest
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.config.model import ProteinsConfig
from giap.model import Interval, Replicon
from giap.model.coords import check_interval
from giap.model.genes import GeneModel

L = 30_000
CIRC = Replicon("chr1", L, "circular", "chromosome")
LIN = Replicon("chr1", L, "linear", "unknown")
PL = Replicon("pl1", 4_000, "circular", "plasmid")


def _tb():  # type: ignore[no-untyped-def]
    from giap.mapping import tblastn

    return tblastn


def _hsp(pid, contig, q0, q1, s_lo, s_hi, strand, qlen, nident, bits):  # type: ignore[no-untyped-def]
    """1-based closed HSP; ``s_lo``/``s_hi`` are 1-based genome bounds, flipped on minus."""
    n = q1 - q0 + 1
    sstart, send = (s_lo, s_hi) if strand == "+" else (s_hi, s_lo)
    return _tb().Hsp(
        qseqid=pid,
        sseqid=contig,
        qstart=q0,
        qend=q1,
        sstart=sstart,
        send=send,
        length=n,
        nident=min(n, nident),
        positive=min(n, nident),
        bitscore=float(bits),
        qlen=qlen,
    )


@st.composite
def random_hsps(draw):  # type: ignore[no-untyped-def]
    """Arbitrary HSPs (not necessarily collinear) on chr1 (circular or linear) and pl1."""
    linear = draw(st.booleans())
    reps = {"chr1": LIN if linear else CIRC, "pl1": PL}
    out = []
    for _ in range(draw(st.integers(0, 12))):
        pid = draw(st.sampled_from(["P1", "P2", "sp|Q1|A_B"]))
        contig = draw(st.sampled_from(["chr1", "chr1", "pl1"]))
        size = reps[contig].length
        qlen = {"P1": 120, "P2": 90, "sp|Q1|A_B": 60}[pid]
        q0 = draw(st.integers(1, qlen))
        q1 = draw(st.integers(q0, qlen))
        nt = 3 * (q1 - q0 + 1)
        s_lo = draw(st.integers(1, size - nt + 1))
        strand = draw(st.sampled_from(["+", "-"]))
        nident = draw(st.integers(0, q1 - q0 + 1))
        bits = draw(st.integers(20, 300))
        out.append(_hsp(pid, contig, q0, q1, s_lo, s_lo + nt - 1, strand, qlen, nident, bits))
    return reps, out


def _key(a):  # type: ignore[no-untyped-def]
    iv = a.interval
    return (a.protein_id, a.rank, iv.contig, iv.start, iv.end, iv.strand, a.q_start, a.q_end)


@settings(max_examples=150, deadline=None)
@given(random_hsps(), st.randoms(use_true_random=False))
def test_input_order_invariance(data, rnd) -> None:  # type: ignore[no-untyped-def]
    reps, hsps = data
    ref = sorted(_tb().merge_hsps(hsps, reps), key=_key)
    shuffled = list(hsps)
    rnd.shuffle(shuffled)
    assert sorted(_tb().merge_hsps(shuffled, reps), key=_key) == ref


@settings(max_examples=150, deadline=None)
@given(random_hsps())
def test_bounds_canonical_and_ranks(data) -> None:  # type: ignore[no-untyped-def]
    reps, hsps = data
    alns = _tb().merge_hsps(hsps, reps)
    assert len(alns) <= len(hsps)
    assert {a.protein_id for a in alns} == {h.qseqid for h in hsps}
    hsp_strands = {(h.qseqid, h.sseqid, "+" if h.sstart <= h.send else "-") for h in hsps}
    for a in alns:
        rep = reps[a.interval.contig]
        check_interval(a.interval, rep)
        assert 0 <= a.interval.start < rep.length
        assert 0.0 < a.query_coverage <= 1.0
        assert 0.0 <= a.identity <= 1.0
        assert 0 <= a.q_start < a.q_end <= a.protein_len_aa
        assert a.route == "tblastn"
        assert a.frameshift == 0 and a.stop_codon == 0
        assert (a.protein_id, a.interval.contig, a.interval.strand) in hsp_strands
    for pid in {a.protein_id for a in alns}:
        ranks = sorted(a.rank for a in alns if a.protein_id == pid)
        assert ranks == list(range(1, len(ranks) + 1))
    # canonical output is accepted by M3-12
    from giap.mapping.protein_assign import assign_alignments

    assign_alignments(alns, GeneModel(reps, (), ()), ProteinsConfig())


@st.composite
def collinear_pieces(draw):  # type: ignore[no-untyped-def]
    """A protein cut into k pieces; gaps in [0, 300] nt; placed anywhere (may cross the origin)."""
    k = draw(st.integers(1, 5))
    sizes = draw(st.lists(st.integers(10, 60), min_size=k, max_size=k))
    gaps = draw(st.lists(st.integers(0, 300), min_size=k - 1, max_size=k - 1))
    strand = draw(st.sampled_from(["+", "-"]))
    span = 3 * sum(sizes) + sum(gaps)
    start0 = draw(st.integers(0, L - 1))  # unwrapped start; may cross the origin of CIRC
    return sizes, gaps, strand, start0, span


def _pieces_to_hsps(sizes, gaps, strand, start0, pid="P1"):  # type: ignore[no-untyped-def]
    qlen = sum(sizes)
    blocks = []  # (q0, q1, g0) with g0 the unwrapped 0-based genome start
    q, g = 1, start0
    for i, n in enumerate(sizes):
        blocks.append((q, q + n - 1, g))
        q += n
        g += 3 * n + (gaps[i] if i < len(gaps) else 0)
    end0 = g
    hsps = []
    for q0, q1, g0 in blocks:
        nt = 3 * (q1 - q0 + 1)
        if strand == "-":  # mirror: N-terminus at the high unwrapped coordinate
            g0 = start0 + (end0 - (g0 + nt))
        hsps.append((q0, q1, g0, nt))
    return qlen, hsps, end0


@settings(max_examples=200, deadline=None)
@given(collinear_pieces(), st.randoms(use_true_random=False))
def test_collinear_pieces_make_one_chain(data, rnd) -> None:  # type: ignore[no-untyped-def]
    sizes, gaps, strand, start0, span = data
    qlen, blocks, _end0 = _pieces_to_hsps(sizes, gaps, strand, start0)
    hsps = []
    for q0, q1, g0, nt in blocks:
        lo = g0 % L
        assume(lo + nt <= L)  # BLAST never reports an HSP spanning the origin
        hsps.append(_hsp("P1", "chr1", q0, q1, lo + 1, lo + nt, strand, qlen, q1 - q0 + 1, 50))
    rnd.shuffle(hsps)
    (a,) = _tb().merge_hsps(hsps, {"chr1": CIRC})
    assert a.interval == Interval("chr1", start0 % L, start0 % L + span, strand)
    assert a.query_coverage == 1.0
    assert a.identity == 1.0
    assert a.rank == 1


@settings(max_examples=150, deadline=None)
@given(
    st.lists(st.integers(10, 60), min_size=2, max_size=5),
    st.integers(301, 2000),
    st.sampled_from(["+", "-"]),
)
def test_far_pieces_on_linear_contig_stay_separate(sizes, gap, strand) -> None:  # type: ignore[no-untyped-def]
    gaps = [gap] * (len(sizes) - 1)
    qlen, blocks, end0 = _pieces_to_hsps(sizes, gaps, strand, 100)
    assume(end0 <= L)
    hsps = [
        _hsp("P1", "chr1", q0, q1, g0 + 1, g0 + nt, strand, qlen, q1 - q0 + 1, 50)
        for q0, q1, g0, nt in blocks
    ]
    alns = _tb().merge_hsps(hsps, {"chr1": LIN})
    assert len(alns) == len(sizes)
    assert sum(a.query_coverage for a in alns) == pytest.approx(1.0)


@settings(max_examples=100, deadline=None)
@given(st.lists(st.integers(20, 400), min_size=1, max_size=6, unique=True))
def test_rank_follows_bitscore(bits) -> None:  # type: ignore[no-untyped-def]
    """Single-HSP loci far apart: rank 1 has the highest bitscore, and so on."""
    hsps = [
        _hsp("P1", "chr1", 1, 100, 1 + 2_000 * i, 300 + 2_000 * i, "+", 100, 100, b)
        for i, b in enumerate(bits)
    ]
    alns = _tb().merge_hsps(hsps, {"chr1": CIRC})
    by_rank = sorted(alns, key=lambda a: a.rank)
    expected_starts = [2_000 * i for i, _ in sorted(enumerate(bits), key=lambda t: -t[1])]
    assert [a.interval.start for a in by_rank] == expected_starts
