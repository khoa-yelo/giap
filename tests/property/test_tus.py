"""Hypothesis properties for operon pairs and TUs (issue M2-12; GIAP_SPEC.md §7.1.3).

Acceptance criteria: TUs partition the genes; every TU is single-strand; consecutive genes in a TU
have p ≥ threshold. Plus: TUs are maximal runs, numbering and gene order are canonical, results do
not depend on input order, ``build_pairs`` merges UniOP and the distance model pair by pair, the
distance model is a decreasing logistic, and ``tu_upstream`` equals a brute-force oracle on linear
contigs (window clipped at the nearest upstream gene, at least ``[S, S+50)``; spec Q2).
"""

from __future__ import annotations

import itertools
import math
import random
from collections.abc import Sequence

from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.annotate.operons import (
    adjacent_pairs,
    build_pairs,
    build_tus,
    distance_model_p,
    intergenic_distance,
    tu_upstream,
)
from giap.annotate.uniop_adapter import PairProb
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

FLAG = "operon_source_distance_model"


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
def genomes(draw: st.DrawFn, max_genes: int = 14) -> tuple[list[Gene], dict[str, Replicon]]:
    """1-2 contigs (circular or linear) with genes placed left to right (gaps or overlaps)."""
    reps: dict[str, Replicon] = {}
    genes: list[Gene] = []
    n_contigs = draw(st.integers(1, 2))
    k = 0
    for c in range(n_contigs):
        cid = f"c{c + 1}"
        pos = draw(st.integers(0, 400))
        n = draw(st.integers(0, max_genes))
        for _ in range(n):
            length = draw(st.integers(60, 1200))
            strand = draw(st.sampled_from("+-"))
            genes.append(_gene(f"G{k:03d}", cid, pos, pos + length, strand))
            k += 1
            pos = max(pos + 1, pos + length + draw(st.integers(-20, 400)))
        tail = draw(st.integers(1, 600))
        topology = draw(st.sampled_from(["circular", "linear"]))
        end = max((g.interval.end for g in genes if g.interval.contig == cid), default=0)
        reps[cid] = Replicon(cid, max(pos, end) + tail, topology)  # type: ignore[arg-type]
    return genes, reps


def _probs(draw: st.DrawFn, n: int) -> list[float]:
    return draw(
        st.lists(
            st.one_of(st.floats(0, 1, allow_nan=False), st.sampled_from([0.0, 0.5, 1.0])),
            min_size=n,
            max_size=n,
        )
    )


def _ring_fully_linked(
    genes: Sequence[Gene], reps: dict[str, Replicon], linked: set[tuple[str, str]]
) -> bool:
    for cid, rep in reps.items():
        on = [g for g in genes if g.interval.contig == cid]
        if rep.topology == "circular" and len(on) >= 2:
            n_links = sum(1 for a, b in linked if a in {g.locus_tag for g in on})
            if n_links >= len(on):
                return True
    return False


@settings(max_examples=150, deadline=None)
@given(g=genomes(), data=st.data(), threshold=st.sampled_from([0.3, 0.5, 0.7]))
def test_tus_partition_single_strand_linked_and_maximal(
    g: tuple[list[Gene], dict[str, Replicon]], data: st.DataObject, threshold: float
) -> None:
    genes, reps = g
    adj = adjacent_pairs(genes, reps)
    ps = _probs(data.draw, len(adj))
    pairs = [
        PairProb(a.locus_tag, b.locus_tag, p, "uniop") for (a, b), p in zip(adj, ps, strict=True)
    ]
    linked = {(p.gene_a, p.gene_b) for p in pairs if p.uniop_p >= threshold}
    assume(not _ring_fully_linked(genes, reps, linked))
    p_of = {(p.gene_a, p.gene_b): p.uniop_p for p in pairs}
    by_tag = {x.locus_tag: x for x in genes}

    tus = build_tus(genes, pairs, threshold, reps)

    # partition
    members = [t for tu in tus for t in tu.genes]
    assert sorted(members) == sorted(by_tag)
    tu_of = {t: tu.tu_id for tu in tus for t in tu.genes}
    # numbering TU_1..TU_n in order
    assert [t.tu_id for t in tus] == [f"TU_{i}" for i in range(1, len(tus) + 1)]
    for tu in tus:
        strands = {by_tag[t].interval.strand for t in tu.genes}
        contigs = {by_tag[t].interval.contig for t in tu.genes}
        assert strands == {tu.strand} and contigs == {tu.contig}
        assert tu.interval.strand == tu.strand and tu.interval.contig == tu.contig
        for x, y in zip(tu.genes, tu.genes[1:], strict=False):
            # transcription order: genome order on +, reverse genome order on -
            assert p_of[(x, y) if tu.strand == "+" else (y, x)] >= threshold
        # interval covers every gene (unwrapped across the origin)
        lo = tu.interval.start
        assert by_tag[tu.genes[0 if tu.strand == "+" else -1]].interval.start == lo
        assert tu.interval.end >= max(by_tag[t].interval.end for t in tu.genes)
        assert tu.interval.length_bp <= reps[tu.contig].length or len(tu.genes) == 1
    # maximal: every linked pair is in one TU
    for a, b in linked:
        assert tu_of[a] == tu_of[b]
    # numbering by contig (FASTA order), then start, then strand
    order = list(reps)
    keys = [(order.index(t.contig), t.interval.start, t.strand == "-") for t in tus]
    assert keys == sorted(keys)


@settings(max_examples=60, deadline=None)
@given(g=genomes(), data=st.data(), seed=st.integers(0, 2**32 - 1))
def test_results_do_not_depend_on_input_order(
    g: tuple[list[Gene], dict[str, Replicon]], data: st.DataObject, seed: int
) -> None:
    genes, reps = g
    shuffled = list(genes)
    random.Random(seed).shuffle(shuffled)
    assert adjacent_pairs(shuffled, reps) == adjacent_pairs(genes, reps)
    adj = adjacent_pairs(genes, reps)
    ps = _probs(data.draw, len(adj))
    pairs = [
        PairProb(a.locus_tag, b.locus_tag, p, "uniop") for (a, b), p in zip(adj, ps, strict=True)
    ]
    assume(
        not _ring_fully_linked(
            genes, reps, {(p.gene_a, p.gene_b) for p in pairs if p.uniop_p >= 0.5}
        )
    )
    pairs_shuffled = list(pairs)
    random.Random(seed + 1).shuffle(pairs_shuffled)
    assert build_tus(shuffled, pairs_shuffled, 0.5, reps) == build_tus(genes, pairs, 0.5, reps)


@settings(max_examples=150, deadline=None)
@given(g=genomes())
def test_adjacent_pairs_invariants(g: tuple[list[Gene], dict[str, Replicon]]) -> None:
    genes, reps = g
    adj = adjacent_pairs(genes, reps)
    as_, bs = [a.locus_tag for a, _ in adj], [b.locus_tag for _, b in adj]
    assert len(set(as_)) == len(as_) and len(set(bs)) == len(bs)
    for a, b in adj:
        assert a.interval.contig == b.interval.contig
        assert a.interval.strand == b.interval.strand
        assert a.locus_tag != b.locus_tag
    for cid, rep in reps.items():
        on = sorted(
            (x for x in genes if x.interval.contig == cid),
            key=lambda x: (x.interval.start, x.interval.end),
        )
        fwd = [
            (x.locus_tag, y.locus_tag)
            for x, y in itertools.pairwise(on)
            if x.interval.strand == y.interval.strand
        ]
        if (
            rep.topology == "circular"
            and len(on) >= 2
            and on[-1].interval.strand == on[0].interval.strand
        ):
            fwd.append((on[-1].locus_tag, on[0].locus_tag))
        got = [(a.locus_tag, b.locus_tag) for a, b in adj if a.interval.contig == cid]
        assert got == fwd


@settings(max_examples=100, deadline=None)
@given(g=genomes(), data=st.data(), fail=st.booleans())
def test_build_pairs_merges_uniop_and_distance_model(
    g: tuple[list[Gene], dict[str, Replicon]], data: st.DataObject, fail: bool
) -> None:
    genes, reps = g
    adj = adjacent_pairs(genes, reps)
    reported = data.draw(st.lists(st.booleans(), min_size=len(adj), max_size=len(adj)))
    ps = _probs(data.draw, len(adj))
    uni = [
        PairProb(a.locus_tag, b.locus_tag, p, "uniop")
        for (a, b), p, r in zip(adj, ps, reported, strict=True)
        if r
    ]
    pairs, flags = build_pairs(genes, reps, None if fail else uni, 50.0, 20.0)
    assert [(p.gene_a, p.gene_b) for p in pairs] == [(a.locus_tag, b.locus_tag) for a, b in adj]
    failed = fail or not uni
    if adj:  # with no adjacent pair an empty uniop.pred is not distinguishable from success
        assert flags == ((FLAG,) if failed else ())
    by_key = {(u.gene_a, u.gene_b): u for u in uni}
    for (a, b), p in zip(adj, pairs, strict=True):
        key = (a.locus_tag, b.locus_tag)
        if not failed and key in by_key:
            assert p == by_key[key]
        else:
            assert p.source == "distance_model"
            d = intergenic_distance(a, b, reps[a.interval.contig])
            assert math.isclose(p.uniop_p, distance_model_p(d, 50.0, 20.0), rel_tol=1e-12)


@settings(max_examples=200, deadline=None)
@given(
    d=st.integers(-(10**6), 10**6),
    x=st.integers(0, 10**4),
    d0=st.floats(-100, 200),
    s=st.floats(1, 100),
)
def test_distance_model_logistic(d: int, x: int, d0: float, s: float) -> None:
    p = distance_model_p(d, d0, s)
    assert 0.0 <= p <= 1.0
    assert distance_model_p(d + x, d0, s) <= p
    lo, hi = distance_model_p(round(d0) - x, d0, s), distance_model_p(round(d0) + x, d0, s)
    assert lo >= hi
    if abs(d - d0) / s < 30:
        assert math.isclose(p, 1.0 / (1.0 + math.exp((d - d0) / s)), rel_tol=1e-9)
        assert math.isclose(p + distance_model_p(2 * d0 - d, d0, s), 1.0, rel_tol=1e-9)


@settings(max_examples=200, deadline=None)
@given(
    length=st.integers(2_000, 20_000),
    data=st.data(),
    bp=st.sampled_from([(-300, 50), (-100, 20), (-500, 0)]),
)
def test_tu_upstream_matches_linear_oracle(
    length: int, data: st.DataObject, bp: tuple[int, int]
) -> None:
    rep = Replicon("lin", length, "linear")
    start = data.draw(st.integers(0, length - 100))
    end = data.draw(st.integers(start + 60, min(length, start + 2000)))
    strand = data.draw(st.sampled_from("+-"))
    first = _gene("F", "lin", start, end, strand)
    others = []
    for i in range(data.draw(st.integers(0, 6))):
        s0 = data.draw(st.integers(0, length - 60))
        e0 = data.draw(st.integers(s0 + 60, min(length, s0 + 1500)))
        others.append(_gene(f"O{i}", "lin", s0, e0, data.draw(st.sampled_from("+-"))))
    up, down = -bp[0], bp[1]
    if strand == "+":
        sc = start
        lo = max([sc - up, 0, *(min(o.interval.end, sc) for o in others if o.interval.start < sc)])
        hi = min(sc + down, length)
    else:
        sc = end
        hi = min(
            [sc + up, length, *(max(o.interval.start, sc) for o in others if o.interval.end > sc)]
        )
        lo = max(sc - down, 0)
    assume(lo < hi)
    assert tu_upstream(first, others, rep, bp=bp) == Interval("lin", lo, hi, strand)  # type: ignore[arg-type]
