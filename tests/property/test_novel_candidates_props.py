"""Properties of novel-ORF reconstruction on random genomes (issue M7-06; GIAP_SPEC.md §7.6.2).

Acceptance properties of docs/specs/M7-06-novel-candidates-and-orfs.md:

- mirroring the genome maps every ORF to its mirror (strand, stop, starts) and its ID to the
  mirrored ID;
- every peptide of an ORF lies between ``start_longest`` (or the first peptide) and the stop, in
  the ORF's frame, with no in-frame stop between start and stop;
- a planted ORF with >= 2 peptides cut from it is reconstructed with the planted stop, and its
  planted start is in ``starts``;
- the second-round windows are M2-13's ``gene_windows`` at every candidate start (one generalised
  function), deduplicated by ``(contig, strand, end3p)``;

plus agreement with the brute-force oracle of tests/helpers/novel_oracle.py (rules 1-6: walks,
grouping by ``(contig, strand, stop)``, starts, open ends, IDs), the OSTIR row span of rule 8,
invariance to the input row order, and that inputs are never mutated. Random circular genomes
(origin-crossing peptides included) and linear contigs (open ends) on both strands.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
from helpers.novel_oracle import (
    START_CODONS,
    STOPS,
    Case,
    cds,
    hit_row,
    il,
    interval_nt,
    mirror_case,
    mirror_codon,
    mirror_interval,
    oracle_orfs,
    orf_view,
    revcomp,
    translate,
    wrapped,
)
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

from giap.annotate.promoters import gene_windows
from giap.config.model import GiapConfig
from giap.model.coords import Replicon
from giap.novel.candidates import (
    Orf,
    candidate_rows,
    ostir_rows,
    promoter_windows,
    reconstruct,
)

PROPS = settings(max_examples=120, deadline=None, suppress_health_check=[HealthCheck.too_slow])
CFG = GiapConfig(seed=1)
CONTIG = "chr1"
ANTI_SD = "ACCTCCTTA"


@st.composite
def scenarios(draw: st.DrawFn) -> Case:
    """A random contig and 1-6 stop-free peptide windows (intergenic, no gene model)."""
    seed = draw(st.integers(0, 2**32 - 1))
    n = draw(st.integers(500, 1200))
    topo = draw(st.sampled_from(["circular", "linear"]))
    k = draw(st.integers(1, 6))
    rng = np.random.default_rng(seed)
    seq = "".join(rng.choice(list("ACGT"), n))
    hits: list[dict[str, Any]] = []
    for _ in range(k):
        for _attempt in range(60):
            strand = "+" if rng.random() < 0.5 else "-"
            width = 3 * int(rng.integers(6, 13))
            if topo == "circular" and rng.random() < 0.3:  # near or across the origin
                s = int(rng.integers(n - 40, n))
            else:
                s = int(rng.integers(0, n if topo == "circular" else n - width + 1))
            aa = translate(interval_nt(seq, s, s + width, strand))
            if "*" in aa or "X" in aa:
                continue
            hits.append(hit_row(aa, CONTIG, s, s + width, strand, "intergenic", n))
            break
    assume(hits)
    return Case({CONTIG: seq}, {CONTIG: Replicon(CONTIG, n, topo)}, hits)  # type: ignore[arg-type]


def _run(case: Case) -> list[Orf]:
    rows = candidate_rows(case.hits_df(), case.peptides_df(), None)
    return reconstruct(rows, case.gm, case.genome, CFG)


def _oracle_or_raise(case: Case) -> list[Any] | None:
    """The oracle's ORFs, or None (skipped) when a walk exceeds (L - 40) / 2 nt, so an ORF
    interval is always shorter than the replicon (the whole-turn boundary and ORFs re-entering
    their own peptides in another phase are left to the unit tests and the implementation)."""
    n = case.replicons[CONTIG].length
    try:
        return oracle_orfs(case.hits, case.genome, case.replicons, max_walk=(n - 40) // 2)
    except ValueError:
        return None


def _tseq(case: Case, strand: str) -> str:
    seq = case.genome[CONTIG]
    return seq if strand == "+" else revcomp(seq)


def _t_codon(x: int, strand: str, n: int) -> int:
    """Transcription coordinate of the codon whose lowest genome coordinate is ``x``."""
    return x if strand == "+" else mirror_codon(x, n)


def _t_interval(h: dict[str, Any], n: int) -> tuple[int, int]:
    return (
        (h["start"], h["end"]) if h["strand"] == "+" else mirror_interval(h["start"], h["end"], n)
    )


# ---------------------------------------------------------------- oracle agreement


@PROPS
@given(scenarios())
def test_reconstruct_matches_the_oracle(case: Case) -> None:
    expected = _oracle_or_raise(case)
    assume(expected is not None)
    got = _run(case)
    assert sorted(map(orf_view, got), key=repr) == sorted(map(orf_view, expected), key=repr)
    assert len({o.orf_id for o in got}) == len(got)


# ---------------------------------------------------------------- mirror


def _five_prime_high(o: Orf, case: Case) -> int:
    """Exclusive genome end of a minus ORF's longest interval (its 5' side)."""
    n = case.replicons[CONTIG].length
    if o.start_longest is not None:
        return o.start_longest + 3
    ends = [h["end"] for h in case.hits if il(h["peptide"]) in o.peptides and h["strand"] == "-"]
    if o.end_3p_open:
        return max(ends)
    return max(ends, key=lambda e: (e - o.stop) % n)


@PROPS
@given(scenarios())
def test_mirror_maps_every_orf_to_its_mirror(case: Case) -> None:
    assume(_oracle_or_raise(case) is not None)
    n = case.replicons[CONTIG].length
    flip = {"+": "-", "-": "+"}
    got = {tuple(o.peptides): o for o in _run(case)}
    mir = {tuple(o.peptides): o for o in _run(mirror_case(case))}
    assert got.keys() == mir.keys()
    for peps, o in got.items():
        m = mir[peps]
        assert m.strand == flip[o.strand]
        assert m.starts == tuple(mirror_codon(x, n) for x in o.starts)
        assert (m.end_5p_open, m.end_3p_open) == (o.end_5p_open, o.end_3p_open)
        if o.end_3p_open:
            continue
        assert m.stop == mirror_codon(o.stop, n)
        plus = o.strand == "+"
        lowest = mirror_codon(o.stop, n) if plus else (n - _five_prime_high(o, case)) % n
        assert m.orf_id == f"NOVEL_{CONTIG}_{lowest + 1}_{'m' if o.strand == '+' else 'p'}"


# ---------------------------------------------------------------- containment and frame


@PROPS
@given(scenarios())
def test_peptides_lie_in_frame_between_start_and_stop(case: Case) -> None:
    assume(_oracle_or_raise(case) is not None)
    n = case.replicons[CONTIG].length
    circ = case.replicons[CONTIG].topology == "circular"
    orfs = _run(case)

    for o in orfs:
        tseq = _tseq(case, o.strand)
        rows = [h for h in case.hits if il(h["peptide"]) in o.peptides]
        assert rows and {il(h["peptide"]) for h in rows} == set(o.peptides)
        assert {h["strand"] for h in rows} == {o.strand}
        tiv = [_t_interval(h, n) for h in rows]
        if o.end_3p_open:  # t3 = end of the last whole codon of the frame (linear only)
            assert not circ
            te0 = tiv[0][1]
            t3 = te0 + 3 * ((n - te0) // 3)
        else:  # t3 = first base of the stop codon
            t3 = _t_codon(o.stop, o.strand, n)
            assert wrapped(tseq, t3, t3 + 3) in STOPS

        def dist(tx: int, t3: int = t3) -> int:
            return (t3 - tx) % n if circ else t3 - tx

        dists = [dist(ts) for ts, _ in tiv]
        for (ts, te), d in zip(tiv, dists, strict=True):
            assert d % 3 == 0  # the ORF's frame
            assert d >= te - ts  # the peptide ends at or before the stop
        if o.start_longest is not None:
            d5 = dist(_t_codon(o.start_longest, o.strand, n))
            assert wrapped(tseq, t3 - d5, t3 - d5 + 3) in START_CODONS
            assert d5 >= max(dists)  # start_longest is upstream of every peptide
            ds = [dist(_t_codon(x, o.strand, n)) for x in o.starts]
            assert ds == sorted(ds, reverse=True)  # 5'->3'
            assert ds[0] == d5
            for dx in ds:
                assert dx % 3 == 0 and dx >= max(dists)
                assert wrapped(tseq, t3 - dx, t3 - dx + 3) in START_CODONS
        else:
            assert o.starts == ()
            d5 = max(dists)
        for p in range(t3 - d5, t3, 3):  # no in-frame stop from the start (or first peptide)
            assert wrapped(tseq, p, p + 3) not in STOPS
    loci: dict[str, set[tuple[Any, ...]]] = {}
    for h in case.hits:
        loci.setdefault(h["peptide"], set()).add((h["start"], h["end"], h["strand"]))
    kept = {il(p) for p, ls in loci.items() if len(ls) == 1}  # a repeated row is one locus
    assert {p for o in orfs for p in o.peptides} == kept  # every retained peptide is in an ORF
    assert sum(len(o.peptides) for o in orfs) == len(kept)  # ... and in exactly one


# ---------------------------------------------------------------- planted ORF


@st.composite
def planted(draw: st.DrawFn) -> tuple[Case, int, int]:
    seed = draw(st.integers(0, 2**32 - 1))
    n = draw(st.integers(600, 1500))
    topo = draw(st.sampled_from(["circular", "linear"]))
    strand = draw(st.sampled_from(["+", "-"]))
    m = draw(st.integers(15, 60))  # codons after the ATG, before the stop
    rng = np.random.default_rng(seed)
    sense = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]
    orf = "ATG" + "".join(rng.choice(sense, m)) + "TAA"
    lo = 0 if topo == "circular" else 3
    hi = n if topo == "circular" else n - len(orf) - 3
    p0 = int(rng.integers(lo, hi))
    t = list("".join(rng.choice(list("ACGT"), n)))
    for k, c in enumerate(orf):
        t[(p0 + k) % n] = c
    tseq = "".join(t)
    seq = tseq if strand == "+" else revcomp(tseq)
    k = draw(st.integers(2, 4))
    hits = []
    for _ in range(k):
        w = int(rng.integers(6, 13))
        i = int(rng.integers(0, m + 1 - w + 1))  # codon index, the stop is codon m + 1
        ts, te = p0 + 3 * i, p0 + 3 * (i + w)
        pep = translate(orf[3 * i : 3 * (i + w)])
        if strand == "+":
            s, e = ts % n, ts % n + (te - ts)
        else:
            s, e = mirror_interval(ts % n, ts % n + (te - ts), n)
        hits.append(hit_row(pep, CONTIG, s, e, strand, "intergenic", n))
    peps = [h["peptide"] for h in hits]
    assume(len({il(p) for p in peps}) >= 2 and len(set(peps)) == len(peps))
    tstop = p0 + 3 * (m + 1)
    if strand == "+":
        stop, start = tstop % n, p0 % n
    else:
        stop, start = mirror_codon(tstop % n, n), mirror_codon(p0 % n, n)
    case = Case({CONTIG: seq}, {CONTIG: Replicon(CONTIG, n, topo)}, hits)  # type: ignore[arg-type]
    return case, stop, start


@PROPS
@given(planted())
def test_planted_orf_is_reconstructed(args: tuple[Case, int, int]) -> None:
    case, stop, start = args
    orfs = _run(case)
    hosts = [o for o in orfs if (o.strand, o.stop) == (case.hits[0]["strand"], stop)]
    assert len(hosts) == 1
    (o,) = hosts
    assert set(o.peptides) == {il(h["peptide"]) for h in case.hits}
    assert start in o.starts
    assert len(o.peptides) >= 2


# ---------------------------------------------------------------- windows and OSTIR rows


@PROPS
@given(scenarios())
def test_windows_are_gene_windows_at_every_start(case: Case) -> None:
    assume(_oracle_or_raise(case) is not None)
    rep = case.replicons[CONTIG]
    seq = case.genome[CONTIG]
    orfs = _run(case)
    got = promoter_windows(orfs, case.gm, case.genome, CFG)
    keys = [(w.contig, w.strand, w.end3p) for w in got]
    assert len(keys) == len(set(keys))
    expected: dict[tuple[str, str, int], str] = {}
    for o in orfs:
        for x in o.starts:
            g = cds("S", CONTIG, x, x + 3, o.strand)
            for w in gene_windows(g, rep, seq, (-200, 20)):
                expected.setdefault((w.contig, w.strand, w.end3p), w.seq)
    assert {k: w.seq for k, w in zip(keys, got, strict=True)} == expected
    owners = {o.orf_id for o in orfs if o.starts}
    assert {w.owner for w in got} <= owners


@PROPS
@given(scenarios())
def test_ostir_rows_span_the_candidate_starts(case: Case) -> None:
    assume(_oracle_or_raise(case) is not None)
    rep = case.replicons[CONTIG]
    n, circ = rep.length, rep.topology == "circular"
    orfs = _run(case)
    rows = {r.id: r for r in ostir_rows(orfs, case.gm, case.genome, ANTI_SD, CFG)}
    assert rows.keys() == {o.orf_id for o in orfs if o.starts}
    for o in orfs:
        if not o.starts:
            continue
        r = rows[o.orf_id]
        tseq = _tseq(case, o.strand)
        t0 = _t_codon(o.starts[0], o.strand, n)
        d = (_t_codon(o.starts[-1], o.strand, n) - t0) % n if circ else (
            _t_codon(o.starts[-1], o.strand, n) - t0)  # fmt: skip
        if circ:
            assert r.seq == wrapped(tseq, t0 - 50, t0 + d + 50)
            assert (r.start, r.end) == (51, 51 + d)
        else:
            lo, hi = max(t0 - 50, 0), min(t0 + d + 50, n)
            assert r.seq == tseq[lo:hi]
            assert (r.start, r.end) == (t0 - lo + 1, t0 + d - lo + 1)
        assert r.seq[r.start - 1 : r.start + 2] in START_CODONS
        assert r.seq[r.end - 1 : r.end + 2] in START_CODONS
        assert r.anti_sd == ANTI_SD


# ---------------------------------------------------------------- order invariance and purity


@PROPS
@given(scenarios(), st.randoms(use_true_random=False))
def test_row_order_does_not_matter_and_inputs_are_untouched(case: Case, rnd: Any) -> None:
    assume(_oracle_or_raise(case) is not None)
    hits, peps = case.hits_df(), case.peptides_df()
    hits0, peps0, genome0 = hits.copy(deep=True), peps.copy(deep=True), dict(case.genome)
    ref = reconstruct(candidate_rows(hits, peps, None), case.gm, case.genome, CFG)
    order = list(range(len(hits)))
    rnd.shuffle(order)
    shuffled = hits.iloc[order].reset_index(drop=True)
    got = reconstruct(candidate_rows(shuffled, peps, None), case.gm, case.genome, CFG)
    assert got == ref
    pd.testing.assert_frame_equal(hits, hits0)
    pd.testing.assert_frame_equal(peps, peps0)
    assert case.genome == genome0
