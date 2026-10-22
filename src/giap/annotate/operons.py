"""Operon pairs, distance-model fallback, transcription units and TU upstream regions
(GIAP_SPEC.md §7.1.3, §9.3, §6.5; M2-12).

Adjacent pairs are consecutive same-strand genes per contig in genome order, plus the wrap pair
(last, first) on circular replicons. UniOP probabilities are used for every pair UniOP reports;
the rest get p(operon | d) = 1 / (1 + exp((d - d₀)/s)). TUs are maximal runs linked by
p ≥ threshold. Coordinates are 0-based half-open; origin-spanning intervals are unwrapped.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from giap.annotate.uniop_adapter import PairProb, genome_order
from giap.io.fasta import write_fasta
from giap.io.tsv import write_tsv
from giap.model.coords import Interval, Replicon, to_one_based_closed
from giap.model.genes import Gene

DISTANCE_FLAG = "operon_source_distance_model"
PAIRS_COLUMNS = ("gene_a", "gene_b", "uniop_p", "source")
TUS_COLUMNS = (
    "tu_id", "contig", "start", "end", "strand", "genes", "upstream_start", "upstream_end",
)  # fmt: skip
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True, slots=True)
class Tu:
    """A transcription unit; ``genes`` in transcription order."""

    tu_id: str
    contig: str
    interval: Interval
    strand: str
    genes: tuple[str, ...]
    upstream: Interval


def adjacent_pairs(
    genes: Sequence[Gene], replicons: Mapping[str, Replicon], same_strand: bool = True
) -> list[tuple[Gene, Gene]]:
    """Consecutive genes per contig (FASTA order), plus the wrap pair if circular; only the
    same-strand ones unless ``same_strand`` is False (the M6-04 pair chain)."""
    out: list[tuple[Gene, Gene]] = []
    for contig, gs in genome_order(genes, list(replicons)).items():
        ring = replicons[contig].topology == "circular" and len(gs) >= 2
        nxt = list(itertools.pairwise(gs)) + ([(gs[-1], gs[0])] if ring else [])
        out += [(a, b) for a, b in nxt if not same_strand or a.interval.strand == b.interval.strand]
    return out


def intergenic_distance(a: Gene, b: Gene, replicon: Replicon) -> int:
    """``b.start - a.end`` in genome order (across the origin on circular); < 0 = overlap."""
    d = b.interval.start - a.interval.end
    if replicon.topology == "circular" and b.interval.start < a.interval.start:
        d += replicon.length
    return d


def distance_model_p(d: int, d0: float, s: float) -> float:
    """p(operon | d) = 1 / (1 + exp((d - d₀)/s)), computed without overflow."""
    z = (d - d0) / s
    if z >= 0:
        e = math.exp(-z)
        return e / (1.0 + e)
    return 1.0 / (1.0 + math.exp(z))


def build_pairs(
    genes: Sequence[Gene],
    replicons: Mapping[str, Replicon],
    uniop: Sequence[PairProb] | None,
    d0: float,
    s: float,
) -> tuple[list[PairProb], tuple[str, ...]]:
    """One probability per adjacent pair; ``uniop`` None or empty = UniOP failure (run flag)."""
    failed = not uniop
    reported = {} if failed else {(p.gene_a, p.gene_b): p for p in uniop or ()}
    out: list[PairProb] = []
    for a, b in adjacent_pairs(genes, replicons):
        hit = reported.get((a.locus_tag, b.locus_tag))
        if hit is None:
            d = intergenic_distance(a, b, replicons[a.interval.contig])
            hit = PairProb(a.locus_tag, b.locus_tag, distance_model_p(d, d0, s), "distance_model")
        out.append(hit)
    return out, ((DISTANCE_FLAG,) if failed else ())


def _window(lo: int, hi: int, replicon: Replicon) -> tuple[int, int]:
    if replicon.topology == "linear":
        return max(lo, 0), min(hi, replicon.length)
    shift = (lo // replicon.length) * replicon.length
    return lo - shift, hi - shift


def tu_upstream(
    first: Gene, others: Sequence[Gene], replicon: Replicon, bp: Sequence[int] = (-300, 50)
) -> Interval:
    """-|bp[0]| to +bp[1] around the start codon, strand-aware, clipped at the nearest upstream
    gene (either strand) but never past the start codon (spec Q2)."""
    iv = first.interval
    shifts = (-replicon.length, 0, replicon.length) if replicon.topology == "circular" else (0,)
    copies = [
        (o.interval.start + k, o.interval.end + k)
        for o in others
        if o.interval.contig == iv.contig
        for k in shifts
    ]
    if replicon.topology == "circular":  # the first gene itself bounds the window on a ring
        copies += [(iv.start + k, iv.end + k) for k in shifts if k]
    if iv.strand == "-":
        sc = iv.end
        hi = min([sc - bp[0], *(max(s, sc) for s, e in copies if e > sc)])
        lo, hi = _window(sc - bp[1], hi, replicon)
    else:
        sc = iv.start
        lo = max([sc + bp[0], *(min(e, sc) for s, e in copies if s < sc)])
        lo, hi = _window(lo, sc + bp[1], replicon)
    return Interval(iv.contig, lo, hi, iv.strand)


def _chains(gs: list[Gene], linked: set[tuple[str, str]], ring: bool) -> list[list[Gene]]:
    """Maximal runs of linked genes (genome order); a fully linked ring breaks at the origin."""
    n = len(gs)
    link = [(gs[i].locus_tag, gs[(i + 1) % n].locus_tag) in linked for i in range(n)]
    if not ring:
        link[-1] = False
    starts = [i for i in range(n) if not link[i - 1]] or [0]
    out = []
    for i in starts:
        run = [gs[i]]
        j = i
        while link[j] and (j + 1) % n != i:
            j = (j + 1) % n
            run.append(gs[j])
        out.append(run)
    return out


def build_tus(
    genes: Sequence[Gene],
    pairs: Sequence[PairProb],
    threshold: float,
    replicons: Mapping[str, Replicon],
    upstream_bp: Sequence[int] = (-300, 50),
) -> list[Tu]:
    """TUs numbered ``TU_<n>`` by contig (FASTA order), start, then strand (+ before -)."""
    linked = {(p.gene_a, p.gene_b) for p in pairs if p.uniop_p >= threshold}
    order = {c: k for k, c in enumerate(replicons)}
    raw: list[tuple[tuple[int, int, int], str, Interval, tuple[str, ...], Interval]] = []
    for contig, gs in genome_order(genes, list(replicons)).items():
        if not gs:
            continue
        rep = replicons[contig]
        ring = rep.topology == "circular" and len(gs) >= 2
        for run in _chains(gs, linked, ring):
            start, strand = run[0].interval.start, run[0].interval.strand
            end = max(g.interval.end + (rep.length if g.interval.start < start else 0) for g in run)
            if len(run) > 1 and rep.topology == "circular":  # never longer than the ring
                end = min(end, start + rep.length)
            first = run[0] if strand != "-" else run[-1]
            tags = tuple(g.locus_tag for g in (run if strand != "-" else run[::-1]))
            up = tu_upstream(first, [g for g in gs if g is not first], rep, upstream_bp)
            key = (order[contig], start, 0 if strand == "+" else 1)
            raw.append((key, contig, Interval(contig, start, end, strand), tags, up))
    raw.sort(key=lambda r: (r[0], r[3]))
    return [
        Tu(f"TU_{n}", contig, iv, iv.strand, tags, up)
        for n, (_, contig, iv, tags, up) in enumerate(raw, start=1)
    ]


def write_operon_pairs_tsv(path: Path, pairs: Sequence[PairProb]) -> None:
    """``annotation/operon_pairs.tsv`` (§9.3)."""
    write_tsv(path, PAIRS_COLUMNS, [(p.gene_a, p.gene_b, p.uniop_p, p.source) for p in pairs])


def write_tus_tsv(path: Path, tus: Sequence[Tu]) -> None:
    """``annotation/tus.tsv`` (§9.3): 1-based closed, unwrapped across the origin."""
    rows = [
        (t.tu_id, t.contig, *to_one_based_closed(t.interval), t.strand, list(t.genes),
         *to_one_based_closed(t.upstream))
        for t in tus
    ]  # fmt: skip
    write_tsv(path, TUS_COLUMNS, rows)


def region_sequence(iv: Interval, seq: str) -> str:
    """Strand-aware sequence of ``iv`` (unwrapped intervals read across the origin)."""
    nt = (seq + seq)[iv.start : iv.end] if iv.end > len(seq) else seq[iv.start : iv.end]
    return nt.translate(_COMP)[::-1] if iv.strand == "-" else nt


def write_tu_regions(outdir: Path, tus: Sequence[Tu], genome: Mapping[str, str]) -> list[Path]:
    """One ``<tu_id>.fa`` per TU, its upstream region on one line; returns the paths in TU order."""
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for t in tus:
        path = outdir / f"{t.tu_id}.fa"
        seq = region_sequence(t.upstream, genome[t.contig])
        # one sequence line (= promoter_calculator.region_fasta_text; §7.1.1 gotcha, M2-14)
        write_fasta(path, [(t.tu_id, seq)], width=max(len(seq), 1))
        paths.append(path)
    return paths
