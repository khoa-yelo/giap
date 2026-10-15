"""Route 2b: tblastn for protein leftovers, collinear HSP chaining (GIAP_SPEC.md §7.2.3; M3-13).

BLAST+ (public domain) is run as a separate program; this module builds the ``makeblastdb`` and
``tblastn`` argument vectors, parses the 15-column tabular output (``OUTFMT``) and merges
collinear HSPs per (protein, contig, strand) into ``ProteinAlignment`` objects that go through
the M3-12 assignment (``assign_alignments``) with the same thresholds.

Chaining (docs/specs/M3-13-tblastn-route.md; lead decisions):
- HSPs are placed in transcription order (minus strand mirrored: ``t = L - x``), so the protein
  N-terminus comes first on both strands.
- ``nxt`` may follow ``prev`` when it starts and ends later on the genome and on the protein,
  the genome gap ``nxt.start - prev.end`` (0-based half-open, may be negative: real frameshift
  HSPs overlap by a few nt) is ``<= max_gap_nt`` and the query overlap is
  ``<= max_q_overlap_aa`` (both limits inclusive).
- On a circular contig the gap is measured on unwrapped coordinates, so an HSP ending near the
  contig end chains with one starting near 0 (origin-crossing locus); a chain never spans more
  than the contig.
- Best chain (review decision 2): a dynamic program over the compatible HSPs of one
  protein x contig x strand finds the chain with the largest summed bitscore; ties -> larger
  union query coverage, then lower canonical start. Its HSPs are removed and the next best
  chain is built from the rest until none is left. Both criteria are additive along a chain
  (query starts and ends increase, so the union grows by ``q1 - max(q0, prev.q1)`` per step);
  sums use exact fractions. The result depends only on the set of HSPs, never on input order.

``max_gap_nt = 300`` and ``max_q_overlap_aa = 10`` are fixed keyword defaults (Khoa,
2026-10-14): 300 nt tolerates a short frameshifted or low-similarity stretch between two HSPs of
one gene while keeping HSPs of neighbouring genes apart; 10 aa absorbs the few residues two
frameshift HSPs typically share at their junction.

Merged statistics: identity = sum nident / sum length, positive = sum positive / sum length,
coverage = |union of HSP query ranges| / qlen (set explicitly on ``ProteinAlignment``), rank per
protein by summed bitscore (1 = best); ``frameshift``, ``stop_codon`` = 0 and ``alt_start``
False (tblastn sees no start codon). Intervals are canonical (``0 <= start < L``, unwrapped
``end <= start + L``) as ``assign_alignments`` requires (M3-12 review decision 4).
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
from typing import TYPE_CHECKING

from giap.mapping.miniprot import ProteinAlignment
from giap.model.coords import Interval, Replicon, Strand

if TYPE_CHECKING:
    import pandas as pd

    from giap.config.model import ProteinsConfig
    from giap.model.genes import GeneModel

OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore "
    "qlen nident positive"
)
MAX_GAP_NT = 300
MAX_Q_OVERLAP_AA = 10
_N_COLS = 15


def makeblastdb_cmd(genome: Path, db: Path) -> list[str]:
    """``makeblastdb -in genome.fa -dbtype nucl -out db``.

    No ``-parse_seqids``: with it BLAST reports accession-style IDs as ``ref|NC_000913.3|``,
    which never match GIAP's contig IDs (review decision 1).
    """
    return ["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)]


def tblastn_cmd(db: Path, proteins: Path, threads: int, args: str) -> list[str]:
    """``tblastn -db D -query P -num_threads N <args> -outfmt OUTFMT``; ``args`` is split."""
    return [
        "tblastn",
        "-db",
        str(db),
        "-query",
        str(proteins),
        "-num_threads",
        str(threads),
        *args.split(),
        "-outfmt",
        OUTFMT,
    ]


@dataclass(frozen=True, slots=True)
class Hsp:
    """One tabular tblastn row; coordinates as BLAST reports them (1-based closed)."""

    qseqid: str
    sseqid: str
    qstart: int
    qend: int
    sstart: int
    send: int
    length: int
    nident: int
    positive: int
    bitscore: float
    qlen: int


def _parse_row(cols: list[str], where: str) -> Hsp:
    if len(cols) != _N_COLS:
        raise ValueError(f"{where}: {len(cols)} columns, expected {_N_COLS} (outfmt {OUTFMT!r})")
    try:
        return Hsp(
            qseqid=cols[0],
            sseqid=cols[1],
            qstart=int(cols[6]),
            qend=int(cols[7]),
            sstart=int(cols[8]),
            send=int(cols[9]),
            length=int(cols[3]),
            nident=int(cols[13]),
            positive=int(cols[14]),
            bitscore=float(cols[11]),
            qlen=int(cols[12]),
        )
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc


def parse_tblastn(path: Path) -> list[Hsp]:
    """All rows of a tabular tblastn file in ``OUTFMT``, in file order (empty file -> ``[]``).

    Raises ``ValueError`` naming the line on a column count other than 15 or a non-numeric field.
    """
    out: list[Hsp] = []
    with path.open(encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.rstrip("\n\r")
            if not line.strip() or line.startswith("#"):
                continue
            out.append(_parse_row(line.split("\t"), f"{path.name} line {lineno}"))
    return out


@dataclass(frozen=True, slots=True)
class _Piece:
    """An HSP in transcription-order coordinates ``[t0, t1)`` and 0-based query ``[q0, q1)``."""

    hsp: Hsp
    t0: int
    t1: int
    q0: int
    q1: int

    @property
    def key(self) -> tuple[float, int, int, int, int, int, int, int]:
        h = self.hsp
        return (-h.bitscore, h.qstart, h.qend, h.sstart, h.send, h.length, h.nident, h.positive)


def _piece(h: Hsp, replicons: Mapping[str, Replicon]) -> tuple[Strand, _Piece]:
    rep = replicons.get(h.sseqid)
    if rep is None:
        raise ValueError(f"tblastn hit of {h.qseqid!r} on unknown contig {h.sseqid!r}")
    if not 1 <= h.qstart <= h.qend <= h.qlen:
        raise ValueError(f"{h.qseqid!r}: query range {h.qstart}-{h.qend} outside 1..{h.qlen}")
    lo, hi = min(h.sstart, h.send) - 1, max(h.sstart, h.send)
    if lo < 0 or hi > rep.length:
        raise ValueError(f"{h.qseqid!r}: hit {h.sstart}-{h.send} outside {h.sseqid!r}")
    strand: Strand = "-" if h.sstart > h.send else "+"
    t0, t1 = (lo, hi) if strand == "+" else (rep.length - hi, rep.length - lo)
    return strand, _Piece(h, t0, t1, h.qstart - 1, h.qend)


@dataclass(frozen=True, slots=True)
class _Limits:
    max_gap: int
    max_ov: int
    length: int
    circular: bool


def _follows(prev: _Piece, ps: int, nxt: _Piece, ns: int, lim: _Limits) -> bool:
    """``nxt`` (shifted by ``ns``) may come right after ``prev`` (shifted by ``ps``)."""
    a0, a1, b0, b1 = prev.t0 + ps, prev.t1 + ps, nxt.t0 + ns, nxt.t1 + ns
    return (
        b0 > a0
        and b1 > a1
        and b0 - a1 <= lim.max_gap
        and nxt.q0 > prev.q0
        and nxt.q1 > prev.q1
        and prev.q1 - nxt.q0 <= lim.max_ov
    )


_Placed = tuple[_Piece, int]  # piece and its unwrapped shift (0 or +-L)


_Score = tuple[Fraction, int]  # summed bitscore, union query coverage (aa)


def _chains_from(
    first: _Piece, placed: list[_Placed], lim: _Limits
) -> list[tuple[_Score, list[_Placed]]]:
    """Best chain (by ``_Score``) starting at ``first`` (shift 0) ending at each placed piece.

    ``placed`` is sorted by unwrapped start; a chain never spans more than the contig.
    """
    best: list[tuple[_Score, list[_Placed]]] = [
        ((Fraction(first.hsp.bitscore), first.q1 - first.q0), [(first, 0)])
    ]
    for p, s in placed:
        if p is first or p.t0 + s <= first.t0 or p.t1 + s - first.t0 > lim.length:
            continue
        top: tuple[_Score, list[_Placed]] | None = None
        for (bits, cov), chain in best:
            last, ls = chain[-1]
            if not _follows(last, ls, p, s, lim):
                continue
            new = (bits + Fraction(p.hsp.bitscore), cov + p.q1 - max(p.q0, last.q1))
            if top is None or new > top[0]:
                top = (new, [*chain, (p, s)])
        if top is not None:
            best.append(top)
    return best


def _span(chain: list[_Placed], strand: Strand, rep: Replicon) -> tuple[int, int]:
    """Canonical ``[lo, hi)`` of a chain (``0 <= lo < L``, unwrapped ``hi``)."""
    t0 = min(p.t0 + s for p, s in chain)
    t1 = max(p.t1 + s for p, s in chain)
    lo, hi = (t0, t1) if strand == "+" else (rep.length - t1, rep.length - t0)
    shift = -rep.length * (lo // rep.length)
    return lo + shift, hi + shift


def _best_chain(pieces: list[_Piece], strand: Strand, rep: Replicon, lim: _Limits) -> list[_Placed]:
    """The best chain among ``pieces`` (non-empty): max bitscore, coverage, then lower start."""
    shifts = (0, lim.length) if lim.circular else (0,)
    placed = sorted(
        ((p, s) for p in pieces for s in shifts), key=lambda x: (x[0].t0 + x[1], x[0].key)
    )
    top: tuple[tuple[Fraction, int, int, int], list[_Placed]] | None = None
    for first in pieces:
        for (bits, cov), chain in _chains_from(first, placed, lim):
            lo, hi = _span(chain, strand, rep)
            k = (-bits, -cov, lo, hi)
            if top is None or k < top[0]:
                top = (k, chain)
    assert top is not None
    return top[1]


def _union_len(ranges: list[tuple[int, int]]) -> int:
    total, end = 0, -1
    for q0, q1 in sorted(ranges):
        if q1 > end:
            total += q1 - max(q0, end)
            end = q1
    return total


def _alignment(
    chain: list[_Placed], strand: Strand, rep: Replicon
) -> tuple[float, ProteinAlignment]:
    pieces = [p for p, _ in chain]
    lo, hi = _span(chain, strand, rep)
    hsps = [p.hsp for p in pieces]
    qlen = hsps[0].qlen
    length = sum(h.length for h in hsps)
    aln = ProteinAlignment(
        protein_id=hsps[0].qseqid,
        protein_len_aa=qlen,
        route="tblastn",
        interval=Interval(rep.id, lo, hi, strand),
        rank=0,
        identity=min(1.0, sum(h.nident for h in hsps) / length) if length else 0.0,
        positive=min(1.0, sum(h.positive for h in hsps) / length) if length else None,
        q_start=min(p.q0 for p in pieces),
        q_end=max(p.q1 for p in pieces),
        frameshift=0,
        stop_codon=0,
        coverage=min(1.0, _union_len([(p.q0, p.q1) for p in pieces]) / qlen),
    )
    return sum(h.bitscore for h in hsps), aln


def _rank_key(item: tuple[float, ProteinAlignment]) -> tuple[float, str, int, int, str, int, int]:
    bits, a = item
    iv = a.interval
    return (-bits, iv.contig, iv.start, iv.end, iv.strand, a.q_start, a.q_end)


def merge_hsps(
    hsps: Sequence[Hsp],
    replicons: Mapping[str, Replicon],
    max_gap_nt: int = MAX_GAP_NT,
    max_q_overlap_aa: int = MAX_Q_OVERLAP_AA,
) -> list[ProteinAlignment]:
    """Chain collinear HSPs per (protein, contig, strand) into ranked ``ProteinAlignment``s.

    Output is sorted by protein, then rank. Raises ``ValueError`` for an HSP on a contig absent
    from ``replicons`` or with coordinates outside the protein or the contig.
    """
    groups: dict[tuple[str, str, Strand], list[_Piece]] = defaultdict(list)
    for h in hsps:
        strand, piece = _piece(h, replicons)
        groups[(h.qseqid, h.sseqid, strand)].append(piece)
    by_protein: dict[str, list[tuple[float, ProteinAlignment]]] = defaultdict(list)
    for (pid, contig, strand), pieces in groups.items():
        rep = replicons[contig]
        lim = _Limits(max_gap_nt, max_q_overlap_aa, rep.length, rep.topology == "circular")
        free = sorted(pieces, key=lambda p: (p.t0, p.key))
        while free:
            chain = _best_chain(free, strand, rep, lim)
            used = {id(p) for p, _ in chain}
            free = [p for p in free if id(p) not in used]
            by_protein[pid].append(_alignment(chain, strand, rep))
    out: list[ProteinAlignment] = []
    for pid in sorted(by_protein):
        ranked = sorted(by_protein[pid], key=_rank_key)
        for rank, (_, aln) in enumerate(ranked, start=1):
            out.append(replace(aln, rank=rank))
    return out


def assign_tblastn(
    hsps: Sequence[Hsp], gm: GeneModel, cfg: ProteinsConfig
) -> tuple[pd.DataFrame, set[str]]:
    """Merge HSPs and assign the chains with the M3-12 rules (``route = tblastn``)."""
    from giap.mapping.protein_assign import assign_alignments

    return assign_alignments(merge_hsps(hsps, gm.replicons), gm, cfg)
