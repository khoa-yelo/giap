"""Brute-force oracle for `giap makedb` six-frame ORFs (issue M3-23; GIAP_SPEC.md §7.9.2).

Independent of ``giap.makedb``: plain-string translation with table 11 and a codon-by-codon walk
over each strand, written from docs/specs/M3-23-makedb.md:

* ORFs are maximal stop-to-stop segments of each of the six frames (stop codons excluded),
  >= ``min_aa`` residues, without ``X``. On a circular replicon the frames run across the origin
  (unwrapped coordinates, ``end > length``); on a linear contig a segment may also end at either
  contig end (truncated ORF, kept only if >= ``min_aa``).
* An ORF that overlaps an annotated CDS or sORF *in frame* (same strand, same phase) is dropped;
  its upstream part becomes an ``nterm_extension`` entry (see :func:`nterm_entries`): upstream
  part + the first 30 codons of the CDS translated codon by codon from the genome (lead decision
  5); header interval = upstream part + those codons (lead decision 9); undefined when the walk
  runs off a linear contig or meets ``X`` (lead decision 4).
* Other ORFs take the class of ``giap.mapping.loci.classify_locus``; ``rna_overlap`` and
  ``pseudogene`` ORFs are omitted.

The oracle only uses ``classify_locus`` (M3-14, merged) for the class of non-dropped ORFs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from giap.mapping.loci import classify_locus
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene, GeneModel

STOPS = frozenset({"TAA", "TAG", "TGA"})
_COMP = str.maketrans("ACGTN", "TGCAN")
_BASES = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11 = {
    a + b + c: _AA[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
ORF_CLASSES = ("intergenic", "cds_outframe", "antisense")
NTERM_PROTEIN_AA = 30  # lead decisions 2 and 5 (M3-23 Q1)
NTERM_MIN_UPSTREAM_AA = 7


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def translate(nt: str) -> str:
    """Table 11, codon by codon (no start-codon rule); non-ACGT codon -> ``X``."""
    return "".join(CODE11.get(nt[i : i + 3], "X") for i in range(0, len(nt) - 2, 3))


def strand_nt(seq: str, start: int, end: int, strand: str) -> str:
    """Nucleotides of unwrapped ``[start, end)`` read on ``strand`` (circular modulo len)."""
    n = len(seq)
    s = "".join(seq[i % n] for i in range(start, end))
    return revcomp(s) if strand == "-" else s


class DegenerateCircle(ValueError):
    """A frame of a circular contig has a stop-free run longer than the contig (not specified)."""


@dataclass(frozen=True)
class Orf:
    contig: str
    start: int  # 0-based, unwrapped
    end: int
    strand: str
    aa: str


def _plus_orfs(seq: str, circular: bool) -> list[tuple[int, int, str]]:
    """Maximal stop-to-stop segments of the three frames of ``seq`` read 5'->3' (stops excluded)."""
    n = len(seq)
    out: list[tuple[int, int, str]] = []
    if circular:
        stop_set = {q for q in range(n) if strand_nt(seq, q, q + 3, "+") in STOPS}
        for q in sorted(stop_set):
            s = q + 3
            p = s
            while (p % n) not in stop_set:
                p += 3
                if p - s > n:  # stop-free stretch longer than the replicon: undefined
                    raise DegenerateCircle(f"stop-free run of > {n} nt on a circular contig")
            if p > s:
                start = s % n
                end = start + (p - s)
                out.append((start, end, translate(strand_nt(seq, start, end, "+"))))
        return out
    for phase in range(3):
        s = phase
        p = phase
        while p + 3 <= n:
            if seq[p : p + 3] in STOPS:
                if p > s:
                    out.append((s, p, translate(seq[s:p])))
                s = p + 3
            p += 3
        end = s + 3 * ((n - s) // 3)
        if end > s:
            out.append((s, end, translate(seq[s:end])))
    return out


def all_orfs(contig: str, seq: str, rep: Replicon) -> list[Orf]:
    """Every maximal stop-to-stop segment of the six frames (no length or X filter)."""
    n = len(seq)
    circ = rep.topology == "circular"
    out = [Orf(contig, a, b, "+", aa) for a, b, aa in _plus_orfs(seq, circ)]
    for a, b, aa in _plus_orfs(revcomp(seq), circ):
        s, e = n - b, n - a
        if s < 0:
            s, e = s + n, e + n
        out.append(Orf(contig, s, e, "-", aa))
    return out


def _shifts(rep: Replicon) -> tuple[int, ...]:
    return (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)


def in_frame_overlap(s: int, e: int, strand: str, g: Gene, rep: Replicon) -> bool:
    gi = g.interval
    if gi.strand != strand:
        return False
    for d in _shifts(rep):
        a, b = s + d, e + d
        if max(a, gi.start) >= min(b, gi.end):
            continue
        if strand == "+" and (a - gi.start) % 3 == 0:
            return True
        if strand == "-" and (gi.end - b) % 3 == 0:
            return True
    return False


def classified_orfs(gm: GeneModel, genome: Mapping[str, str], min_aa: int) -> list[tuple[str, Orf]]:
    """(class, ORF) of every written non-nterm ORF, sorted by contig, start, end, strand."""
    coding = [g for g in gm.genes if g.feature_type in ("cds", "sorf")]
    rows: list[tuple[str, Orf]] = []
    order = {c: k for k, c in enumerate(gm.replicons)}
    for contig, seq in genome.items():
        rep = gm.replicons[contig]
        for orf in all_orfs(contig, seq, rep):
            if len(orf.aa) < min_aa or "X" in orf.aa:
                continue
            if any(
                in_frame_overlap(orf.start, orf.end, orf.strand, g, rep)
                for g in coding
                if g.interval.contig == contig
            ):
                continue
            cls, _ = classify_locus(
                Interval(contig, orf.start, orf.end, orf.strand),  # type: ignore[arg-type]
                gm,
                genome,
            )
            if cls in ("rna_overlap", "pseudogene"):
                continue
            if cls not in ORF_CLASSES:
                raise AssertionError(f"oracle: unexpected class {cls} for {orf}")
            rows.append((cls, orf))
    rows.sort(key=lambda r: (order[r[1].contig], r[1].start, r[1].end, r[1].strand == "-"))
    return rows


@dataclass(frozen=True)
class Nterm:
    locus_tag: str
    contig: str
    strand: str
    up_start: int  # 0-based unwrapped interval of the upstream part [up_start, up_end)
    up_end: int
    aa: str  # upstream translation + the first 30 codons of the CDS (translated from the genome)
    n_protein: int  # residues taken from the annotated CDS (the junction part)


def _upstream(seq: str, g: Gene, rep: Replicon) -> tuple[int, int] | None:
    """Upstream part (unwrapped plus coordinates), or None when it is undefined: the in-frame walk
    runs off a linear contig, meets a non-ACGT codon (``X``), or finds no stop within one turn of
    a circular contig (lead decision 4; the last case is degenerate and pinned by a unit test)."""
    n = rep.length
    circ = rep.topology == "circular"
    gi = g.interval
    if gi.strand == "+":
        q = gi.start - 3
        while True:
            if not circ and q < 0:
                return None
            cod = strand_nt(seq, q, q + 3, "+")
            if cod in STOPS:
                break
            if cod not in CODE11:
                return None
            q -= 3
            if gi.start - q > n:
                return None
        a, b = q + 3, gi.start
    else:
        q = gi.end
        while True:
            if not circ and q + 3 > n:
                return None
            cod = strand_nt(seq, q, q + 3, "-")
            if cod in STOPS:
                break
            if cod not in CODE11:
                return None
            q += 3
            if q - gi.end > n:
                return None
        a, b = gi.end, q
    if a < 0:
        a, b = a + n, b + n
    elif a >= n:
        a, b = a - n, b - n
    return a, b


def junction(seq: str, g: Gene) -> str:
    """First ``min(30, protein length)`` codons of the CDS, translated codon by codon from the
    genome (GTG -> V, TTG -> L; lead decision 5), not copied from Bakta's ``aa``."""
    gi = g.interval
    k = min(NTERM_PROTEIN_AA, (gi.end - gi.start) // 3 - 1)
    if gi.strand == "+":
        return translate(strand_nt(seq, gi.start, gi.start + 3 * k, "+"))
    return translate(strand_nt(seq, gi.end - 3 * k, gi.end, "-"))


def nterm_entries(gm: GeneModel, genome: Mapping[str, str]) -> tuple[list[Nterm], set[str]]:
    """Expected ``nterm_extension`` entries and the locus tags whose expectation is undefined.

    Undefined (not asserted either way; lead decision 4): the in-frame upstream walk runs off a
    linear contig or meets an ``X`` (or, degenerate, finds no stop within one circular turn).
    """
    out: list[Nterm] = []
    undefined: set[str] = set()
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf"):
            continue
        rep = gm.replicons[g.interval.contig]
        seq = genome[g.interval.contig]
        up = _upstream(seq, g, rep)
        if up is None:
            undefined.add(g.locus_tag)
            continue
        a, b = up
        up_aa = translate(strand_nt(seq, a, b, g.interval.strand))
        if len(up_aa) < NTERM_MIN_UPSTREAM_AA:
            continue
        prot = junction(seq, g)
        out.append(
            Nterm(g.locus_tag, g.interval.contig, g.interval.strand, a, b, up_aa + prot, len(prot))
        )
    return out, undefined


def orf_header(cls: str, contig: str, start: int, end: int, strand: str) -> str:
    return f"GIAPORF_{cls}_{contig}_{start + 1}_{end}_{'p' if strand == '+' else 'm'}"


def nterm_header(nt: Nterm, gm: GeneModel) -> str:
    """The pinned header (lead decision 9): the interval covers the upstream part plus the
    junction codons (upstream + first 30 codons of the CDS), unwrapped across the origin."""
    span = 3 * nt.n_protein
    if nt.strand == "+":
        a, b = nt.up_start, nt.up_end + span
    else:
        n = gm.replicons[nt.contig].length
        lo = nt.up_start - span
        a, b = (lo + n, nt.up_end + n) if lo < 0 else (lo, nt.up_end)
    return orf_header("nterm_extension", nt.contig, a, b, nt.strand)
