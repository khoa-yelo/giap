"""Transcription and translation signals (issue M2-05; GIPPA_SPEC.md §11.3).

Each transcription unit (TU) laid out by the v2 simulator carries a ``LEADER`` upstream of its
first start codon and a ``TRAILER`` downstream of its last stop codon (both background sequence,
part of the TU span used by the layout). Signals are planted with probability = strength:

* RBS: ``AGGAGG`` with a 5-10 nt spacer (nt between its last ``G`` and the start codon). For the
  first gene it lies in the leader; for later genes in the intra-TU gap or, when the gap is too
  short, inside the upstream gene's coding sequence (codons rewritten without creating a stop).
  Intra-TU gaps of 8-10 bp fit no spacer (the core would overlap the upstream stop codon), so the
  v2 layout widens them to 11-13 bp. Genes with an N-terminal extension (annotated start moved
  downstream of the true start) get no planted RBS.
* One Shine-Dalgarno core per gene: ``AGGAGG`` overlaps itself only at a shift of 3, so a second
  core within the 5-10 nt spacer window can only come from an ``AGG`` right after the planted
  core (spacer >= 8) or right before it (spacer <= 7). That ``AGG`` is broken by a ``C`` (its
  first, resp. last base), in the gene's reading orientation on either strand.
* sigma-70 promoter: ``TTGACA`` + 17 nt + ``TATAAT`` ending 3-9 nt before the TSS; the TSS lies
  2-10 nt upstream of the (possibly unplanted) RBS.
* Terminator: a 7-bp GC stem, 4-nt loop and an 8-nt U-tract, 5-15 nt after the TU's last stop.

Plans are made in TU-local transcription coordinates (0 = 5' end of the TU span) and applied to
the genome once it is written. ``put``/``get`` address a sequence relative to an anchor
``[s, e)`` on a strand: offset 0 is the anchor's 5' base, negative offsets lie upstream; both
read and write in the strand's reading (transcription) orientation. Signals are planted at their
strengths also when ``write_precomputed`` is off (then no truth tables are written).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

import numpy as np

from gippasim.codons import SENSE_USAGE, codons_of, revcomp, sense_usage, stops

LEADER = 70
TRAILER = 50
RBS = "AGGAGG"
MINUS35, MINUS10 = "TTGACA", "TATAAT"
U_TRACT = "TTTTTTTT"
SPACERS = range(5, 11)
_SENSE = sorted(SENSE_USAGE)
_COMP = {ord("A"): ord("T"), ord("C"): ord("G"), ord("G"): ord("C"), ord("T"): ord("A")}
_COMP_BYTES = bytes.maketrans(b"ACGTN", b"TGCAN")


@dataclass
class TuPlan:
    """Signals of one TU; offsets are TU-local, in transcription orientation."""

    genes: list[Any]  # gippasim.simulate._Gene, transcription order
    length: int
    strand: str = "+"
    pos: int = 0  # genome start of the TU span (set by the layout)
    writes: list[tuple[int, str]] = field(default_factory=list)
    fixes: list[tuple[int, int]] = field(default_factory=list)  # (RBS offset, spacer) to clean
    rbs: list[tuple[Any, int]] = field(default_factory=list)  # (gene, spacer)
    tss: int | None = None
    term: tuple[int, int] | None = None

    def span(self, x0: int, x1: int) -> tuple[int, int]:
        """Genome ``[start, end)`` of the TU-local ``[x0, x1)``."""
        if self.strand == "+":
            return self.pos + x0, self.pos + x1
        return self.pos + self.length - x1, self.pos + self.length - x0


def _index(s: int, e: int, strand: str, x: int, n: int) -> int:
    return (s + x if strand == "+" else e - 1 - x) % n


def put(seq: bytearray, s: int, e: int, strand: str, x: int, w: str) -> None:
    """Write ``w`` (transcription orientation) at offset ``x`` of the anchor ``[s, e)``."""
    n = len(seq)
    for i, ch in enumerate(w.encode()):
        seq[_index(s, e, strand, x + i, n)] = ch if strand == "+" else _COMP[ch]


def get(seq: bytearray, s: int, e: int, strand: str, x: int, k: int) -> str:
    """``k`` nt at offset ``x`` of the anchor ``[s, e)``, in ``strand``'s reading orientation."""
    n = len(seq)
    out = bytes(seq[_index(s, e, strand, x + i, n)] for i in range(k))
    return (out if strand == "+" else out.translate(_COMP_BYTES)).decode()


def rbs_fixes(x: int, sp: int) -> list[tuple[int, int]]:
    """(``AGG`` triplet offset, offset of the base set to ``C``) guarding the RBS at offset ``x``
    with spacer ``sp`` against a second core: the triplet after it when ``sp >= 8`` (first base),
    the one before it when ``sp <= 7`` (last base, never part of an upstream stop codon)."""
    return [(x + 6, x + 6)] if sp >= 8 else [(x - 3, x - 1)]


def clean_rbs(seq: bytearray, s: int, e: int, strand: str, x: int, sp: int) -> None:
    """Break the ``AGG`` that would put a second core in the spacer window of the RBS at ``x``."""
    for t, b in rbs_fixes(x, sp):
        if get(seq, s, e, strand, t, 3) == "AGG":
            put(seq, s, e, strand, b, "C")


def hairpin(rng: np.random.Generator) -> str:
    stem = "".join(rng.choice(list("GC"), size=7).tolist())
    loop = "".join(rng.choice(list("ACGT"), size=4).tolist())
    return stem + loop + revcomp(stem) + U_TRACT


def _sense_codons(table: int) -> list[str]:
    """Sorted sense codons of genetic code ``table`` (TGA included at 4/25)."""
    return _SENSE if table == 11 else sorted(sense_usage(table))


def _body_rbs(rng: np.random.Generator, a: Any, rel: int, sp: int, table: int = 11) -> None:
    """Write the RBS into ``a``'s codons at nt offset ``rel`` without creating a stop codon of
    genetic code ``table`` (at 4/25 a TGA it creates is sense and stays)."""
    nt = a.nt[:rel] + RBS + a.nt[rel + 6 :]
    c0 = rel - rel % 3
    if nt[c0 : c0 + 3] in stops(table):
        tail = nt[rel : c0 + 3]
        options = [c for c in _sense_codons(table) if c.endswith(tail)]
        nt = nt[:c0] + options[int(rng.integers(len(options)))] + nt[c0 + 3 :]
    for t, b in rbs_fixes(rel, sp):  # a C never completes a stop codon
        if t >= 0 and nt[t : t + 3] == "AGG":
            nt = nt[:b] + "C" + nt[b + 1 :]
    a.codons = codons_of(nt)


def _intra_rbs(
    rng: np.random.Generator, plan: TuPlan, a: Any, b: Any, table: int = 11
) -> int | None:
    """Plant the RBS of ``b`` (after ``a`` in the same TU); returns its spacer."""
    ok = [sp for sp in SPACERS if b.start - sp <= a.end - 3 or b.start - sp - 6 >= a.end]
    if not ok:
        return None
    sp = ok[int(rng.integers(len(ok)))]
    x0 = b.start - sp - 6
    if x0 >= a.end:
        plan.writes.append((x0, RBS))
        plan.fixes.append((x0, sp))
    else:
        _body_rbs(rng, a, x0 - a.start, sp, table)
    return sp


def plan_tu(rng: np.random.Generator, cfg: Any, genes: list[Any], length: int) -> TuPlan:
    """Shift ``genes`` by ``LEADER`` and plan the TU's signals (strengths and
    ``genetic_code`` of ``cfg``)."""
    table = int(getattr(cfg, "genetic_code", 11))
    for g in genes:
        g.start += LEADER
        g.end += LEADER
    plan = TuPlan(genes, length + LEADER + TRAILER, strand=genes[0].strand)
    sp0, utr = int(rng.integers(5, 11)), int(rng.integers(2, 11))
    if rng.random() < cfg.rbs_strength:
        x = LEADER - sp0 - 6
        plan.writes.append((x, RBS))
        plan.rbs.append((genes[0], sp0))
        plan.fixes.append((x, sp0))
    tss = LEADER - sp0 - 6 - utr
    d = int(rng.integers(3, 10))
    if rng.random() < cfg.promoter_strength:
        plan.writes += [(tss - d - 29, MINUS35), (tss - d - 6, MINUS10)]
        plan.tss = tss
    for a, b in pairwise(genes):
        if rng.random() < cfg.rbs_strength:
            sp = _intra_rbs(rng, plan, a, b, table)
            if sp is not None:
                plan.rbs.append((b, sp))
    end, off = genes[-1].end, int(rng.integers(5, 16))
    if rng.random() < cfg.terminator_strength:
        term = hairpin(rng)
        plan.writes.append((end + off, term))
        plan.term = (end + off, end + off + len(term))
    return plan


def apply_plan(seq: bytearray, plan: TuPlan) -> None:
    s, e = plan.pos, plan.pos + plan.length
    for x, w in plan.writes:
        put(seq, s, e, plan.strand, x, w)
    for x, sp in plan.fixes:
        clean_rbs(seq, s, e, plan.strand, x, sp)


def plant_gene_rbs(
    rng: np.random.Generator, strength: float, seq: bytearray, s: int, e: int, strand: str
) -> int | None:
    """RBS upstream of a gene outside the TU layout (planted blocks); returns its spacer."""
    if rng.random() >= strength:
        return None
    sp = int(rng.integers(5, 11))
    put(seq, s, e, strand, -sp - 6, RBS)
    clean_rbs(seq, s, e, strand, -sp - 6, sp)
    return sp
