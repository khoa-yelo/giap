"""Brute-force oracle for ±1-nt frameshift splits (issue M3-08; GIPPA_SPEC.md §7.2.4, §11.2).

Written by the test author from docs/specs/M3-08-frameshift-splits.md only; it does not import
``gippa.mapping.frameshift``, ``gippa.mapping.tolerant_dp`` or ``gippa.mapping.diffops``. Side costs
come from the M3-06/M3-07 test oracles (``tests/property/oracle_tolerant.py``,
``tests/property/brute_tolerant.py``), which encode the §7.2.4 operations from the spec text.

Model (normative text of M3-08, "Semantics"), in one strand's nucleotide coordinates ``x`` of the
*extended* sequence (circular replicons carry their wrap):

- a split at query position ``p`` (``min_side <= p <= L - min_side``) aligns the prefix
  ``q[0:p]`` to a window of frame ``F`` that *ends* at the junction (free start) and the suffix
  ``q[p:L]`` to a window of frame ``G`` (same contig and strand) that *starts* at the junction
  (free end); the first suffix codon starts ``delta`` nt (``±1``) after the last prefix codon ends;
- cost = 1 + M3-06 cost of the prefix window + M3-06 cost of the suffix window (<= k);
- decoys (test-author reading of "searched identically" / "mirrored target intervals satisfy the
  same adjacency rule"): a decoy frame is read in its own direction. Every decoy residue sits on
  the codon of its mirrored target residue (M3-03 ``mirror_aa``); for ``shuffle`` that is the same
  codon, for ``reverse`` the frame runs 3'->5' over the target codons. The adjacency is applied in
  the decoy's reading direction: the suffix's first mirrored codon lies ``3 + delta`` nt *after*
  (``shuffle``) or *before* (``reverse``) the prefix's last mirrored codon. Hence the reverse-decoy
  split placements of ``Q`` are exactly the target split placements of ``reversed(Q)``;
- precedence at one locus (M3-08 edge case "at the same locus, the non-frameshift explanation
  wins"): a split placement is dropped when the whole query has a non-frameshift alignment (M3-06
  operations, any frame of the same contig and strand, same search) of cost <= the split's cost
  whose (mirrored) interval overlaps the split's locus (modulo the replicon length when circular);
- per peptide and search the result is the minimal cost among the remaining placements and every
  remaining placement at that cost.

Genome coordinates: the locus start is reduced into ``[0, length)``, the end and the two blocks
are shifted by the same amount (unwrapped, §6.5), blocks ascending. ``frame`` is the TSV frame label
of the prefix block (for decoys: of the prefix's mirrored block).
"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from dataclasses import dataclass
from functools import cache
from pathlib import Path

import numpy as np

_PROPERTY = str(Path(__file__).resolve().parents[1] / "property")
if _PROPERTY not in sys.path:
    sys.path.insert(0, _PROPERTY)

from brute_tolerant import (  # noqa: E402
    INF,
    NRES,
    ONE_TWO,
    TWO_ONE,
    TWO_TWO,
    overlaps_mod,
    semi_global_costs,
)
from oracle_tolerant import PRECEDENCE, align_window  # noqa: E402

from gippa.mapping.sixframe import Frame, SixFrameIndex, decode, encode_query  # noqa: E402
from gippa.model import Interval  # noqa: E402

TYPE_ORDER: tuple[str, ...] = (*PRECEDENCE, "frameshift")  # M3-05 DIFF_TYPES order


def semi_global_matrix(q: np.ndarray, t: np.ndarray, kmax: int) -> np.ndarray:
    """``D[i, j]``: minimal cost of ``q[:i]`` against any window of ``t`` ending at ``j``.

    Same recurrences as ``brute_tolerant.semi_global_costs`` (free start), all rows kept; values
    above ``kmax`` are ``INF``. ``D[0, j] = 0``.
    """
    q = np.asarray(q, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    L, n = len(q), len(t)
    D = np.full((L + 1, n + 1), INF, dtype=np.int64)
    D[0, :] = 0
    tok = t < NRES
    pair_comp = np.full(n + 1, -1, dtype=np.int64)
    if n >= 2:
        lo = np.minimum(t[:-1], t[1:])
        hi = np.maximum(t[:-1], t[1:])
        ok = tok[:-1] & tok[1:]
        pair_comp[2:] = np.where(ok, lo * 32 + hi, -1)
    t1 = np.full(n + 1, -1, dtype=np.int64)
    t1[1:] = t
    t0 = np.full(n + 1, -1, dtype=np.int64)
    t0[2:] = t[:-1]
    ok1 = np.zeros(n + 1, dtype=bool)
    ok1[1:] = tok
    for i in range(1, L + 1):
        row = np.full(n + 1, INF, dtype=np.int64)
        b = int(q[i - 1])
        prev = np.full(n + 1, INF, dtype=np.int64)
        prev[1:] = D[i - 1, :-1]
        step = np.where(t1 == b, 0, 1)
        row = np.where(ok1, np.minimum(row, prev + step), row)
        if b < NRES and b in ONE_TWO and n >= 2:
            prev2 = np.full(n + 1, INF, dtype=np.int64)
            prev2[2:] = D[i - 1, :-2]
            mask = np.isin(pair_comp, list(ONE_TWO[b]))
            row = np.where(mask, np.minimum(row, prev2 + 1), row)
        if i >= 2:
            a = int(q[i - 2])
            if a < NRES and b < NRES:
                qc = min(a, b) * 32 + max(a, b)
                if n >= 2:
                    prev22 = np.full(n + 1, INF, dtype=np.int64)
                    prev22[2:] = D[i - 2, :-2]
                    mask = np.zeros(n + 1, dtype=bool)
                    if a != b:
                        mask |= (t0 == b) & (t1 == a)
                    if qc in TWO_TWO:
                        mask |= np.isin(pair_comp, list(TWO_TWO[qc]))
                    row = np.where(mask, np.minimum(row, prev22 + 1), row)
                if qc in TWO_ONE:
                    prev21 = np.full(n + 1, INF, dtype=np.int64)
                    prev21[1:] = D[i - 2, :-1]
                    mask = ok1 & np.isin(t1, list(TWO_ONE[qc]))
                    row = np.where(mask, np.minimum(row, prev21 + 1), row)
        row[0] = INF
        D[i] = np.where(row > kmax, INF, row)
    return D


@dataclass(frozen=True, order=True)
class Placement:
    contig: str
    strand: str
    start: int  # genome, in [0, length)
    end: int  # unwrapped
    blocks: tuple[tuple[int, int], tuple[int, int]]  # ascending, same shift as start
    prefix_block: tuple[int, int]
    cost: int
    types: tuple[str, ...]  # sorted names, "frameshift" included
    frame: int

    @property
    def blocks_str(self) -> str:
        return ";".join(f"{a}-{b}" for a, b in self.blocks)

    @property
    def diff_types(self) -> str:
        return ";".join(sorted(self.types))


@dataclass(frozen=True)
class Result:
    cost: int | None
    placements: tuple[Placement, ...]


class _FrameView:
    """One (possibly decoy) frame read in its own direction with mirrored codon positions."""

    def __init__(self, fid: int, frame: Frame, kind: str | None) -> None:
        self.fid = fid
        self.frame = frame
        self.n = len(frame.aa)
        self.o = frame.offset
        self.rev = kind == "reverse"
        self.text = decode(frame.aa)

    def cs(self, i: np.ndarray | int) -> np.ndarray | int:
        """Strand-coordinate start of the codon carrying residue ``i``."""
        return self.o + 3 * (self.n - 1 - i) if self.rev else self.o + 3 * i

    def inv(self, v: np.ndarray) -> np.ndarray:
        """Residue index whose codon starts at ``v`` (float; integral when valid)."""
        if self.rev:
            return self.n - 1 - (v - self.o) / 3
        return (v - self.o) / 3

    def strand_iv(self, s: int, e: int) -> tuple[int, int]:
        """Strand-coordinate interval of residues ``[s, e)``."""
        a, b = int(self.cs(s)), int(self.cs(e - 1))
        return min(a, b), max(a, b) + 3


def _to_genome(frame: Frame, a: int, b: int) -> tuple[int, int]:
    if frame.strand == "+":
        return a, b
    return frame.ext_nt - b, frame.ext_nt - a


def _rank_key(types: Sequence[str]) -> tuple[int, ...]:
    return tuple(sorted(TYPE_ORDER.index(t) for t in types))


def frameshift_oracle(
    peptide: str,
    index: SixFrameIndex,
    frames: Sequence[Frame],
    kind: str | None,
    k: int,
    min_side: int,
) -> Result:
    """All split placements of ``peptide`` at the minimal cost (see the module docstring)."""
    q = encode_query(peptide)
    letters = decode(q)
    L = len(q)
    if k < 1 or 2 * min_side > L:
        return Result(None, ())
    views = [_FrameView(fid, f, kind) for fid, f in enumerate(frames)]
    groups: dict[tuple[str, str], list[_FrameView]] = {}
    for v in views:
        groups.setdefault((v.frame.contig, v.frame.strand), []).append(v)
    rem = k - 1
    pre = {v.fid: semi_global_matrix(q, v.frame.aa, rem) for v in views}
    suf = {v.fid: semi_global_matrix(q[::-1], v.frame.aa[::-1], rem) for v in views}

    # (contig, strand, prefix strand-iv, suffix strand-iv) -> best (cost, ranks, prefix fid)
    cands: dict[tuple[str, str, tuple[int, int], tuple[int, int]], tuple[int, tuple[int, ...]]] = {}
    fid_of: dict[tuple[str, str, tuple[int, int], tuple[int, int]], int] = {}
    side_cache: dict[tuple[int, str, int, int], tuple[int, tuple[int, ...]] | None] = {}

    def side(text: str, fid: int, qs: str, s: int, e: int) -> tuple[int, tuple[int, ...]] | None:
        key = (fid, qs, s, e)
        if key not in side_cache:
            side_cache[key] = align_window(qs, text, s, e, None, rem) if s < e else None
        return side_cache[key]

    for (contig, strand), vs in groups.items():
        sign = -1 if kind == "reverse" else 1
        for F in vs:
            for G in vs:
                for delta in (1, -1):
                    J = np.arange(1, F.n + 1)
                    want = np.asarray(F.cs(J - 1)) + sign * (3 + delta)
                    jp = G.inv(want)
                    ok = (jp == np.round(jp)) & (jp >= 0) & (jp < G.n)
                    for p in range(min_side, L - min_side + 1):
                        pc = pre[F.fid][p, J]
                        jj = np.where(ok, np.round(jp), 0).astype(np.int64)
                        sc = suf[G.fid][L - p, G.n - jj]
                        tot = 1 + pc + sc
                        hit = np.nonzero(ok & (pc < INF) & (sc < INF) & (tot <= k))[0]
                        for h in hit:
                            j, j2 = int(J[h]), int(jj[h])
                            qa, qb = letters[:p], letters[p:]
                            for s in range(max(0, j - p - 2 * rem), j - p + 2 * rem + 1):
                                ga = side(F.text, F.fid, qa, s, j)
                                if ga is None:
                                    continue
                                for e in range(j2 + (L - p) - 2 * rem, j2 + (L - p) + 2 * rem + 1):
                                    if e <= j2 or e > G.n:
                                        continue
                                    gb = side(G.text, G.fid, qb, j2, e)
                                    if gb is None or 1 + ga[0] + gb[0] > k:
                                        continue
                                    cost = 1 + ga[0] + gb[0]
                                    ranks = tuple(
                                        sorted((*ga[1], *gb[1], TYPE_ORDER.index("frameshift")))
                                    )
                                    key = (contig, strand, F.strand_iv(s, j), G.strand_iv(j2, e))
                                    best = (cost, ranks)
                                    if key not in cands or best < cands[key]:
                                        cands[key] = best
                                    fid_of[key] = F.fid
    return _finish(index, frames, groups, cands, fid_of, letters)


def _finish(
    index: SixFrameIndex,
    frames: Sequence[Frame],
    groups: dict[tuple[str, str], list[_FrameView]],
    cands: dict[tuple[str, str, tuple[int, int], tuple[int, int]], tuple[int, tuple[int, ...]]],
    fid_of: dict[tuple[str, str, tuple[int, int], tuple[int, int]], int],
    letters: str,
) -> Result:
    placements: dict[tuple[str, str, int, int, str], Placement] = {}
    for key, (cost, ranks) in cands.items():
        contig, strand, piv, siv = key
        rep = index.replicons[contig]
        Ln = rep.length
        frame = frames[fid_of[key]]
        a = min(piv[0], siv[0])
        b = max(piv[1], siv[1])
        if b - a > Ln:
            continue
        if _suppressed(letters, groups[(contig, strand)], (a, b), cost, Ln, rep.topology):
            continue
        gp = _to_genome(frame, *piv)
        gs = _to_genome(frame, *siv)
        g0 = min(gp[0], gs[0])
        shift = (g0 % Ln) - g0 if rep.topology == "circular" else 0
        gp = (gp[0] + shift, gp[1] + shift)
        gs = (gs[0] + shift, gs[1] + shift)
        blocks = tuple(sorted((gp, gs)))
        start, end = blocks[0][0], blocks[1][1]
        pstart = gp[0] % Ln
        label = index.frame_label(Interval(contig, pstart, pstart + gp[1] - gp[0], strand))  # type: ignore[arg-type]
        types = tuple(sorted(TYPE_ORDER[r] for r in ranks))
        pl = Placement(contig, strand, start, end, blocks, gp, cost, types, label)  # type: ignore[arg-type]
        pkey = (contig, strand, start, end, pl.blocks_str)
        old = placements.get(pkey)
        if old is None or (pl.cost, _rank_key(pl.types)) < (old.cost, _rank_key(old.types)):
            placements[pkey] = pl
    if not placements:
        return Result(None, ())
    best = min(p.cost for p in placements.values())
    keep = tuple(sorted(p for p in placements.values() if p.cost == best))
    return Result(best, keep)


@cache
def _nonfs_windows(letters: str, text: str, kmax: int) -> tuple[tuple[int, int, int], ...]:
    """Every ``(s, e, cost)`` window of ``text`` aligning the whole query at cost <= ``kmax``."""
    from gippa.mapping.sixframe import encode_target

    q = encode_query(letters)
    t = encode_target(text)
    L = len(q)
    fwd = semi_global_costs(q, t, kmax)
    rev = semi_global_costs(q[::-1], t[::-1], kmax)[::-1]
    E = np.nonzero(fwd <= kmax)[0]
    S = np.nonzero(rev <= kmax)[0]
    out = []
    for s in S:
        for e in E:
            if not L - 2 * kmax <= e - s <= L + 2 * kmax:
                continue
            got = align_window(letters, text, int(s), int(e), None, kmax)
            if got is not None:
                out.append((int(s), int(e), got[0]))
    return tuple(out)


def _suppressed(
    letters: str,
    views: list[_FrameView],
    locus: tuple[int, int],
    cost: int,
    length: int,
    topology: str,
) -> bool:
    a, b = locus
    for v in views:
        for s, e, c in _nonfs_windows(letters, v.text, cost):
            if c > cost:
                continue
            x, y = v.strand_iv(s, e)
            if topology == "circular":
                if overlaps_mod(a, b - a, x, y - x, length):
                    return True
            elif a < y and x < b:
                return True
    return False
