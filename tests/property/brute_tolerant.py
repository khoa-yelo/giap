"""Brute-force tolerant matcher over whole frames (issue M3-07; GIPPA_SPEC.md §7.2.4, §11.2).

Written by the test author from the spec text only (it does not import ``gippa.mapping.diffops``
or ``gippa.mapping.tolerant_dp``). Two parts:

- :func:`semi_global_costs`: an unbanded numpy DP of the query against *every* window of one
  frame (free start and end in the frame). Row ``i`` holds, for every frame position ``j``, the
  minimal cost of aligning ``query[:i]`` to some ``frame[s:j]``. The operations are those of
  §7.2.4 (match, substitution, adjacent swap, mass-equivalent block, near-equal block; each
  costs 1, a target ``*``/``X`` is never consumed and a query NOMATCH only by a substitution).
- :func:`brute_force`: the minimal cost ``C`` over all frames of one search and every window
  ``[s, e)`` reaching it (start set from the DP on the reversed query and frame, end set from the
  forward DP, each pair confirmed with the M3-06 test oracle ``align_window`` unbanded), grouped
  into clusters of overlapping windows per frame.

A tolerant result is *complete and sound* when its minimal cost equals ``C``, every cluster is
covered by a reported locus (same contig, strand and frame phase, overlapping modulo the replicon
length), and every reported window really costs ``C``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from oracle_tolerant import align_window

from gippa.mapping.decoys import mirror_aa
from gippa.mapping.sixframe import ALPHABET, Frame, SixFrameIndex, decode, encode_query

INF = 1 << 20
NRES = 19  # residue codes 0..18; 19 '*', 20 'X', 21 query NOMATCH

# §7.2.4, as written (compositions, both orders of each pair)
_BLOCK_RAW = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_RAW = (("W", "AD"), ("W", "EG"), ("R", "GV"))


def _code(letter: str) -> int:
    return ALPHABET.index("L" if letter == "I" else letter)


def _comp(seg: str) -> int:
    a, b = sorted(_code(c) for c in seg)
    return a * 32 + b


ONE_TWO: dict[int, set[int]] = {}  # single residue -> compositions of the 2-residue partner
TWO_ONE: dict[int, set[int]] = {}  # 2-residue composition -> single residues
TWO_TWO: dict[int, set[int]] = {}  # 2-residue composition -> 2-residue compositions
for _a, _b in (*_BLOCK_RAW, *_NEAR_RAW):
    if len(_a) == 1:
        ONE_TWO.setdefault(_code(_a), set()).add(_comp(_b))
        TWO_ONE.setdefault(_comp(_b), set()).add(_code(_a))
    else:
        TWO_TWO.setdefault(_comp(_a), set()).add(_comp(_b))
        TWO_TWO.setdefault(_comp(_b), set()).add(_comp(_a))


def semi_global_costs(q: np.ndarray, t: np.ndarray, kmax: int) -> np.ndarray:
    """``D[L, j]`` for ``j = 0..n``: minimal cost of ``q`` against any window ending at ``j``.

    Values above ``kmax`` are reported as ``INF``.
    """
    q = np.asarray(q, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    L, n = len(q), len(t)
    D = np.full((L + 1, n + 1), INF, dtype=np.int64)
    D[0, :] = 0
    tok = t < NRES
    pair_ok = np.zeros(n + 1, dtype=bool)  # target pair t[j-2], t[j-1] usable, indexed by j
    pair_comp = np.full(n + 1, -1, dtype=np.int64)
    if n >= 2:
        lo = np.minimum(t[:-1], t[1:])
        hi = np.maximum(t[:-1], t[1:])
        pair_ok[2:] = tok[:-1] & tok[1:]
        pair_comp[2:] = np.where(pair_ok[2:], lo * 32 + hi, -1)
    t1 = np.full(n + 1, -1, dtype=np.int64)  # t[j-1]
    t1[1:] = t
    t0 = np.full(n + 1, -1, dtype=np.int64)  # t[j-2]
    t0[2:] = t[:-1]
    ok1 = np.zeros(n + 1, dtype=bool)
    ok1[1:] = tok
    for i in range(1, L + 1):
        row = np.full(n + 1, INF, dtype=np.int64)
        b = int(q[i - 1])
        # match / substitution from (i-1, j-1)
        prev = np.full(n + 1, INF, dtype=np.int64)
        prev[1:] = D[i - 1, :-1]
        step = np.where(t1 == b, 0, 1)
        row = np.where(ok1, np.minimum(row, prev + step), row)
        # 1 -> 2 exchange from (i-1, j-2)
        if b < NRES and b in ONE_TWO and n >= 2:
            prev2 = np.full(n + 1, INF, dtype=np.int64)
            prev2[2:] = D[i - 1, :-2]
            mask = np.isin(pair_comp, list(ONE_TWO[b]))
            row = np.where(mask, np.minimum(row, prev2 + 1), row)
        if i >= 2:
            a = int(q[i - 2])
            if a < NRES and b < NRES:
                qc = min(a, b) * 32 + max(a, b)
                # swap and 2 <-> 2 exchange from (i-2, j-2)
                if n >= 2:
                    prev22 = np.full(n + 1, INF, dtype=np.int64)
                    prev22[2:] = D[i - 2, :-2]
                    mask = np.zeros(n + 1, dtype=bool)
                    if a != b:
                        mask |= (t0 == b) & (t1 == a)
                    if qc in TWO_TWO:
                        mask |= np.isin(pair_comp, list(TWO_TWO[qc]))
                    row = np.where(mask, np.minimum(row, prev22 + 1), row)
                # 2 -> 1 exchange from (i-2, j-1)
                if qc in TWO_ONE:
                    prev21 = np.full(n + 1, INF, dtype=np.int64)
                    prev21[1:] = D[i - 2, :-1]
                    mask = ok1 & np.isin(t1, list(TWO_ONE[qc]))
                    row = np.where(mask, np.minimum(row, prev21 + 1), row)
        row[0] = INF
        D[i] = np.where(row > kmax, INF, row)
    out: np.ndarray = D[L]
    return out


@dataclass(frozen=True)
class Brute:
    cost: int | None  # minimal cost over all frames (None when above k)
    clusters: tuple[tuple[tuple[int, int, int], ...], ...]  # overlapping (frame_id, s, e) windows


def brute_force(peptide: str, frames: Sequence[Frame], k: int) -> Brute:
    """Minimal cost of ``peptide`` over every window of ``frames`` and the windows reaching it."""
    q = encode_query(peptide)
    L = len(q)
    qletters = decode(q)
    ends: list[np.ndarray] = []
    starts: list[np.ndarray] = []
    best = INF
    for f in frames:
        fwd = semi_global_costs(q, f.aa, k)
        rev = semi_global_costs(q[::-1], f.aa[::-1], k)[::-1]  # rev[s]: best window starting at s
        ends.append(fwd)
        starts.append(rev)
        best = min(best, int(fwd.min()))
    if best >= INF:
        return Brute(None, ())
    clusters: list[tuple[tuple[int, int, int], ...]] = []
    for fid, f in enumerate(frames):
        E = np.nonzero(ends[fid] == best)[0]
        S = np.nonzero(starts[fid] == best)[0]
        if len(E) == 0:
            continue
        text = decode(f.aa)
        wins: list[tuple[int, int]] = []
        for s in S:
            for e in E:
                if not L - best <= e - s <= L + best:
                    continue
                got = align_window(qletters, text, int(s), int(e), None, best)
                if got is not None and got[0] == best:
                    wins.append((int(s), int(e)))
        assert wins, "brute force: an end without a start (oracle disagreement)"
        wins.sort()
        cur: list[tuple[int, int, int]] = []
        reach = -1
        for s, e in wins:
            if cur and s >= reach:
                clusters.append(tuple(cur))
                cur = []
                reach = -1
            cur.append((fid, s, e))
            reach = max(reach, e)
        clusters.append(tuple(cur))
    return Brute(best, tuple(clusters))


def window_interval(
    index: SixFrameIndex, frames: Sequence[Frame], kind: str | None, fid: int, s: int, e: int
) -> tuple[str, int, int, str, int]:
    """(contig, wrapped start, length, strand, frame label) of a window (mirrored for decoys)."""
    a, b = (s, e) if kind is None else mirror_aa(kind, len(frames[fid].aa), s, e)  # type: ignore[arg-type]
    iv = index.to_interval(fid, a, b)
    length = index.replicons[iv.contig].length
    return iv.contig, iv.start % length, iv.end - iv.start, iv.strand, index.frame_label(iv)


def overlaps_mod(a_start: int, a_len: int, b_start: int, b_len: int, length: int) -> bool:
    """Do ``[a, a+a_len)`` and ``[b, b+b_len)`` overlap on a circle of ``length`` (or a line)?"""
    for shift in (-length, 0, length):
        b0 = b_start + shift
        if a_start < b0 + b_len and b0 < a_start + a_len:
            return True
    return False
