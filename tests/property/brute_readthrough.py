"""Brute-force tolerant matcher with consumable stop cells (M7-16 rule 4; ADR-0058 item 1).

Written by the test author from docs/specs/M7-16-stop-codon-readthrough.md only. It does not
import ``giap.mapping.readthrough``, ``giap.mapping.tolerant_dp`` or ``giap.mapping.diffops``,
and the M3-06/M3-07 oracles (``oracle_tolerant.py``, ``brute_tolerant.py``) are not edited.

It is the unbanded M3-07 brute force (``brute_tolerant.semi_global_costs``) with a boolean mask of
**consumable cells** per frame. At a marked cell the frame holds ``*`` (code 19); the only
operation that may consume it is a substitution (one query residue, any code, ``NOMATCH``
included) at cost 1. A marked cell never matches, and no swap, block or near-block uses it.
Every other ``*`` and every ``X`` stays unconsumable, exactly as in §7.2.4.

- :func:`costs_rt`: row ``D[L]`` of the DP of the query against every window of one frame (free
  start, or a fixed start ``s``), values above ``kmax`` reported as ``INF``.
- :func:`window_cost_rt`: the minimal cost of the query against exactly ``frame[s:e]``.
- :func:`brute_force_rt`: the minimal cost over all frames of one search and every window
  reaching it, grouped into clusters of overlapping windows per frame (as ``brute_force``).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from giap.mapping.sixframe import ALPHABET, Frame, encode_query

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


ONE_TWO: dict[int, set[int]] = {}
TWO_ONE: dict[int, set[int]] = {}
TWO_TWO: dict[int, set[int]] = {}
for _a, _b in (*_BLOCK_RAW, *_NEAR_RAW):
    if len(_a) == 1:
        ONE_TWO.setdefault(_code(_a), set()).add(_comp(_b))
        TWO_ONE.setdefault(_comp(_b), set()).add(_code(_a))
    else:
        TWO_TWO.setdefault(_comp(_a), set()).add(_comp(_b))
        TWO_TWO.setdefault(_comp(_b), set()).add(_comp(_a))


def _shift(row: np.ndarray, by: int) -> np.ndarray:
    out = np.full(row.shape[0], INF, dtype=np.int64)
    out[by:] = row[: row.shape[0] - by]
    return out


def costs_rt(
    q: np.ndarray, t: np.ndarray, cells: np.ndarray, kmax: int, start: int | None = None
) -> np.ndarray:
    """``D[L, j]`` for ``j = 0..n``: minimal cost of ``q`` against a window ending at ``j``.

    ``cells`` is a boolean mask over ``t`` (consumable cells; each must hold ``*``). With
    ``start`` the window must start at ``start``; otherwise the start is free.
    """
    q = np.asarray(q, dtype=np.int64)
    t = np.asarray(t, dtype=np.int64)
    cells = np.asarray(cells, dtype=bool)
    L, n = len(q), len(t)
    if cells.shape != (n,):
        raise ValueError("cell mask must have the frame's length")
    if np.any(t[cells] != NRES):
        raise ValueError("a consumable cell must hold '*'")
    D = np.full((L + 1, n + 1), INF, dtype=np.int64)
    if start is None:
        D[0, :] = 0
    else:
        D[0, start] = 0
    tok = t < NRES
    pair_comp = np.full(n + 1, -1, dtype=np.int64)
    if n >= 2:
        lo = np.minimum(t[:-1], t[1:])
        hi = np.maximum(t[:-1], t[1:])
        ok2 = tok[:-1] & tok[1:]
        pair_comp[2:] = np.where(ok2, lo * 32 + hi, -1)
    t1 = np.full(n + 1, -1, dtype=np.int64)
    t1[1:] = t
    t0 = np.full(n + 1, -1, dtype=np.int64)
    t0[2:] = t[:-1]
    res1 = np.zeros(n + 1, dtype=bool)  # t[j-1] is a residue
    res1[1:] = tok
    cell1 = np.zeros(n + 1, dtype=bool)  # t[j-1] is a consumable cell
    cell1[1:] = cells
    for i in range(1, L + 1):
        row = np.full(n + 1, INF, dtype=np.int64)
        b = int(q[i - 1])
        prev = _shift(D[i - 1], 1)
        # match / substitution over a residue
        row = np.where(res1, np.minimum(row, prev + np.where(t1 == b, 0, 1)), row)
        # the read-through case: any query code over a marked cell, cost 1
        row = np.where(cell1, np.minimum(row, prev + 1), row)
        if b < NRES and b in ONE_TWO and n >= 2:
            mask = np.isin(pair_comp, list(ONE_TWO[b]))
            row = np.where(mask, np.minimum(row, _shift(D[i - 1], 2) + 1), row)
        if i >= 2:
            a = int(q[i - 2])
            if a < NRES and b < NRES:
                qc = min(a, b) * 32 + max(a, b)
                if n >= 2:
                    mask = np.zeros(n + 1, dtype=bool)
                    if a != b:
                        mask |= (t0 == b) & (t1 == a)
                    if qc in TWO_TWO:
                        mask |= np.isin(pair_comp, list(TWO_TWO[qc]))
                    row = np.where(mask, np.minimum(row, _shift(D[i - 2], 2) + 1), row)
                if qc in TWO_ONE:
                    mask = res1 & np.isin(t1, list(TWO_ONE[qc]))
                    row = np.where(mask, np.minimum(row, _shift(D[i - 2], 1) + 1), row)
        row[0] = INF
        D[i] = np.where(row > kmax, INF, row)
    out: np.ndarray = D[L]
    return out


def window_cost_rt(
    q: np.ndarray, t: np.ndarray, cells: np.ndarray, s: int, e: int, kmax: int
) -> int | None:
    """Minimal cost (<= ``kmax``) of ``q`` against exactly ``t[s:e]``, else None."""
    if not 0 <= s < e <= len(t):
        return None
    c = int(costs_rt(q, t, cells, kmax, start=s)[e])
    return None if c >= INF else c


def consumes_cell(
    q: np.ndarray, t: np.ndarray, cells: np.ndarray, s: int, e: int, cost: int
) -> bool:
    """Does every alignment of ``q`` to ``t[s:e]`` at ``cost`` need a cell?"""
    plain = np.zeros(len(t), dtype=bool)
    c = window_cost_rt(q, t, plain, s, e, cost)
    return c is None or c > cost


@dataclass(frozen=True)
class Brute:
    cost: int | None
    clusters: tuple[tuple[tuple[int, int, int], ...], ...]


def brute_force_rt(
    peptide: str, frames: Sequence[Frame], cells: Sequence[np.ndarray], k: int
) -> Brute:
    """Minimal cost of ``peptide`` over every window of ``frames`` with consumable ``cells``
    (one boolean mask per frame) and the windows reaching it, clustered per frame."""
    q = encode_query(peptide)
    L = len(q)
    ends: list[np.ndarray] = []
    starts: list[np.ndarray] = []
    best = INF
    for f, m in zip(frames, cells, strict=True):
        fwd = costs_rt(q, f.aa, m, k)
        rev = costs_rt(q[::-1], f.aa[::-1], m[::-1], k)[::-1]
        ends.append(fwd)
        starts.append(rev)
        best = min(best, int(fwd.min()))
    if best >= INF:
        return Brute(None, ())
    clusters: list[tuple[tuple[int, int, int], ...]] = []
    for fid, (f, m) in enumerate(zip(frames, cells, strict=True)):
        E = np.nonzero(ends[fid] == best)[0]
        S = np.nonzero(starts[fid] == best)[0]
        if len(E) == 0:
            continue
        wins: list[tuple[int, int]] = []
        for s in S:
            for e in E:
                if not L - 2 * best <= e - s <= L + 2 * best:
                    continue
                if window_cost_rt(q, f.aa, m, int(s), int(e), best) == best:
                    wins.append((int(s), int(e)))
        assert wins, "brute force: an end without a start"
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


def masks(frames: Sequence[Frame], cells: Sequence[np.ndarray]) -> list[np.ndarray]:
    """Boolean masks from per-frame position arrays."""
    out = []
    for f, c in zip(frames, cells, strict=True):
        m = np.zeros(len(f.aa), dtype=bool)
        m[np.asarray(c, dtype=np.int64)] = True
        out.append(m)
    return out
