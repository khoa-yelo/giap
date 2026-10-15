"""Banded tolerant verification of candidate diagonals (GIAP_SPEC.md §7.2.4, G4).

The whole query ``q[0:L]`` is aligned globally to a target substring ``t[s:e]`` (free ends) of
one frame inside a band of ``±2`` around the candidate diagonal ``d0``: every DP cell ``(i, j)``
satisfies ``|(j - i) - d0| <= 2``. The only operations are a match (cost 0; codes already fold
I into L), a substitution, an adjacent swap ``XY <-> YX`` (``X != Y``), a mass-equivalent block
(1↔2, 2↔1, 2↔2 from ``BLOCKS``) and a near-equal block (1↔2, 2↔1 from ``NEAR_BLOCKS``), each
costing one difference. A target ``*`` (``STOP``) or ``X`` (``AMBIG``) is never consumed, and a
query ``NOMATCH`` residue can only be consumed by a substitution.

Tie-breaks (docs/specs/M3-06-tolerant-verifier.md): lowest cost, then the lexicographically
smallest diff-type multiset under the ``DIFF_TYPES`` precedence (substitution < swap < block <
near_block), then the smallest ``t_start`` (then the smallest ``t_end``). The multiset order is
encoded in one integer key per DP cell, which is additive over operations, so the DP keeps the
exact optimum without a traceback: with base ``B`` the key is
``cost·B³ + (n_swap + n_block + n_near)·B² + (n_block + n_near)·B + n_near``.

Every window inside one ``±2`` band overlaps the others for queries of more than four residues,
so a call reports at most one alignment (one locus per candidate diagonal). Queries of at most
``2·BAND`` (4) residues are rejected: the band can then hold two disjoint best windows (``AK`` in
``AKGAK``); callers with such short pieces need their own anchored check (M3-08).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit

from giap.mapping.diffops import DIFF_TYPES, BlockTables, DiffType
from giap.mapping.sixframe import MAX_DIFF, STOP

__all__ = ["BAND", "Alignment", "decode_key", "exchange_lookup", "verify", "verify_batch"]

BAND: int = MAX_DIFF
_NDIAG: int = 2 * BAND + 1
_NRES: int = STOP  # residue codes 0..18; every code >= 19 is a stop, X or NOMATCH
_BASE: int = 1 << 15  # per-component base of the cell key (cost < 2**15)
_RANK_SUB, _RANK_SWAP, _RANK_BLOCK, _RANK_NEAR = 0, 1, 2, 3
# offsets of the 1->2, 2->1 and 2->2 exchange tables inside the flat lookup array
_OFF12: int = 0
_OFF21: int = _NRES**3
_OFF22: int = 2 * _NRES**3
_EX_SIZE: int = 2 * _NRES**3 + _NRES**4


@dataclass(frozen=True, slots=True)
class Alignment:
    """One tolerant hit: ``frame[t_start:t_end]`` (0-based half-open) at ``cost`` differences.

    ``diff_types`` lists one entry per difference, sorted by the ``DIFF_TYPES`` precedence.
    """

    cost: int
    t_start: int
    t_end: int
    diff_types: tuple[DiffType, ...]


def exchange_lookup(tables: BlockTables) -> np.ndarray:
    """Flat int8 lookup ``(query segment, target segment) -> rank`` (-1 when no exchange).

    ``BlockTables`` rows are directional, so each row is entered in both orientations. This is
    the ``ex`` argument of :func:`verify_batch`.
    """
    ex = np.full(_EX_SIZE, -1, dtype=np.int8)
    n = _NRES
    for r in range(tables.kind.shape[0]):
        rank = _RANK_BLOCK if int(tables.kind[r]) == 0 else _RANK_NEAR
        left = [int(c) for c in tables.left[r, : int(tables.left_len[r])]]
        right = [int(c) for c in tables.right[r, : int(tables.right_len[r])]]
        for x, y in ((left, right), (right, left)):
            if any(c >= n for c in x + y):
                continue
            if len(x) == 1 and len(y) == 2:
                idx = _OFF12 + (x[0] * n + y[0]) * n + y[1]
            elif len(x) == 2 and len(y) == 1:
                idx = _OFF21 + (x[0] * n + x[1]) * n + y[0]
            elif len(x) == 2 and len(y) == 2:
                idx = _OFF22 + ((x[0] * n + x[1]) * n + y[0]) * n + y[1]
            else:
                continue
            if ex[idx] < 0 or rank < ex[idx]:
                ex[idx] = rank
    return ex


@njit(cache=True, nogil=True)
def _rank_weight(rank: int, b: int) -> int:  # pragma: no cover - compiled by numba; via verify
    """Key increment of one difference of ``rank`` (0 substitution .. 3 near_block)."""
    b2 = b * b
    w = b2 * b
    if rank >= 1:
        w += b2
    if rank >= 2:
        w += b
    if rank >= 3:
        w += 1
    return w


@njit(cache=True, nogil=True)
def _verify_kernel(
    q: np.ndarray, frame: np.ndarray, ex: np.ndarray, d0: int, k: int
) -> tuple[int, int, int]:  # pragma: no cover - compiled by numba; exercised through verify
    """Best ``(key, t_start, t_end)`` in the band; ``(-1, -1, -1)`` when cost > ``k``."""
    size = (q.shape[0] + 1) * 5
    key = np.empty(size, dtype=np.int64)
    st = np.empty(size, dtype=np.int64)
    return _verify_core(q, frame, ex, d0, k, key, st)


@njit(cache=True, nogil=True)
def _verify_core(
    q: np.ndarray,
    frame: np.ndarray,
    ex: np.ndarray,
    d0: int,
    k: int,
    key: np.ndarray,
    st: np.ndarray,
) -> tuple[int, int, int]:  # pragma: no cover - compiled by numba; exercised through verify
    """:func:`_verify_kernel` on caller-owned DP buffers of ``(len(q) + 1) * 5`` cells each.

    Every cell read is written first, so the buffers may be reused between calls.
    """
    nres = 19
    nd = 5
    b = 32768
    big = 1 << 62
    L = q.shape[0]
    n = frame.shape[0]
    limit = (k + 1) * b * b * b  # key < limit <=> cost <= k
    w_sub = _rank_weight(0, b)
    w_swap = _rank_weight(1, b)
    off21 = nres * nres * nres
    off22 = 2 * off21
    lo = d0 - 2
    for dd in range(nd):
        j = lo + dd
        if 0 <= j <= n:
            key[dd] = 0
            st[dd] = j
        else:
            key[dd] = big
            st[dd] = 0
    for i in range(1, L + 1):
        row = i * nd
        prow = row - nd
        prow2 = row - 2 * nd
        q1 = int(q[i - 1])  # Python int: uint8 arithmetic overflows without numba
        q0 = int(q[i - 2]) if i >= 2 else nres
        for dd in range(nd):
            j = i + lo + dd
            best = big
            bs = 0
            if 1 <= j <= n:
                t1 = int(frame[j - 1])
                t0 = int(frame[j - 2]) if j >= 2 else nres
                if t1 < nres:
                    # match / substitution from (i-1, j-1)
                    pk = key[prow + dd]
                    if pk < big:
                        c = pk + (0 if q1 == t1 else w_sub)
                        if c < limit:
                            best = c
                            bs = st[prow + dd]
                    # 2->2 swap or block from (i-2, j-2)
                    if i >= 2 and t0 < nres and q0 < nres and q1 < nres:
                        pk = key[prow2 + dd]
                        if pk < big:
                            w = big
                            if q0 != q1 and q0 == t1 and q1 == t0:
                                w = w_swap
                            r = int(ex[off22 + ((q0 * nres + q1) * nres + t0) * nres + t1])
                            if r >= 0:
                                rw = _rank_weight(r, b)
                                if rw < w:
                                    w = rw
                            if w < big:
                                c = pk + w
                                s = st[prow2 + dd]
                                if c < limit and (c < best or (c == best and s < bs)):
                                    best = c
                                    bs = s
                    # 1->2 block from (i-1, j-2): diagonal dd-1
                    if dd >= 1 and t0 < nres and q1 < nres:
                        pk = key[prow + dd - 1]
                        if pk < big:
                            r = int(ex[(q1 * nres + t0) * nres + t1])
                            if r >= 0:
                                c = pk + _rank_weight(r, b)
                                s = st[prow + dd - 1]
                                if c < limit and (c < best or (c == best and s < bs)):
                                    best = c
                                    bs = s
                    # 2->1 block from (i-2, j-1): diagonal dd+1
                    if i >= 2 and dd + 1 < nd and q0 < nres and q1 < nres:
                        pk = key[prow2 + dd + 1]
                        if pk < big:
                            r = int(ex[off21 + (q0 * nres + q1) * nres + t1])
                            if r >= 0:
                                c = pk + _rank_weight(r, b)
                                s = st[prow2 + dd + 1]
                                if c < limit and (c < best or (c == best and s < bs)):
                                    best = c
                                    bs = s
            key[row + dd] = best
            st[row + dd] = bs
    row = L * nd
    best = big
    bs = -1
    be = -1
    for dd in range(nd):
        c = key[row + dd]
        if c < best or (c == best and c < big and st[row + dd] < bs):
            best = c
            bs = st[row + dd]
            be = L + lo + dd
    if best >= big:
        return -1, -1, -1
    return best, bs, be


_TABLE_CACHE: dict[int, tuple[BlockTables, np.ndarray]] = {}
_TYPES_CACHE: dict[int, tuple[int, tuple[DiffType, ...]]] = {}


def _lookup_for(tables: BlockTables) -> np.ndarray:
    hit = _TABLE_CACHE.get(id(tables))
    if hit is not None and hit[0] is tables:
        return hit[1]
    ex = exchange_lookup(tables)
    _TABLE_CACHE[id(tables)] = (tables, ex)  # keeps ``tables`` alive, so the id stays unique
    return ex


def decode_key(key: int) -> tuple[int, tuple[DiffType, ...]]:
    """``(cost, diff_types)`` from a key of :func:`verify_batch` (counts per precedence rank)."""
    hit = _TYPES_CACHE.get(key)
    if hit is not None:
        return hit
    b = _BASE
    cost, rest = divmod(key, b**3)
    ge1, rest = divmod(rest, b * b)
    ge2, ge3 = divmod(rest, b)
    counts = (cost - ge1, ge1 - ge2, ge2 - ge3, ge3)
    types: tuple[DiffType, ...] = tuple(
        name for name, c in zip(DIFF_TYPES[:4], counts, strict=True) for _ in range(c)
    )
    out = (cost, types)
    _TYPES_CACHE[key] = out
    return out


def verify(
    query: np.ndarray, frame: np.ndarray, d0: int, k: int, tables: BlockTables
) -> list[Alignment]:
    """Minimal-cost alignment (cost ≤ ``k``) of ``query`` to ``frame`` in the band around ``d0``.

    ``query`` and ``frame`` are M3-02 codes (:func:`~giap.mapping.sixframe.encode_query`,
    :func:`~giap.mapping.sixframe.encode_target`); ``d0`` is the candidate diagonal (target
    position of query residue 0). Returns ``[]`` when the minimum exceeds ``k``, otherwise one
    :class:`Alignment` (see the module docstring for the tie-breaks).

    Inputs of any integer dtype or stride are normalised to C-contiguous ``uint8`` so that a
    single compiled kernel specialisation serves every call.

    Raises:
        ValueError: if ``len(query) <= 2 * BAND`` (4 residues): the ``±2`` band can then hold
            two disjoint equally good windows, so one diagonal no longer identifies one locus.
    """
    L = query.shape[0]
    if L <= 2 * BAND:
        raise ValueError(f"query of {L} residues is too short: verify needs more than {2 * BAND}")
    if k < 0:
        return []
    if k > L:
        k = L  # every operation consumes a query residue
    if k >= _BASE:  # pragma: no cover - queries are far shorter than 2**15
        raise ValueError("k must be below 32768")
    q = np.ascontiguousarray(query, dtype=np.uint8)
    f = np.ascontiguousarray(frame, dtype=np.uint8)
    key, s, e = _verify_kernel(q, f, _lookup_for(tables), int(d0), int(k))
    if key < 0:
        return []
    cost, types = decode_key(int(key))
    return [Alignment(cost, int(s), int(e), types)]  # plain ints also when run without numba


@njit(cache=True, nogil=True)
def verify_batch(
    query: np.ndarray, frame: np.ndarray, d0s: np.ndarray, k: int, ex: np.ndarray
) -> np.ndarray:  # pragma: no cover - compiled by numba; exercised through the tests and M3-07
    """:func:`verify` for every candidate diagonal of ``d0s``, callable from njit code.

    ``query`` and ``frame`` are C-contiguous ``uint8`` codes, ``d0s`` ``int64`` and ``ex`` the
    :func:`exchange_lookup` of the block tables. Returns an ``int64`` array of shape
    ``(len(d0s), 3)``: rows ``(key, t_start, t_end)`` with ``key`` decoded by :func:`decode_key`,
    or ``(-1, -1, -1)`` when the band holds nothing at cost ``<= k``. Row ``r`` equals
    ``verify(query, frame, d0s[r], k, tables)``.

    Raises:
        ValueError: if ``len(query) <= 2 * BAND`` (4 residues), as :func:`verify`, or ``k`` is
            not below 32768.
    """
    L = query.shape[0]
    if L <= 4:
        raise ValueError("verify_batch: queries of at most 2 * BAND (4) residues are rejected")
    if k >= 32768:
        raise ValueError("verify_batch: k must be below 32768")
    n = d0s.shape[0]
    out = np.full((n, 3), -1, dtype=np.int64)
    if k < 0:
        return out
    kk = L if k > L else k  # every operation consumes a query residue
    size = (L + 1) * 5
    key = np.empty(size, dtype=np.int64)
    st = np.empty(size, dtype=np.int64)
    for r in range(n):
        c, s, e = _verify_core(query, frame, ex, int(d0s[r]), kk, key, st)
        if c >= 0:
            out[r, 0] = c
            out[r, 1] = s
            out[r, 2] = e
    return out
