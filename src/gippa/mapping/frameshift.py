"""Tolerant matching, ±1-nt frameshift splits (GIPPA_SPEC.md §7.2.4; issue M3-08).

A split at query position ``p`` places the prefix ``q[0:p]`` on the codons of a window of frame
``F`` that ends at the junction and the suffix ``q[p:L]`` on a window of frame ``G`` (same contig
and strand) whose first codon starts ``δ = ±1`` nt after the last prefix codon ends. Cost =
1 + the M3-06 cost of each side; the remaining budget ``k - 1`` (``k <= 2``) may be spent on
either side, and both sides hold at least ``min_side_aa`` query residues.

**Seeds (lossless).** In one strand's nucleotide coordinates a query q-gram at offset ``i`` that
hits the codon starting at ``x`` lies on the *nucleotide diagonal* ``x - 3 i``. The q-grams of
an aligned prefix starting at ``P0`` lie on ``P0`` (``P0 ± 3`` after a 1↔2 block), those of the
suffix on ``P0 + δ`` (± 3), so all surviving q-grams of one split lie within five adjacent
nucleotide diagonals (one difference at most besides the shift). A split loses at most the
``q - 1`` q-grams across the junction and ``q + 1`` per difference on its side; the threshold
``t`` is the minimum over the allowed split points. Every window of five diagonals starting at a
hit diagonal ``u`` with ``>= t`` distinct offsets yields the anchor range ``P0 in [u-4, u+4]``.
When no ``q`` gives ``t >= 1`` every anchor of the strand is verified.

**Verification (anchored).** For every anchor ``P0`` the prefix frame and start are fixed; a
banded DP (``±(k-1)`` diagonals, the M3-06 operations and key layout of
:mod:`gippa.mapping.tolerant_dp`) runs over the prefix, switches into the shifted frame at every
allowed split point (entry key = prefix key, junction recorded) and continues over the suffix.
Each end column gives the best split of one locus. Short sides (even ≤ 4 aa) are part of this
fixed-start DP, so :func:`~gippa.mapping.tolerant_dp.verify` is never called on them.

**Precedence.** A split is dropped when the whole query aligns without a frameshift (any frame of
the same contig and strand and search, cost ≤ the split's cost) on an interval overlapping the
split locus (modulo the replicon length). Per peptide and search all remaining split loci at the
minimal cost are reported, one row per interval.

**Decoys.** Searched identically, each frame read in its own direction: a ``reverse`` decoy frame
reversed again lies on the target codons, so the reverse-decoy search of ``Q`` is the split search
of ``reversed(Q)`` on the reversed decoy frames (prefix of ``Q`` = suffix of the reversed query);
``shuffle`` decoys keep the target codons.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from functools import cache
from math import comb
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numba import njit, prange

from gippa.config.model import TolerantConfig
from gippa.mapping.decoys import DecoyType
from gippa.mapping.diffops import block_tables
from gippa.mapping.hits import hits_frame
from gippa.mapping.qgram import MAX_QUERY
from gippa.mapping.sixframe import Frame, SixFrameIndex, decode, encode_query
from gippa.mapping.tolerant import Search, budget, decoy_marks, parse_max_diff
from gippa.mapping.tolerant_dp import decode_key, exchange_lookup, rank_weight
from gippa.model.coords import Interval, Replicon

if TYPE_CHECKING:
    from gippa.mapping.readthrough import StopCells

__all__ = ["band_row", "match_frameshift", "nonfs_overlap"]

Arr = np.ndarray
# search space: text, frame starts, frame lengths, frame offsets, frame groups, group frames
# (by offset), group extended lengths, replicon lengths, circular flags
Space = tuple[Arr, Arr, Arr, Arr, Arr, Arr, Arr, Arr, Arr]

_SEARCH_ORDER = {"target": 0, "decoy": 1}
_CAP = 16  # placements per query kept by the parallel pass; more are recomputed serially
# placement columns (7): key, P0, prefix end, suffix start, suffix end, group, delta
_QS = (4, 3)  # a q = 5 index costs more to build than its sharper seeds save
_VERIFY_COST = 20000.0  # relative cost of verifying one anchor range vs one q-gram hit


# --- banded split DP ---------------------------------------------------------------------------


@njit(cache=True, nogil=True)
def band_row(
    q: np.ndarray,
    i: int,
    base: int,
    text: np.ndarray,
    fa: int,
    fn: int,
    ex: np.ndarray,
    key: np.ndarray,
    jn: np.ndarray,
    nd: int,
    limit: int,
    lo_row: int,
) -> None:  # pragma: no cover - compiled by numba; exercised through match_frameshift
    """Fill DP row ``i`` (columns ``base + i + dd``) from rows ``>= lo_row`` (M3-06 operations).

    The frame is ``text[fa : fa + fn]``; ``jn`` carries the junction of the best predecessor.
    """
    nres = 19
    big = 1 << 62
    b = 32768
    w_sub = rank_weight(0, b)
    w_swap = rank_weight(1, b)
    off21 = nres * nres * nres
    off22 = 2 * off21
    row = i * nd
    prow = row - nd
    prow2 = row - 2 * nd
    one = i - 1 >= lo_row
    two = i - 2 >= lo_row
    q1 = int(q[i - 1])
    q0 = int(q[i - 2]) if two else nres
    for dd in range(nd):
        j = base + i + dd
        best = big
        bj = -1
        if 1 <= j <= fn:
            t1 = int(text[fa + j - 1])
            t0 = int(text[fa + j - 2]) if j >= 2 else nres
            if t1 < nres:
                if one:
                    pk = key[prow + dd]
                    if pk < big:
                        c = pk + (0 if q1 == t1 else w_sub)
                        if c < limit:
                            best = c
                            bj = jn[prow + dd]
                if two and t0 < nres and q0 < nres and q1 < nres:
                    pk = key[prow2 + dd]
                    if pk < big:
                        w = big
                        if q0 != q1 and q0 == t1 and q1 == t0:
                            w = w_swap
                        r = int(ex[off22 + ((q0 * nres + q1) * nres + t0) * nres + t1])
                        if r >= 0:
                            rw = rank_weight(r, b)
                            if rw < w:
                                w = rw
                        if w < big and pk + w < limit and pk + w < best:
                            best = pk + w
                            bj = jn[prow2 + dd]
                if one and dd >= 1 and t0 < nres and q1 < nres:
                    pk = key[prow + dd - 1]
                    r = int(ex[(q1 * nres + t0) * nres + t1])
                    if pk < big and r >= 0:
                        c = pk + rank_weight(r, b)
                        if c < limit and c < best:
                            best = c
                            bj = jn[prow + dd - 1]
                if two and dd + 1 < nd and q0 < nres and q1 < nres:
                    pk = key[prow2 + dd + 1]
                    r = int(ex[off21 + (q0 * nres + q1) * nres + t1])
                    if pk < big and r >= 0:
                        c = pk + rank_weight(r, b)
                        if c < limit and c < best:
                            best = c
                            bj = jn[prow2 + dd + 1]
        key[row + dd] = best
        jn[row + dd] = bj


@njit(cache=True, nogil=True)
def _verify_anchor(
    q: np.ndarray,
    p0: int,
    g: int,
    rem: int,
    pmin: int,
    text: np.ndarray,
    fstart: np.ndarray,
    flen: np.ndarray,
    gfr: np.ndarray,
    ex: np.ndarray,
    k1: np.ndarray,
    j1: np.ndarray,
    k2: np.ndarray,
    j2: np.ndarray,
    out: np.ndarray,
    n: int,
) -> int:  # pragma: no cover - compiled by numba; exercised through match_frameshift
    """Append the best split of every locus whose prefix starts at strand nt ``p0`` to ``out``."""
    big = 1 << 62
    L = q.shape[0]
    nd = 2 * rem + 1
    limit = (rem + 1) * (1 << 45)
    pmax = L - pmin
    of = p0 % 3
    s = p0 // 3
    f = gfr[g, of]
    fa = int(fstart[f])
    fn = int(flen[f])
    if s >= fn:
        return n
    base1 = s - rem
    for dd in range(nd):
        k1[dd] = 0 if dd == rem else big
        j1[dd] = -1
    for i in range(1, pmax + 1):
        band_row(q, i, base1, text, fa, fn, ex, k1, j1, nd, limit, 0)
    for delta in (1, -1):
        og = (of + delta) % 3
        cd = (of + delta - og) // 3
        gf = gfr[g, og]
        ga = int(fstart[gf])
        gn = int(flen[gf])
        base2 = base1 + cd
        for i in range(pmin, L + 1):
            row = i * nd
            if i > pmin:
                band_row(q, i, base2, text, ga, gn, ex, k2, j2, nd, limit, pmin)
            else:
                for dd in range(nd):
                    k2[row + dd] = big
                    j2[row + dd] = -1
            if i <= pmax:
                for dd in range(nd):
                    c = k1[row + dd]
                    jj = base2 + i + dd
                    if c < k2[row + dd] and 0 <= jj <= gn:
                        k2[row + dd] = c
                        j2[row + dd] = row + dd
        row = L * nd
        for dd in range(nd):
            c = k2[row + dd]
            if c < limit and n < out.shape[0]:
                jun = j2[row + dd]
                jp = base1 + jun // nd + jun % nd  # prefix end column in F
                out[n, 0] = c
                out[n, 1] = p0
                out[n, 2] = of + 3 * jp
                out[n, 3] = og + 3 * (jp + cd)
                out[n, 4] = og + 3 * (base2 + L + dd)
                out[n, 5] = g
                out[n, 6] = delta
                n += 1
    return n


# --- precedence: whole-query alignments without a frameshift -----------------------------------


@njit(cache=True, nogil=True)
def nonfs_overlap(
    q: np.ndarray, text: np.ndarray, fa: int, fn: int, o: int, a: int, b: int, c: int, ex: Arr
) -> bool:  # pragma: no cover - compiled by numba; exercised through match_frameshift
    """Does ``q`` align (cost ``<= c``, no frameshift) to a window of the frame overlapping
    strand nt ``[a, b)``? Start columns are restricted to codons starting before ``b``."""
    nres = 19
    big = 1 << 30
    off21 = nres * nres * nres
    off22 = 2 * off21
    L = q.shape[0]
    smax = min((b - o - 1) // 3, fn - 1)
    lo = max(0, (a - o) // 3 - L - 2 * c - 2)
    if smax < lo:
        return False
    hi = min(fn, smax + L + 2 * c + 2)
    emin = (a - o) // 3 + 1
    m = hi - lo + 1
    d = np.full((3, m), big, dtype=np.int64)
    for x in range(lo, smax + 1):
        d[0, x - lo] = 0
    for i in range(1, L + 1):
        cur = d[i % 3]
        p1 = d[(i - 1) % 3]
        p2 = d[(i - 2) % 3]
        q1 = int(q[i - 1])
        q0 = int(q[i - 2]) if i >= 2 else nres
        cur[0] = big
        for x in range(1, m):
            j = lo + x
            best = big
            t1 = int(text[fa + j - 1])
            if t1 == 22:  # READTHROUGH: a substitution, the only move over a consumable stop
                best = p1[x - 1] + 1
            elif t1 < nres:
                v = p1[x - 1] + (0 if q1 == t1 else 1)
                best = min(best, v)
                if x >= 2:
                    t0 = int(text[fa + j - 2])
                    if t0 < nres and q1 < nres:
                        if int(ex[(q1 * nres + t0) * nres + t1]) >= 0:
                            best = min(best, p1[x - 2] + 1)
                        if i >= 2 and q0 < nres:
                            sw = q0 != q1 and q0 == t1 and q1 == t0
                            r = int(ex[off22 + ((q0 * nres + q1) * nres + t0) * nres + t1])
                            if sw or r >= 0:
                                best = min(best, p2[x - 2] + 1)
                ok21 = i >= 2 and q0 < nres and q1 < nres
                if ok21 and int(ex[off21 + (q0 * nres + q1) * nres + t1]) >= 0:
                    best = min(best, p2[x - 1] + 1)
            cur[x] = best if best <= c else big
    return bool((d[L % 3][max(emin - lo, 0) :] <= c).any())


@njit(cache=True, nogil=True)
def _suppressed(
    q: Arr, pl: Arr, cost: int, sp: Space, ex: Arr
) -> bool:  # pragma: no cover - compiled by numba; exercised through match_frameshift
    """Precedence: a non-frameshift alignment (cost <= ``cost``) overlaps the split locus."""
    text, fstart, flen, _, _, gfr, _, g_len, g_circ = sp
    g = pl[5]
    a = pl[1]
    b = pl[4]
    ln = g_len[g]
    for sh in (-1, 0, 1):
        if sh != 0 and not g_circ[g]:
            continue
        for o in range(3):
            f = gfr[g, o]
            fa = int(fstart[f])
            fn = int(flen[f])
            if nonfs_overlap(q, text, fa, fn, o, a + sh * ln, b + sh * ln, cost, ex):
                return True
    return False


# --- seeds -------------------------------------------------------------------------------------


@njit(cache=True, nogil=True)
def _anchor_ranges(
    qc: Arr, q: int, offsets: Arr, positions: Arr, t: int, sp: Space, gspan: int, pmh: int
) -> Arr:  # pragma: no cover - compiled by numba; exercised through match_frameshift
    """Merged anchor ranges ``(group, P0_lo, P0_hi)`` (inclusive); prefilter above ``pmh`` hits."""
    _, fstart, _, f_off, f_grp, _, _, _, _ = sp
    L = qc.shape[0]
    nq = L - q + 1
    codes = np.full(max(nq, 0), -1, dtype=np.int64)
    total = 0
    valid = 0
    for i in range(nq):
        c = 0
        ok = True
        for j in range(q):
            r = int(qc[i + j])
            if r >= 19:
                ok = False
                break
            c = c * 19 + r
        if ok:
            codes[i] = c
            m = offsets[c + 1] - offsets[c]
            if m > 0:
                total += m
                valid += 1
    if valid < t:
        return np.empty((0, 3), dtype=np.int64)
    keys = np.empty(total, dtype=np.int64)
    m = 0
    nfr = fstart.shape[0] - 1
    for i in range(nq):
        c = codes[i]
        if c < 0:
            continue
        for x in range(offsets[c], offsets[c + 1]):
            p = int(positions[x])
            fid = np.searchsorted(fstart, p, side="right") - 1
            if fid >= nfr:
                fid = nfr - 1
            g = f_off[fid] + 3 * (p - fstart[fid]) - 3 * i + 3 * L + 8
            keys[m] = ((f_grp[fid] * gspan + g) << 10) | i
            m += 1
    if t >= 2 and total > pmh:  # drop hits whose two 8-diagonal bins hold < t hits (hashed)
        cnt = np.zeros(1 << 18, dtype=np.uint8)
        for x in range(total):
            h = ((keys[x] >> 13) * 0x9E3779B1) & 0x3FFFF
            cnt[h] = min(cnt[h] + 1, 255)
        m = 0
        for x in range(total):
            bn = keys[x] >> 13
            c0 = int(cnt[((bn - 1) * 0x9E3779B1) & 0x3FFFF])
            c1 = int(cnt[(bn * 0x9E3779B1) & 0x3FFFF])
            c2 = int(cnt[((bn + 1) * 0x9E3779B1) & 0x3FFFF])
            if c0 + c1 >= t or c1 + c2 >= t:
                keys[m] = keys[x]
                m += 1
        total = m
        keys = keys[:m]
    keys.sort()
    ud = np.empty(total, dtype=np.int64)
    us = np.empty(total + 1, dtype=np.int64)
    nu = 0
    for x in range(total):
        d = keys[x] >> 10
        if nu == 0 or d != ud[nu - 1]:
            ud[nu] = d
            us[nu] = x
            nu += 1
    us[nu] = total
    seen = np.full(max(nq, 1), -1, dtype=np.int64)
    out = np.empty((nu, 3), dtype=np.int64)
    no = 0
    hi = 0
    for u in range(nu):
        while hi < nu and ud[hi] <= ud[u] + 4:
            hi += 1
        if us[hi] - us[u] < t:
            continue
        nd = 0
        for x in range(us[u], us[hi]):
            i = keys[x] & 1023
            if seen[i] != u:
                seen[i] = u
                nd += 1
        if nd < t:
            continue
        grp = ud[u] // gspan
        p_lo = ud[u] % gspan - 3 * L - 8 - 4
        p_hi = p_lo + 8
        if no > 0 and out[no - 1, 0] == grp and p_lo <= out[no - 1, 2] + 1:
            out[no - 1, 2] = max(out[no - 1, 2], p_hi)
        else:
            out[no, 0] = grp
            out[no, 1] = p_lo
            out[no, 2] = p_hi
            no += 1
    return out[:no].copy()


# --- one query ---------------------------------------------------------------------------------


@njit(cache=True, nogil=True)
def _search_one(
    qc: Arr, k: int, t: int, q: int, pmin: int, offsets: Arr, positions: Arr, sp: Space,
    gspan: int, ex: Arr, pmh: int,
) -> Arr:  # pragma: no cover - compiled by numba; exercised through match_frameshift  # fmt: skip
    """Placements (7 columns) of one query at the minimal unsuppressed split cost."""
    text, fstart, flen, _, _, gfr, g_ext, g_len, _ = sp
    L = qc.shape[0]
    rem = k - 1
    if t >= 1:
        rng = _anchor_ranges(qc, q, offsets, positions, t, sp, gspan, pmh)
    else:
        ng = gfr.shape[0]
        rng = np.empty((ng, 3), dtype=np.int64)
        for g in range(ng):
            rng[g, 0] = g
            rng[g, 1] = 0
            rng[g, 2] = g_ext[g] - 1
    nd = 2 * rem + 1
    size = (L + 1) * nd
    k1 = np.empty(size, dtype=np.int64)
    j1 = np.empty(size, dtype=np.int64)
    k2 = np.empty(size, dtype=np.int64)
    j2 = np.empty(size, dtype=np.int64)
    cap = 64
    pls = np.empty((cap, 7), dtype=np.int64)
    n = 0
    for r in range(rng.shape[0]):
        g = rng[r, 0]
        for p0 in range(max(rng[r, 1], 0), min(rng[r, 2], g_ext[g] - 1) + 1):
            if pls.shape[0] - n < 2 * nd:
                grown = np.empty((2 * pls.shape[0], 7), dtype=np.int64)
                grown[:n] = pls[:n]
                pls = grown
            n = _verify_anchor(
                qc, p0, g, rem, pmin, text, fstart, flen, gfr, ex, k1, j1, k2, j2, pls, n
            )
    keep = np.zeros(n, dtype=np.bool_)
    for cost in range(1, k + 1):
        found = False
        for x in range(n):
            if pls[x, 0] // (1 << 45) + 1 != cost or pls[x, 4] - pls[x, 1] > g_len[pls[x, 5]]:
                continue
            if not _suppressed(qc, pls[x], cost, sp, ex):
                keep[x] = True
                found = True
        if found:
            break
    m = 0
    for x in range(n):
        if keep[x]:
            pls[m] = pls[x]
            m += 1
    return pls[:m].copy()


@njit(cache=True, parallel=True)
def _search_many(
    qcat: Arr, qoff: Arr, ks: Arr, ts: Arr, q: int, pmin: int, offsets: Arr, positions: Arr,
    sp: Space, gspan: int, ex: Arr, pmh: int, cap: int,
) -> tuple[Arr, Arr]:  # pragma: no cover - compiled by numba; via match_frameshift  # fmt: skip
    """:func:`_search_one` for every query (placement counts and the first ``cap`` rows)."""
    nq = qoff.shape[0] - 1
    counts = np.zeros(nq, dtype=np.int64)
    rows = np.full((nq, cap, 7), -1, dtype=np.int64)
    for j in prange(nq):
        qc = qcat[qoff[j] : qoff[j + 1]]
        res = _search_one(qc, ks[j], ts[j], q, pmin, offsets, positions, sp, gspan, ex, pmh)
        counts[j] = res.shape[0]
        for r in range(min(res.shape[0], cap)):
            rows[j, r] = res[r]
    return counts, rows


# --- Python driver -----------------------------------------------------------------------------


def _threshold(L: int, k: int, q: int, pmin: int) -> int:
    """Lossless q-gram count of a split (minimum over the allowed split points)."""
    rem = k - 1
    best = 1 << 30
    for p in range(pmin, L - pmin + 1):
        a, b = max(0, p - q + 1), max(0, L - p - q + 1)
        g = a + b if rem == 0 else min(max(0, a - q - 1) + b, a + max(0, b - q - 1))
        best = min(best, g)
    return best


@cache
def _choose_q(L: int, k: int, pmin: int, n_text: int) -> tuple[int, int]:
    """``(q, t)`` minimising estimated hits + verified anchor ranges; ``t = 0``: no filter."""
    best: tuple[float, int, int] | None = None
    for q in _QS:
        t = _threshold(L, k, q, pmin)
        if t < 1:
            continue
        nq = L - q + 1
        hits = nq * n_text / 19.0**q
        cands = 3.0 * n_text * comb(nq, t) * (5.0 / 19.0**q) ** t
        est = hits + _VERIFY_COST * cands
        if best is None or est < best[0]:
            best = (est, q, t)
    return (3, 0) if best is None else (best[1], best[2])


class _FSearch:
    """One search space: frames (read in the search direction), q-gram indexes, strand groups."""

    def __init__(
        self,
        frames: tuple[Frame, ...],
        replicons: Mapping[str, Replicon],
        marked: Sequence[np.ndarray] | None = None,
    ) -> None:
        self.base = Search(frames, marked)
        self.frames = frames
        groups: dict[tuple[str, str], int] = {}
        gfr: list[list[int]] = []
        ext: list[int] = []
        length: list[int] = []
        circ: list[bool] = []
        f_grp = np.zeros(len(frames), dtype=np.int64)
        for fid, f in enumerate(frames):
            key = (f.contig, f.strand)
            if key not in groups:
                groups[key] = len(gfr)
                gfr.append([-1, -1, -1])
                rep = replicons[f.contig]
                ext.append(f.ext_nt)
                length.append(rep.length)
                circ.append(rep.topology == "circular")
            gfr[groups[key]][f.offset] = fid
            f_grp[fid] = groups[key]
        if any(-1 in row for row in gfr):
            raise ValueError("every contig strand needs frames at offsets 0, 1 and 2")
        self.f_grp = f_grp
        self.f_off = np.array([f.offset for f in frames], dtype=np.int64)
        self.flen = np.array([len(f.aa) for f in frames], dtype=np.int64)
        self.gfr = np.array(gfr, dtype=np.int64).reshape(-1, 3)
        self.g_ext = np.array(ext, dtype=np.int64)
        self.g_len = np.array(length, dtype=np.int64)
        self.g_circ = np.array(circ, dtype=np.bool_)
        self.space: Space = (
            self.base.text, self.base.starts, self.flen, self.f_off, self.f_grp, self.gfr,
            self.g_ext, self.g_len, self.g_circ,
        )  # fmt: skip
        self.gspan = (max(ext) if ext else 0) + 3 * MAX_QUERY + 64
        if len(gfr) * self.gspan >= 1 << 52:
            raise ValueError("genome too large for the frameshift seed keys")

    def run(
        self, codes: Mapping[str, Arr], budgets: Mapping[str, int], pmin: int, ex: Arr, pmh: int
    ) -> dict[str, np.ndarray]:
        """Placements of every query (keyed by its text); ``pmh``: prefilter hit threshold."""
        b = self.base
        n_text = len(b.text)
        groups: dict[int, list[tuple[str, int, int]]] = {}
        for text in sorted(codes):
            k = budgets[text]
            q, t = _choose_q(len(codes[text]), k, pmin, n_text)
            groups.setdefault(q, []).append((text, k, t))
        out: dict[str, np.ndarray] = {}
        for q, items in sorted(groups.items()):
            idx = b.index(q)
            cs = [np.ascontiguousarray(codes[text], dtype=np.uint8) for text, _, _ in items]
            qoff = np.zeros(len(cs) + 1, dtype=np.int64)
            np.cumsum([len(c) for c in cs], out=qoff[1:])
            qcat = np.ascontiguousarray(np.concatenate(cs), dtype=np.uint8)
            ks = np.array([k for _, k, _ in items], dtype=np.int64)
            ts = np.array([t for _, _, t in items], dtype=np.int64)
            args = (idx.offsets, idx.positions, self.space, self.gspan, ex, pmh)
            counts, rows = _search_many(qcat, qoff, ks, ts, q, pmin, *args, _CAP)
            for j, (text, k, t) in enumerate(items):
                m = int(counts[j])
                if m > _CAP:
                    out[text] = _search_one(cs[j], k, t, q, pmin, *args)
                else:
                    out[text] = rows[j, :m]
        return out


def _like(f: Frame, aa: np.ndarray) -> Frame:
    return Frame(f.contig, f.strand, f.offset, np.ascontiguousarray(aa, dtype=np.uint8), f.ext_nt)


def _strand_to_genome(frame: Frame, a: int, b: int) -> tuple[int, int]:
    return (a, b) if frame.strand == "+" else (frame.ext_nt - b, frame.ext_nt - a)


def _rows_of(
    pep: str,
    text: str,
    pls: np.ndarray,
    search: _FSearch,
    index: SixFrameIndex,
    kind: DecoyType | None,
) -> list[dict[str, object]]:
    """One row per genome interval (lowest key first, then the first placement)."""
    best: dict[tuple[str, int, int, str], tuple[int, dict[str, object]]] = {}
    for pl in sorted(tuple(int(v) for v in r) for r in pls):
        key, p0, pe, ss, se, g, _ = pl
        blocks = [(p0, pe), (ss, se)]
        fprefix = int(search.gfr[g, p0 % 3])
        pre_block = blocks[0]
        if kind == "reverse":  # prefix of Q = suffix of the reversed query
            fprefix = int(search.gfr[g, ss % 3])
            pre_block = blocks[1]
        frame = search.frames[fprefix]
        rep = index.replicons[frame.contig]
        gb = [_strand_to_genome(frame, a, b) for a, b in blocks]
        gp = _strand_to_genome(frame, *pre_block)
        g0 = min(a for a, _ in gb)
        shift = (g0 % rep.length) - g0 if rep.topology == "circular" else 0
        gb = sorted((a + shift, b + shift) for a, b in gb)
        start, end = gb[0][0], gb[1][1]
        pstart = (gp[0] + shift) % rep.length
        piv = Interval(frame.contig, pstart, pstart + gp[1] - gp[0], frame.strand)
        label = index.frame_label(piv)
        a_aa, b_aa = (pre_block[0] - frame.offset) // 3, (pre_block[1] - frame.offset) // 3
        if kind == "reverse":
            n_aa = len(frame.aa)
            a_aa, b_aa = n_aa - b_aa, n_aa - a_aa
        cost, types = decode_key(key)
        iv = (frame.contig, start, end, frame.strand)
        if iv in best and best[iv][0] <= key:
            continue
        best[iv] = (
            key,
            {
                "peptide": pep,
                "query": text if kind != "reverse" else text[::-1],
                "search": "target" if kind is None else "decoy",
                "decoy_type": kind,
                "contig": frame.contig,
                "frame_id": fprefix,
                "aa_start": a_aa,
                "aa_end": b_aa,
                "start": start,
                "end": end,
                "strand": frame.strand,
                "frame": label,
                "route": "tolerant",
                "n_diff": cost + 1,
                "diff_types": ";".join(sorted([*types, "frameshift"])),
                "blocks": ";".join(f"{a}-{b}" for a, b in gb),
            },
        )
    return [row for _, row in best.values()]


def match_frameshift(
    queries: Sequence[str],
    index: SixFrameIndex,
    decoys: Mapping[DecoyType, tuple[Frame, ...]],
    decoy_of: Mapping[str, DecoyType],
    target_queries: set[str],
    cfg: TolerantConfig,
    *,
    prefilter_min_hits: int = 256,
    cells: StopCells | None = None,
) -> pd.DataFrame:
    """Frameshift-split target and decoy hits (``HIT_COLUMNS``, ``route = "tolerant"``).

    Peptides are searched as by :func:`~gippa.mapping.tolerant.match_tolerant` (same length
    classes and budgets); a split needs ``k >= 1`` and ``min_side_aa`` residues on each side.
    Rows carry ``diff_types`` with ``frameshift`` and ``blocks = "s1-e1;s2-e2"`` (genome,
    ascending, unwrapped like ``start``/``end``); ``frame`` labels the prefix block. Rows are
    sorted by peptide, search (target first), decoy type, contig order, start, end, strand.
    ``prefilter_min_hits``: the hashed seed-bin prefilter runs for a query with more than this
    many q-gram hits (it never changes the output, only the work done). ``cells``: consumable
    annotated stops (ADR-0058 item 3): split sides never consume them, but the precedence check
    sees read-through alignments; each search marks its own cells (target, or the mirrored cells
    of its decoy type in the search direction).

    Raises:
        ValueError: as :func:`~gippa.mapping.tolerant.match_tolerant` (a ``decoy_of`` type without
            decoy frames, a searched peptide without a decoy type, an invalid ``max_diff``).
    """
    for pep, kind in sorted(decoy_of.items()):
        if kind not in decoys:
            raise ValueError(f"decoy type {kind!r} of peptide {pep!r} has no decoy frames")
    ranges = parse_max_diff(cfg.max_diff)
    if not cfg.frameshift.enabled:
        return hits_frame()
    pmin = int(cfg.frameshift.min_side_aa)
    norm: dict[str, str] = {}
    budgets: dict[str, int] = {}
    for pep in sorted(set(queries)):
        k = budget(len(pep), cfg, ranges)
        if k is None:
            continue
        if pep not in decoy_of:
            raise ValueError(f"peptide {pep!r} has no decoy type in decoy_of")
        if k < 1 or 2 * pmin > len(pep):
            continue
        text = decode(encode_query(pep))
        norm[pep] = text
        budgets[text] = k
        budgets[text[::-1]] = k
    if not norm:
        return hits_frame()

    ex = exchange_lookup(block_tables())
    # search spaces: (frames read in the search direction, kind, {peptide: searched text}, marks)
    Plan = tuple[tuple[Frame, ...], DecoyType | None, dict[str, str], tuple[Arr, ...] | None]
    plans: list[Plan] = []
    tmarks = None if cells is None else cells.target
    targets = {p: t for p, t in norm.items() if p in target_queries}
    plans.append((index.frames, None, targets, tmarks))
    for kind in sorted(decoys):
        frames = decoys[kind]
        marks = None if cells is None else decoy_marks(cells, kind)
        if kind == "reverse":
            if marks is not None:  # into the re-reversed orientation
                marks = tuple(
                    np.sort(len(f.aa) - 1 - m) for f, m in zip(frames, marks, strict=True)
                )
            frames = tuple(_like(f, f.aa[::-1]) for f in frames)
            if all(np.array_equal(a.aa, b.aa) for a, b in zip(frames, index.frames, strict=True)):
                frames = index.frames  # reversed reverse decoys are the target frames
        rev = kind == "reverse"
        texts = {p: t[::-1] if rev else t for p, t in norm.items() if decoy_of[p] == kind}
        plans.append((frames, kind, texts, marks))

    def space_key(frames: tuple[Frame, ...], marks: tuple[Arr, ...] | None) -> tuple[object, ...]:
        return (id(frames), None if marks is None else tuple(m.tobytes() for m in marks))

    found: dict[tuple[object, ...], dict[str, np.ndarray]] = {}
    spaces: dict[tuple[object, ...], _FSearch] = {}
    for frames, _, _, marks in plans:
        sk = space_key(frames, marks)
        if sk in found:
            continue
        wanted = {
            t for fr, _, texts, mk in plans if space_key(fr, mk) == sk for t in texts.values()
        }
        if wanted:
            spaces[sk] = _FSearch(frames, index.replicons, marks)
            codes = {t: encode_query(t) for t in sorted(wanted)}
            found[sk] = spaces[sk].run(codes, budgets, pmin, ex, prefilter_min_hits)
    rows: list[dict[str, object]] = []
    for frames, kd, texts, marks in plans:
        sk = space_key(frames, marks)
        for pep, t in texts.items():
            rows.extend(_rows_of(pep, t, found[sk][t], spaces[sk], index, kd))

    contig_order = {c: i for i, c in enumerate(index.replicons)}
    rows.sort(
        key=lambda r: (
            r["peptide"],
            _SEARCH_ORDER[str(r["search"])],
            str(r["decoy_type"] or ""),
            contig_order[str(r["contig"])],
            r["start"],
            r["end"],
            r["strand"],
        )
    )
    return hits_frame(rows)
