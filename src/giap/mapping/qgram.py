"""q-gram seeds for tolerant matching (GIAP_SPEC.md §7.2.4, G4; issue M3-07).

The frames of one search are concatenated into one text with ``GAP`` stop codes between them
(``frame_starts[f]`` is the global position of residue 0 of frame ``f``), so no q-gram crosses a
frame boundary and a ``±2`` band around any diagonal of a query of at most ``MAX_QUERY`` residues
touches one frame only. :class:`QGramIndex` is a CSR over the ``19**q`` q-gram codes; q-grams
holding ``*`` or ``X`` are skipped.

**Lossless threshold.** A difference destroys at most ``q + 1`` of the ``L - q + 1`` query
q-grams, so a window within ``k`` differences shares ``t = (L - q + 1) - k (q + 1)`` q-grams
with it, all on the diagonals of its alignment path (at most ``k + 1`` adjacent diagonals; three
for ``k <= 2``). A centre diagonal ``c`` is a *raw* candidate when the number of distinct query
q-gram offsets hit on ``c - 1, c, c + 1`` is ``>= t``. Raw candidates are merged greedily: the
smallest uncovered raw ``r`` gives the reported centre ``r + 1``, which covers the raw candidates
``r .. r + 2`` (its ``±2`` verification band contains each one's ``±1`` band). Reported centres
are therefore at least three diagonals apart and every alignment path is inside the band of one.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numba import njit

from giap.mapping.sixframe import STOP, Frame

__all__ = [
    "GAP",
    "MAX_QUERY",
    "QGramIndex",
    "build_qgram_index",
    "candidate_diagonals",
    "choose_q",
    "concat_frames",
]

GAP: int = 1024  # stop codes between frames in the concatenated text
MAX_QUERY: int = 512  # longest query: bands stay inside one frame, offsets fit 10 bits
_NRES: int = 19
_QS: tuple[int, ...] = (5, 4, 3)


def choose_q(L: int, k: int) -> tuple[int, int]:
    """Largest ``q`` in {5, 4, 3} with a lossless threshold ``t >= 1``, as ``(q, t)``.

    Raises:
        ValueError: when no ``q`` gives ``t >= 1``.
    """
    for q in _QS:
        t = (L - q + 1) - k * (q + 1)
        if t >= 1:
            return q, t
    raise ValueError(f"no lossless q-gram threshold for a query of {L} residues with k = {k}")


@dataclass(frozen=True)
class QGramIndex:
    """CSR over q-gram codes: positions of code ``c`` are ``positions[offsets[c]:offsets[c+1]]``.

    ``positions`` are global ``int32`` positions in the concatenated text (ascending within a
    code); ``frame_starts`` has one entry per frame plus the text length.
    """

    q: int
    offsets: np.ndarray
    positions: np.ndarray
    frame_starts: np.ndarray


def concat_frames(frames: Sequence[Frame]) -> tuple[np.ndarray, np.ndarray]:
    """``(text, frame_starts)``: frames joined with ``GAP`` stop codes after each one."""
    starts = np.zeros(len(frames) + 1, dtype=np.int64)
    np.cumsum([len(f.aa) + GAP for f in frames], out=starts[1:])
    total = int(starts[-1])
    if total >= 2**31:
        raise ValueError(f"concatenated frames of {total} residues exceed the int32 positions")
    text = np.full(total, STOP, dtype=np.uint8)
    for f, frame in enumerate(frames):
        a = int(starts[f])
        text[a : a + len(frame.aa)] = frame.aa
    return text, starts


@njit(cache=True, nogil=True)
def _build_kernel(
    text: np.ndarray, q: int
) -> tuple[np.ndarray, np.ndarray]:  # pragma: no cover - compiled by numba; via build_qgram_index
    """CSR ``(offsets, positions)`` of every q-gram of residue codes ``< 19`` in ``text``."""
    n = text.shape[0]
    nb = 1
    for _ in range(q):
        nb *= 19
    mod = nb // 19
    m = n - q + 1 if n >= q else 0
    codes = np.full(m, -1, dtype=np.int32)
    run = 0
    code = 0
    for p in range(n):
        c = int(text[p])
        if c < 19:
            code = (code % mod) * 19 + c
            run += 1
            if run >= q:
                codes[p - q + 1] = code
        else:
            run = 0
            code = 0
    offsets = np.zeros(nb + 1, dtype=np.int64)
    for i in range(m):
        if codes[i] >= 0:
            offsets[codes[i] + 1] += 1
    for c in range(nb):
        offsets[c + 1] += offsets[c]
    positions = np.empty(offsets[nb], dtype=np.int32)
    fill = offsets[:nb].copy()
    for i in range(m):
        c = codes[i]
        if c >= 0:
            positions[fill[c]] = i
            fill[c] += 1
    return offsets, positions


def _index_from_text(text: np.ndarray, starts: np.ndarray, q: int) -> QGramIndex:
    if q not in _QS:
        raise ValueError(f"q must be one of {_QS}, got {q}")
    offsets, positions = _build_kernel(np.ascontiguousarray(text, dtype=np.uint8), q)
    return QGramIndex(q=q, offsets=offsets, positions=positions, frame_starts=starts)


def build_qgram_index(frames: Sequence[Frame], q: int) -> QGramIndex:
    """q-gram index of ``frames`` (q-grams with ``*`` or ``X`` skipped)."""
    text, starts = concat_frames(frames)
    return _index_from_text(text, starts, q)


@njit(cache=True, nogil=True)
def _candidates_kernel(
    qc: np.ndarray, q: int, offsets: np.ndarray, positions: np.ndarray, t: int
) -> np.ndarray:  # pragma: no cover - compiled by numba; exercised through candidate_diagonals
    """Merged centre diagonals (global, ascending) of query codes ``qc`` (see module docstring)."""
    L = qc.shape[0]
    nq = L - q + 1
    if nq < t:
        return np.empty(0, dtype=np.int64)
    codes = np.full(nq, -1, dtype=np.int64)
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
            n = offsets[c + 1] - offsets[c]
            if n > 0:
                total += n
                valid += 1
    if valid < t:
        return np.empty(0, dtype=np.int64)
    # one key per hit: shifted diagonal (p - i + L >= 1) in the high bits, offset i in the low 10
    keys = np.empty(total, dtype=np.int64)
    m = 0
    for i in range(nq):
        c = codes[i]
        if c >= 0:
            for x in range(offsets[c], offsets[c + 1]):
                keys[m] = (int(positions[x]) - i + L) * 1024 + i
                m += 1
    keys.sort()
    ud = np.empty(total, dtype=np.int64)  # distinct shifted diagonals
    us = np.empty(total + 1, dtype=np.int64)  # first hit of each
    nu = 0
    for x in range(total):
        d = keys[x] >> 10
        if nu == 0 or d != ud[nu - 1]:
            ud[nu] = d
            us[nu] = x
            nu += 1
    us[nu] = total
    seen = np.full(nq, -2, dtype=np.int64)
    out = np.empty(nu, dtype=np.int64)
    no = 0
    last = -(1 << 62)
    cover = -(1 << 62)
    lo = 0
    hi = 0
    for u in range(nu):
        for c in (ud[u] - 1, ud[u], ud[u] + 1):
            if c <= last:
                continue
            last = c
            if c <= cover:  # inside the band of a reported centre
                continue
            while ud[lo] < c - 1:
                lo += 1
            while hi < nu and ud[hi] <= c + 1:
                hi += 1
            if us[hi] - us[lo] < t:
                continue
            cnt = 0
            for x in range(us[lo], us[hi]):
                i = keys[x] & 1023
                if seen[i] != c:
                    seen[i] = c
                    cnt += 1
            if cnt >= t:
                out[no] = c + 1 - L
                no += 1
                cover = c + 2
    return out[:no].copy()


def candidate_diagonals(query: np.ndarray, idx: QGramIndex, t: int) -> np.ndarray:
    """Merged candidate centre diagonals (global ``p - i``, ``int64``, ascending, >= 3 apart).

    ``query`` holds M3-02 query codes of any integer dtype; q-grams with a ``NOMATCH`` residue
    are skipped (``t`` still counts them, so such queries may find nothing).

    Raises:
        ValueError: if ``t < 1`` or the query is longer than ``MAX_QUERY``.
    """
    if t < 1:
        raise ValueError(f"threshold t must be >= 1, got {t}")
    qc = np.ascontiguousarray(query, dtype=np.uint8)
    if qc.shape[0] > MAX_QUERY:
        raise ValueError(f"query of {qc.shape[0]} residues exceeds MAX_QUERY = {MAX_QUERY}")
    out: np.ndarray = _candidates_kernel(qc, idx.q, idx.offsets, idx.positions, int(t))
    return out
