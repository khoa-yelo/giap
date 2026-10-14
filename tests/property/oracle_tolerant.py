"""Independent brute-force oracle for the tolerant verifier (issue M3-06; GIAP_SPEC.md §7.2.4).

Written by the test author from the spec text only (it does not import ``giap.mapping.diffops``):
the operation set, the block lists and the tie-break order are copied from §7.2.4 and the M3-06
semantics table. It works on *letters* (``decode`` output: residues with I folded into L, ``*``
for a stop, ``X`` for an ambiguous codon, ``#`` for a query NOMATCH).

The query is aligned globally to a target window ``frame[s:e]`` (free ends) by exhaustive
recursion over operation sequences (memoised on the cell). ``oracle_verify`` enumerates the
windows of length ``L - 2k .. L + 2k``; with ``band=True`` every visited cell ``(i, j)`` must
satisfy ``|(j - i) - d0| <= 2`` (the DP's band), with ``band=False`` there is no band at all and
every window of the frame is tried.
"""

from __future__ import annotations

from functools import cache

PRECEDENCE: tuple[str, ...] = ("substitution", "swap", "block", "near_block")
RANK = {name: i for i, name in enumerate(PRECEDENCE)}

# §7.2.4, as compositions (sorted letters); "both orders of each pair".
_BLOCK_PAIRS = (("N", "GG"), ("Q", "AG"), ("GT", "AS"), ("AV", "GL"), ("AN", "GQ"))
_NEAR_PAIRS = (("W", "AD"), ("W", "EG"), ("R", "GV"))
BAND = 2
UNUSABLE_TARGET = frozenset("*X")


def _comp(s: str) -> str:
    return "".join(sorted(s))


def _pairs(raw: tuple[tuple[str, str], ...]) -> frozenset[frozenset[str]]:
    return frozenset(frozenset({_comp(a), _comp(b)}) for a, b in raw)


BLOCK_SET = _pairs(_BLOCK_PAIRS)
NEAR_SET = _pairs(_NEAR_PAIRS)

Best = tuple[int, tuple[int, ...]]  # (cost, sorted diff-type ranks)


def _is_exchange(qseg: str, tseg: str, pairs: frozenset[frozenset[str]]) -> bool:
    if set(tseg) & UNUSABLE_TARGET or "#" in qseg:
        return False
    key = frozenset({_comp(qseg), _comp(tseg)})
    return len(key) == 2 and key in pairs


def _ops(q: str, t: str, i: int, j: int, e: int) -> list[tuple[int, int, int | None]]:
    """All operations applicable at cell (i, j): (a, b, diff rank or None for a match)."""
    out: list[tuple[int, int, int | None]] = []
    L = len(q)
    if i < L and j < e and t[j] not in UNUSABLE_TARGET:
        if q[i] == t[j] and q[i] != "#":
            out.append((1, 1, None))
        else:
            out.append((1, 1, RANK["substitution"]))
    if i + 2 <= L and j + 2 <= e:
        a, b = q[i], q[i + 1]
        x, y = t[j], t[j + 1]
        if a != b and a == y and b == x and not ({x, y} & UNUSABLE_TARGET) and "#" not in (a, b):
            out.append((2, 2, RANK["swap"]))
    for qa, tb in ((1, 2), (2, 1), (2, 2)):
        if i + qa <= L and j + tb <= e:
            qseg, tseg = q[i : i + qa], t[j : j + tb]
            if _is_exchange(qseg, tseg, BLOCK_SET):
                out.append((qa, tb, RANK["block"]))
            if (qa, tb) != (2, 2) and _is_exchange(qseg, tseg, NEAR_SET):
                out.append((qa, tb, RANK["near_block"]))
    return out


def align_window(q: str, t: str, s: int, e: int, d0: int | None, kmax: int) -> Best | None:
    """Minimal (cost, sorted ranks) aligning all of ``q`` to exactly ``t[s:e]``, cost <= kmax.

    ``d0`` not None enforces the band on every visited cell (frame coordinates).
    """
    L = len(q)

    def in_band(i: int, j: int) -> bool:
        return d0 is None or abs((j - i) - d0) <= BAND

    @cache
    def rec(i: int, j: int) -> Best | None:
        if not in_band(i, j):
            return None
        if i == L and j == e:
            return (0, ())
        best: Best | None = None
        for a, b, rank in _ops(q, t, i, j, e):
            sub = rec(i + a, j + b)
            if sub is None:
                continue
            cost = sub[0] + (0 if rank is None else 1)
            if cost > kmax:
                continue
            types = sub[1] if rank is None else tuple(sorted((*sub[1], rank)))
            cand = (cost, types)
            if best is None or cand < best:
                best = cand
        return best

    return rec(0, s)


def oracle_verify(
    q: str, frame: str, d0: int, k: int, *, band: bool = True
) -> tuple[int, tuple[str, ...], int] | None:
    """Best (cost, diff_types, t_start) over windows, or None when the minimum exceeds ``k``.

    Order: lowest cost, then the lexicographically smallest diff-type multiset, then the smallest
    ``t_start``. Window lengths are ``L - 2k .. L + 2k``; with ``band`` the window ends must lie on
    band cells as well.
    """
    L, n = len(q), len(frame)
    found: list[tuple[int, tuple[int, ...], int]] = []
    starts = range(max(0, d0 - BAND), min(n, d0 + BAND) + 1) if band else range(n + 1)
    for s in starts:
        for length in range(max(1, L - 2 * k), L + 2 * k + 1):
            e = s + length
            if e > n:
                break
            if band and abs((e - L) - d0) > BAND:
                continue
            got = align_window(q, frame, s, e, d0 if band else None, k)
            if got is not None:
                found.append((got[0], got[1], s))
    if not found:
        return None
    cost, ranks, s = min(found)
    return cost, tuple(PRECEDENCE[r] for r in ranks), s


def window_cost(
    q: str, frame: str, s: int, e: int, d0: int, k: int
) -> tuple[int, tuple[str, ...]] | None:
    """Banded best (cost, diff_types) of one given window (for checking a reported window)."""
    got = align_window(q, frame, s, e, d0, k)
    return None if got is None else (got[0], tuple(PRECEDENCE[r] for r in got[1]))
