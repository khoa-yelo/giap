"""Residue-level alignment of every accepted placement (GIPPA_SPEC.md §7.5.2 step 1; issue M7-01).

The merged mapper keeps no traceback (M3-06 encodes the optimum in one additive key), so Module 5
recomputes it here: exact rows are a zip (I/L only), tolerant rows a global banded DP (``±2``
around diagonal 0) with traceback over the M3-06 operations and alphabet (I folded into L for the
cost only), and frameshift rows (M3-08) align the prefix on the 5' block and the suffix on the
other block at every feasible peptide junction. The recomputed cost must equal ``n_diff``.

Co-optimal alignments are ordered by the M3-06 key (cost, then the diff-type multiset under the
``DIFF_TYPES`` precedence), then by the lexicographically smallest tuple of the differences'
peptide start positions (M7-01 lead decision 1). Coordinates are 0-based half-open and unwrapped
(§6.5); ``codons`` hold the lowest coordinate of each consumed codon, in peptide order. An exact
row whose residue 0 is ``M`` over an annotated non-ATG start codon (M11-20 rule 6) has the op
``initiator`` there.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache
from typing import Literal

import numpy as np
import pandas as pd
from numba import njit

from gippa.mapping.diffops import block_tables
from gippa.mapping.sixframe import encode_query, encode_target, translate_text
from gippa.mapping.tolerant_dp import BAND, decode_key, exchange_lookup, rank_weight
from gippa.model import Replicon

__all__ = [
    "FrameshiftEvent",
    "HitAlignment",
    "OpKind",
    "ResidueOp",
    "align_hit",
    "align_hits",
    "translate_locus",
]

OpKind = Literal["match", "il", "substitution", "swap", "block", "near_block", "initiator"]
_DIFF_KINDS: tuple[OpKind, ...] = ("substitution", "swap", "block", "near_block")
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_BIG: int = 1 << 62
_B3: int = (1 << 15) ** 3  # cost unit of the M3-06 key


@dataclass(frozen=True, slots=True)
class ResidueOp:
    """One operation: peptide residues ``pep`` against the genome codons ``codons``."""

    kind: OpKind
    pep: tuple[int, int]
    codons: tuple[int, ...]
    ref: str
    alt: str


@dataclass(frozen=True, slots=True)
class FrameshiftEvent:
    """The ±1-nt junction of a frameshift split (M3-08)."""

    pep_junction: int
    last_prefix_codon: int
    first_suffix_codon: int
    shift: Literal[-1, 1]


@dataclass(frozen=True, slots=True)
class HitAlignment:
    """Residue correspondence of one ``peptide_hits.tsv`` row."""

    peptide: str
    contig: str
    strand: Literal["+", "-"]
    frame: int
    ops: tuple[ResidueOp, ...]
    frameshift: FrameshiftEvent | None


def translate_locus(
    seq: str, rep: Replicon, start: int, end: int, strand: str, *, table: int = 11
) -> str:
    """Genetic-code ``table`` residues of ``[start, end)`` on ``strand`` (no I -> L; ``X`` for
    ambiguous).

    Reads across the origin of a circular replicon (``end > length``). Raises ValueError on an
    empty or over-long interval, or past the end of a linear contig.
    """
    n = rep.length
    if not 0 <= start < end or end - start > n:
        raise ValueError(f"{rep.id}: invalid locus [{start}, {end})")
    if end > n and rep.topology != "circular":
        raise ValueError(f"{rep.id}: locus [{start}, {end}) beyond the end of a linear contig")
    shift = start - start % n
    s, e = start - shift, end - shift
    nt = seq[s:e] if e <= n else seq[s:] + seq[: e - n]
    if strand == "-":
        nt = nt.translate(_COMP)[::-1]
    return translate_text(nt, table)


@cache
def _exchange() -> np.ndarray:
    return exchange_lookup(block_tables())


@cache
def _pos_weights(n: int) -> np.ndarray:
    """Weights of the position tuple: slot ``s`` (the ``s``-th difference) has ``P**(K-1-s)``."""
    p, k = n + 1, 1
    while p ** (k + 1) <= 1 << 62:
        k += 1
    return np.array([p ** (k - 1 - s) for s in range(k)], dtype=np.int64)


@njit(cache=True, nogil=True)
def _dp(
    q: np.ndarray, t: np.ndarray, ex: np.ndarray, pw: np.ndarray
) -> tuple[int, int, np.ndarray, np.ndarray]:  # pragma: no cover - compiled by numba
    """Global ``±2``-band DP: best ``(key, pos)`` at ``(len(q), len(t))``, moves and ranks."""
    nres, nd, b = 19, 5, 32768
    b3 = b * b * b
    big = 1 << 62
    m, n = q.shape[0], t.shape[0]
    key = np.full((m + 1) * nd, big, dtype=np.int64)
    pos = np.zeros((m + 1) * nd, dtype=np.int64)
    mv = np.full((m + 1) * nd, -1, dtype=np.int8)
    rk = np.full((m + 1) * nd, -1, dtype=np.int8)
    key[2] = 0
    nw = pw.shape[0]
    for i in range(1, m + 1):
        q1 = int(q[i - 1])
        q0 = int(q[i - 2]) if i >= 2 else nres
        for dd in range(nd):
            j = i + dd - 2
            if j < 1 or j > n:
                continue
            t1 = int(t[j - 1])
            t0 = int(t[j - 2]) if j >= 2 else nres
            rt = t1 == 22  # a read-through stop (M7-16 rule 10): only move 0, cost 1
            if t1 >= nres and not rt:
                continue
            bk, bp, bm, br = big, 0, -1, -1
            for move in range(4):
                if move == 0:
                    src, pi, r = (i - 1) * nd + dd, i - 1, -1 if q1 == t1 else 0
                elif rt:
                    continue
                elif move == 1:
                    if i < 2 or j < 2 or t0 >= nres or q0 >= nres or q1 >= nres:
                        continue
                    src, pi = (i - 2) * nd + dd, i - 2
                    r = int(ex[2 * nres**3 + ((q0 * nres + q1) * nres + t0) * nres + t1])
                    if q0 != q1 and q0 == t1 and q1 == t0:
                        r = 1
                elif move == 2:
                    if dd < 1 or j < 2 or t0 >= nres or q1 >= nres:
                        continue
                    src, pi = (i - 1) * nd + dd - 1, i - 1
                    r = int(ex[(q1 * nres + t0) * nres + t1])
                else:
                    if i < 2 or dd > 3 or q0 >= nres or q1 >= nres:
                        continue
                    src, pi = (i - 2) * nd + dd + 1, i - 2
                    r = int(ex[nres**3 + (q0 * nres + q1) * nres + t1])
                if key[src] >= big or (move > 0 and r < 0):
                    continue
                c, cp = key[src], pos[src]
                if r >= 0:
                    slot = c // b3
                    c += rank_weight(r, b)
                    if slot < nw:
                        cp += pi * pw[slot]
                if c < bk or (c == bk and cp < bp):
                    bk, bp, bm, br = c, cp, move, r
            key[i * nd + dd], pos[i * nd + dd] = bk, bp
            mv[i * nd + dd], rk[i * nd + dd] = bm, br
    end = m * nd + n - m + 2
    return int(key[end]), int(pos[end]), mv, rk


_STEP = ((1, 1, 0), (2, 2, 0), (1, 2, -1), (2, 1, 1))  # (peptide, genome, Δdd) of each move


def _ops(pep: str, g: str, codons: list[int], off: int) -> tuple[int, list[ResidueOp]]:
    """M3-06 key and ops of the optimal banded alignment of ``pep`` to residues ``g``; every
    ``*`` of ``g`` (an accepted stop site, :func:`_residues`) is consumable as code 22."""
    m, n = len(pep), len(g)
    if abs(n - m) > BAND:
        return _BIG, []
    t = encode_target(g)
    t[np.frombuffer(g.encode(), dtype=np.uint8) == ord("*")] = 22
    key, _, mv, rk = _dp(encode_query(pep), t, _exchange(), _pos_weights(m))
    if key >= _BIG:
        return key, []
    out: list[ResidueOp] = []
    i, dd = m, n - m + 2
    while i > 0:
        a, b, d = _STEP[mv[i * 5 + dd]]
        r, j = int(rk[i * 5 + dd]), i + dd - 2
        ref, alt = g[j - b : j], pep[i - a : i]
        kind: OpKind = _DIFF_KINDS[r] if r >= 0 else "match" if ref == alt else "il"
        out.append(ResidueOp(kind, (i - a + off, i + off), tuple(codons[j - b : j]), ref, alt))
        i, dd = i - a, dd + d
    return key, out[::-1]


def _codons(s: int, e: int, strand: str) -> list[int]:
    return list(range(s, e - 2, 3)) if strand == "+" else list(range(e - 3, s - 1, -3))


def _wrap(codon: int, rep: Replicon) -> int:
    return codon % rep.length if rep.topology == "circular" else codon


def _frame(s: int, e: int, strand: str, n: int) -> int:
    return 1 + (s % n) % 3 if strand == "+" else -(1 + ((n - e) % n) % 3)


def _residues(
    peptide: str,
    seq: str,
    rep: Replicon,
    s: int,
    e: int,
    strand: str,
    stops: frozenset[tuple[str, str, int]] = frozenset(),
    table: int = 11,
) -> str:
    """Genome residues of the locus; a ``*`` is allowed only at a codon in ``stops``."""
    g = translate_locus(seq, rep, s, e, strand, table=table)
    cods = _codons(s, e, strand)
    bad = "X" in g or any(
        a == "*" and (rep.id, strand, c % rep.length) not in stops
        for a, c in zip(g, cods, strict=True)
    )
    if bad:
        raise ValueError(f"peptide {peptide}: genome residue X or * in its locus [{s}, {e})")
    return g


def _parse_blocks(peptide: str, blocks: str) -> tuple[tuple[int, int], tuple[int, int]]:
    try:
        lo, hi = sorted(tuple(int(v) for v in b.split("-")) for b in blocks.split(";"))
    except ValueError as exc:
        raise ValueError(f"peptide {peptide}: malformed blocks {blocks!r}") from exc
    if len(lo) != 2 or len(hi) != 2 or hi[0] - lo[1] not in (-1, 1):
        raise ValueError(f"peptide {peptide}: blocks {blocks!r} are not a ±1-nt split")
    return (lo[0], lo[1]), (hi[0], hi[1])


def align_hit(
    peptide: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    n_diff: int,
    blocks: str | None,
    seq: str,
    rep: Replicon,
    stop_sites: frozenset[tuple[str, str, int]] = frozenset(),
    *,
    table: int = 11,
    starts: frozenset[tuple[str, str, int]] = frozenset(),
) -> HitAlignment:
    """Residue alignment of one hit ``[start, end)``; raises ValueError if it disagrees with
    ``n_diff`` (or, with ``n_diff = 0``, if the letters differ beyond I/L, except a residue-0
    ``M`` over a ``(contig, strand, codon_start)`` of ``starts``: op ``initiator``, M11-20). A
    genome ``*`` is a substitution only at a ``(contig, strand, codon_start)`` of ``stop_sites``
    and only on a non-frameshift row (M7-16 rule 10). ``table`` is the genetic code of the genome
    residues."""
    if strand not in ("+", "-"):
        raise ValueError(f"peptide {peptide}: strand {strand!r}")
    std: Literal["+", "-"] = "+" if strand == "+" else "-"
    if blocks is None:
        g = _residues(peptide, seq, rep, start, end, strand, stop_sites, table)
        cods = _codons(start, end, strand)
        frame = _frame(start, end, strand, rep.length)
        if n_diff == 0:
            fold_p, fold_g = peptide.replace("I", "L"), g.replace("I", "L")
            init = (
                len(g) == len(peptide) > 0
                and fold_p != fold_g
                and peptide[0] == "M"
                and fold_p[1:] == fold_g[1:]
                and (contig, strand, _wrap(cods[0], rep)) in starts
            )
            if not init and (len(g) != len(peptide) or fold_g != fold_p):
                raise ValueError(f"peptide {peptide}: not exact on {g} (n_diff 0)")
            ops = [
                ResidueOp("match" if a == r else "il", (k, k + 1), (cods[k],), r, a)
                for k, (a, r) in enumerate(zip(peptide, g, strict=True))
            ]
            if init:
                ops[0] = ResidueOp("initiator", (0, 1), (cods[0],), g[0], "M")
            return HitAlignment(peptide, contig, std, frame, tuple(ops), None)
        key, ops = _ops(peptide, g, cods, 0)
        _check(peptide, n_diff, key)
        return HitAlignment(peptide, contig, std, frame, tuple(ops), None)
    lo, hi = _parse_blocks(peptide, blocks)
    pre, suf = (lo, hi) if strand == "+" else (hi, lo)
    gp = _residues(peptide, seq, rep, *pre, strand, table=table)
    gs = _residues(peptide, seq, rep, *suf, strand, table=table)
    cp, cs = _codons(*pre, strand), _codons(*suf, strand)
    best: tuple[int, tuple[int, ...], int, list[ResidueOp]] | None = None
    for p in range(max(1, len(gp) - BAND), min(len(peptide), len(gp) + BAND + 1)):
        ka, oa = _ops(peptide[:p], gp, cp, 0)
        kb, ob = _ops(peptide[p:], gs, cs, p)
        if ka >= _BIG or kb >= _BIG:
            continue
        ops = oa + ob
        diffs = tuple(o.pep[0] for o in ops if o.kind not in ("match", "il"))
        cand = (ka + kb, diffs, p, ops)
        if best is None or cand[:3] < best[:3]:
            best = cand
    _check(peptide, n_diff, _BIG if best is None else best[0], frameshift=True)
    assert best is not None
    shift: Literal[-1, 1] = 1 if hi[0] - lo[1] == 1 else -1
    event = FrameshiftEvent(best[2], cp[-1], cs[0], shift)
    frame = _frame(*pre, strand, rep.length)
    return HitAlignment(peptide, contig, std, frame, tuple(best[3]), event)


def _check(peptide: str, n_diff: int, key: int, frameshift: bool = False) -> None:
    got = "none within the band" if key >= _BIG else str(decode_key(key)[0] + int(frameshift))
    if got != str(n_diff):
        raise ValueError(f"peptide {peptide}: n_diff {n_diff}, recomputed {got}")


def align_hits(
    peptide_hits: pd.DataFrame,
    hits_tolerant: pd.DataFrame | None,
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    stop_sites: frozenset[tuple[str, str, int]] = frozenset(),
    *,
    table: int = 11,
    starts: frozenset[tuple[str, str, int]] = frozenset(),
) -> list[HitAlignment]:
    """One :class:`HitAlignment` per ``peptide_hits`` row (0-based half-open), in row order;
    ``starts`` holds the initiator start-site keys (M11-20, :func:`align_hit`).

    Frameshift rows take their ``blocks`` from the matching target row of ``hits_tolerant``
    (peptide, contig, start, end, strand); a missing row or missing blocks raise ValueError.
    """
    work: dict[tuple[str, str, int, int, str], object] = {}
    if hits_tolerant is not None:
        for r in hits_tolerant[hits_tolerant["search"] == "target"].to_dict("records"):
            k = (r["peptide"], r["contig"], int(r["start"]), int(r["end"]), r["strand"])
            work.setdefault(k, r["blocks"])
    out: list[HitAlignment] = []
    for row in peptide_hits.to_dict("records"):
        pep, contig, s, e, strand = (row["peptide"], row["contig"], int(row["start"]),
                                     int(row["end"]), row["strand"])  # fmt: skip
        blocks: str | None = None
        dt = row["diff_types"]
        if isinstance(dt, str) and "frameshift" in dt:
            if (pep, contig, s, e, strand) not in work:
                raise ValueError(f"peptide {pep}: frameshift row without a hits_tolerant row")
            b = work[(pep, contig, s, e, strand)]
            if not isinstance(b, str) or not b:
                raise ValueError(f"peptide {pep}: frameshift row without blocks")
            blocks = b
        n_diff = int(row["n_diff"])
        out.append(
            align_hit(pep, contig, s, e, strand, n_diff, blocks, genome[contig], replicons[contig],
                      stop_sites, table=table, starts=starts)
        )  # fmt: skip
    return out
