"""Route 2: tolerant matching, seed and verify (GIPPA_SPEC.md §7.2.4, G4; issue M3-07).

Every peptide of ``min_len`` .. 40 aa (more than ``2 * BAND`` residues; ``k`` from
``max_diff``) is seeded with the lossless q-gram threshold of :func:`~gippa.mapping.qgram.
choose_q` and each merged candidate diagonal is verified with the banded DP of M3-06
(:func:`~gippa.mapping.tolerant_dp.verify_batch`). Per peptide and search all loci at the minimal
cost (``<= k``) are kept. The target search runs only for ``target_queries`` (exact leftovers,
M3-10 Q1); the decoy search runs for every searched peptide against the decoy frames of its own
type. Target windows map with ``to_interval`` and are de-duplicated modulo the replicon length
(circular extension); decoy windows keep their mirrored target interval (M3-03).

The frames of one search are concatenated with stop separators (:func:`~gippa.mapping.qgram.
concat_frames`), so verification runs on one text and never crosses frames. Queries run in a
``numba.prange``; each query's result is independent and rows are sorted afterwards, so the
output does not depend on the thread count (§6.6). q-gram indexes are built lazily, only for
the q values needed.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from numba import njit, prange

from gippa.config.model import TolerantConfig
from gippa.mapping.decoys import DecoyType, mirror_aa
from gippa.mapping.diffops import block_tables
from gippa.mapping.hits import hits_frame
from gippa.mapping.qgram import (
    MAX_QUERY,
    QGramIndex,
    _candidates_kernel,
    _index_from_text,
    choose_q,
    concat_frames,
)
from gippa.mapping.sixframe import READTHROUGH, STOP, Frame, SixFrameIndex, decode, encode_query
from gippa.mapping.tolerant_dp import BAND, decode_key, exchange_lookup, verify_batch
from gippa.model.coords import Interval

if TYPE_CHECKING:
    from gippa.mapping.readthrough import StopCells

__all__ = ["Search", "budget", "decoy_marks", "match_tolerant", "parse_max_diff"]

_SEARCH_ORDER = {"target": 0, "decoy": 1}
_CAP = 16  # rows per query kept by the parallel pass; larger results are recomputed serially
_COST_UNIT = 1 << 45  # 32768**3: key // _COST_UNIT is the cost (tolerant_dp key layout)
_Row = tuple[int, int, int, int]  # (key, frame_id, aa_start, aa_end)


@njit(cache=True, nogil=True)
def _search_one(
    qc: np.ndarray,
    k: int,
    t: int,
    q: int,
    offsets: np.ndarray,
    positions: np.ndarray,
    text: np.ndarray,
    ex: np.ndarray,
) -> np.ndarray:  # pragma: no cover - compiled by numba; exercised through match_tolerant
    """Rows ``(key, t_start, t_end)`` (global) of every minimal-cost candidate of one query."""
    cands = _candidates_kernel(qc, q, offsets, positions, t)
    if cands.shape[0] == 0:
        return np.empty((0, 3), dtype=np.int64)
    res = verify_batch(qc, text, cands, k, ex)
    unit = 1 << 45
    best = 1 << 62
    for r in range(res.shape[0]):
        if res[r, 0] >= 0:
            c = res[r, 0] // unit
            if c < best:
                best = c
    n = 0
    for r in range(res.shape[0]):
        if res[r, 0] >= 0 and res[r, 0] // unit == best:
            n += 1
    out = np.empty((n, 3), dtype=np.int64)
    n = 0
    for r in range(res.shape[0]):
        if res[r, 0] >= 0 and res[r, 0] // unit == best:
            out[n] = res[r]
            n += 1
    return out


@njit(cache=True, parallel=True)
def _search_many(
    qcat: np.ndarray,
    qoff: np.ndarray,
    ks: np.ndarray,
    t: np.ndarray,
    q: int,
    offsets: np.ndarray,
    positions: np.ndarray,
    text: np.ndarray,
    ex: np.ndarray,
    cap: int,
) -> tuple[np.ndarray, np.ndarray]:  # pragma: no cover - compiled by numba; via match_tolerant
    """:func:`_search_one` for every query of ``qcat`` (row counts and the first ``cap`` rows)."""
    nq = qoff.shape[0] - 1
    counts = np.zeros(nq, dtype=np.int64)
    rows = np.full((nq, cap, 3), -1, dtype=np.int64)
    for j in prange(nq):
        res = _search_one(qcat[qoff[j] : qoff[j + 1]], ks[j], t[j], q, offsets, positions, text, ex)
        m = res.shape[0]
        counts[j] = m
        for r in range(min(m, cap)):
            rows[j, r, 0] = res[r, 0]
            rows[j, r, 1] = res[r, 1]
            rows[j, r, 2] = res[r, 2]
    return counts, rows


def parse_max_diff(max_diff: Mapping[str, int]) -> list[tuple[int, int, int]]:
    """``[(lo, hi, k)]`` from config keys ``"lo-hi"``.

    Raises:
        ValueError: if a key is not ``"lo-hi"``, ``hi`` exceeds ``MAX_QUERY`` (512) or ``k`` is
            outside ``0 .. BAND`` (2): the q-gram/band search is lossless only for ``k <= 2``.
    """
    out = []
    for key, k in max_diff.items():
        lo, sep, hi = key.partition("-")
        if not sep or not lo.strip().isdigit() or not hi.strip().isdigit():
            raise ValueError(f"mapping.tolerant.max_diff key {key!r} is not 'lo-hi'")
        if int(hi) > MAX_QUERY:
            raise ValueError(
                f"mapping.tolerant.max_diff class {key!r} exceeds MAX_QUERY = {MAX_QUERY} residues"
            )
        if not 0 <= int(k) <= BAND:
            raise ValueError(
                f"mapping.tolerant.max_diff[{key!r}] = {k}: k must lie in 0..{BAND} (lossless band)"
            )
        out.append((int(lo), int(hi), int(k)))
    return sorted(out)


def budget(L: int, cfg: TolerantConfig, ranges: list[tuple[int, int, int]]) -> int | None:
    """Difference budget ``k`` of a query of ``L`` residues, or ``None`` when not searched."""
    if not cfg.min_len <= L <= MAX_QUERY or L <= 2 * BAND:
        return None
    for lo, hi, k in ranges:
        if lo <= L <= hi:
            return k
    return None


class Search:
    """One search space (target frames or the decoy frames of one type) with lazy indexes."""

    def __init__(
        self, frames: tuple[Frame, ...], marked: Sequence[np.ndarray] | None = None
    ) -> None:
        """``marked``: per frame, the cells written as ``READTHROUGH`` into the search text."""
        self.frames = frames
        self.text, self.starts = concat_frames(frames)
        if marked is not None:
            if len(marked) != len(frames):
                raise ValueError(f"{len(marked)} cell arrays for {len(frames)} frames")
            for fid, cells in enumerate(marked):
                pos = int(self.starts[fid]) + np.asarray(cells, dtype=np.int64)
                if np.any(self.text[pos] != STOP):
                    raise ValueError(f"a marked cell of frame {fid} does not hold '*'")
                self.text[pos] = READTHROUGH
        self._indexes: dict[int, QGramIndex] = {}

    def index(self, q: int) -> QGramIndex:
        """The q-gram index of the concatenated frames (built on first use)."""
        if q not in self._indexes:
            self._indexes[q] = _index_from_text(self.text, self.starts, q)
        return self._indexes[q]

    def run(
        self, queries: Sequence[str], budget: Mapping[str, int], ex: np.ndarray
    ) -> dict[str, list[_Row]]:
        """Minimal-cost windows ``(key, frame_id, aa_start, aa_end)`` of every query."""
        groups: dict[int, list[tuple[str, int, int]]] = {}
        for text in queries:
            k = budget[text]
            q, t = choose_q(len(text), k)
            groups.setdefault(q, []).append((text, k, t))
        out: dict[str, list[_Row]] = {}
        for q, items in sorted(groups.items()):
            idx = self.index(q)
            codes = [encode_query(text) for text, _, _ in items]
            qoff = np.zeros(len(codes) + 1, dtype=np.int64)
            np.cumsum([len(c) for c in codes], out=qoff[1:])
            qcat = np.ascontiguousarray(np.concatenate(codes), dtype=np.uint8)
            ks = np.array([k for _, k, _ in items], dtype=np.int64)
            ts = np.array([t for _, _, t in items], dtype=np.int64)
            counts, rows = _search_many(
                qcat, qoff, ks, ts, q, idx.offsets, idx.positions, self.text, ex, _CAP
            )
            for j, (text, k, t) in enumerate(items):
                m = int(counts[j])
                if m > _CAP:
                    res = _search_one(
                        np.ascontiguousarray(codes[j], dtype=np.uint8),
                        k,
                        t,
                        q,
                        idx.offsets,
                        idx.positions,
                        self.text,
                        ex,
                    )
                else:
                    res = rows[j, :m]
                out[text] = [self._local(int(r[0]), int(r[1]), int(r[2])) for r in res]
        return out

    def _local(self, key: int, s: int, e: int) -> _Row:
        fid = int(np.searchsorted(self.starts, s, side="right")) - 1
        a = int(self.starts[fid])
        return key, fid, s - a, e - a


def decoy_marks(cells: StopCells, kind: DecoyType) -> tuple[np.ndarray, ...]:
    """The mirrored cells of decoy type ``kind`` (never the target cells).

    Raises:
        ValueError: if ``cells`` carry no cells for ``kind``.
    """
    if kind not in cells.decoys:
        raise ValueError(f"stop cells have no {kind} decoy cells")
    return cells.decoys[kind]


def match_tolerant(
    queries: Sequence[str],
    index: SixFrameIndex,
    decoys: Mapping[DecoyType, tuple[Frame, ...]],
    decoy_of: Mapping[str, DecoyType],
    target_queries: set[str],
    cfg: TolerantConfig,
    *,
    cells: StopCells | None = None,
) -> pd.DataFrame:
    """Tolerant target and decoy hits (``HIT_COLUMNS``, ``route = "tolerant"``).

    ``cells``: consumable annotated stops (ADR-0058); the target search marks ``cells.target``,
    each decoy search the mirrored cells of its type. ``None`` or empty cells change nothing.

    Peptides shorter than ``cfg.min_len``, of at most ``2 * BAND`` (4) residues or outside every
    ``max_diff`` length class are not searched. Rows are sorted by peptide, search (target
    first), decoy type, contig order, start, end, strand; ``diff_types`` is ``;``-joined and
    sorted.

    Raises:
        ValueError: if a ``decoy_of`` type has no decoy frames, a searched peptide has no decoy
            type, a ``max_diff`` class ends above ``MAX_QUERY`` (512) residues, or a configured
            ``k`` lies outside ``0 .. BAND`` (2) or leaves no lossless q-gram threshold.
    """
    for pep, kind in sorted(decoy_of.items()):
        if kind not in decoys:
            raise ValueError(f"decoy type {kind!r} of peptide {pep!r} has no decoy frames")
    ranges = parse_max_diff(cfg.max_diff)
    norm: dict[str, str] = {}
    budgets: dict[str, int] = {}
    for pep in sorted(set(queries)):
        k = budget(len(pep), cfg, ranges)
        if k is None:
            continue
        if pep not in decoy_of:
            raise ValueError(f"peptide {pep!r} has no decoy type in decoy_of")
        text = decode(encode_query(pep))
        norm[pep] = text
        budgets[text] = k
    if not norm:
        return hits_frame()

    ex = exchange_lookup(block_tables())
    contig_order = {c: i for i, c in enumerate(index.replicons)}
    rows: list[dict[str, object]] = []

    def emit(
        found: Mapping[str, list[_Row]],
        frames: tuple[Frame, ...],
        kind: DecoyType | None,
    ) -> None:
        for pep, text in norm.items():
            if kind is None and pep not in target_queries:
                continue
            if kind is not None and decoy_of[pep] != kind:
                continue
            loci: dict[Interval, tuple[int, int, int, int]] = {}
            for key, fid, a, b in sorted(found.get(text, ())):
                rep_len = index.replicons[frames[fid].contig].length
                if 3 * (b - a) > rep_len:
                    continue
                ta, tb = (a, b) if kind is None else mirror_aa(kind, len(frames[fid].aa), a, b)
                loci.setdefault(index.to_interval(fid, ta, tb), (key, fid, a, b))
            for iv, (key, fid, a, b) in loci.items():
                cost, types = decode_key(key)
                rows.append(
                    {
                        "peptide": pep,
                        "query": text,
                        "search": "target" if kind is None else "decoy",
                        "decoy_type": kind,
                        "contig": iv.contig,
                        "frame_id": fid,
                        "aa_start": a,
                        "aa_end": b,
                        "start": iv.start,
                        "end": iv.end,
                        "strand": iv.strand,
                        "frame": index.frame_label(iv),
                        "route": "tolerant",
                        "n_diff": cost,
                        "diff_types": ";".join(sorted(types)),
                        "blocks": None,
                    }
                )

    targets = sorted({norm[p] for p in norm if p in target_queries})
    if targets:
        marks = None if cells is None else cells.target
        emit(Search(index.frames, marks).run(targets, budgets, ex), index.frames, None)
    for kind in sorted(decoys):
        wanted = sorted({text for p, text in norm.items() if decoy_of[p] == kind})
        if wanted:
            frames = decoys[kind]
            marks = None if cells is None else decoy_marks(cells, kind)
            emit(Search(frames, marks).run(wanted, budgets, ex), frames, kind)

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
