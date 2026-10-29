"""Route 1: exact peptide matching with Aho-Corasick (GIPPA_SPEC.md §7.2.2, decision 5).

All unique queries (I -> L, letters outside the residue alphabet never match) of
``min_len..max_len`` aa go into one automaton that scans every target frame and, per peptide,
the decoy frames of its own source (``decoy_of``). Every locus is kept; origin copies on circular
replicons are mapped once (decoy hits on their mirrored target interval). Rows follow
``HIT_COLUMNS`` and are sorted by peptide, search (target first), contig order, start, end,
strand.

With ``starts`` (M11-20 rule 3) an ``M``-initial query is also credited at the annotated non-ATG
start cells: window ``[p, p + k)`` in the target and ``shuffle`` frames, the mirrored
``[c - k + 1, c + 1)`` in ``reverse`` frames; credited rows are ordinary exact rows.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import ahocorasick
import numpy as np
import pandas as pd

from gippa.mapping.decoys import DecoyType, mirror_aa
from gippa.mapping.hits import hits_frame
from gippa.mapping.initiator import StartCells, credited_window
from gippa.mapping.sixframe import Frame, SixFrameIndex, decode, encode_query
from gippa.model.coords import Interval

_SEARCH_ORDER = {"target": 0, "decoy": 1}
# (frame_id, aa_start, aa_end) of the first occurrence of each mapped interval, per key.
_Loci = dict[int, dict[Interval, tuple[int, int, int]]]


def _normalise(peptide: str) -> str | None:
    """I -> L query text, or ``None`` when a letter can never match (decodes to ``#``)."""
    q = decode(encode_query(peptide))
    return None if "#" in q else q


def _scan(
    automaton: ahocorasick.Automaton,
    frames: Sequence[Frame],
    texts: Sequence[str],
    index: SixFrameIndex,
    lengths: Sequence[int],
    wanted: set[int] | None,
    kind: DecoyType | None,
) -> _Loci:
    """Every hit of the automaton in ``frames``, deduplicated on its (mirrored) interval."""
    out: _Loci = {}
    for fid, (frame, text) in enumerate(zip(frames, texts, strict=True)):
        n_aa = len(text)
        rep_len = index.replicons[frame.contig].length
        for end, key in automaton.iter(text):
            if wanted is not None and key not in wanted:
                continue
            k = lengths[key]
            if 3 * k > rep_len:
                continue
            a, b = end + 1 - k, end + 1
            ta, tb = (a, b) if kind is None else mirror_aa(kind, n_aa, a, b)
            iv = index.to_interval(fid, ta, tb)
            out.setdefault(key, {}).setdefault(iv, (fid, a, b))
    return out


def _credit(
    out: _Loci,
    by_len: Mapping[int, Mapping[str, int]],
    frames: Sequence[Frame],
    texts: Sequence[str],
    cells: Sequence[np.ndarray],
    index: SixFrameIndex,
    wanted: set[int] | None,
    kind: DecoyType | None,
) -> None:
    """Add the credited windows of ``M``-initial queries at ``cells`` to ``out`` (rule 3)."""
    rule = kind or "target"
    for fid, (frame, text, cl) in enumerate(zip(frames, texts, cells, strict=True)):
        if not len(cl):
            continue
        n_aa = len(text)
        rep_len = index.replicons[frame.contig].length
        for c in (int(x) for x in cl):
            for k, rests in by_len.items():
                w = credited_window(rule, n_aa, c, k)
                if w is None or 3 * k > rep_len:
                    continue
                a, b = w
                key = rests.get(text[a : b - 1] if rule == "reverse" else text[a + 1 : b])
                if key is None or (wanted is not None and key not in wanted):
                    continue
                ta, tb = (a, b) if kind is None else mirror_aa(kind, n_aa, a, b)
                iv = index.to_interval(fid, ta, tb)
                out.setdefault(key, {}).setdefault(iv, (fid, a, b))


def match_exact(
    queries: Sequence[str],
    index: SixFrameIndex,
    decoys: Mapping[DecoyType, tuple[Frame, ...]],
    decoy_of: Mapping[str, DecoyType],
    *,
    min_len: int = 7,
    max_len: int = 40,
    starts: StartCells | None = None,
) -> pd.DataFrame:
    """Exact target and decoy hits of every query of ``min_len..max_len`` aa (``HIT_COLUMNS``);
    with ``starts`` also the initiator-credited windows (M11-20 rule 3).

    Raises:
        ValueError: on a peptide without a decoy type, a decoy type without frames, or ``starts``
            with sites but no cells for a decoy type the queries use.
    """
    for pep, kind in sorted(decoy_of.items()):
        if kind not in decoys:
            raise ValueError(f"decoy type {kind!r} of peptide {pep!r} has no decoy frames")
    keys: dict[str, int] = {}
    peptides_of: list[list[str]] = []
    for pep in sorted(set(queries)):
        if not min_len <= len(pep) <= max_len:
            continue
        if pep not in decoy_of:
            raise ValueError(f"peptide {pep!r} has no decoy type in decoy_of")
        q = _normalise(pep)
        if q is None:
            continue
        key = keys.setdefault(q, len(keys))
        if key == len(peptides_of):
            peptides_of.append([])
        peptides_of[key].append(pep)
    if not keys:
        return hits_frame()

    automaton = ahocorasick.Automaton(ahocorasick.STORE_INTS)
    for q, key in keys.items():
        automaton.add_word(q, key)
    automaton.make_automaton()
    texts = list(keys)
    lengths = [len(q) for q in texts]
    by_len: dict[int, dict[str, int]] = {}
    if starts is not None and starts.sites:
        for q, key in keys.items():
            if q[0] == "M":
                by_len.setdefault(len(q), {})[q[1:]] = key

    contig_order = {c: i for i, c in enumerate(index.replicons)}
    rows: list[dict[str, object]] = []

    def emit(loci: _Loci, search: str, kind: DecoyType | None) -> None:
        for key, ivs in loci.items():
            for pep in peptides_of[key]:
                if kind is not None and decoy_of[pep] != kind:
                    continue
                for iv, (fid, a, b) in ivs.items():
                    rows.append(
                        {
                            "peptide": pep,
                            "query": texts[key],
                            "search": search,
                            "decoy_type": kind,
                            "contig": iv.contig,
                            "frame_id": fid,
                            "aa_start": a,
                            "aa_end": b,
                            "start": iv.start,
                            "end": iv.end,
                            "strand": iv.strand,
                            "frame": index.frame_label(iv),
                            "route": "exact",
                            "n_diff": 0,
                            "diff_types": "",
                            "blocks": None,
                        }
                    )

    target_texts = [decode(f.aa) for f in index.frames]
    target = _scan(automaton, index.frames, target_texts, index, lengths, None, None)
    if by_len and starts is not None:
        _credit(target, by_len, index.frames, target_texts, starts.target, index, None, None)
    emit(target, "target", None)
    for kind, frames in decoys.items():
        wanted = {
            key for key, peps in enumerate(peptides_of) if any(decoy_of[p] == kind for p in peps)
        }
        if wanted:
            if starts is not None and starts.sites and kind not in starts.decoys:
                raise ValueError(f"start cells have no {kind} decoy cells")
            decoy_texts = [decode(f.aa) for f in frames]
            loci = _scan(automaton, frames, decoy_texts, index, lengths, wanted, kind)
            if by_len and starts is not None:
                _credit(loci, by_len, frames, decoy_texts, starts.decoys[kind], index, wanted, kind)
            emit(loci, "decoy", kind)

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
