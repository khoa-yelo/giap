"""Annotated stop codons the tolerant matcher may consume (ADR-0058 items 1-2; M7-16 rules 1-3).

A *stop site* is the last codon of a CDS or sORF (pseudogenes never) when that codon is a table-11
stop. A *cell* is the frame position (aa index) of a site's ``*`` in the six-frame index, one per
copy in the origin-extended frame of a circular replicon. Decoy cells mirror the target cells:
``n - 1 - p`` in a ``reverse`` decoy frame of ``n`` residues, ``p`` in a ``shuffle`` decoy (stops
keep their positions). Everything here is pure and deterministic; nothing draws randomness (§6.6).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from giap.mapping.decoys import DecoyType
from giap.mapping.sixframe import STOP, Frame, SixFrameIndex
from giap.model.genes import GeneModel

__all__ = [
    "StopCells",
    "StopSite",
    "annotated_stops",
    "consumed_sites",
    "frame_cells",
    "mirror_cells",
    "site_keys",
]

STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


@dataclass(frozen=True, slots=True)
class StopSite:
    """One consumable stop codon: ``[codon_start, codon_start + 3)`` (start wrapped) on ``strand``.

    ``codon`` is read on the gene strand; ``locus_tags`` (sorted) share the codon.
    """

    contig: str
    strand: Literal["+", "-"]
    codon_start: int
    locus_tags: tuple[str, ...]
    codon: Literal["TAA", "TAG", "TGA"]
    kind: Literal["terminal", "recoded"]


@dataclass(frozen=True, slots=True)
class StopCells:
    """Sites and their cells: per frame id, sorted ``int64`` aa positions (target and decoys)."""

    sites: tuple[StopSite, ...]
    target: tuple[np.ndarray, ...]
    decoys: Mapping[DecoyType, tuple[np.ndarray, ...]]


def _codon(seq: str, start: int, n: int) -> str:
    return "".join(seq[(start + i) % n] for i in range(3)).upper()


def annotated_stops(gm: GeneModel, genome: Mapping[str, str]) -> tuple[tuple[StopSite, ...], int]:
    """Terminal stop sites of every CDS/sORF and the number of features whose last codon is not a
    stop (or runs off a linear contig). Sorted by (contig order, strand, codon_start)."""
    tags: dict[tuple[str, str, int], tuple[str, list[str]]] = {}
    n_not = 0
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf"):
            continue
        iv = g.interval
        rep = gm.replicons[iv.contig]
        n, circular = rep.length, rep.topology == "circular"
        start = iv.end - 3 if iv.strand == "+" else iv.start
        if not circular and (start < 0 or start + 3 > n):
            n_not += 1
            continue
        codon = _codon(genome[iv.contig], start, n)
        if iv.strand == "-":
            codon = codon.translate(_COMP)[::-1]
        if codon not in STOP_CODONS:
            n_not += 1
            continue
        key = (iv.contig, iv.strand, start % n)
        tags.setdefault(key, (codon, []))[1].append(g.locus_tag)
    order = {c: i for i, c in enumerate(gm.replicons)}
    sites = tuple(
        StopSite(c, s, p, tuple(sorted(lt)), codon, "terminal")  # type: ignore[arg-type]
        for (c, s, p), (codon, lt) in sorted(
            tags.items(), key=lambda kv: (order[kv[0][0]], *kv[0][1:])
        )
    )
    return sites, n_not


def frame_cells(index: SixFrameIndex, sites: Sequence[StopSite]) -> tuple[np.ndarray, ...]:
    """Cells of ``sites`` per frame id (rule 2), origin copies included.

    Raises:
        ValueError: if a cell does not map back to its site's codon or does not hold ``*``.
    """
    fid_of = {(f.contig, f.strand, f.offset): i for i, f in enumerate(index.frames)}
    out: list[list[int]] = [[] for _ in index.frames]
    for site in sites:
        rep = index.replicons[site.contig]
        ext = index.frames[fid_of[(site.contig, site.strand, 0)]].ext_nt
        xs = [site.codon_start]
        if rep.topology == "circular" and site.codon_start + rep.length + 3 <= ext:
            xs.append(site.codon_start + rep.length)
        for x in xs:
            if site.strand == "+":
                off, p = x % 3, x // 3
            else:
                off = (ext - 3 - x) % 3
                p = (ext - off - 3 - x) // 3
            fid = fid_of[(site.contig, site.strand, off)]
            aa = index.frames[fid].aa
            ok = 0 <= p < len(aa) and int(aa[p]) == STOP
            if ok:
                iv = index.to_interval(fid, p, p + 1)
                ok = iv.start == site.codon_start and iv.strand == site.strand
            if not ok:
                raise ValueError(f"stop site {site} is not a stop in frame {fid} at aa {p}")
            out[fid].append(p)
    return tuple(np.array(sorted(c), dtype=np.int64) for c in out)


def mirror_cells(
    kind: DecoyType, frames: tuple[Frame, ...], cells: tuple[np.ndarray, ...]
) -> tuple[np.ndarray, ...]:
    """Decoy cells of ``kind`` in ``frames`` (that decoy's frames): ``n - 1 - p`` (``reverse``) or
    ``p`` (``shuffle``). Draws no randomness.

    Raises:
        ValueError: if the frame counts differ or a mirrored cell does not hold ``*``.
    """
    if len(frames) != len(cells):
        raise ValueError(f"{len(cells)} cell arrays for {len(frames)} {kind} decoy frames")
    out = []
    for fid, (f, c) in enumerate(zip(frames, cells, strict=True)):
        c = np.asarray(c, dtype=np.int64)
        m = np.sort(len(f.aa) - 1 - c) if kind == "reverse" else c.copy()
        if len(m) and (m.min() < 0 or m.max() >= len(f.aa) or np.any(f.aa[m] != STOP)):
            raise ValueError(f"a mirrored {kind} cell of frame {fid} does not hold '*'")
        out.append(m)
    return tuple(out)


def site_keys(sites: Sequence[StopSite]) -> frozenset[tuple[str, str, int]]:
    """``(contig, strand, codon_start)`` of every site (the M7-01 re-alignment set)."""
    return frozenset((s.contig, s.strand, s.codon_start) for s in sites)


def consumed_sites(
    rows: pd.DataFrame,
    sites: Sequence[StopSite],
    lengths: Mapping[str, int],
    circular: Mapping[str, bool],
) -> list[tuple[StopSite, ...]]:
    """Per row, the sites its window consumes, sorted by ``codon_start`` (rule 15).

    Only target rows without ``blocks`` consume; a site is consumed when a copy ``x`` of its codon
    (``codon_start``, plus ``+ n`` on a circular replicon) lies in frame inside ``[start, end)``.
    """
    by: dict[tuple[str, str], list[StopSite]] = {}
    for s in sorted(sites, key=lambda s: s.codon_start):
        by.setdefault((s.contig, s.strand), []).append(s)
    blocks = rows["blocks"] if "blocks" in rows.columns else [None] * len(rows)
    out: list[tuple[StopSite, ...]] = []
    for search, contig, strand, a, b, blk in zip(
        rows["search"],
        rows["contig"],
        rows["strand"],
        rows["start"],
        rows["end"],
        blocks,
        strict=True,
    ):
        cand = by.get((str(contig), str(strand)), [])
        if search != "target" or (isinstance(blk, str) and blk) or not cand:
            out.append(())
            continue
        a, b, n = int(a), int(b), int(lengths[str(contig)])
        shifts = (0, n) if circular.get(str(contig), False) else (0,)
        out.append(
            tuple(
                s
                for s in cand
                if any(
                    a <= x and x + 3 <= b and (x - a) % 3 == 0
                    for x in (s.codon_start + d for d in shifts)
                )
            )
        )
    return out
