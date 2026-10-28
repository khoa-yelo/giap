"""Annotated stop codons the tolerant matcher may consume (ADR-0058 items 1-2, 6; M7-16 rules 1-3,
21-22).

A *stop site* is the last codon of a CDS or sORF (pseudogenes never) when that codon is a stop of
the run's genetic code (kind ``terminal``), or, with ``recoded=True``, an internal codon the
annotated protein reads as Sec (``U`` over TGA, table 11 only) or Pyl (``O`` over TAG) (kind
``recoded``; ``recoded_codons``, ADR-0067). A *cell* is the frame
position (aa index) of a site's ``*`` in the six-frame index, one per copy in the origin-extended
frame of a circular replicon. Decoy cells mirror the target cells:
``n - 1 - p`` in a ``reverse`` decoy frame of ``n`` residues, ``p`` in a ``shuffle`` decoy (stops
keep their positions). Everything here is pure and deterministic; nothing draws randomness (§6.6).
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from gippa.mapping.decoys import DecoyType
from gippa.mapping.sixframe import STOP, Frame, SixFrameIndex
from gippa.model.genes import GeneModel
from gippa.model.genetic_code import recoded_codons, stop_codons

__all__ = [
    "StopCells",
    "StopSite",
    "annotated_stops",
    "consumed_sites",
    "frame_cells",
    "mirror_cells",
    "recoded_skips",
    "site_keys",
]

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_RECODED_LETTERS = frozenset("UO")  # Bakta Sec / Pyl; the codon each needs comes from the table


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


def _read(seq: str, start: int, n: int, strand: str) -> str:
    codon = _codon(seq, start, n)
    return codon.translate(_COMP)[::-1] if strand == "-" else codon


def _scan(
    gm: GeneModel, genome: Mapping[str, str], recoded: bool, table: int
) -> tuple[dict[tuple[str, str, int], tuple[str, set[str], bool]], int, dict[str, int]]:
    """Rule 1 (terminal) and rule 21 (recoded, when ``recoded``): merged sites keyed by
    ``(contig, strand, codon_start)`` -> (codon, locus tags, recoded?), ``n_not`` and skips.
    A ``U`` / ``O`` with no recoding codon in ``table`` is a ``codon_mismatch`` (M11-16 D6)."""
    stops, recoding = stop_codons(table), recoded_codons(table)
    tags: dict[tuple[str, str, int], tuple[str, set[str], bool]] = {}
    skips = {"length_mismatch": 0, "codon_mismatch": 0}
    n_not = 0

    def add(key: tuple[str, str, int], codon: str, tag: str, rec: bool) -> None:
        old = tags.get(key)
        tags[key] = (codon, ({tag} | old[1]) if old else {tag}, rec or (old[2] if old else False))

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
        seq = genome[iv.contig]
        codon = _read(seq, start, n, iv.strand)
        if codon not in stops:
            n_not += 1
            continue
        add((iv.contig, iv.strand, start % n), codon, g.locus_tag, False)
        if not recoded or not g.aa:  # empty protein: not scanned, not counted (lead decision)
            continue
        span = iv.end - iv.start
        if span % 3 or len(g.aa) != span // 3 - 1:
            skips["length_mismatch"] += 1
            continue
        for i in range(1, len(g.aa)):
            if g.aa[i] not in _RECODED_LETTERS:
                continue
            want = recoding.get(g.aa[i])
            cs = iv.start + 3 * i if iv.strand == "+" else iv.end - 3 * i - 3
            if want is None or _read(seq, cs, n, iv.strand) != want:
                skips["codon_mismatch"] += 1
                continue
            add((iv.contig, iv.strand, cs % n), want, g.locus_tag, True)
    return tags, n_not, skips


def annotated_stops(
    gm: GeneModel, genome: Mapping[str, str], recoded: bool = False, *, table: int = 11
) -> tuple[tuple[StopSite, ...], int]:
    """Terminal stop sites of every CDS/sORF and the number of features whose last codon is not a
    stop of ``table`` (or runs off a linear contig). With ``recoded``, also every annotated Sec
    (``U`` over TGA, table 11) and Pyl (``O`` over TAG) codon of those features (rule 21); a key
    shared with a terminal site is one ``recoded`` site. Sorted by (contig order, strand,
    codon_start)."""
    tags, n_not, _ = _scan(gm, genome, recoded, table)
    order = {c: i for i, c in enumerate(gm.replicons)}
    sites = tuple(
        StopSite(c, s, p, tuple(sorted(lt)), codon, "recoded" if rec else "terminal")  # type: ignore[arg-type]
        for (c, s, p), (codon, lt, rec) in sorted(
            tags.items(), key=lambda kv: (order[kv[0][0]], *kv[0][1:])
        )
    )
    return sites, n_not


def recoded_skips(gm: GeneModel, genome: Mapping[str, str], *, table: int = 11) -> dict[str, int]:
    """Rule 21 (b)/(e) counts: ``length_mismatch`` per feature, ``codon_mismatch`` per residue."""
    return _scan(gm, genome, True, table)[2]


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
    # starts flattened once per (contig, strand), so each row bisects instead of scanning
    starts = {k: [s.codon_start for s in v] for k, v in by.items()}
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
        key = (str(contig), str(strand))
        cand = by.get(key, [])
        if search != "target" or (isinstance(blk, str) and blk) or not cand:
            out.append(())
            continue
        xs = starts[key]
        a, b, n = int(a), int(b), int(lengths[str(contig)])
        shifts = (0, n) if circular.get(str(contig), False) else (0,)
        hits: set[int] = set()
        for d in shifts:
            # x = codon_start + d with a <= x and x + 3 <= b
            lo = bisect_left(xs, a - d)
            hi = bisect_right(xs, b - 3 - d)
            hits.update(k for k in range(lo, hi) if (xs[k] + d - a) % 3 == 0)
        out.append(tuple(cand[k] for k in sorted(hits)))
    return out
