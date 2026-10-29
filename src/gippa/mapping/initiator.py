"""Initiator M at annotated non-ATG starts, exact route (M11-20 option A; ADR-0068).

A *start site* is the first codon of a CDS or sORF (pseudogenes never) whose Bakta protein begins
with ``M`` and whose genome codon, read on the gene strand, is in ``annotation.start_codons`` and
is not ``ATG`` (rules 1-3). A *cell* is the frame position (aa index) of a site's codon in the
six-frame index, one per copy in the origin-extended frame of a circular replicon. Decoy cells:
``n - 1 - p`` in a ``reverse`` decoy frame of ``n`` residues (the cell is the *last* position of
its credited window, so the window mirrors onto the target's), ``p`` in a ``shuffle`` decoy.
Everything here is pure and deterministic; nothing draws randomness (§6.6).
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd

from gippa.mapping.decoys import DecoyType
from gippa.mapping.sixframe import STOP, Frame, SixFrameIndex, encode_target, translate_text
from gippa.model.genes import GeneModel

__all__ = [
    "SKIP_KEYS",
    "StartCells",
    "StartSite",
    "annotated_starts",
    "initiator_rows",
    "mirror_start_cells",
    "start_frame_cells",
    "start_keys",
    "usable_windows",
]

SKIP_KEYS: tuple[str, ...] = ("atg", "not_a_start", "aa_not_m", "off_contig")
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_Key = tuple[str, str, int]


@dataclass(frozen=True, slots=True)
class StartSite:
    """One annotated non-ATG start codon ``[codon_start, codon_start + 3)`` (wrapped) on
    ``strand``; ``codon`` is read on the gene strand; ``locus_tags`` (sorted) share the codon."""

    contig: str
    strand: Literal["+", "-"]
    codon_start: int
    locus_tags: tuple[str, ...]
    codon: str


@dataclass(frozen=True, slots=True)
class StartCells:
    """Sites and their cells: per frame id, sorted ``int64`` aa positions (target and decoys)."""

    sites: tuple[StartSite, ...]
    target: tuple[np.ndarray, ...]
    decoys: Mapping[DecoyType, tuple[np.ndarray, ...]]


def _read(seq: str, start: int, n: int, strand: str) -> str:
    codon = "".join(seq[(start + i) % n] for i in range(3)).upper()
    return codon.translate(_COMP)[::-1] if strand == "-" else codon


def annotated_starts(
    gm: GeneModel, genome: Mapping[str, str], start_codons: Collection[str], *, table: int
) -> tuple[tuple[StartSite, ...], dict[str, int]]:
    """Start sites of every CDS/sORF (rules 1-3) and the skip counts ``atg``, ``not_a_start``,
    ``aa_not_m`` and ``off_contig``; features sharing a codon are one site. Sorted by (contig
    order, strand, codon_start). ``table`` is accepted for the interface: rule 2 reads codons."""
    del table
    want = {c.upper() for c in start_codons}
    skipped = dict.fromkeys(SKIP_KEYS, 0)
    found: dict[_Key, tuple[str, set[str]]] = {}
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf"):
            continue
        if not g.aa or g.aa[0] != "M":
            skipped["aa_not_m"] += 1
            continue
        iv = g.interval
        rep = gm.replicons[iv.contig]
        n = rep.length
        x = iv.start if iv.strand == "+" else iv.end - 3
        if rep.topology != "circular" and not 0 <= x <= n - 3:
            skipped["off_contig"] += 1
            continue
        codon = _read(genome[iv.contig], x, n, iv.strand)
        if codon == "ATG":
            skipped["atg"] += 1
        elif codon not in want:
            skipped["not_a_start"] += 1
        else:
            key = (iv.contig, iv.strand, x % n)
            found.setdefault(key, (codon, set()))[1].add(g.locus_tag)
    order = {c: i for i, c in enumerate(gm.replicons)}
    sites = tuple(
        StartSite(c, s, p, tuple(sorted(tags)), codon)  # type: ignore[arg-type]
        for (c, s, p), (codon, tags) in sorted(
            found.items(), key=lambda kv: (order[kv[0][0]], *kv[0][1:])
        )
    )
    return sites, skipped


def start_keys(sites: Sequence[StartSite]) -> frozenset[_Key]:
    """``(contig, strand, codon_start)`` of every site (the M7-01 initiator set)."""
    return frozenset((s.contig, s.strand, s.codon_start) for s in sites)


def start_frame_cells(index: SixFrameIndex, sites: Sequence[StartSite]) -> tuple[np.ndarray, ...]:
    """Cells of ``sites`` per frame id (rule 2, as ``readthrough.frame_cells``), origin copies
    included.

    Raises:
        ValueError: if a cell does not map back to its site's codon or does not hold the table's
            translation of the site's codon.
    """
    fid_of = {(f.contig, f.strand, f.offset): i for i, f in enumerate(index.frames)}
    out: list[list[int]] = [[] for _ in index.frames]
    for site in sites:
        rep = index.replicons[site.contig]
        ext = index.frames[fid_of[(site.contig, site.strand, 0)]].ext_nt
        letter = translate_text(site.codon, index.translation_table)
        want = int(encode_target(letter)[0]) if len(letter) == 1 else STOP
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
            ok = want < STOP and 0 <= p < len(aa) and int(aa[p]) == want
            if ok:
                iv = index.to_interval(fid, p, p + 1)
                ok = iv.start == site.codon_start and iv.strand == site.strand
            if not ok:
                raise ValueError(
                    f"start site {site} does not read {letter!r} in frame {fid} at {p}"
                )
            out[fid].append(p)
    return tuple(np.array(sorted(c), dtype=np.int64) for c in out)


def mirror_start_cells(
    kind: DecoyType, frames: tuple[Frame, ...], cells: tuple[np.ndarray, ...]
) -> tuple[np.ndarray, ...]:
    """Decoy cells of ``kind`` in ``frames`` (that decoy's frames): ``n - 1 - p`` (``reverse``) or
    ``p`` (``shuffle``; the residue the shuffle put there does not matter). Draws no randomness.

    Raises:
        ValueError: if the frame counts differ, a cell is outside its frame or ``kind`` is unknown.
    """
    if kind not in ("reverse", "shuffle"):
        raise ValueError(f"unknown decoy type {kind!r}")
    if len(frames) != len(cells):
        raise ValueError(f"{len(cells)} cell arrays for {len(frames)} {kind} decoy frames")
    out = []
    for fid, (f, c) in enumerate(zip(frames, cells, strict=True)):
        c = np.asarray(c, dtype=np.int64)
        if len(c) and (c.min() < 0 or c.max() >= len(f.aa)):
            raise ValueError(f"a start cell of frame {fid} lies outside the frame")
        out.append(np.sort(len(f.aa) - 1 - c) if kind == "reverse" else c.copy())
    return tuple(out)


def credited_window(kind: str, n: int, cell: int, k: int) -> tuple[int, int] | None:
    """Rule 3: ``[c - k + 1, c + 1)`` in a ``reverse`` decoy, else ``[p, p + k)``; None when the
    window leaves the frame of ``n`` residues."""
    a, b = (cell - k + 1, cell + 1) if kind == "reverse" else (cell, cell + k)
    return (a, b) if a >= 0 and b <= n else None


def usable_windows(
    kind: str,
    index: SixFrameIndex,
    frames: tuple[Frame, ...],
    cells: tuple[np.ndarray, ...],
    k: int,
) -> int:
    """Credited windows of length ``k`` (``3k`` within the replicon) whose residues other than
    the cell are all stop-free and unambiguous; ``kind`` is ``target``, ``reverse`` or
    ``shuffle``."""
    total = 0
    for f, cl in zip(frames, cells, strict=True):
        if 3 * k > index.replicons[f.contig].length:
            continue
        for c in (int(x) for x in cl):
            w = credited_window(kind, len(f.aa), c, k)
            if w is None:
                continue
            a, b = w
            rest = f.aa[a : b - 1] if kind == "reverse" else f.aa[a + 1 : b]
            total += bool(np.all(rest < STOP))
    return total


def initiator_rows(
    rows: pd.DataFrame,
    sites: Sequence[StartSite],
    lengths: Mapping[str, int],
    circular: Mapping[str, bool],
) -> list[bool]:
    """Per row, whether it is a credited initiator row (rule 5): ``search = target``,
    ``route = exact``, the raw peptide begins with ``M`` and the row's 5' codon (``start`` on
    ``+``, ``end - 3`` on ``-``, wrapped on a circular replicon) is a start site."""
    keys = start_keys(sites)
    if not keys or not len(rows):
        return [False] * len(rows)
    out: list[bool] = []
    for pep, search, route, contig, a, b, strand in zip(
        rows["peptide"], rows["search"], rows["route"], rows["contig"], rows["start"],
        rows["end"], rows["strand"], strict=True,
    ):  # fmt: skip
        if search != "target" or route != "exact" or not str(pep).startswith("M"):
            out.append(False)
            continue
        c = str(contig)
        x = int(a) if strand == "+" else int(b) - 3
        if circular.get(c, False):
            x %= int(lengths[c])
        out.append((c, str(strand), x) in keys)
    return out
