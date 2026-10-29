"""Toy genomes and an independent model of the initiator-M credit rule (docs/specs/
M11-20-initiator-methionine.md § Definition rules 1-3, § Examples). Test tooling only: nothing here
imports ``gippa.mapping.initiator``, so the helpers build inputs before the implementation exists.

- ``c1`` (example 1): circular 120 nt, ``C x 12`` + G1 ``+ [12, 39)`` = ``GTG GCT GAA CTG TGG CAT
  AAA GGC TAA`` (aa ``MAELWHKG``) + ``C x 81``; frame 0 = ``PPPPVAELWHKG*PPP...`` twice (80 aa).
- ``c2`` (example 6): linear 60 nt, G2 ``- [6, 30)`` reading ``TTG AAA CAG TGG TAT TTT CGT TAG``.
- ``c3`` (example 7): circular 90 nt, G3 ``+ [87, 111)`` ``TTG GCT GAA CTG TGG CAT AAA TAA``.
- ``c4`` (example 7b): circular 90 nt, G4 ``+ [89, 113)``, the GTG at bases 89, 0, 1.
- ``c5`` (example 11): linear 60 nt, ``GTG GCT GAA TGA CTG CAT AAA GGC TAA`` at ``[6, 33)``.

The model: :func:`ref_cells` finds every frame position whose one-residue window maps
(``to_interval``) onto a start key, origin copies included; :func:`credited_rows` is the brute force
of rule 3 (target and ``shuffle``: cell ``p``, window ``[p, p + k)``; ``reverse``: cell ``c``,
window ``[c - k + 1, c + 1)`` mirrored onto the target); :func:`usable` counts credited windows.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import replace

import numpy as np

from gippa.mapping.decoys import build_decoys, mirror_aa
from gippa.mapping.sixframe import STOP, Frame, SixFrameIndex, build_sixframe, decode, encode_query
from gippa.model import Replicon
from gippa.model.genes import GeneModel
from gippa.rng import for_step
from helpers.readthrough_toys import gene, model, rc

Key = tuple[str, str, int]  # (contig, strand, codon_start)
Row = tuple[str, str, str | None, str, int, int, str]  # peptide, search, decoy_type, iv
STARTS = ("ATG", "GTG", "TTG")
G1_NT = "GTGGCTGAACTGTGGCATAAAGGCTAA"  # MAELWHKG*
G2_NT = "TTGAAACAGTGGTATTTTCGTTAG"  # MKQWYFR* (gene strand)
G3_NT = "TTGGCTGAACTGTGGCATAAATAA"  # MAELWHK*
G4_NT = "GTGGCTGAACTGTGGCATAAATAA"  # MAELWHK*
G5_NT = "GTGGCTGAATGACTGCATAAAGGCTAA"  # M A E (TGA) L H K G *


def cds(
    tag: str,
    contig: str,
    start: int,
    end: int,
    strand: str,
    aa: str,
    codon: str,
    feature_type: str = "cds",
):  # type: ignore[no-untyped-def]
    """A Gene with its Bakta ``aa`` and ``start_codon``."""
    return replace(gene(tag, contig, start, end, strand, feature_type), aa=aa, start_codon=codon)


def toy_c1(first: str = "GTG", copy_at_60: bool = False):  # type: ignore[no-untyped-def]
    """Example 1 (``first`` = G1's first codon; example 3: ``ATG``); example 2 with the copy."""
    seq = list("C" * 12 + first + G1_NT[3:] + "C" * 81)
    if copy_at_60:  # plus a frame-0 TAA at [90, 93), so the copy is not G1's N-terminal extension
        seq[60:81], seq[90:93] = "ATGGCTGAACTGTGGCATAAA", "TAA"
    reps = {"c1": Replicon("c1", 120, "circular")}
    gm = model(reps, [cds("G1", "c1", 12, 39, "+", "MAELWHKG", first)])
    return {"c1": "".join(seq)}, reps, gm


def toy_c2():  # type: ignore[no-untyped-def]
    nt = "C" * 6 + rc(G2_NT) + "C" * 30
    assert nt[27:30] == "CAA"
    reps = {"c2": Replicon("c2", 60, "linear")}
    return {"c2": nt}, reps, model(reps, [cds("G2", "c2", 6, 30, "-", "MKQWYFR", "TTG")])


def toy_c3():  # type: ignore[no-untyped-def]
    seq = list("C" * 90)
    for i, b in enumerate(G3_NT):
        seq[(87 + i) % 90] = b
    reps = {"c3": Replicon("c3", 90, "circular")}
    return (
        {"c3": "".join(seq)},
        reps,
        model(reps, [cds("G3", "c3", 87, 111, "+", "MAELWHK", "TTG")]),
    )


def toy_c4():  # type: ignore[no-untyped-def]
    seq = list("C" * 90)
    for i, b in enumerate(G4_NT):
        seq[(89 + i) % 90] = b
    assert seq[89] + seq[0] + seq[1] == "GTG"
    reps = {"c4": Replicon("c4", 90, "circular")}
    return (
        {"c4": "".join(seq)},
        reps,
        model(reps, [cds("G4", "c4", 89, 113, "+", "MAELWHK", "GTG")]),
    )


def toy_c5(table: int):  # type: ignore[no-untyped-def]
    aa = {4: "MAEWLHKG", 25: "MAEGLHKG"}[table]
    reps = {"c5": Replicon("c5", 60, "linear")}
    nt = "C" * 6 + G5_NT + "C" * 27
    return {"c5": nt}, reps, model(reps, [cds("G5", "c5", 6, 33, "+", aa, "GTG")])


def with_decoys(
    genome: Mapping[str, str], reps: Mapping[str, Replicon], table: int = 11, seed: int = 7
) -> tuple[SixFrameIndex, dict[str, tuple[Frame, ...]]]:
    """The index and both decoy types (the step's ``for_step("sixframe_index", seed)`` stream)."""
    index = build_sixframe(genome, reps, table=table)
    rng = for_step("sixframe_index", seed)
    return index, {k: build_decoys(index, k, 50, rng) for k in ("reverse", "shuffle")}  # type: ignore[arg-type]


def start_keys(
    gm: GeneModel, genome: Mapping[str, str], codons: Iterable[str] = STARTS, table: int = 11
) -> set[Key]:
    """Rule 1 from the gene model, independently: cds/sorf, ``aa`` starting with M, the gene-strand
    first codon (wrapped) in ``codons`` and not ATG, inside a linear contig."""
    del table  # rule 2 reads codons, not residues
    want = {c.upper() for c in codons} - {"ATG"}
    out: set[Key] = set()
    for g in gm.genes:
        if g.feature_type not in ("cds", "sorf") or not g.aa or g.aa[0] != "M":
            continue
        iv, rep = g.interval, gm.replicons[g.interval.contig]
        x = iv.start if iv.strand == "+" else iv.end - 3
        if rep.topology != "circular" and not 0 <= x <= rep.length - 3:
            continue
        seq = genome[iv.contig]
        nt = "".join(seq[(x + i) % rep.length] for i in range(3)).upper()
        nt = rc(nt) if iv.strand == "-" else nt
        if nt in want:
            out.add((iv.contig, iv.strand, x % rep.length))
    return out


def ref_cells(index: SixFrameIndex, keys: Iterable[Key]) -> tuple[np.ndarray, ...]:
    """Every frame position whose residue maps onto a start key (origin copies included)."""
    ks = set(keys)
    out = []
    for fid, f in enumerate(index.frames):
        cells = [p for p in range(len(f.aa)) if _key(index, fid, p) in ks] if ks else []
        out.append(np.array(cells, dtype=np.int64))
    return tuple(out)


def _key(index: SixFrameIndex, fid: int, p: int) -> Key:
    iv = index.to_interval(fid, p, p + 1)
    return (iv.contig, iv.strand, iv.start)


def mirrored(
    kind: str, frames: Sequence[Frame], cells: Sequence[np.ndarray]
) -> tuple[np.ndarray, ...]:
    if kind == "shuffle":
        return tuple(np.asarray(c, dtype=np.int64) for c in cells)
    return tuple(
        np.sort(len(f.aa) - 1 - np.asarray(c, dtype=np.int64))
        for f, c in zip(frames, cells, strict=True)
    )


def windows(kind: str, n: int, cells: np.ndarray, k: int) -> list[tuple[int, int, int]]:
    """(cell, aa_start, aa_end) of the credited windows of one frame (rule 3), in range only."""
    out = []
    for c in (int(x) for x in cells):
        a, b = (c - k + 1, c + 1) if kind == "reverse" else (c, c + k)
        if a >= 0 and b <= n:
            out.append((c, a, b))
    return out


def rest_of(kind: str, aa: np.ndarray, a: int, b: int) -> np.ndarray:
    """The frame residues the query's positions 1..k-1 must equal (in query order)."""
    return aa[a : b - 1] if kind == "reverse" else aa[a + 1 : b]


def usable(kind: str, frame: Frame, cells: np.ndarray, k: int, rep_len: int) -> int:
    if 3 * k > rep_len:
        return 0
    return sum(
        bool(np.all(rest_of(kind, frame.aa, a, b) < STOP))
        for _, a, b in windows(kind, len(frame.aa), cells, k)
    )


def credited_rows(
    queries: Iterable[str],
    index: SixFrameIndex,
    decoys: Mapping[str, Sequence[Frame]],
    decoy_of: Mapping[str, str],
    cells: Mapping[str, Sequence[np.ndarray]],
    *,
    min_len: int = 7,
    max_len: int = 40,
) -> set[Row]:
    """Brute force of rule 3: credited rows ``(peptide, search, decoy_type, contig, start, end,
    strand)``; ``cells`` maps ``target`` / ``reverse`` / ``shuffle`` to per-frame cells."""
    out: set[Row] = set()
    for pep in set(queries):
        q = decode(encode_query(pep))
        k = len(q)
        if "#" in q or not min_len <= k <= max_len or q[0] != "M":
            continue
        spaces: list[tuple[str, str | None, Sequence[Frame]]] = [("target", None, index.frames)]
        if decoy_of[pep] in decoys:
            spaces.append(("decoy", decoy_of[pep], decoys[decoy_of[pep]]))
        for search, kind, frames in spaces:
            rule = kind or "target"
            for fid, f in enumerate(frames):
                if 3 * k > index.replicons[f.contig].length:
                    continue
                for _, a, b in windows(rule, len(f.aa), cells[rule][fid], k):
                    if decode(rest_of(rule, f.aa, a, b)) != q[1:]:
                        continue
                    ta, tb = (a, b) if kind is None else mirror_aa(kind, len(f.aa), a, b)  # type: ignore[arg-type]
                    iv = index.to_interval(fid, ta, tb)
                    out.add((pep, search, kind, iv.contig, iv.start, iv.end, iv.strand))
    return out


def row_set(df) -> set[Row]:  # type: ignore[no-untyped-def]
    return {
        (
            str(r.peptide),
            str(r.search),
            None if r.decoy_type is None or r.decoy_type != r.decoy_type else str(r.decoy_type),
            str(r.contig),
            int(r.start),
            int(r.end),
            str(r.strand),
        )
        for r in df.itertuples(index=False)
    }
