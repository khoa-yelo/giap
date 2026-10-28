"""``--annotation`` GFF3/GBFF: Bakta ``--regions`` input and ``annotation/id_map.tsv`` (M2-11).

GIPPA_SPEC.md §5.3 ("Existing annotation"): with a GFF3 or GBFF from another annotator, Bakta runs
with ``--regions`` on that file's CDS so the gene coordinates match, and ``annotation/id_map.tsv``
maps the source IDs to Bakta locus tags by identical CDS coordinates. §7.1.1: ``--regions`` takes
CDS only. Intervals are 0-based half-open internally; origin-crossing CDS on circular replicons are
stored unwrapped (``end > length``) and written once (§6.5).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from pathlib import Path
from typing import Literal, cast

from gippa.io._open import open_binary, open_text, open_text_write
from gippa.io.gff3 import escape_attr, read_gff3
from gippa.model import Gene, Interval, Replicon
from gippa.model.coords import Strand, check_interval, interval_sort_key

AnnotationFormat = Literal["bakta_json", "gff3", "gbff"]
MatchKind = Literal["identical", "unmatched"]

ID_MAP_COLUMNS = ("source_id", "locus_tag", "contig", "start", "end", "strand", "match")
_SNIFF_BYTES = 4096
# Multi-row CDS parts may overlap or be at most this far apart (frameshifts, origin split).
_MAX_PART_GAP = 3


@dataclass(frozen=True, slots=True)
class UserCds:
    """One CDS of the user's annotation, on the safe contig ID, 0-based half-open."""

    source_id: str
    interval: Interval
    pseudo: bool


@dataclass(frozen=True, slots=True)
class IdMapRow:
    """One ``id_map.tsv`` row; coordinates are 0-based half-open (the writer converts)."""

    source_id: str
    locus_tag: str | None
    contig: str
    start: int
    end: int
    strand: str
    match: MatchKind


# --------------------------------------------------------------------------- format detection
def detect_format(path: Path) -> AnnotationFormat:
    """Detect the annotation format by content (gzip aware): JSON, GFF3 or GenBank flat file."""
    with open_binary(path) as fh:
        head = fh.read(_SNIFF_BYTES)
    text = head.decode("utf-8", errors="replace").lstrip("﻿ \t\r\n")
    if text.startswith("{"):
        return "bakta_json"
    if text.startswith("##gff-version"):
        return "gff3"
    if text.startswith("LOCUS"):
        return "gbff"
    raise ValueError(f"{path}: unrecognised annotation format (expected Bakta JSON, GFF3 or GBFF)")


# --------------------------------------------------------------------------- shared helpers
def _resolve_contig(
    name: str, replicons: Mapping[str, Replicon], contig_map: Mapping[str, str]
) -> Replicon:
    """Translate an original contig name to its replicon (via ``contig_map``) or raise."""
    safe = contig_map.get(name)
    if safe is None and name in replicons:
        safe = name
    rep = replicons.get(safe) if safe is not None else None
    if rep is None:
        raise ValueError(f"annotation contig {name} not in genome")
    return rep


def _fallback_id(contig: str, start: int, strand: str) -> str:
    """``<contig>_<start1>_<strand>`` for a CDS without any identifier."""
    return f"{contig}_{start + 1}_{strand}"


def _merge_parts(
    parts: Sequence[tuple[int, int]], rep: Replicon, where: str, max_gap: int | None = None
) -> tuple[int, int]:
    """Merge CDS parts (any order) into one span; parts split at the origin are unwrapped.

    On a circular replicon the span is cut at the largest gap between parts, so rows split at
    the origin (``[s, length)`` and ``[0, e)``) become ``[s, length + e)``. With ``max_gap``,
    consecutive parts (after unwrapping) must overlap or leave a gap of at most ``max_gap`` nt,
    else ``ValueError``; on a linear replicon this also rejects rows split at the origin.
    """
    ordered = sorted(parts)
    if rep.topology == "circular" and len(ordered) > 1:
        best_gap = ordered[0][0] + rep.length - max(e for _, e in ordered)
        cut = 0
        reach = ordered[0][1]
        for i in range(1, len(ordered)):
            gap = ordered[i][0] - reach
            if gap > best_gap:
                best_gap, cut = gap, i
            reach = max(reach, ordered[i][1])
        if cut:
            ordered = ordered[cut:] + [(s + rep.length, e + rep.length) for s, e in ordered[:cut]]
    if max_gap is not None:
        reach = ordered[0][1]
        for s, e in ordered[1:]:
            if s - reach > max_gap:
                raise ValueError(
                    f"{where}: parts {reach}-{s} apart by {s - reach} nt (more than {max_gap} nt"
                    f" on {rep.topology} contig {rep.id}); not one CDS"
                )
            reach = max(reach, e)
    start, end = min(s for s, _ in ordered), max(e for _, e in ordered)
    if start >= end:
        raise ValueError(f"{where}: empty CDS span {start}-{end}")
    return start, end


def _make_cds(
    source_id: str, rep: Replicon, start: int, end: int, strand: str, pseudo: bool, where: str
) -> UserCds:
    iv = Interval(rep.id, start, end, cast("Strand", strand))
    try:
        check_interval(iv, rep)
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    return UserCds(source_id, iv, pseudo)


def _sorted_cds(cds: Iterable[UserCds], order: Mapping[str, int]) -> list[UserCds]:
    return sorted(cds, key=lambda c: interval_sort_key(c.interval, order, c.source_id))


# --------------------------------------------------------------------------- GFF3
def _gff3_seqids(path: Path) -> list[str]:
    """Seqids of all feature lines (before ``##FASTA``), in file order."""
    seen: dict[str, None] = {}
    with open_text(path, newline="") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith(("##FASTA", ">")):
                break
            if line.strip() and not line.startswith("#"):
                seen.setdefault(line.split("\t", 1)[0], None)
    return list(seen)


def _attr(attrs: Mapping[str, object], key: str) -> str | None:
    value = attrs.get(key)
    if value is None:
        return None
    first = value if isinstance(value, str) else next(iter(cast("Sequence[str]", value)), "")
    return first or None


def _gff3_pseudo(attrs: Mapping[str, object]) -> bool:
    flag = _attr(attrs, "pseudo")
    return (flag is not None and flag.lower() == "true") or "pseudogene" in attrs


def _read_gff3_cds(
    path: Path, replicons: Mapping[str, Replicon], contig_map: Mapping[str, str]
) -> list[UserCds]:
    # The M1-05 reader validates coordinates against replicons keyed (and named) by the file's
    # own seqids, so build that view first; translation to safe IDs happens afterwards.
    by_name = {}
    for name in _gff3_seqids(path):
        rep = _resolve_contig(name, replicons, contig_map)
        by_name[name] = Replicon(name, rep.length, rep.topology, rep.kind)
    groups: dict[tuple[str, str], list[tuple[int, int, str, Mapping[str, object]]]] = {}
    singles: list[tuple[str, int, int, str, Mapping[str, object]]] = []
    for f in read_gff3(path, by_name):
        if f.type != "CDS":
            continue
        iv = f.interval
        if iv.strand not in ("+", "-"):
            raise ValueError(f"{path}: CDS at {iv.contig}:{iv.start + 1} has no strand")
        fid = _attr(f.attributes, "ID")
        if fid is None:
            singles.append((iv.contig, iv.start, iv.end, iv.strand, f.attributes))
        else:
            groups.setdefault((iv.contig, fid), []).append(
                (iv.start, iv.end, iv.strand, f.attributes)
            )
    out: list[UserCds] = []
    for name, start, end, strand, attrs in singles:
        rep = _resolve_contig(name, replicons, contig_map)
        sid = (
            _attr(attrs, "locus_tag") or _attr(attrs, "Name") or _fallback_id(rep.id, start, strand)
        )
        out.append(_make_cds(sid, rep, start, end, strand, _gff3_pseudo(attrs), str(path)))
    for (name, fid), rows in groups.items():
        rep = _resolve_contig(name, replicons, contig_map)
        strands = {r[2] for r in rows}
        if len(strands) != 1:
            raise ValueError(f"{path}: CDS {fid} has rows on both strands")
        start, end = _merge_parts(
            [(r[0], r[1]) for r in rows], rep, f"{path}: CDS {fid}", _MAX_PART_GAP
        )
        attrs = rows[0][3]
        sid = _attr(attrs, "locus_tag") or fid
        pseudo = any(_gff3_pseudo(r[3]) for r in rows)
        out.append(_make_cds(sid, rep, start, end, rows[0][2], pseudo, f"{path}: CDS {fid}"))
    return out


# --------------------------------------------------------------------------- GBFF
def _qual(quals: Mapping[str, Sequence[str]], key: str) -> str | None:
    values = quals.get(key)
    if not values:
        return None
    return str(values[0]) or None


def _gbff_span(
    parts: Sequence[tuple[int, int]], strand: str, rep: Replicon, where: str
) -> tuple[int, int]:
    """Span of a (possibly joined) location; parts are in biological order (INSDC)."""
    if len(parts) == 1:
        return parts[0]
    step = 1 if strand == "+" else -1
    wraps = any(step * (b[0] - a[0]) < 0 for a, b in pairwise(parts))
    if wraps and rep.topology != "circular":
        raise ValueError(f"{where}: join across the origin on linear contig {rep.id}")
    return _merge_parts(parts, rep, where)


def _read_gbff_cds(
    path: Path, replicons: Mapping[str, Replicon], contig_map: Mapping[str, str]
) -> list[UserCds]:
    from Bio import SeqIO  # BSD-3, in-process (spec: allowed dependency)

    out: list[UserCds] = []
    with open_text(path) as fh:
        for rec in SeqIO.parse(fh, "genbank"):  # type: ignore[no-untyped-call]
            name = rec.name if rec.name in contig_map or rec.name in replicons else rec.id
            rep = _resolve_contig(str(name), replicons, contig_map)
            for feat in rec.features:
                if feat.type != "CDS" or feat.location is None:
                    continue
                loc = feat.location
                where = f"{path}: {rep.id} CDS {loc}"
                if loc.strand not in (1, -1):
                    raise ValueError(f"{where}: CDS has no strand")
                strand = "+" if loc.strand == 1 else "-"
                parts = [(int(p.start), int(p.end)) for p in loc.parts]
                start, end = _gbff_span(parts, strand, rep, where)
                q = feat.qualifiers
                sid = (
                    _qual(q, "locus_tag")
                    or _qual(q, "protein_id")
                    or _qual(q, "gene")
                    or _fallback_id(rep.id, start, strand)
                )
                pseudo = "pseudo" in q or "pseudogene" in q
                out.append(_make_cds(sid, rep, start, end, strand, pseudo, where))
    return out


# --------------------------------------------------------------------------- public readers
def read_user_cds(
    path: Path,
    fmt: AnnotationFormat,
    replicons: Mapping[str, Replicon],
    contig_map: Mapping[str, str],
) -> list[UserCds]:
    """Read the CDS of a user GFF3/GBFF (gzip aware) on safe contig IDs, sorted per §6.5.

    ``contig_map`` maps original contig names to safe IDs; a contig that is not in the genome
    raises ``ValueError("annotation contig <name> not in genome")``.
    """
    if fmt == "gff3":
        cds = _read_gff3_cds(path, replicons, contig_map)
    elif fmt == "gbff":
        cds = _read_gbff_cds(path, replicons, contig_map)
    else:
        raise ValueError(f"read_user_cds: format {fmt!r} is used as is (Bakta is not run)")
    order = {rid: i for i, rid in enumerate(replicons)}
    return _sorted_cds(cds, order)


# --------------------------------------------------------------------------- writers
def write_bakta_regions(path: Path, cds: Sequence[UserCds], replicons: Sequence[Replicon]) -> None:
    """Write Bakta's ``--regions`` GFF3: ``##sequence-region`` pragmas and CDS lines only.

    Origin-crossing CDS are written once, unwrapped (``end > length``).
    """
    by_id = {r.id: r for r in replicons}
    order = {r.id: i for i, r in enumerate(replicons)}
    lines = ["##gff-version 3"]
    lines.extend(f"##sequence-region {r.id} 1 {r.length}" for r in replicons)
    for c in sorted(cds, key=lambda c: interval_sort_key(c.interval, order, c.source_id)):
        iv = c.interval
        check_interval(iv, by_id[iv.contig])
        cols = [iv.contig, "gippa", "CDS", str(iv.start + 1), str(iv.end), ".", iv.strand, "0"]
        lines.append("\t".join([*cols, f"ID={escape_attr(c.source_id)}"]))
    with open_text_write(path) as fh:
        fh.write("\n".join(lines) + "\n")


def build_id_map(
    cds: Sequence[UserCds],
    genes: Sequence[Gene],
    contig_order: Sequence[str],
) -> list[IdMapRow]:
    """Map each user CDS to the gene with identical ``(contig, start, end, strand)``.

    One row per user CDS; without a match ``locus_tag`` is ``None`` and ``match`` is
    ``unmatched``. Rows are sorted per §6.5; ``contig_order`` lists the contig IDs in input
    FASTA order (required, no alphabetical fallback).
    """
    index: dict[tuple[str, int, int, str], str] = {}
    for g in sorted(genes, key=lambda g: (g.feature_type != "cds", g.locus_tag)):
        iv = g.interval
        index.setdefault((iv.contig, iv.start, iv.end, iv.strand), g.locus_tag)
    order = {cid: i for i, cid in enumerate(contig_order)}
    rows = []
    for c in sorted(cds, key=lambda c: interval_sort_key(c.interval, order, c.source_id)):
        iv = c.interval
        tag = index.get((iv.contig, iv.start, iv.end, iv.strand))
        match: MatchKind = "unmatched" if tag is None else "identical"
        rows.append(IdMapRow(c.source_id, tag, iv.contig, iv.start, iv.end, iv.strand, match))
    return rows


def write_id_map(path: Path, rows: Sequence[IdMapRow]) -> None:
    """Write ``annotation/id_map.tsv`` (1-based closed, ``NA`` for a missing locus tag)."""
    lines = ["\t".join(ID_MAP_COLUMNS)]
    for r in rows:
        tag = "NA" if r.locus_tag is None else r.locus_tag
        cells = [r.source_id, tag, r.contig, str(r.start + 1), str(r.end), r.strand, r.match]
        lines.append("\t".join(cells))
    with open_text_write(path) as fh:
        fh.write("\n".join(lines) + "\n")
