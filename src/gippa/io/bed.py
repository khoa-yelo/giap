"""BED6/BED12 writer (GIPPA_SPEC.md §6.5, decision 16; ADR-0005 item 1).

BED is 0-based half-open like the internal model, so coordinates are written unchanged. An
origin-crossing feature on a circular replicon becomes two records ``<name>_part1`` (``[start,
length)``) and ``<name>_part2`` (``[0, end - length)``); blocks and thick range are split with it.
Records are validated, never clamped (ADR-0006 item 3).
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from gippa.io._open import open_text_write
from gippa.model import Interval, Replicon, check_interval, interval_sort_key, split_at_origin

Block = tuple[int, int]  # (size, start relative to the record start)

_WHITESPACE = re.compile(r"\s")


@dataclass(frozen=True, slots=True)
class BedRecord:
    """One BED feature; ``blocks == ()`` writes BED6, otherwise BED12."""

    interval: Interval
    name: str
    score: int = 0
    thick: Interval | None = None
    rgb: str = "0"
    blocks: tuple[Block, ...] = ()


def _validate(rec: BedRecord, replicon: Replicon) -> None:
    check_interval(rec.interval, replicon)
    if _WHITESPACE.search(rec.name):
        raise ValueError(f"BED name may not contain whitespace: {rec.name!r}")
    if type(rec.score) is not int or not 0 <= rec.score <= 1000:
        raise ValueError(f"BED score must be an int in [0, 1000], got {rec.score!r} ({rec.name})")
    if not rec.blocks:
        if rec.thick is not None or rec.rgb != "0":
            raise ValueError(f"BED6 record {rec.name} may not set thick or rgb (give blocks)")
        return
    pos = 0  # blocks must tile [0, length) in order without overlap (gaps are introns)
    for i, (size, rel) in enumerate(rec.blocks):
        if size < 1 or rel < pos or (i == 0 and rel != 0):
            raise ValueError(
                f"{rec.name}: BED12 blocks must start at 0, be sorted and non-overlapping, "
                f"got {rec.blocks}"
            )
        pos = rel + size
    if pos != rec.interval.length_bp:
        raise ValueError(
            f"{rec.name}: last BED12 block ends at {pos}, record length is {rec.interval.length_bp}"
        )


def _thick_range(rec: BedRecord, replicon: Replicon) -> tuple[int, int]:
    """Thick range in the record's unwrapped frame; it must lie inside the interval."""
    iv = rec.interval
    if rec.thick is None:
        return iv.start, iv.end
    t0, t1 = rec.thick.start, rec.thick.end
    if t0 < iv.start and iv.end > replicon.length:  # thick given in wrapped coordinates
        t0, t1 = t0 + replicon.length, t1 + replicon.length
    if not iv.start <= t0 <= t1 <= iv.end:
        raise ValueError(f"{rec.name}: thick range {rec.thick} lies outside {iv}")
    return t0, t1


def bed_parts(rec: BedRecord, replicon: Replicon) -> list[tuple[Interval, list[str]]]:
    """Return ``(part interval, fields)`` per output line: one, or two if origin-crossing."""
    _validate(rec, replicon)
    iv = rec.interval
    parts = split_at_origin(iv, replicon)
    shift = replicon.length
    t0, t1 = _thick_range(rec, replicon)
    out: list[tuple[Interval, list[str]]] = []
    for k, part in enumerate(parts):
        # Work in unwrapped coordinates; part 2 lies at [length, end) before re-basing.
        off = shift if k == 1 else 0
        lo, hi = part.start + off, part.end + off
        blocks: list[tuple[int, int]] = []  # absolute unwrapped [b0, b1)
        for size, rel in rec.blocks:
            b0, b1 = max(iv.start + rel, lo), min(iv.start + rel + size, hi)
            if b0 < b1:
                blocks.append((b0, b1))
        if rec.blocks:
            assert blocks, "validated blocks tile the record, so every part has one"
            lo, hi = blocks[0][0], blocks[-1][1]
        name = rec.name if len(parts) == 1 else f"{rec.name}_part{k + 1}"
        piv = Interval(iv.contig, lo - off, hi - off, iv.strand)
        fields = [iv.contig, str(piv.start), str(piv.end), name, str(rec.score), iv.strand]
        if rec.blocks:
            ts, te = max(t0, lo), min(t1, hi)
            if ts >= te:
                ts = te = lo
            fields += [
                str(ts - off),
                str(te - off),
                rec.rgb,
                str(len(blocks)),
                ",".join(str(b1 - b0) for b0, b1 in blocks),
                ",".join(str(b0 - lo) for b0, _ in blocks),
            ]
        out.append((piv, fields))
    return out


def bed_lines(rec: BedRecord, replicon: Replicon) -> list[str]:
    """BED lines (no newline) for ``rec``: 1 line, or 2 if it crosses the origin."""
    return ["\t".join(f) for _, f in bed_parts(rec, replicon)]


def sorted_lines(
    parts: Iterable[tuple[Interval, list[str]]], replicons: Mapping[str, Replicon]
) -> list[str]:
    """Sort split records by contig order, start, end, strand, name (§6.5) and join fields."""
    order = {cid: i for i, cid in enumerate(replicons)}
    keyed = sorted(parts, key=lambda p: interval_sort_key(p[0], order, p[1][3]))
    return ["\t".join(f) for _, f in keyed]


def write_bed(path: Path, records: Iterable[BedRecord], replicons: Mapping[str, Replicon]) -> None:
    """Write sorted BED records; origin-crossing records are split (``_part2`` sorts first)."""
    parts = [p for rec in records for p in bed_parts(rec, replicons[rec.interval.contig])]
    with open_text_write(path) as fh:
        fh.writelines(line + "\n" for line in sorted_lines(parts, replicons))
