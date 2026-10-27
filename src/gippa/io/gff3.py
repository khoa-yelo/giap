"""GFF3 writer and minimal reader (GIPPA_SPEC.md §6.5, §9.5, decision 16; ADR-0005 item 1).

Columns are 1-based closed (``start + 1``, ``end``); origin-crossing features on circular
replicons keep ``end > length``. Every replicon gets a ``##sequence-region`` pragma and circular
ones also a ``region`` feature with ``Is_circular=true``, which :func:`read_gff3` skips.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import cast

from gippa.io._open import open_text, open_text_write
from gippa.io.tsv import format_value
from gippa.model import Interval, Replicon, Strand, check_interval, interval_sort_key

AttrValue = str | Sequence[str]

_RESERVED = frozenset(";=&,%")
_FASTA_WIDTH = 80
_SAFE_ID = re.compile(r"[A-Za-z0-9._-]+")
_DIGITS = re.compile(r"[0-9]+")
_ESCAPE = re.compile(r"%([0-9A-Fa-f]{2})")


@dataclass(frozen=True, slots=True)
class GffFeature:
    """One GFF3 feature line; ``attributes`` keep insertion order."""

    interval: Interval
    source: str
    type: str
    score: float | None = None
    phase: int | None = None
    attributes: Mapping[str, AttrValue] = field(default_factory=dict)


def escape_attr(value: str) -> str:
    """Percent-encode ``; = & , %``, TAB, LF, CR and other control characters (GFF3 column 9)."""
    return "".join(
        f"%{ord(c):02X}" if c in _RESERVED or ord(c) < 0x20 or ord(c) == 0x7F else c for c in value
    )


def _check_field(name: str, value: str) -> str:
    if not value or any(ord(c) < 0x20 or ord(c) == 0x7F for c in value):
        raise ValueError(f"GFF3 {name} must be non-empty without control characters: {value!r}")
    return value


def _check_seqid(seqid: str) -> str:
    if not _SAFE_ID.fullmatch(seqid):
        raise ValueError(f"GFF3 seqid must match [A-Za-z0-9._-]+, got {seqid!r}")
    return seqid


def unescape_attr(text: str) -> str:
    """Decode ``%XX`` escapes as strict UTF-8; a ``%`` not followed by two hex digits raises."""
    out = bytearray()
    pos = 0
    for m in _ESCAPE.finditer(text):
        out += text[pos : m.start()].encode()
        out.append(int(m.group(1), 16))
        pos = m.end()
    out += text[pos:].encode()
    if "%" in _ESCAPE.sub("", text):
        raise ValueError(f"invalid percent-escape in {text!r}")
    try:
        return out.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise ValueError(f"percent-escapes in {text!r} are not valid UTF-8") from exc


def _attr_column(attrs: Mapping[str, AttrValue]) -> str:
    parts = []
    for key, value in attrs.items():
        if not key or escape_attr(key) != key:
            raise ValueError(f"invalid GFF3 attribute key {key!r}")
        items = [value] if isinstance(value, str) else list(value)
        if not items or not all(isinstance(x, str) and x for x in items):
            raise ValueError(f"GFF3 attribute {key}: values must be non-empty strings: {value!r}")
        parts.append(f"{key}=" + ",".join(escape_attr(x) for x in items))
    return ";".join(parts) if parts else "."


def _feature_line(f: GffFeature, replicon: Replicon) -> str:
    check_interval(f.interval, replicon)
    if f.score is not None and not math.isfinite(f.score):
        raise ValueError(f"GFF3 score must be finite, got {f.score!r}")
    if f.phase is not None and f.phase not in (0, 1, 2):
        raise ValueError(f"GFF3 phase must be 0, 1 or 2, got {f.phase!r}")
    if f.type == "CDS" and f.phase is None:
        raise ValueError(f"{f.interval.contig}:{f.interval.start}: CDS needs a phase (GFF3 spec)")
    iv = f.interval
    return "\t".join(
        [
            iv.contig,
            _check_field("source", f.source),
            _check_field("type", f.type),
            str(iv.start + 1),
            str(iv.end),
            "." if f.score is None else format_value(f.score),
            iv.strand,
            "." if f.phase is None else str(f.phase),
            _attr_column(f.attributes),
        ]
    )


def _region_line(r: Replicon) -> str:
    return f"{r.id}\tgippa\tregion\t1\t{r.length}\t.\t.\t.\tID={r.id};Is_circular=true"


def write_gff3(
    path: Path,
    features: Iterable[GffFeature],
    replicons: Sequence[Replicon],
    fasta: Mapping[str, str] | None = None,
) -> None:
    """Write sorted GFF3 (§6.5 order) with pragmas, region lines and an optional ``##FASTA``."""
    by_id = {_check_seqid(r.id): r for r in replicons}
    order = {r.id: i for i, r in enumerate(replicons)}
    keyed = []
    for f in features:
        rep = by_id.get(f.interval.contig)
        if rep is None:
            raise ValueError(f"feature on unknown contig {f.interval.contig!r}")
        line = _feature_line(f, rep)
        key = interval_sort_key(f.interval, order, str(f.attributes.get("ID", "")))
        keyed.append((key, line))
    keyed.sort()
    lines = ["##gff-version 3"]
    for r in replicons:
        lines.append(f"##sequence-region {r.id} 1 {r.length}")
        if r.topology == "circular":
            lines.append(_region_line(r))
    lines.extend(line for _, line in keyed)
    if fasta is not None:
        lines.append("##FASTA")
        for rid, seq in fasta.items():
            lines.append(f">{rid}")
            lines.extend(seq[i : i + _FASTA_WIDTH] for i in range(0, len(seq), _FASTA_WIDTH))
    with open_text_write(path) as fh:
        fh.write("\n".join(lines) + "\n")


def _parse_attrs(col: str, where: str) -> dict[str, AttrValue]:
    attrs: dict[str, AttrValue] = {}
    if col in ("", "."):
        return attrs
    for part in col.split(";"):
        if not part:
            continue
        key, sep, raw = part.partition("=")
        if not sep or not key:
            raise ValueError(f"{where}: malformed attribute {part!r}")
        try:
            values = [unescape_attr(v) for v in raw.split(",")]
            name = unescape_attr(key)
        except ValueError as exc:
            raise ValueError(f"{where}: {exc}") from exc
        if name in attrs:
            raise ValueError(f"{where}: duplicate attribute key {name!r}")
        attrs[name] = values[0] if len(values) == 1 else tuple(values)
    return attrs


def _parse_line(line: str, replicons: Mapping[str, Replicon], where: str) -> GffFeature:
    cols = line.split("\t")
    if len(cols) != 9:
        raise ValueError(f"{where}: expected 9 tab-separated columns, got {len(cols)}")
    seqid, source, ftype, start_s, end_s, score_s, strand, phase_s, attr_s = cols
    rep = replicons.get(seqid)
    if rep is None:
        raise ValueError(f"{where}: unknown seqid {seqid!r}")
    if not (_DIGITS.fullmatch(start_s) and _DIGITS.fullmatch(end_s)):
        raise ValueError(f"{where}: start/end must be plain digits, got {start_s!r}, {end_s!r}")
    try:
        start, end = int(start_s), int(end_s)
        score = None if score_s == "." else float(score_s)
        phase = None if phase_s == "." else int(phase_s)
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    if start < 1 or start > end:
        raise ValueError(f"{where}: need 1 <= start <= end, got {start}-{end}")
    if strand == "?":
        strand = "."
    if strand not in ("+", "-", "."):
        raise ValueError(f"{where}: invalid strand {strand!r}")
    if phase is not None and phase not in (0, 1, 2):
        raise ValueError(f"{where}: phase must be 0, 1 or 2, got {phase}")
    iv = Interval(seqid, start - 1, end, cast("Strand", strand))
    try:
        check_interval(iv, rep)
    except ValueError as exc:
        raise ValueError(f"{where}: {exc}") from exc
    return GffFeature(iv, source, ftype, score, phase, _parse_attrs(attr_s, where))


def read_gff3(path: Path, replicons: Mapping[str, Replicon]) -> list[GffFeature]:
    """Read GFF3 feature lines (attributes unescaped); stops at ``##FASTA``.

    Comments, other pragmas, blank lines and GIPPA's own region lines are skipped.
    """
    with open_text(path, newline="") as fh:
        text = fh.read()
    region_lines = {_region_line(r) for r in replicons.values() if r.topology == "circular"}
    out: list[GffFeature] = []
    for lineno, raw in enumerate(text.split("\n"), start=1):
        line = raw.rstrip("\r")
        if line.startswith(("##FASTA", ">")):
            break
        if not line.strip() or line.startswith("#") or line in region_lines:
            continue
        out.append(_parse_line(line, replicons, f"{path}: line {lineno}"))
    return out
