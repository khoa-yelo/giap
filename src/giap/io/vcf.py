"""VCF 4.3 writer (GIAP_SPEC.md §6.5, §7.5.5, §9.5, Appendix B.7, G14; ADR-0005 item 1).

POS is ``wrap_position(pos0) + 1`` within ``[1, length]``; REF/ALT are forward-strand bases given
by the caller. A REF allele crossing the contig end raises ``ValueError``.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from giap.io._open import open_text_write
from giap.io.tsv import format_value
from giap.model import Replicon, wrap_position

InfoValue = str | int | float | Sequence[str] | Sequence[int] | Sequence[float]

_SHA256 = re.compile(r"[0-9a-f]{64}")
_INFO_RESERVED = frozenset(":;=%,")
_FIELD_TEXT = re.compile(r"[^\t\n\r]+")


@dataclass(frozen=True, slots=True)
class VcfRecord:
    """One VCF record at 0-based ``pos0`` (may be unwrapped on a circular contig)."""

    contig: str
    pos0: int
    ref: str
    alts: tuple[str, ...]
    info: Mapping[str, InfoValue]
    gt: str
    ad: tuple[int, ...]
    dp: int
    id: str = "."
    qual: str = "."
    filter: str = "PASS"


@dataclass(frozen=True, slots=True)
class VcfHeaderField:
    """An ``##INFO`` or ``##FORMAT`` definition."""

    id: str
    number: str
    type: str
    description: str


def _escape_info(value: str) -> str:
    return "".join(
        "".join(f"%{b:02X}" for b in c.encode()) if c in _INFO_RESERVED or c.isspace() else c
        for c in value
    )


def _info_item(v: object) -> str:
    return _escape_info(v) if isinstance(v, str) else format_value(v)


def _info_column(info: Mapping[str, InfoValue]) -> str:
    parts = []
    for key, value in info.items():
        if isinstance(value, str | int | float):
            text = _info_item(value)
        else:
            text = ("|" if key == "FLAGS" else ",").join(_info_item(x) for x in value)
        parts.append(f"{key}={text}")
    return ";".join(parts) if parts else "."


def _header_line(kind: str, f: VcfHeaderField) -> str:
    desc = f.description.replace("\\", "\\\\").replace('"', '\\"')
    return f'##{kind}=<ID={f.id},Number={f.number},Type={f.type},Description="{desc}">'


def _record_key(rec: VcfRecord, by_id: Mapping[str, tuple[int, Replicon]]) -> tuple[int, int, str]:
    if rec.contig not in by_id:
        raise ValueError(f"VCF record on unknown contig {rec.contig!r}")
    idx, rep = by_id[rec.contig]
    pos0 = wrap_position(rec.pos0, rep)
    if not rec.ref or pos0 + len(rec.ref) > rep.length:
        raise ValueError(
            f"{rec.contig}:{rec.pos0}: REF {rec.ref!r} must be non-empty and end within the "
            f"contig (length {rep.length}); origin-crossing alleles are not representable"
        )
    if len(rec.ad) != 1 + len(rec.alts):
        raise ValueError(
            f"{rec.contig}:{pos0 + 1}: AD has {len(rec.ad)} values, need {1 + len(rec.alts)}"
        )
    cols = [
        rec.contig,
        str(pos0 + 1),
        rec.id,
        rec.ref,
        ",".join(rec.alts) if rec.alts else ".",
        rec.qual,
        rec.filter,
        _info_column(rec.info),
        "GT:AD:DP",
        f"{rec.gt}:{','.join(str(a) for a in rec.ad)}:{rec.dp}",
    ]
    for c in cols:
        if not _FIELD_TEXT.fullmatch(c):
            raise ValueError(f"{rec.contig}:{pos0 + 1}: invalid VCF field {c!r}")
    return idx, pos0, "\t".join(cols)


def write_vcf(
    path: Path,
    records: Iterable[VcfRecord],
    replicons: Sequence[Replicon],
    *,
    giap_version: str,
    genome_sha256: str,
    sample: str,
    info_fields: Sequence[VcfHeaderField],
    format_fields: Sequence[VcfHeaderField],
) -> None:
    """Write a VCF 4.3 file; records are sorted by contig order, then POS."""
    if not _SHA256.fullmatch(genome_sha256):
        raise ValueError(f"genome_sha256 must be 64 lowercase hex digits, got {genome_sha256!r}")
    if not _FIELD_TEXT.fullmatch(sample):
        raise ValueError(f"invalid VCF sample name {sample!r}")
    by_id = {r.id: (i, r) for i, r in enumerate(replicons)}
    body = sorted(_record_key(rec, by_id) for rec in records)
    lines = [
        "##fileformat=VCFv4.3",
        f"##source=giap-{giap_version}",
        f"##reference=sha256:{genome_sha256}",
        *(f"##contig=<ID={r.id},length={r.length}>" for r in replicons),
        *(_header_line("INFO", f) for f in info_fields),
        *(_header_line("FORMAT", f) for f in format_fields),
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
        ),
        *(line for _, _, line in body),
    ]
    with open_text_write(path) as fh:
        fh.write("\n".join(lines) + "\n")
