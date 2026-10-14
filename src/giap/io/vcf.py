"""VCF 4.3 writer (GIAP_SPEC.md §6.5, §7.5.5, §9.5, Appendix B.7, G14; ADR-0005 item 1).

POS is ``wrap_position(pos0) + 1`` within ``[1, length]``; REF/ALT are forward-strand bases given
by the caller. A REF allele crossing the contig end raises ``ValueError``. Header and records are
validated strictly (ADR-0005 item 1, after the M1-05 review) so htslib reads them without warnings.
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from giap.io._open import open_text_write
from giap.io.tsv import format_value
from giap.model import Replicon, wrap_position

InfoValue = str | int | float | bool | Sequence[str] | Sequence[int] | Sequence[float]

_SHA256 = re.compile(r"[0-9a-f]{64}")
_INFO_RESERVED = frozenset(":;=%,")
_FIELD_TEXT = re.compile(r"[^\t\n\r]+")
_SAFE_ID = re.compile(r"[A-Za-z0-9._-]+")
_KEY = re.compile(r"[A-Za-z_][0-9A-Za-z_.]*")
_NUMBER = re.compile(r"[0-9]+|[ARG.]")
_TYPES = frozenset({"Integer", "Float", "Flag", "Character", "String"})
_BASES = re.compile(r"[ACGTN]+")
_GT = re.compile(r"\.|[1-9][0-9]*")
_ID_COL = re.compile(r"[^\s;]+")
_FILTER_ID = re.compile(r"[^\s;,=<>\"]+")


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


def _check_text(name: str, value: str) -> str:
    if any(ord(c) < 0x20 or ord(c) == 0x7F for c in value):
        raise ValueError(f"VCF {name} must not contain control characters: {value!r}")
    return value


def _escape_info(value: str) -> str:
    return "".join(
        "".join(f"%{b:02X}" for b in c.encode()) if c in _INFO_RESERVED or c.isspace() else c
        for c in value
    )


def _info_item(key: str, v: object, vtype: str) -> str:
    if isinstance(v, bool):
        raise ValueError(f"INFO {key}: bool is only allowed for Type=Flag")
    if isinstance(v, str):
        if vtype in ("Integer", "Float") or not v:
            raise ValueError(f"INFO {key}: invalid {vtype} value {v!r}")
        return _escape_info(v)
    if isinstance(v, float) and math.isnan(v):
        return "."
    if not isinstance(v, int | float) or (vtype == "Integer" and not isinstance(v, int)):
        raise ValueError(f"INFO {key}: invalid {vtype} value {v!r}")
    return format_value(v)


def _info_column(info: Mapping[str, InfoValue], declared: Mapping[str, VcfHeaderField]) -> str:
    parts = []
    for key, value in info.items():
        field = declared.get(key)
        if field is None:
            raise ValueError(f"INFO key {key!r} is not declared in the header")
        if field.type == "Flag" or isinstance(value, bool):
            if field.type != "Flag" or not isinstance(value, bool):
                raise ValueError(f"INFO {key}: Type=Flag needs a bool value, got {value!r}")
            if value:
                parts.append(key)
            continue
        if isinstance(value, str | int | float):
            text = _info_item(key, value, field.type)
        else:
            items = list(value)
            if not items:
                raise ValueError(f"INFO {key}: empty list")
            if key == "FLAGS" and any(isinstance(x, str) and "|" in x for x in items):
                raise ValueError(f"INFO FLAGS items must not contain '|': {items!r}")
            sep = "|" if key == "FLAGS" else ","
            text = sep.join(_info_item(key, x, field.type) for x in items)
        parts.append(f"{key}={text}")
    return ";".join(parts) if parts else "."


def _header_line(kind: str, f: VcfHeaderField) -> str:
    if not _KEY.fullmatch(f.id) or not _NUMBER.fullmatch(f.number) or f.type not in _TYPES:
        raise ValueError(f"invalid ##{kind} definition {f!r}")
    if kind == "FORMAT" and f.type == "Flag":
        raise ValueError(f"##FORMAT {f.id}: Type=Flag is not allowed")
    if f.type == "Flag" and f.number != "0":
        raise ValueError(f"##{kind} {f.id}: Type=Flag needs Number=0")
    _check_text(f"##{kind} {f.id} description", f.description)
    desc = f.description.replace("\\", "\\\\").replace('"', '\\"')
    return f'##{kind}=<ID={f.id},Number={f.number},Type={f.type},Description="{desc}">'


def _check_alleles(where: str, ref: str, alts: tuple[str, ...], gt: str) -> None:
    if not _BASES.fullmatch(ref):
        raise ValueError(f"{where}: REF must be uppercase [ACGTN]+, got {ref!r}")
    for alt in alts:
        if not _BASES.fullmatch(alt) or alt == ref or len(alt) != len(ref):
            raise ValueError(
                f"{where}: ALT {alt!r} must be uppercase [ACGTN]+, differ from REF {ref!r} "
                "and have the same length"
            )
    if not _GT.fullmatch(gt) or (gt != "." and int(gt) > len(alts)):
        raise ValueError(f"{where}: GT must be '.' or 1..{len(alts)}, got {gt!r}")


def _record_key(
    rec: VcfRecord,
    by_id: Mapping[str, tuple[int, Replicon]],
    info_decl: Mapping[str, VcfHeaderField],
    sequences: Mapping[str, str] | None,
) -> tuple[int, int, str]:
    if rec.contig not in by_id:
        raise ValueError(f"VCF record on unknown contig {rec.contig!r}")
    idx, rep = by_id[rec.contig]
    pos0 = wrap_position(rec.pos0, rep)
    where = f"{rec.contig}:{pos0 + 1}"
    if not rec.ref or pos0 + len(rec.ref) > rep.length:
        raise ValueError(
            f"{rec.contig}:{rec.pos0}: REF {rec.ref!r} must be non-empty and end within the "
            f"contig (length {rep.length}); origin-crossing alleles are not representable"
        )
    _check_alleles(where, rec.ref, rec.alts, rec.gt)
    if sequences is not None:
        seq = sequences.get(rec.contig)
        if seq is None:
            raise ValueError(f"{where}: no sequence given for contig {rec.contig!r}")
        found = seq[pos0 : pos0 + len(rec.ref)].upper()
        if found != rec.ref:
            raise ValueError(f"{where}: REF {rec.ref!r} does not match the genome ({found!r})")
    if len(rec.ad) != 1 + len(rec.alts):
        raise ValueError(f"{where}: AD has {len(rec.ad)} values, need {1 + len(rec.alts)}")
    if any(a < 0 for a in rec.ad) or rec.dp < 0:
        raise ValueError(f"{where}: AD and DP must be non-negative")
    if not _ID_COL.fullmatch(rec.id):
        raise ValueError(f"{where}: invalid ID {rec.id!r}")
    if rec.qual != ".":
        try:
            float(rec.qual)
        except ValueError as exc:
            raise ValueError(f"{where}: QUAL must be '.' or a number, got {rec.qual!r}") from exc
    cols = [
        rec.contig,
        str(pos0 + 1),
        rec.id,
        rec.ref,
        ",".join(rec.alts) if rec.alts else ".",
        rec.qual,
        rec.filter,
        _info_column(rec.info, info_decl),
        "GT:AD:DP",
        f"{rec.gt}:{','.join(str(a) for a in rec.ad)}:{rec.dp}",
    ]
    for c in cols:
        if not _FIELD_TEXT.fullmatch(c):
            raise ValueError(f"{where}: invalid VCF field {c!r}")
    return idx, pos0, "\t".join(cols)


def _filter_ids(text: str) -> list[str]:
    if text in ("PASS", "."):
        return []
    ids = text.split(";")
    for f in ids:
        if not _FILTER_ID.fullmatch(f) or f == "PASS":
            raise ValueError(f"invalid FILTER value {text!r}")
    return ids


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
    sequences: Mapping[str, str] | None = None,
) -> None:
    """Write a VCF 4.3 file; records are sorted by contig order, then POS.

    Every INFO key must be declared in ``info_fields`` and GT, AD and DP in ``format_fields``.
    FILTER values other than ``PASS`` get a generic ``##FILTER`` line. When ``sequences``
    (contig ID -> forward-strand sequence) is given, each REF must match the genome at POS.
    """
    if not _SHA256.fullmatch(genome_sha256):
        raise ValueError(f"genome_sha256 must be 64 lowercase hex digits, got {genome_sha256!r}")
    if not sample or _check_text("sample name", sample) != sample.strip():
        raise ValueError(f"invalid VCF sample name {sample!r}")
    _check_text("giap_version", giap_version)
    for r in replicons:
        if not _SAFE_ID.fullmatch(r.id):
            raise ValueError(f"VCF contig ID must match [A-Za-z0-9._-]+, got {r.id!r}")
    info_decl = {f.id: f for f in info_fields}
    format_ids = {f.id for f in format_fields}
    if len(info_decl) != len(info_fields) or len(format_ids) != len(format_fields):
        raise ValueError("duplicate ##INFO or ##FORMAT ID")
    missing = {"GT", "AD", "DP"} - format_ids
    if missing:
        raise ValueError(f"FORMAT fields {sorted(missing)} must be declared in the header")
    info_lines = [_header_line("INFO", f) for f in info_fields]
    format_lines = [_header_line("FORMAT", f) for f in format_fields]
    by_id = {r.id: (i, r) for i, r in enumerate(replicons)}
    recs = list(records)
    body = sorted(_record_key(rec, by_id, info_decl, sequences) for rec in recs)
    filters = sorted({f for rec in recs for f in _filter_ids(rec.filter)})
    lines = [
        "##fileformat=VCFv4.3",
        f"##source=giap-{giap_version}",
        f"##reference=sha256:{genome_sha256}",
        *(f"##contig=<ID={r.id},length={r.length}>" for r in replicons),
        *(f'##FILTER=<ID={f},Description="GIAP filter {f}">' for f in filters),
        *info_lines,
        *format_lines,
        "\t".join(
            ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample]
        ),
        *(line for _, _, line in body),
    ]
    with open_text_write(path) as fh:
        fh.write("\n".join(lines) + "\n")
