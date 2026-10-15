"""Module 2 outputs (GIAP_SPEC.md §7.2.9, §9.2, §9.3, §6.5, §6.6; docs/specs/M3-16).

``write_mapping_outputs`` writes ``peptide_hits.tsv``, ``protein_hits.tsv``, ``mappings.probed``,
``mappings.bed``, ``mappings.gff3``, ``unmapped.tsv`` and ``mapping_fdr.tsv``. Inputs carry 0-based
half-open, unwrapped coordinates; the shift to 1-based closed happens here, in the TSV rows, and
in the M1-04/M1-05 writers (BED/proBed split origin-crossing hits, GFF3 keeps ``end > length``;
``match_part`` children starting beyond the origin are shifted by -L).
The row preparation is pure; file I/O is delegated to ``giap.io``.
"""

from __future__ import annotations

import dataclasses
import math
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any, cast

import pandas as pd

from giap.io.bed import BedRecord, Block, write_bed
from giap.io.gff3 import AttrValue, GffFeature, write_gff3
from giap.io.probed import ProBedRecord, write_probed
from giap.io.tsv import format_value, write_tsv
from giap.mapping.contaminants import UnmappedRow
from giap.mapping.protein_assign import PROTEIN_HITS_COLUMNS
from giap.model import Interval, Replicon, interval_sort_key

__all__ = [
    "MAPPING_FDR_COLUMNS",
    "PEPTIDE_HITS_COLUMNS",
    "PROTEIN_HITS_COLUMNS",
    "UNMAPPED_COLUMNS",
    "write_mapping_outputs",
]

PEPTIDE_HITS_COLUMNS: tuple[str, ...] = (
    "peptide",
    "peptide_len",
    "source",
    "best_q_input",
    "pep_input",
    "score_input",
    "contig",
    "start",
    "end",
    "strand",
    "frame",
    "crosses_origin",
    "route",
    "n_diff",
    "diff_types",
    "placement_q",
    "length_class",
    "locus_class",
    "uniqueness",
    "gene",
    "novel_candidate_id",
    "flags",
)
UNMAPPED_COLUMNS: tuple[str, ...] = tuple(f.name for f in dataclasses.fields(UnmappedRow))
MAPPING_FDR_COLUMNS: tuple[str, ...] = (
    "length_class",
    "n_diff",
    "targets",
    "decoys",
    "pi_wrong",
    "fdr_stratum",
    "n_accepted",
    "decoy_type",
    "section",
    "score_bin",
    "precision",
    "ci_low",
    "ci_high",
    "c_value",
    "flags",
)

PROBED_UNIQUENESS = {
    "unique": "unique",
    "multi_identical": "not-unique[same-set]",
    "shared": "not-unique[conflict]",
}
_INT_COLUMNS = frozenset(
    {"peptide_len", "start", "end", "frame", "n_diff", "protein_len_aa", "rank"}
    | {"targets", "decoys", "n_accepted", "score_bin"}
)
_LEADING_INT = re.compile(r"\d+")
Row = dict[str, Any]


def _missing(v: object) -> bool:
    if v is None or v is pd.NA or (isinstance(v, str) and v == ""):
        return True
    if isinstance(v, float) and math.isnan(v):
        return True
    return isinstance(v, list | tuple) and len(v) == 0


def _cell(row: Mapping[str, Any], col: str) -> Any:
    """One normalised cell: missing/empty values become ``None``, integral floats ``int``."""
    v = row.get(col)
    if _missing(v):
        return None
    if col in _INT_COLUMNS and isinstance(v, float) and v.is_integer():
        return int(v)
    return v


def _records(df: pd.DataFrame) -> list[Row]:
    return cast("list[Row]", df.to_dict("records"))


def _interval(row: Mapping[str, Any]) -> Interval:
    return Interval(str(row["contig"]), int(row["start"]), int(row["end"]), row["strand"])


def _sorted_hits(rows: list[Row], order: Mapping[str, int], id_col: str) -> list[Row]:
    """Sort by contig (FASTA order), start, end, strand, ID; ties by the formatted row."""

    def key(r: Row) -> tuple[Any, ...]:
        full = tuple(format_value(_cell(r, c)) for c in sorted(r) if c != "blocks")
        return (interval_sort_key(_interval(r), order, str(r[id_col])), full)

    return sorted(rows, key=key)


def _tsv_row(row: Row, columns: Sequence[str], replicon: Replicon) -> list[object]:
    """§6.5: 1-based closed, ``end > length`` kept for origin-crossing hits."""
    cells = {c: _cell(row, c) for c in columns}
    cells["start"] = int(row["start"]) + 1
    cells["end"] = int(row["end"])
    if "crosses_origin" in cells:
        cells["crosses_origin"] = int(row["end"]) > replicon.length
    return [cells[c] for c in columns]


def _blocks(spec: object, start: int, end: int) -> tuple[Block, ...]:
    """BED12 ``(size, relative start)`` blocks from M3-08's ``"s1-e1;s2-e2"`` (absolute).

    A block overlapping the next one (``delta = -1``) is shortened so blocks do not overlap.
    """
    if _missing(spec):
        return ((end - start, 0),)
    spans = sorted((int(s), int(e)) for s, e in (p.split("-") for p in str(spec).split(";") if p))
    out: list[Block] = []
    for i, (s, e) in enumerate(spans):
        if i + 1 < len(spans):
            e = min(e, spans[i + 1][0])
        out.append((e - s, s - start))
    return tuple(out)


def _bed_record(row: Row, name: str, score: int) -> BedRecord:
    iv = _interval(row)
    return BedRecord(iv, name, score, blocks=_blocks(row.get("blocks"), iv.start, iv.end))


def _attrs(pairs: Iterable[tuple[str, Any]]) -> dict[str, AttrValue]:
    """GFF3 attributes in order; missing values are omitted, ``;`` lists become value lists."""
    out: dict[str, AttrValue] = {}
    for key, v in pairs:
        if _missing(v):
            continue
        if isinstance(v, list | tuple):
            out[key] = [format_value(x) for x in v]
        elif key in ("flags", "diff_types") and isinstance(v, str):
            out[key] = [x for x in v.split(";") if x]
        else:
            out[key] = format_value(v)
    return out


def _wrapped_parts(s: int, e: int, length: int) -> list[tuple[int, int]]:
    """Unwrapped ``[s, e)`` with ``s >= length`` shifted by -L.

    A block spanning the origin keeps ``end > length`` like every origin-crossing GFF3 feature
    (§6.5; pinned by tests/unit/test_mapping_outputs_review.py).
    """
    if s >= length:
        return [(s - length, e - length)]
    return [(s, e)]


def _frameshift_base(spec: str, length: int) -> tuple[str, int] | None:
    """The odd base of a +-1 frameshift as ``(attribute, 1-based position wrapped into [1, L])``.

    ``delta = -1`` (second block starts 1 nt before the first ends): the base read twice;
    ``delta = +1`` (a 1-nt gap): the base skipped (M3-16 lead decision 2, M3-08 edge cases).
    """
    spans = sorted((int(s), int(e)) for s, e in (p.split("-") for p in spec.split(";") if p))
    if len(spans) != 2:
        return None
    (_, e1), (s2, _) = spans
    if s2 == e1 - 1:
        return "frameshift_reread", s2 % length + 1
    if s2 == e1 + 1:
        return "frameshift_skipped", e1 % length + 1
    return None


def _gff_features(
    row: Row, fid: str, name: str, extra: Iterable[tuple[str, Any]], replicon: Replicon
) -> list[GffFeature]:
    """``nucleotide_to_protein_match`` (SO:0000568), plus ``match_part`` children per block.

    The parent keeps ``end > length`` when it crosses the origin; a child starting at or beyond
    ``length`` is shifted by -L, a child spanning the origin keeps ``end > length``.
    """
    iv = _interval(row)
    spec = row.get("blocks")
    odd = None if _missing(spec) else _frameshift_base(str(spec), replicon.length)
    attrs = _attrs([("ID", fid), ("Name", name), *extra, *([odd] if odd else [])])
    feats = [GffFeature(iv, "GIAP", "nucleotide_to_protein_match", attributes=attrs)]
    if not _missing(spec):
        k = 0
        for part in str(spec).split(";"):
            if not part:
                continue
            s, e = (int(x) for x in part.split("-"))
            for ps, pe in _wrapped_parts(s, e, replicon.length):
                k += 1
                child = Interval(iv.contig, ps, pe, iv.strand)
                part_attrs = {"ID": f"{fid}.{k}", "Parent": fid, "Name": name}
                feats.append(GffFeature(child, "GIAP", "match_part", attributes=part_attrs))
    return feats


def _length_class_key(label: object) -> tuple[int, int, str]:
    """Sort key for M3-10 length-class labels: ``lo-hi`` (numeric ``lo``), then ``all``, then NA.

    Expects M3-10's ``lo-hi`` labels (e.g. ``7-8``, ``9-30``) and the pooled ``all``; any other
    label sorts by its first integer (0 if none), ties broken by the text.
    """
    if _missing(label):
        return (2, 0, "")
    text = str(label)
    if text == "all":
        return (1, 0, text)
    m = _LEADING_INT.search(text)
    return (0, int(m.group()) if m else 0, text)


def _fdr_sort_key(row: Row) -> tuple[Any, ...]:
    """Section, decoy type, n_diff, length class (numeric lower bound, ``all`` last), bin."""

    def opt(v: Any) -> tuple[int, Any]:
        return (1, 0) if v is None else (0, v)

    return (
        str(_cell(row, "section") or ""),
        opt(_cell(row, "decoy_type")),
        opt(_cell(row, "n_diff")),
        _length_class_key(_cell(row, "length_class")),
        opt(_cell(row, "score_bin")),
        tuple(format_value(_cell(row, c)) for c in MAPPING_FDR_COLUMNS),
    )


def _peptide_extra(r: Row) -> list[tuple[str, Any]]:
    keys = ("route", "n_diff", "diff_types", "placement_q", "locus_class", "uniqueness")
    return [(k, r.get(k)) for k in (*keys, "gene", "flags")]


def _protein_extra(r: Row) -> list[tuple[str, Any]]:
    return [(k, r.get(k)) for k in ("identity", "query_coverage", "class", "rank", "gene", "flags")]


def _features(
    rows: list[Row],
    prefix: str,
    id_col: str,
    extra: Callable[[Row], list[tuple[str, Any]]],
    replicons: Mapping[str, Replicon],
) -> list[GffFeature]:
    out: list[GffFeature] = []
    for i, r in enumerate(rows, start=1):
        rep = replicons[str(r["contig"])]
        out += _gff_features(r, f"{prefix}{i}", str(r[id_col]), extra(r), rep)
    return out


def write_mapping_outputs(
    outdir: Path,
    peptide_hits: pd.DataFrame,
    protein_hits: pd.DataFrame,
    unmapped: Sequence[UnmappedRow],
    fdr_rows: pd.DataFrame,
    replicons: Mapping[str, Replicon],
    genome_sha256: str,
) -> None:
    """Write the seven Module 2 output files into ``outdir`` (deterministic, no paths/times)."""
    order = {cid: i for i, cid in enumerate(replicons)}
    peps = _sorted_hits(_records(peptide_hits), order, "peptide")
    prots = _sorted_hits(_records(protein_hits), order, "protein_id")

    write_tsv(
        outdir / "peptide_hits.tsv",
        PEPTIDE_HITS_COLUMNS,
        (_tsv_row(r, PEPTIDE_HITS_COLUMNS, replicons[r["contig"]]) for r in peps),
    )
    write_tsv(
        outdir / "protein_hits.tsv",
        PROTEIN_HITS_COLUMNS,
        (_tsv_row(r, PROTEIN_HITS_COLUMNS, replicons[r["contig"]]) for r in prots),
    )
    unmapped_rows = sorted(
        (dataclasses.astuple(u) for u in unmapped), key=lambda t: (t[2], t[0], t[1:])
    )
    write_tsv(
        outdir / "unmapped.tsv",
        UNMAPPED_COLUMNS,
        ([None if _missing(v) else v for v in t] for t in unmapped_rows),
    )
    fdr = sorted(_records(fdr_rows), key=_fdr_sort_key)
    write_tsv(
        outdir / "mapping_fdr.tsv",
        MAPPING_FDR_COLUMNS,
        ([_cell(r, c) for c in MAPPING_FDR_COLUMNS] for r in fdr),
    )

    pep_beds = [
        _bed_record(r, str(r["peptide"]), 1000 if r["route"] == "exact" else 500) for r in peps
    ]
    prot_beds = [
        _bed_record(r, str(r["protein_id"]), 1000 if r["class"] == "identical" else 500)
        for r in prots
    ]
    write_bed(outdir / "mappings.bed", [*pep_beds, *prot_beds], replicons)
    probed = [
        ProBedRecord(
            bed,
            {
                "proteinAccession": _cell(r, "gene"),
                "peptideSequence": str(r["peptide"]),
                "uniqueness": PROBED_UNIQUENESS.get(str(r["uniqueness"])),
                "psmScore": _cell(r, "score_input"),
                "fdr": _cell(r, "placement_q"),
            },
        )
        for r, bed in zip(peps, pep_beds, strict=True)
    ]
    write_probed(outdir / "mappings.probed", probed, replicons, f"sha256:{genome_sha256}")
    features = _features(peps, "peptide_hit_", "peptide", _peptide_extra, replicons)
    features += _features(prots, "protein_hit_", "protein_id", _protein_extra, replicons)
    write_gff3(outdir / "mappings.gff3", features, list(replicons.values()))
