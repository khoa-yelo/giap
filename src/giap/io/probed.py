"""proBed 1.0 writer (GIAP_SPEC.md §7.2.9, §6.5; ADR-0005 item 1).

Always BED12 + 13 proBed columns (25 fields) under a ``# proBed-version 1.0`` header; a record
without blocks gets one block spanning it. Missing proBed values are ``.`` (not ``NA``).
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

from giap.io.bed import BedRecord, bed_parts, sorted_lines
from giap.io.tsv import NA, Cell, format_value, write_text
from giap.model import Interval, Replicon

PROBED_EXTRA = (
    "proteinAccession",
    "peptideSequence",
    "uniqueness",
    "genomeReferenceVersion",
    "psmScore",
    "fdr",
    "modifications",
    "charge",
    "expMassToCharge",
    "calcMassToCharge",
    "psmRank",
    "datasetID",
    "uri",
)
HEADER = "# proBed-version 1.0"


@dataclass(frozen=True, slots=True)
class ProBedRecord:
    """A BED feature plus proBed columns; keys of ``extra`` must be in ``PROBED_EXTRA``."""

    bed: BedRecord
    extra: Mapping[str, Cell]


def _extra_fields(extra: Mapping[str, Cell], genome_version: str) -> list[str]:
    unknown = set(extra) - set(PROBED_EXTRA)
    if unknown:
        raise ValueError(f"unknown proBed columns: {sorted(unknown)}")
    values = {"genomeReferenceVersion": genome_version, **extra}
    out = []
    for key in PROBED_EXTRA:
        s = format_value(values.get(key))
        out.append("." if s in (NA, "") else s)
    return out


def write_probed(
    path: Path,
    records: Iterable[ProBedRecord],
    replicons: Mapping[str, Replicon],
    genome_version: str,
) -> None:
    """Write sorted proBed 1.0 records; origin-crossing records are split like BED."""
    parts: list[tuple[Interval, list[str]]] = []
    for rec in records:
        bed = rec.bed
        if not bed.blocks:
            bed = dataclasses.replace(bed, blocks=((bed.interval.length_bp, 0),))
        extra = _extra_fields(rec.extra, genome_version)
        for piv, fields in bed_parts(bed, replicons[bed.interval.contig]):
            parts.append((piv, fields + extra))
    lines = [HEADER, *sorted_lines(parts, replicons)]
    write_text(path, "".join(line + "\n" for line in lines))
