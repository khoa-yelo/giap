"""geNOMAD 1.12.0 summary parser (GIAP_SPEC.md §7.1.1, §7.1.7; M2-18).

Reads ``<outdir>/<prefix>_summary/<prefix>_{virus_summary,plasmid_summary,virus_genes}.tsv``.
``<contig>|provirus_<a>_<b>`` (1-based closed) → ``Interval(contig, a-1, b, ".")``; a ``seq_name``
without ``|provirus_`` is a whole-contig virus (``interval=None``). Virus gene coordinates are
absolute 1-based closed contig positions (M2-18 lead decision 1); ``NA`` → ``None`` / ``()``.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from giap.model.coords import Interval, Replicon, Strand

_PROVIRUS = re.compile(r"^(?P<contig>.+)\|provirus_(?P<a>\d+)_(?P<b>\d+)$")
_STRAND: dict[str, Strand] = {"1": "+", "-1": "-"}


@dataclass(frozen=True, slots=True)
class Provirus:
    contig: str
    interval: Interval | None  # None = whole contig classified viral
    virus_score: float
    fdr: float | None
    n_hallmarks: int
    taxonomy: str


@dataclass(frozen=True, slots=True)
class PlasmidCall:
    contig: str
    plasmid_score: float
    fdr: float | None
    n_hallmarks: int
    conjugation_genes: tuple[str, ...]
    amr_genes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class VirusGene:
    contig: str
    interval: Interval
    marker: str | None
    description: str | None
    virus_hallmark: bool


@dataclass(frozen=True, slots=True)
class GenomadResult:
    proviruses: tuple[Provirus, ...]
    plasmids: tuple[PlasmidCall, ...]
    virus_genes: tuple[VirusGene, ...]


def _na(value: str) -> str | None:
    return None if value == "NA" else value


def _opt_float(value: str) -> float | None:
    v = _na(value)
    return None if v is None else float(v)


def _gene_list(value: str) -> tuple[str, ...]:
    v = _na(value)
    return () if not v else tuple(v.split(";"))


def _replicon(contig: str, replicons: Mapping[str, Replicon]) -> Replicon:
    if contig not in replicons:
        raise ValueError(f"geNOMAD: contig {contig!r} not in the genome")
    return replicons[contig]


def _interval(
    contig: str, start1: int, end1: int, strand: Strand, replicons: Mapping[str, Replicon]
) -> Interval:
    rep = _replicon(contig, replicons)
    if not 1 <= start1 <= end1 <= rep.length:
        raise ValueError(f"geNOMAD: {contig}:{start1}-{end1} outside contig (length {rep.length})")
    return Interval(contig, start1 - 1, end1, strand)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _provirus(row: dict[str, str], replicons: Mapping[str, Replicon]) -> Provirus:
    name = row["seq_name"]
    m = _PROVIRUS.match(name)
    if m:
        contig = m["contig"]
        iv: Interval | None = _interval(contig, int(m["a"]), int(m["b"]), ".", replicons)
    else:
        contig, iv = _replicon(name, replicons).id, None
    return Provirus(
        contig,
        iv,
        float(row["virus_score"]),
        _opt_float(row["fdr"]),
        int(row["n_hallmarks"]),
        row["taxonomy"],
    )


def _plasmid(row: dict[str, str], replicons: Mapping[str, Replicon]) -> PlasmidCall:
    return PlasmidCall(
        _replicon(row["seq_name"], replicons).id,
        float(row["plasmid_score"]),
        _opt_float(row["fdr"]),
        int(row["n_hallmarks"]),
        _gene_list(row["conjugation_genes"]),
        _gene_list(row["amr_genes"]),
    )


def _virus_gene(row: dict[str, str], replicons: Mapping[str, Replicon]) -> VirusGene:
    seq_name = row["gene"].rsplit("_", 1)[0]
    contig = seq_name.split("|", 1)[0]
    if row["strand"] not in _STRAND:
        raise ValueError(f"geNOMAD: gene {row['gene']!r}: invalid strand {row['strand']!r}")
    iv = _interval(contig, int(row["start"]), int(row["end"]), _STRAND[row["strand"]], replicons)
    return VirusGene(
        contig,
        iv,
        _na(row["marker"]),
        _na(row["annotation_description"]),
        row["virus_hallmark"] == "1",
    )


def parse_genomad(outdir: Path, prefix: str, replicons: Mapping[str, Replicon]) -> GenomadResult:
    """Parse the three geNOMAD summary TSVs of ``prefix`` under ``outdir``."""
    d = outdir / f"{prefix}_summary"
    return GenomadResult(
        tuple(_provirus(r, replicons) for r in _read(d / f"{prefix}_virus_summary.tsv")),
        tuple(_plasmid(r, replicons) for r in _read(d / f"{prefix}_plasmid_summary.tsv")),
        tuple(_virus_gene(r, replicons) for r in _read(d / f"{prefix}_virus_genes.tsv")),
    )
