"""UniOP adapter: Prodigal-style GFF in, ``uniop.pred`` out (GIAP_SPEC.md §7.1.1, §7.1.3; M2-12).

Bakta CDS, sORFs and pseudogenes are written in genome order as a Prodigal-style GFF with IDs
``<contig>_<k>`` (``k`` 1-based per contig); the returned map leads back to locus tags. UniOP is
never imported (ADR-0008); it runs as an external process (``giap.annotate.uniop_cmd``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from giap.io._open import open_text, open_text_write
from giap.model.coords import Replicon, to_one_based_closed
from giap.model.genes import Gene

PairSource = Literal["uniop", "distance_model"]


@dataclass(frozen=True, slots=True)
class PairProb:
    """Probability that adjacent same-strand genes ``gene_a`` → ``gene_b`` share an operon."""

    gene_a: str
    gene_b: str
    uniop_p: float
    source: PairSource


def genome_order(genes: Sequence[Gene], contigs: Sequence[str]) -> dict[str, list[Gene]]:
    """Genes per contig (``contigs`` order) sorted by start, end, strand and locus tag."""
    out: dict[str, list[Gene]] = {c: [] for c in contigs}
    for g in genes:
        if g.interval.contig not in out:
            raise ValueError(f"gene {g.locus_tag}: unknown contig {g.interval.contig!r}")
        out[g.interval.contig].append(g)
    for gs in out.values():
        gs.sort(key=lambda g: (g.interval.start, g.interval.end, g.interval.strand, g.locus_tag))
    return out


def write_prodigal_gff(
    path: Path, genes: Sequence[Gene], replicons: Sequence[Replicon]
) -> list[tuple[str, str]]:
    """Write the UniOP input GFF; returns ``[(prodigal_id, locus_tag)]`` in file order."""
    id_map: list[tuple[str, str]] = []
    lines = ["##gff-version 3"]
    for contig, gs in genome_order(genes, [r.id for r in replicons]).items():
        for k, g in enumerate(gs, start=1):
            pid = f"{contig}_{k}"
            start1, end1 = to_one_based_closed(g.interval)
            cols = [contig, "GIAP", "CDS", str(start1), str(end1), ".", g.interval.strand, "0"]
            lines.append("\t".join([*cols, f"ID={pid}"]))
            id_map.append((pid, g.locus_tag))
    with open_text_write(path, compress=False) as fh:
        fh.write("\n".join(lines) + "\n")
    return id_map


def parse_uniop_pred(path: Path, id_map: Mapping[str, str]) -> list[PairProb]:
    """Parse ``uniop.pred`` (header row; first two columns gene IDs, probability column ``p`` or
    the last column) into locus-tag pairs in file order. An unknown ID raises ``ValueError``."""
    with open_text(path) as fh:
        rows = [ln.split() for ln in fh.read().splitlines() if ln.strip()]
    if not rows:
        return []
    header, body = rows[0], rows[1:]
    col = header.index("p") if "p" in header else len(header) - 1
    out: list[PairProb] = []
    for r in body:
        missing = [x for x in r[:2] if x not in id_map]
        if missing:
            raise ValueError(f"{path}: UniOP gene ID {missing[0]!r} not in the input map")
        out.append(PairProb(id_map[r[0]], id_map[r[1]], float(r[col]), "uniop"))
    return out
