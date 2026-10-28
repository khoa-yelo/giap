"""UniOP adapter: Prodigal-style GFF in, ``uniop.pred`` out (GIPPA_SPEC.md §7.1.1, §7.1.3; M2-12).

Bakta CDS, sORFs and pseudogenes are written in genome order as a Prodigal-style GFF with IDs
``<contig>_<k>`` (``k`` 1-based per contig); the returned map leads back to locus tags. UniOP is
never imported (ADR-0008); it runs as an external process (``gippa.annotate.uniop_cmd``).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from gippa.io._open import open_text, open_text_write
from gippa.model.coords import Replicon, to_one_based_closed
from gippa.model.genes import Gene

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
            cols = [contig, "GIPPA", "CDS", str(start1), str(end1), ".", g.interval.strand, "0"]
            lines.append("\t".join([*cols, f"ID={pid}"]))
            id_map.append((pid, g.locus_tag))
    with open_text_write(path, compress=False) as fh:
        fh.write("\n".join(lines) + "\n")
    return id_map


UNIOP_TITLE = "operon prediction scores"
UNIOP_HEADER = ["Gene A", "Gene B", "Prediction"]


def uniop_gene_order(gff: Path) -> list[tuple[str, str, str]]:
    """``(ID, contig, strand)`` of every GFF gene in UniOP's order: stable sort by start."""
    rows: list[tuple[int, str, str, str]] = []
    with open_text(gff) as fh:
        for ln in fh.read().splitlines():
            if ln and not ln.startswith("#"):
                f = ln.split("\t")
                rows.append((int(f[3]), f[8].split(";")[0].removeprefix("ID="), f[0], f[6]))
    return [(i, c, s) for _, i, c, s in sorted(rows, key=lambda r: r[0])]


def parse_uniop_pred(path: Path, gff: Path, id_map: Mapping[str, str]) -> list[PairProb]:
    """Parse UniOP's ``uniop.pred`` (pin 3ed14dfc) into locus-tag pairs in row order.

    Title ``operon prediction scores``, header ``Gene A/Gene B/Prediction``, rows ``k  k+1  p`` of
    1-based row numbers into UniOP's gene order (:func:`uniop_gene_order` of ``gff``, the exact
    input UniOP ran on). Row ``k`` pairs sorted gene ``k`` with the next gene on its contig; empty
    probabilities are skipped. An empty file or title + header only gives ``[]``. ``ValueError``
    on a missing title/header, a row number outside the gene list, a probability without a
    same-contig same-strand partner, or a GFF ID absent from ``id_map``.
    """
    with open_text(path) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    if not lines:
        return []
    if lines[0].strip() != UNIOP_TITLE or len(lines) < 2:
        raise ValueError(f"{path}: not a UniOP prediction file (title {UNIOP_TITLE!r} missing)")
    if [c.strip() for c in lines[1].split("\t")] != UNIOP_HEADER:
        raise ValueError(f"{path}: unexpected UniOP header {lines[1]!r}")
    order = uniop_gene_order(gff)
    out: list[PairProb] = []
    for ln in lines[2:]:
        cells = ln.split("\t")
        if len(cells) < 3 or not cells[2].strip():
            continue
        k = int(cells[0])
        if not 1 <= k <= len(order):
            raise ValueError(f"{path}: row number {k} outside the {len(order)} input genes")
        a_id, contig, strand = order[k - 1]
        nxt = next((g for g in order[k:] if g[1] == contig), None)
        if nxt is None or nxt[2] != strand:
            raise ValueError(f"{path}: row {k} ({a_id}) has no same-strand partner on {contig}")
        missing = [x for x in (a_id, nxt[0]) if x not in id_map]
        if missing:
            raise ValueError(f"{path}: UniOP gene ID {missing[0]!r} not in the input map")
        out.append(PairProb(id_map[a_id], id_map[nxt[0]], float(cells[2]), "uniop"))
    return out
