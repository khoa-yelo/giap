"""Pure layers derived from a parsed Bakta annotation (GIPPA_SPEC.md §7.1.2, §9.3; M2-04).

Ribosomal-protein set (CAI reference), 16S rRNA selection and anti-Shine-Dalgarno 9-mer, tRNA
copy numbers per anticodon (tAI), db_xref splitting, ``annotation_source`` and origin-aware
nearest genes (TE6). No file I/O; the ``parse_bakta`` step in ``gippa.steps.annotation`` writes.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Sequence
from typing import Literal

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, NcFeature

RIBOSOMAL_PRODUCT: re.Pattern[str] = re.compile(r"^(30S|50S) ribosomal protein")
RIBOSOMAL_GENE: re.Pattern[str] = re.compile(r"^(rps|rpl|rpm)[A-Z]")
DEFAULT_ANTI_SD = "ACCTCCTTA"
ANTI_SD_DEFAULT_FLAG = "anti_sd_default"

_DNA9 = re.compile(r"^[ACGT]{9}$")

AnnotationSource = Literal["bakta_full", "bakta_light", "user"]


def is_ribosomal_protein(g: Gene) -> bool:
    """§7.1.2: product ``^(30S|50S) ribosomal protein`` or gene ``^(rps|rpl|rpm)[A-Z]``."""
    return bool(RIBOSOMAL_PRODUCT.match(g.product or "")) or bool(
        g.gene and RIBOSOMAL_GENE.match(g.gene)
    )


def rrna16s(features: Sequence[NcFeature]) -> list[NcFeature]:
    """The 16S rRNA copies (product contains ``16S`` or gene ``rrs``), in input order."""
    return [f for f in features if f.type == "rrna" and ("16S" in f.product or f.gene == "rrs")]


def anti_sd(rrna16s_nts: Sequence[str], configured: str = "auto") -> tuple[str, tuple[str, ...]]:
    """Anti-SD 9-mer and run flags (§7.1.2, §5.2).

    ``configured`` other than ``auto`` must match ``^[ACGT]{9}$`` and is used as is. Otherwise
    the most common 3' 9-mer across 16S copies (ties: first copy in input order; copies shorter
    than 9 nt are skipped) is used if it contains ``CCTCC``, else the E. coli default and the
    ``anti_sd_default`` flag.
    """
    if configured != "auto":
        if not _DNA9.match(configured):
            raise ValueError(
                f"annotation.ostir.anti_sd must be 'auto' or 9 nt ACGT: {configured!r}"
            )
        return configured, ()
    tails = [s[-9:].upper() for s in rrna16s_nts if len(s) >= 9]
    if tails:
        counts = Counter(tails)  # insertion-ordered: most_common keeps first-seen order on ties
        best = counts.most_common(1)[0][0]
        if "CCTCC" in best:
            return best, ()
    return DEFAULT_ANTI_SD, (ANTI_SD_DEFAULT_FLAG,)


def trna_counts(features: Sequence[NcFeature]) -> dict[str, int]:
    """tRNA gene copies per anticodon (§7.1.2), keys sorted."""
    counts = Counter(f.anticodon for f in features if f.type == "trna" and f.anticodon)
    return {k: counts[k] for k in sorted(counts)}


def split_xrefs(genes: Sequence[Gene]) -> list[tuple[str, str, str]]:
    """``(locus_tag, db, id)`` per db_xref; db is text before the first ``:``; GO keeps it."""
    rows: list[tuple[str, str, str]] = []
    for g in genes:
        for x in g.db_xrefs:
            db, _, ident = x.partition(":")
            rows.append((g.locus_tag, db, x if db == "GO" else ident))
    return rows


def annotation_source(db_type: str) -> AnnotationSource:
    """§9.3 enum from the Bakta DB type (M2 lead decision 7): ``user`` when unknown."""
    if db_type == "full":
        return "bakta_full"
    if db_type == "light":
        return "bakta_light"
    return "user"


def _gap(a: Interval, b: Interval, rep: Replicon) -> int:
    """Gap in bp between two intervals (0 when overlapping or adjacent); wraps if circular."""
    if rep.topology != "circular":
        return max(0, b.start - a.end, a.start - b.end)
    n = rep.length
    la, lb = a.end - a.start, b.end - b.start
    if (b.start - a.start) % n < la or (a.start - b.start) % n < lb:
        return 0
    return min((b.start - a.end) % n, (a.start - b.end) % n)


def nearest_genes(
    feature: NcFeature, genes: Sequence[Gene], replicon: Replicon, n: int = 2
) -> list[tuple[str, int]]:
    """The ``n`` genes on ``feature``'s contig closest to it: ``(locus_tag, gap_bp)`` ascending.

    Ties keep the input order. Circular replicons are origin-aware; linear contigs never wrap.
    """
    on = [g for g in genes if g.interval.contig == feature.interval.contig]
    dist = [(g.locus_tag, _gap(feature.interval, g.interval, replicon)) for g in on]
    return sorted(dist, key=lambda x: x[1])[:n]
