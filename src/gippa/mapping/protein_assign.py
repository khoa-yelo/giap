"""Assignment of protein alignments to Bakta CDS/sORF (GIPPA_SPEC.md §7.2.3; issue M3-12).

An alignment is assigned to a same-strand CDS or sORF when identity >= ``min_identity``, query
coverage >= ``min_coverage`` and reciprocal overlap >= 0.80 (fixed by §7.2.3). Among passing
genes the largest overlap in bp wins, ties go to the lower locus tag. Every alignment passing
identity and coverage becomes a ``protein_hits`` row (``gene`` NA when no gene passes). A protein
without an assigned alignment is a leftover for tblastn (M3-13). Frameshifts and internal stops
are flags only (no identity penalty).
"""

from __future__ import annotations

import bisect
from collections.abc import Sequence
from typing import Literal

import pandas as pd

from gippa.config.model import ProteinsConfig
from gippa.mapping.miniprot import ProteinAlignment
from gippa.model.coords import Interval, Replicon, check_interval
from gippa.model.genes import Gene, GeneModel

SimilarityClass = Literal["identical", "near", "divergent"]

MIN_RECIPROCAL_OVERLAP = 0.80  # §7.2.3, not configurable

PROTEIN_HITS_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "protein_len_aa",
    "route",
    "contig",
    "start",
    "end",
    "strand",
    "rank",
    "identity",
    "query_coverage",
    "reciprocal_overlap",
    "class",
    "gene",
    "flags",
)

_TARGET_TYPES = ("cds", "sorf")


def _shifts(replicon: Replicon) -> tuple[int, ...]:
    n = replicon.length
    return (-n, 0, n) if replicon.topology == "circular" else (0,)


def _overlap_bp(a: Interval, b: Interval, replicon: Replicon) -> int:
    """Overlap of two (possibly unwrapped) intervals, positions taken modulo the replicon."""
    total = 0
    for shift in _shifts(replicon):
        lo = max(a.start + shift, b.start)
        hi = min(a.end + shift, b.end)
        if hi > lo:
            total += hi - lo
    return total


def reciprocal_overlap(aln: Interval, cds: Interval, replicon: Replicon) -> float:
    """``min(overlap / len(cds), overlap / len(aln))``; origin copies count when circular."""
    ov = _overlap_bp(aln, cds, replicon)
    if ov == 0:
        return 0.0
    return min(ov / cds.length_bp, ov / aln.length_bp)


def similarity_class(identity: float, coverage: float, near: float) -> SimilarityClass:
    """``identical`` (identity, coverage 1.0), ``near`` (identity >= near), else ``divergent``."""
    if identity >= 1.0 and coverage >= 1.0:
        return "identical"
    if identity >= near:
        return "near"
    return "divergent"


class _GeneIndex:
    """Same-strand CDS/sORF per (contig, strand), sorted by start, for interval queries."""

    def __init__(self, gm: GeneModel) -> None:
        buckets: dict[tuple[str, str], list[Gene]] = {}
        for g in gm.genes:
            if g.feature_type in _TARGET_TYPES:
                buckets.setdefault((g.interval.contig, g.interval.strand), []).append(g)
        self._genes: dict[tuple[str, str], list[Gene]] = {}
        self._starts: dict[tuple[str, str], list[int]] = {}
        self._max_len: dict[tuple[str, str], int] = {}
        for key, genes in buckets.items():
            genes.sort(key=lambda g: (g.interval.start, g.interval.end, g.locus_tag))
            self._genes[key] = genes
            self._starts[key] = [g.interval.start for g in genes]
            self._max_len[key] = max(g.interval.length_bp for g in genes)

    def candidates(self, iv: Interval, replicon: Replicon) -> list[Gene]:
        key = (iv.contig, iv.strand)
        genes = self._genes.get(key)
        if not genes:
            return []
        starts, max_len = self._starts[key], self._max_len[key]
        seen: dict[str, Gene] = {}
        for shift in _shifts(replicon):
            lo, hi = iv.start + shift, iv.end + shift
            i = bisect.bisect_left(starts, lo - max_len + 1)
            j = bisect.bisect_left(starts, hi)
            for g in genes[i:j]:
                if g.interval.end > lo:
                    seen[g.locus_tag] = g
        return list(seen.values())


def _best_gene(
    aln: ProteinAlignment, index: _GeneIndex, replicon: Replicon
) -> tuple[str, float] | None:
    best: tuple[int, str, float] | None = None
    for g in index.candidates(aln.interval, replicon):
        ro = reciprocal_overlap(aln.interval, g.interval, replicon)
        if ro < MIN_RECIPROCAL_OVERLAP:
            continue
        ov = _overlap_bp(aln.interval, g.interval, replicon)
        if best is None or (-ov, g.locus_tag) < (-best[0], best[1]):
            best = (ov, g.locus_tag, ro)
    return None if best is None else (best[1], best[2])


def _flags(aln: ProteinAlignment) -> str | None:
    flags = [f for f, n in (("fs", aln.frameshift), ("st", aln.stop_codon)) if n > 0]
    return ";".join(flags) if flags else None


def _sort_key(
    aln: ProteinAlignment,
) -> tuple[str, int, str, int, int, str, float, bool, float, int, int, int, int, int, str]:
    iv = aln.interval
    return (
        aln.protein_id,
        aln.rank,
        iv.contig,
        iv.start,
        iv.end,
        iv.strand,
        aln.identity,
        aln.positive is None,
        aln.positive if aln.positive is not None else 0.0,
        aln.q_start,
        aln.q_end,
        aln.frameshift,
        aln.stop_codon,
        aln.protein_len_aa,
        aln.route,
    )


def _empty_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_id": pd.Series(dtype=object),
            "protein_len_aa": pd.Series(dtype="int64"),
            "route": pd.Series(dtype=object),
            "contig": pd.Series(dtype=object),
            "start": pd.Series(dtype="int64"),
            "end": pd.Series(dtype="int64"),
            "strand": pd.Series(dtype=object),
            "rank": pd.Series(dtype="int64"),
            "identity": pd.Series(dtype="float64"),
            "query_coverage": pd.Series(dtype="float64"),
            "reciprocal_overlap": pd.Series(dtype="float64"),
            "class": pd.Series(dtype=object),
            "gene": pd.Series(dtype=object),
            "flags": pd.Series(dtype=object),
        }
    )


def _check_canonical(aln: ProteinAlignment, gm: GeneModel) -> None:
    iv = aln.interval
    replicon = gm.replicons.get(iv.contig)
    if replicon is None:
        raise ValueError(f"alignment of {aln.protein_id!r} is on unknown contig {iv.contig!r}")
    if iv.start < 0:
        raise ValueError(f"alignment of {aln.protein_id!r} has negative start {iv.start}")
    try:
        check_interval(iv, replicon)
    except ValueError as exc:
        raise ValueError(f"alignment of {aln.protein_id!r}: {exc}") from exc


def assign_alignments(
    alns: Sequence[ProteinAlignment], gm: GeneModel, cfg: ProteinsConfig
) -> tuple[pd.DataFrame, set[str]]:
    """Return ``(protein_hits rows, leftover protein IDs)``.

    Rows (0-based half-open, unwrapped ``start``/``end``) are every alignment passing identity and
    coverage, sorted by protein, rank, coordinates; ``gene`` and ``reciprocal_overlap`` are NA
    when no CDS/sORF passes. Leftovers are the aligned proteins without an assigned row.

    Raises ``ValueError`` unless every interval is canonical (``0 <= start < L``, unwrapped
    ``end <= start + L``, linear contigs never wrap); callers normalise chains first.
    """
    for aln in alns:
        _check_canonical(aln, gm)
    index = _GeneIndex(gm)
    rows: list[dict[str, object]] = []
    assigned: set[str] = set()
    for aln in sorted(alns, key=_sort_key):
        cov = aln.query_coverage
        if aln.identity < cfg.min_identity or cov < cfg.min_coverage:
            continue
        iv = aln.interval
        replicon = gm.replicons[iv.contig]
        best = _best_gene(aln, index, replicon)
        if best is not None:
            assigned.add(aln.protein_id)
        rows.append(
            {
                "protein_id": aln.protein_id,
                "protein_len_aa": aln.protein_len_aa,
                "route": aln.route,
                "contig": iv.contig,
                "start": iv.start,
                "end": iv.end,
                "strand": iv.strand,
                "rank": aln.rank,
                "identity": aln.identity,
                "query_coverage": cov,
                "reciprocal_overlap": float("nan") if best is None else best[1],
                "class": similarity_class(aln.identity, cov, cfg.near_identity),
                "gene": None if best is None else best[0],
                "flags": _flags(aln),
            }
        )
    df = _empty_frame()
    if rows:
        df = pd.DataFrame(rows, columns=list(PROTEIN_HITS_COLUMNS)).astype(df.dtypes.to_dict())
    leftovers = {a.protein_id for a in alns} - assigned
    return df, leftovers
