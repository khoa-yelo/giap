"""Gene x run detection matrix and exact (UpSet) intersections (M9-22; GIPPA_SPEC.md §7.9.3).

Both functions are pure and return frames with the ``compare_detection`` columns (ADR-0046: one
long table, ``section`` in {matrix, intersection}).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

__all__ = ["COLUMNS", "detection_matrix", "gene_order", "intersections"]

COLUMNS = (
    "section", "locus_tag", "run", "detected", "detection_status", "intersection", "n_genes",
)  # fmt: skip
NONE = "none"


def gene_order(genes: pd.DataFrame) -> list[str]:
    """Locus tags in §6.5 order: contig as first seen, start, end, strand, locus tag."""
    contig_rank = {c: i for i, c in enumerate(pd.unique(genes["contig"]))}
    keyed = genes.assign(_c=genes["contig"].map(contig_rank))
    keyed = keyed.sort_values(["_c", "start", "end", "strand", "locus_tag"], kind="stable")
    return [str(t) for t in keyed["locus_tag"]]


def _typed(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.astype({"detected": "boolean", "n_genes": "Int64"})[list(COLUMNS)]


def detection_matrix(runs: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """The matrix section: one row per gene x run, genes in §6.5 order of the first run, runs in
    label order; ``detected`` and ``detection_status`` copied; ``intersection``/``n_genes`` NA."""
    labels = list(runs)
    order = gene_order(runs[labels[0]])
    by_run = {lab: runs[lab].set_index("locus_tag") for lab in labels}
    rows = [
        ("matrix", tag, lab, by_run[lab].at[tag, "detected"],
         by_run[lab].at[tag, "detection_status"], None, None)
        for tag in order for lab in labels
    ]  # fmt: skip
    return _typed(pd.DataFrame(rows, columns=list(COLUMNS)))


def intersections(matrix: pd.DataFrame, labels: Sequence[str]) -> pd.DataFrame:
    """The intersection section: every non-empty exact intersection of the detected sets (runs
    joined by ``&`` in label order) and ``none`` (genes detected in no run); sorted by
    ``n_genes`` descending, then name, ``none`` last among equal counts."""
    hit = matrix.assign(_hit=matrix["detected"].fillna(False).astype(bool))
    sets = hit[hit["_hit"]].groupby("locus_tag", sort=False)["run"].agg(set)
    names = [
        "&".join(lab for lab in labels if lab in sets.get(tag, ())) or NONE
        for tag in pd.unique(matrix["locus_tag"])
    ]
    counts = pd.Series(names, dtype=object).value_counts()
    if NONE not in counts.index:
        counts[NONE] = 0
    keyed = sorted(counts.items(), key=lambda kv: (-int(kv[1]), kv[0] == NONE, str(kv[0])))
    rows = [("intersection", None, None, None, None, name, int(n)) for name, n in keyed]
    return _typed(pd.DataFrame(rows, columns=list(COLUMNS)))
