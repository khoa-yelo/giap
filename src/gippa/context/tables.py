"""Row assembly of the two M5 context tables and of the genes.tsv context group.

docs/specs/M5-09-context-tables-step.md "Row assembly"; GIPPA_SPEC.md §9.3.

* ``autocorrelation`` = the M5-05 rows (:func:`gippa.context.autocorrelation.autocorrelation_rows`)
  and the M5-08 rows (:func:`gippa.context.oriter.oriter_rows`), sorted by contig in C1 order,
  then ``statistic`` and ``null_model`` in §9.3 enum order, then ``lag`` (NA last), then ``stat``.
* ``islands`` = the M5-06 frame (:func:`gippa.context.islands.islands_frame`) sorted by contig in
  C1 order, then ``start``, ``end`` and the ``ISL_<n>`` number.
* ``gene_context`` = the M5-06 gene→island map and the M5-08 replichore columns, left-joined
  onto the gene table, one row per gene in gene-table order.

Pure functions: no randomness, no file I/O, no configuration.
"""

from __future__ import annotations

__all__ = [
    "GENE_CONTEXT_COLUMNS",
    "NULL_MODEL_ORDER",
    "STATISTIC_ORDER",
    "autocorrelation_table",
    "gene_context_table",
    "islands_table",
]

import math
from collections.abc import Sequence

import pandas as pd

from gippa.context.autocorrelation import AUTOCORR_COLUMNS
from gippa.context.islands import ISLAND_COLUMNS, gene_island_map
from gippa.context.oriter import gene_replichores

#: §9.3 ``statistic`` enum order (the sort key, not the producers' row order).
STATISTIC_ORDER: tuple[str, ...] = ("morans_i", "join_count_bb", "oriter_slope")
#: §9.3 ``null_model`` enum order.
NULL_MODEL_ORDER: tuple[str, ...] = ("gene_label", "tu_block", "bootstrap")
#: §9.3 genes.tsv context group, as ``work/context/gene_context.parquet`` holds it.
GENE_CONTEXT_COLUMNS: tuple[str, ...] = (
    "locus_tag", "island_id", "replichore", "replichore_position",
)  # fmt: skip


def _contig_rank(replicons: pd.DataFrame) -> dict[str, int]:
    """Position of each replicon in C1 order; an unknown contig sorts last, by name."""
    return {str(c): i for i, c in enumerate(replicons["contig"])}


def _rank(order: Sequence[str], value: object) -> int:
    """Index of ``value`` in an enum order; an unknown value sorts last."""
    text = str(value)
    return order.index(text) if text in order else len(order)


def _island_number(island_id: object) -> int:
    """``n`` of an ``ISL_<n>`` ID; a non-conforming ID sorts last."""
    tail = str(island_id).rsplit("_", 1)[-1]
    return int(tail) if tail.isdigit() else 2**62


def _sorted(out: pd.DataFrame, key: dict[str, list[object]]) -> pd.DataFrame:
    """``out`` reordered by a stable sort on ``key`` (ties keep the producer's order)."""
    order = pd.DataFrame(key).sort_values(list(key), kind="mergesort").index
    return out.loc[order].reset_index(drop=True)


def autocorrelation_table(parts: Sequence[pd.DataFrame], replicons: pd.DataFrame) -> pd.DataFrame:
    """The ``AUTOCORR_COLUMNS`` frames of ``parts`` (M5-05 rows, M5-08 ori-ter rows)
    concatenated and sorted as §9.3 orders the table: contig in ``replicons`` order, then
    ``statistic`` (:data:`STATISTIC_ORDER`), ``null_model`` (:data:`NULL_MODEL_ORDER`),
    ``lag`` (NA last) and ``stat``."""
    frames = [p for p in parts if len(p)]
    if not frames:
        return pd.DataFrame(columns=list(AUTOCORR_COLUMNS))
    out = pd.concat(frames, ignore_index=True)[list(AUTOCORR_COLUMNS)]
    rank = _contig_rank(replicons)
    return _sorted(
        out,
        {
            "contig": [rank.get(str(v), len(rank)) for v in out["contig"]],
            "statistic": [_rank(STATISTIC_ORDER, v) for v in out["statistic"]],
            "null_model": [_rank(NULL_MODEL_ORDER, v) for v in out["null_model"]],
            "lag": [math.inf if pd.isna(v) else float(v) for v in out["lag"]],
            "stat": [str(v) for v in out["stat"]],
        },
    )


def islands_table(islands: pd.DataFrame, replicons: pd.DataFrame) -> pd.DataFrame:
    """The 0-based :func:`gippa.context.islands.islands_frame` frame sorted as §9.3 orders the
    table: contig in ``replicons`` order, then ``start``, ``end`` and the ``ISL_<n>`` number.
    The producer's IDs are kept, not renumbered."""
    out = islands[list(ISLAND_COLUMNS)]
    if not len(out):
        return out.reset_index(drop=True)
    rank = _contig_rank(replicons)
    return _sorted(
        out,
        {
            "contig": [rank.get(str(v), len(rank)) for v in out["contig"]],
            "start": [math.inf if pd.isna(v) else float(v) for v in out["start"]],
            "end": [math.inf if pd.isna(v) else float(v) for v in out["end"]],
            "island_id": [_island_number(v) for v in out["island_id"]],
        },
    )


def gene_context_table(
    genes: pd.DataFrame,
    islands: pd.DataFrame,
    replicons: pd.DataFrame,
    features: pd.DataFrame | None,
) -> pd.DataFrame:
    """:data:`GENE_CONTEXT_COLUMNS`, one row per gene of ``genes`` in gene-table order: the
    M5-06 gene->island map (step 9) of ``islands`` and the M5-08 replichore columns of
    ``replicons`` / the C15 ``features``, left-joined and NA where a gene has neither.
    ``features`` None leaves both replichore columns NA."""
    tags = genes["locus_tag"].astype(str)
    island_of = gene_island_map(islands).set_index("locus_tag")["island_id"]
    replichores = gene_replichores(genes, replicons, features).set_index("locus_tag")
    return pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(),
            "island_id": tags.map(island_of).to_numpy(object),
            "replichore": tags.map(replichores["replichore"]).to_numpy(object),
            "replichore_position": tags.map(replichores["replichore_position"]).to_numpy(float),
        },
        columns=list(GENE_CONTEXT_COLUMNS),
    )
