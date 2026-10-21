"""Row assembly of the two M5 context tables and of the genes.tsv context group.

docs/specs/M5-09-context-tables-step.md "Row assembly"; GIAP_SPEC.md §9.3.

* ``autocorrelation`` = the M5-05 rows (:func:`giap.context.autocorrelation.autocorrelation_rows`)
  and the M5-08 rows (:func:`giap.context.oriter.oriter_rows`), sorted by contig in C1 order,
  then ``statistic`` and ``null_model`` in §9.3 enum order, then ``lag`` (NA last), then ``stat``.
* ``islands`` = the M5-06 frame (:func:`giap.context.islands.islands_frame`) sorted by contig in
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

from giap.context.autocorrelation import AUTOCORR_COLUMNS
from giap.context.islands import ISLAND_COLUMNS, gene_island_map
from giap.context.oriter import gene_replichores

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


def autocorrelation_table(parts: Sequence[pd.DataFrame], replicons: pd.DataFrame) -> pd.DataFrame:
    """Concatenate the autocorrelation row sources and sort them as §9.3 orders the table.

    Parameters
    ----------
    parts : sequence of pandas.DataFrame
        Frames with the columns ``giap.context.autocorrelation.AUTOCORR_COLUMNS`` (the M5-05
        rows and the M5-08 ori-ter rows).
    replicons : pandas.DataFrame
        C1 replicons in table order; their order is the ``contig`` sort key.

    Returns
    -------
    pandas.DataFrame
        ``AUTOCORR_COLUMNS``, sorted by contig (C1 order), ``statistic``
        (:data:`STATISTIC_ORDER`), ``null_model`` (:data:`NULL_MODEL_ORDER`), ``lag`` (NA last)
        and ``stat``, with a fresh range index. The sort is stable, so rows that tie keep the
        order their producer gave them.
    """
    frames = [p for p in parts if len(p)]
    if not frames:
        return pd.DataFrame(columns=list(AUTOCORR_COLUMNS))
    out = pd.concat(frames, ignore_index=True)[list(AUTOCORR_COLUMNS)]
    rank = _contig_rank(replicons)
    unknown = len(rank)
    key = pd.DataFrame(
        {
            "contig": [rank.get(str(v), unknown) for v in out["contig"]],
            "statistic": [_rank(STATISTIC_ORDER, v) for v in out["statistic"]],
            "null_model": [_rank(NULL_MODEL_ORDER, v) for v in out["null_model"]],
            "lag": [math.inf if pd.isna(v) else float(v) for v in out["lag"]],
            "stat": [str(v) for v in out["stat"]],
        }
    )
    order = key.sort_values(list(key.columns), kind="mergesort").index
    return out.loc[order].reset_index(drop=True)


def islands_table(islands: pd.DataFrame, replicons: pd.DataFrame) -> pd.DataFrame:
    """Sort the M5-06 island frame as §9.3 orders the table.

    Parameters
    ----------
    islands : pandas.DataFrame
        The frame of :func:`giap.context.islands.islands_frame`
        (``giap.context.islands.ISLAND_COLUMNS``, 0-based half-open).
    replicons : pandas.DataFrame
        C1 replicons in table order.

    Returns
    -------
    pandas.DataFrame
        ``ISLAND_COLUMNS``, sorted by contig (C1 order), ``start``, ``end`` and the ``ISL_<n>``
        number, with a fresh range index. The island IDs are the producer's and are not
        renumbered.
    """
    out = islands[list(ISLAND_COLUMNS)]
    if not len(out):
        return out.reset_index(drop=True)
    rank = _contig_rank(replicons)
    unknown = len(rank)
    key = pd.DataFrame(
        {
            "contig": [rank.get(str(v), unknown) for v in out["contig"]],
            "start": [math.inf if pd.isna(v) else float(v) for v in out["start"]],
            "end": [math.inf if pd.isna(v) else float(v) for v in out["end"]],
            "island_id": [_island_number(v) for v in out["island_id"]],
        }
    )
    order = key.sort_values(list(key.columns), kind="mergesort").index
    return out.loc[order].reset_index(drop=True)


def gene_context_table(
    genes: pd.DataFrame,
    islands: pd.DataFrame,
    replicons: pd.DataFrame,
    features: pd.DataFrame | None,
) -> pd.DataFrame:
    """The §9.3 genes.tsv context group, one row per gene in gene-table order.

    Parameters
    ----------
    genes : pandas.DataFrame
        The gene table (``locus_tag``, ``contig``, ``start``, ``end``).
    islands : pandas.DataFrame
        The island table; its ``genes`` lists give the gene→island map (M5-06 step 9).
    replicons : pandas.DataFrame
        C1 replicons.
    features : pandas.DataFrame or None
        C15 features; its ``oric`` rows give each replicon's origin (M5-08). None leaves both
        replichore columns NA.

    Returns
    -------
    pandas.DataFrame
        :data:`GENE_CONTEXT_COLUMNS`: ``locus_tag``, ``island_id`` (None off every island),
        ``replichore`` (None off the ori-ter analysis) and ``replichore_position`` (NaN there).
    """
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
