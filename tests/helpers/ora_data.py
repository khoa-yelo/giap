"""Gene frames, gene sets and exact expectations for the M8-04b ORA tests (docs/specs/M8-04).

A gene frame carries the gene-table columns ``detection_enrichment`` / ``category_summary`` read:
``locus_tag``, ``tu_id``, ``tu_position``, ``detected``, ``p_detect``, ``partial``,
``abundance_log2`` (plus the M8-05 columns ``frame`` writes). The toy universe T6 is the spec's
(§Examples): ten eligible genes in six TUs.

``cpb_raw`` is the expected conditional law computed from the RAW TU q of the frame with the
merged primitive and no tilt (the law is tilt-invariant; spec §1 example 5), so a test compares
the ORA layer to the primitive on the TU sets the spec prescribes.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from giap.stats.poisson_binomial import CpbTails, conditional_poisson_binomial
from helpers.gsea_data import Gene, as_list, frame, gene_sets

__all__ = [
    "ORA_COLUMNS",
    "SUMMARY_COLUMNS",
    "T6",
    "Gene",
    "as_list",
    "cpb_raw",
    "frame",
    "gene_sets",
    "kegg",
    "ora",
    "raw_q",
    "row",
    "singletons",
    "summary",
    "t6_frame",
]

#: §9.3 enrichment_detection.tsv + the §9.2 test columns + the M8-04 additions (spec §4 table)
ORA_COLUMNS = frozenset(
    {
        "collection",
        "term",
        "name",
        "direction",
        "n_universe",
        "n_in_set",
        "n_detected",
        "k_detected_in_set",
        "n_in_set_genes",
        "k_detected_in_set_genes",
        "odds_weight",
        "p_classic",
        "p_noelim",
        "parent_group",
        "coverage",
        "stat",
        "effect",
        "ci_low",
        "ci_high",
        "p_value",
        "q_value",
        "method",
        "null_model",
        "n_perm",
        "flags",
    }
)
#: §9.3 category_summary.tsv + ``flags`` (ci_few_units, M4-16)
SUMMARY_COLUMNS = frozenset(
    {
        "collection",
        "category",
        "name",
        "n_genes",
        "n_detected",
        "expected",
        "o_minus_e",
        "ci_low",
        "ci_high",
        "median_abundance_log2",
        "coverage",
        "flags",
    }
)

#: Toy universe T6 (spec §Examples): TU -> [(gene, p_detect)]; detected g2, g4, g6, g9
T6: dict[str, list[tuple[str, float]]] = {
    "TU_1": [("g1", 0.9), ("g2", 0.5)],
    "TU_2": [("g3", 0.5)],
    "TU_3": [("g4", 0.2), ("g5", 0.2)],
    "TU_4": [("g6", 0.5)],
    "TU_5": [("g7", 0.3), ("g8", 0.3), ("g9", 0.3)],
    "TU_6": [("g10", 0.1)],
}
T6_DETECTED = frozenset({"g2", "g4", "g6", "g9"})


def t6_frame(detected: Iterable[str] = T6_DETECTED) -> pd.DataFrame:
    det = set(detected)
    return frame(
        {
            tu: [Gene(g, 10.0 if g in det else None, g in det, p_detect=p) for g, p in gs]
            for tu, gs in T6.items()
        }
    )


def singletons(
    p: Sequence[float], detected: Iterable[int], prefix: str = "g", values: Sequence[float] = ()
) -> pd.DataFrame:
    """TU_<i+1> = {<prefix><i+1>} with p_detect p[i]; ``detected`` lists 0-based indices."""
    det = set(detected)
    vals = list(values) or [10.0] * len(p)
    return frame(
        {
            f"TU_{i + 1}": [
                Gene(f"{prefix}{i + 1}", vals[i] if i in det else None, i in det, p_detect=float(x))
            ]
            for i, x in enumerate(p)
        }
    )


def ora():  # type: ignore[no-untyped-def]
    from giap.function import ora as mod

    return mod


def summary():  # type: ignore[no-untyped-def]
    from giap.function import summary as mod

    return mod


def raw_q(genes: pd.DataFrame) -> tuple[list[str], np.ndarray, np.ndarray]:
    """(TU ids in first-seen order, raw q = 1 - prod(1 - p), detected) over non-partial genes."""
    el = genes[~genes["partial"].astype(bool)]
    ids: list[str] = []
    q: dict[str, float] = {}
    d: dict[str, bool] = {}
    for tu, p, det in zip(el["tu_id"], el["p_detect"], el["detected"], strict=True):
        if tu not in q:
            ids.append(tu)
            q[tu], d[tu] = 0.0, False
        q[tu] += math.log1p(-float(p))
        d[tu] = d[tu] or bool(det)
    return ids, np.array([-math.expm1(q[t]) for t in ids]), np.array([d[t] for t in ids])


def cpb_raw(genes: pd.DataFrame, tus: Iterable[str]) -> CpbTails:
    """The conditional law of K over the TU set ``tus`` given N, from the raw q (no tilt)."""
    ids, q, d = raw_q(genes)
    m = np.isin(np.array(ids), list(tus))
    return conditional_poisson_binomial(q[m], q[~m], int(d[m].sum()), int(d.sum()))


def row(out: pd.DataFrame, term: str, direction: str, collection: str | None = None) -> pd.Series:
    sel = (out["term"] == term) & (out["direction"] == direction)
    if collection is not None:
        sel &= out["collection"] == collection
    rows = out[sel]
    assert len(rows) == 1, f"expected one {direction} row for {term}, got {len(rows)}"
    return rows.iloc[0]


def kegg(sets: Mapping[str, Iterable[str]], coverage: float = 0.5):  # type: ignore[no-untyped-def]
    return gene_sets("kegg", sets, coverage=coverage)
