"""Category summary: observed vs expected detection per COG, KEGG and GO-slim category (spec
§7.7.2, §9.3; issue M8-04 §3).

Pure and descriptive (no test). Gene-level: one row per category, whatever its size, with the
eligible members, detected members, expected = sum of p_detect, O - E with a 95 % TU-block
bootstrap CI of the per-TU sums of d - p (M4-16 ``basic_rows``; NA and ``ci_few_units`` below
``CI_MIN_UNITS`` TUs), the median abundance of detected members and the collection's coverage.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from giap.function.collections import GeneSets
from giap.function.ora import eligible_genes
from giap.stats.bootstrap import tu_block_bootstrap

__all__ = ["COLUMNS", "category_summary"]

SUMMARY_COLLECTIONS = ("cog", "kegg", "go_slim")
CI_MIN_UNITS = 10
COLUMNS = (
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
)
_FLOAT = ("expected", "o_minus_e", "ci_low", "ci_high", "median_abundance_log2", "coverage")


def _sum_v(df: pd.DataFrame) -> float:
    return float(df["v"].sum())


def category_summary(
    genes: pd.DataFrame, colls: Sequence[GeneSets], B: int, rng: np.random.Generator
) -> pd.DataFrame:
    """The category_summary rows (§9.3) of the cog, kegg and go_slim collections, sorted by
    collection in that order, then category; ``attrs["coverage"]`` maps every such collection
    passed (empty or not) to its coverage. Bootstraps draw from ``rng`` in row order.

    Raises:
        ValueError: as ``giap.function.ora.eligible_genes``.
    """
    el = eligible_genes(genes)  # natural (TU, locus tag) order, so per-TU sums are in TU order
    pos = {g: i for i, g in enumerate(el["locus_tag"].astype(str))}
    det = el["detected"].astype("boolean").fillna(False).to_numpy(bool)
    p = pd.to_numeric(el["p_detect"], errors="coerce").to_numpy(float)
    tu = el["tu_id"].astype(str).to_numpy()
    abundance = np.full(len(el), math.nan)
    if "abundance_log2" in el:
        abundance = pd.to_numeric(el["abundance_log2"], errors="coerce").to_numpy(float)

    by_name = {c.collection: c for c in colls}
    rows: list[dict[str, Any]] = []
    for coll in (by_name[n] for n in SUMMARY_COLLECTIONS if n in by_name):
        for cat in sorted(coll.sets):
            idx = np.array(sorted(pos[g] for g in coll.sets[cat] if g in pos), dtype=np.int64)
            o, e = int(det[idx].sum()), float(p[idx].sum())
            per_tu = pd.Series(det[idx] - p[idx]).groupby(tu[idx], sort=False).sum()
            lo = hi = math.nan
            flags = []
            if per_tu.size < CI_MIN_UNITS:
                flags.append("ci_few_units")
            else:
                units = pd.DataFrame({"tu": per_tu.index, "v": per_tu.to_numpy()})
                ci = tu_block_bootstrap(units, "tu", _sum_v, B, rng)
                lo, hi = ci.ci_low, ci.ci_high
            a = abundance[idx][det[idx]]
            a = a[np.isfinite(a)]
            rows.append(
                {
                    "collection": coll.collection,
                    "category": cat,
                    "name": coll.names.get(cat, "NA"),
                    "n_genes": int(idx.size),
                    "n_detected": o,
                    "expected": e,
                    "o_minus_e": o - e,
                    "ci_low": lo,
                    "ci_high": hi,
                    "median_abundance_log2": float(np.median(a)) if a.size else math.nan,
                    "coverage": float(coll.coverage),
                    "flags": flags,
                }
            )
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    for c in _FLOAT:
        df[c] = df[c].astype("float64")
    for c in ("n_genes", "n_detected"):
        df[c] = df[c].astype("int64")
    df.attrs["coverage"] = {
        n: float(by_name[n].coverage) for n in SUMMARY_COLLECTIONS if n in by_name
    }
    return df
