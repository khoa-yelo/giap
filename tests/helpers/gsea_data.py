"""Gene frames and gene-set collections for the M8-05 GSEA tests (M8-05 spec).

A gene frame carries the gene-table columns ``gsea_table`` reads: ``locus_tag``, ``tu_id``,
``tu_position``, ``detected``, ``abundance_log2``, ``abundance_imputed``, ``abundance_scale``,
``partial``, ``p_detect`` (eligibility, §7.3: non-partial with a finite ``p_detect``).
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from gippa.function.collections import GeneSets

#: the §9.3 columns of enrichment_gsea.tsv (+ the §9.2 test columns and the ADR-A ``flags``)
GSEA_COLUMNS = (
    "collection",
    "term",
    "name",
    "size_tus",
    "es",
    "nes",
    "log2err",
    "leading_edge",
    "redundant",
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
)


@dataclass(frozen=True)
class Gene:
    tag: str
    value: float | None = None  # abundance_log2 (None = NA)
    detected: bool = True
    imputed: bool = False
    partial: bool = False
    p_detect: float = 0.5


def frame(tus: Mapping[str, Sequence[Gene]], scale: str | None = "molar") -> pd.DataFrame:
    """Gene frame; genes of a TU get ``tu_position`` 1, 2, ... in the given order."""
    rows = []
    for tu, genes in tus.items():
        for pos, g in enumerate(genes, start=1):
            rows.append(
                {
                    "locus_tag": g.tag,
                    "tu_id": tu,
                    "tu_position": pos,
                    "detected": g.detected,
                    "abundance_log2": math.nan if g.value is None else float(g.value),
                    "abundance_imputed": g.imputed,
                    "abundance_scale": scale,
                    "partial": g.partial,
                    "p_detect": g.p_detect,
                }
            )
    df = pd.DataFrame(rows)
    df["abundance_log2"] = df["abundance_log2"].astype("float64")
    df["abundance_scale"] = df["abundance_scale"].astype(object)
    return df


def singleton_frame(values: Sequence[float], prefix: str = "g") -> pd.DataFrame:
    """One detected gene per TU: TU_<i+1> holds ``<prefix><i+1>`` with abundance ``values[i]``."""
    return frame(
        {f"TU_{i + 1}": [Gene(f"{prefix}{i + 1}", float(v))] for i, v in enumerate(values)}
    )


def gene_sets(
    collection: str,
    sets: Mapping[str, Iterable[str]],
    *,
    parents: Mapping[str, Iterable[str]] | None = None,
    coverage: float = 0.5,
) -> GeneSets:
    s = {t: frozenset(g) for t, g in sorted(sets.items())}
    names = {t: f"name of {t}" for t in s}
    par = {t: frozenset((parents or {}).get(t, ())) for t in s} if parents is not None else {}
    return GeneSets(collection, s, names, par, coverage, 0)  # type: ignore[arg-type]


def spread(n: int, rng: np.random.Generator, loc: float = 20.0) -> np.ndarray:
    """``n`` distinct abundances (no ties): a random permutation of loc + 0.01 * (0 .. n-1)."""
    return loc + 0.01 * rng.permutation(n).astype(float)


def as_list(v: object) -> list[str]:
    """A list cell (list, tuple or array) as a list; NA / None as []."""
    if v is None:
        return []
    if isinstance(v, float) and math.isnan(v):
        return []
    if isinstance(v, str):
        return [x for x in v.split(";") if x]
    return [str(x) for x in v]  # type: ignore[attr-defined]
