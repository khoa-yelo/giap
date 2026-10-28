"""GSEA running sums of the top sets, for V14(a) (docs/specs/M9-17-plot-support-analysis.md).

Pure. The top ``GSEA_TOP`` rows of ``enrichment_gsea`` with a p-value (smallest ``q_value``,
then largest |``nes``|, then ``collection``, ``term``) are re-walked on the ranking the engine
saw (``gsea_table``'s ``attrs["rank_labels"]``, ranked as ``gippa.stats.fgsea.prepare`` does) by
``gippa.stats.fgsea.running_es``; a walk whose max deviation is not the row's ``es`` (1e-9)
raises.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from gippa.function.collections import GeneSets
from gippa.function.gsea import tu_ranking
from gippa.stats.fgsea import running_es

__all__ = ["COLUMNS", "GSEA_TOP", "gsea_running"]

GSEA_TOP = 6
COLUMNS = ("collection", "term", "rank_position", "tu_id", "score", "in_set", "running_es")
ES_TOL = 1e-9


def _top(gsea: pd.DataFrame, k: int) -> list[tuple[str, str]]:
    t = gsea[gsea["p_value"].notna()]
    key = pd.DataFrame({"q": t["q_value"].astype(float), "a": -t["nes"].astype(float).abs(),
                        "c": t["collection"].astype(str), "t": t["term"].astype(str)})  # fmt: skip
    key = key.sort_values(["q", "a", "c", "t"], kind="stable").head(k)
    return list(zip(key["c"], key["t"], strict=True))


def _max_dev(walk: np.ndarray) -> tuple[float, bool]:
    """(the walk's max deviation from 0, whether it sits on the ``maxP == -minP`` tie)."""
    hi, lo = max(0.0, float(walk.max())), min(0.0, float(walk.min()))
    return (hi if hi > -lo else lo), abs(hi + lo) <= ES_TOL


def gsea_running(
    genes: pd.DataFrame, colls: Sequence[GeneSets], gsea: pd.DataFrame, k: int = GSEA_TOP
) -> pd.DataFrame | None:
    """``gsea_running`` rows (``COLUMNS``) of the top ``k`` sets; None with no tested set.

    Raises:
        ValueError: a running sum whose max deviation is not the table's ``es``.
    """
    labels = gsea.attrs.get("rank_labels")
    top = _top(gsea, k) if labels is not None else []
    if not top:
        return None
    rank = tu_ranking(genes)
    score = dict(zip(rank.tu_ids, rank.scores.tolist(), strict=True))
    stats = np.array([score[t] for t in labels], dtype=np.float64)
    order = np.argsort(-stats, kind="stable")
    ranked = [str(labels[int(i)]) for i in order]
    scores = stats[order]
    pos = {t: i + 1 for i, t in enumerate(ranked)}
    tu_of = {g: tu for tu, gs in rank.tu_genes.items() for g in gs}
    sets: dict[str, Mapping[str, frozenset[str]]] = {str(c.collection): c.sets for c in colls}
    es_of = {(str(c), str(t)): float(e) for c, t, e in zip(gsea["collection"], gsea["term"],
                                                            gsea["es"], strict=True)}  # fmt: skip
    parts = []
    for coll, term in top:
        hits = np.array(sorted({pos[tu_of[g]] for g in sets[coll][term] if g in tu_of}))
        walk = running_es(scores, hits.astype(np.int64))
        dev, tie = _max_dev(walk)
        es = es_of[(coll, term)]
        near_tie = tie and (es == 0.0 or abs(abs(dev) - abs(es)) <= ES_TOL)
        if abs(dev - es) > ES_TOL and not near_tie:
            raise ValueError(f"gsea_running {coll} {term}: running sum deviates {dev}, es {es}")
        in_set = np.zeros(len(ranked), dtype=bool)
        in_set[hits - 1] = True
        parts.append(pd.DataFrame({
            "collection": coll, "term": term, "rank_position": np.arange(1, len(ranked) + 1),
            "tu_id": ranked, "score": scores, "in_set": in_set, "running_es": walk,
        }))  # fmt: skip
    return pd.concat(parts, ignore_index=True)[list(COLUMNS)]
