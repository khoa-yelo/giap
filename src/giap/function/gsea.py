"""Abundance GSEA on a TU-collapsed ranking (spec §7.7.4, §8.1, §8.2; issue M8-05).

Pure. The ranked units are TUs: a TU's score is the mean ``abundance_log2`` of its detected,
eligible genes (imputed values included), median-centred (ADR-B), and a TU is in a set when any
of its scoring genes is. Sets are sized in TUs. Per collection, M8-02's ``fgsea_simple`` runs on
the TUs in a seeded random order (M8-02 lead decision 7: ties are never broken by TU ID); q is BH
within the collection. A significant GO parent is ``redundant`` when, re-tested without the genes
of its significant descendants, its p substituted into the collection's BH vector gives q > alpha
(lead decisions 3, 6 and 7: a re-test with 0 TUs is redundant, an unestimated re-test leaves
``redundant`` NA). The writer joins the list cells ``leading_edge`` and ``flags``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config.model import GseaConfig
from giap.function.collections import GeneSets
from giap.stats.fgsea import Prepared, SimpleResult, fgsea_simple, prepare
from giap.stats.multitest import bh

__all__ = ["COLUMNS", "TuRanking", "gsea_table", "tu_ranking", "tu_sets"]

F64 = npt.NDArray[np.float64]
GO_NAMESPACES = frozenset({"go_bp", "go_mf", "go_cc"})
METHOD = "fgsea_simple"
NULL_MODEL = "tu_set_permutation"
COLUMNS = (
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
_FLOAT = ("es", "nes", "log2err", "coverage", "stat", "effect", "ci_low", "ci_high",
          "p_value", "q_value")  # fmt: skip
_INT = ("size_tus", "n_perm")


@dataclass(frozen=True, slots=True)
class TuRanking:
    """Ranked TUs in natural TU-ID order, their median-centred scores and scoring genes."""

    tu_ids: tuple[str, ...]
    scores: F64
    tu_genes: Mapping[str, tuple[str, ...]]  # scoring genes per ranked TU, tu_position order


def _natural_key(s: str) -> tuple[tuple[int, int | str], ...]:
    return tuple((0, int(p)) if p.isdigit() else (1, p) for p in re.split(r"(\d+)", s) if p)


def _scoring(genes: pd.DataFrame) -> pd.Series:
    """Detected, eligible (non-partial, finite ``p_detect``) genes with a non-NA abundance."""
    ok = genes["detected"].fillna(False).astype(bool) & genes["abundance_log2"].notna()
    if "partial" in genes:
        ok &= ~genes["partial"].fillna(False).astype(bool)
    if "p_detect" in genes:
        ok &= np.isfinite(pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float))
    return ok


def tu_ranking(genes: pd.DataFrame) -> TuRanking:
    """TU scores = mean ``abundance_log2`` of the scoring genes, minus the median over TUs."""
    sub = genes.loc[_scoring(genes), ["tu_id", "tu_position", "locus_tag", "abundance_log2"]]
    rows = sorted(
        zip(sub["tu_id"].astype(str), sub["tu_position"], sub["locus_tag"].astype(str),
            sub["abundance_log2"].astype(float), strict=True),
        key=lambda r: (_natural_key(r[0]), r[1], r[2]),
    )  # fmt: skip
    tags: dict[str, list[str]] = {}
    vals: dict[str, list[float]] = {}
    for tu, _, tag, v in rows:
        tags.setdefault(tu, []).append(tag)
        vals.setdefault(tu, []).append(v)
    ids = tuple(tags)
    means = np.array([np.mean(vals[t]) for t in ids], dtype=np.float64)
    scores = means - np.median(means) if ids else means
    return TuRanking(ids, scores, {t: tuple(g) for t, g in tags.items()})


def _tus_of(rank: TuRanking) -> dict[str, str]:
    return {g: tu for tu, gs in rank.tu_genes.items() for g in gs}


def _members(tu_of: Mapping[str, str], genes: Iterable[str]) -> frozenset[str]:
    return frozenset(tu_of[g] for g in genes if g in tu_of)


def tu_sets(rank: TuRanking, coll: GeneSets) -> dict[str, frozenset[str]]:
    """Per term, the ranked TUs with at least one scoring gene in the term's gene set."""
    tu_of = _tus_of(rank)
    return {t: _members(tu_of, gs) for t, gs in coll.sets.items()}


def _no_abundance(genes: pd.DataFrame) -> bool:
    if "abundance_scale" not in genes or "abundance_log2" not in genes:
        return True
    return bool(genes["abundance_scale"].isna().all() and genes["abundance_log2"].isna().all())


def _frame(rows: list[dict[str, Any]], run_flags: list[str]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(COLUMNS))
    for c in _FLOAT:
        df[c] = df[c].astype("float64")
    for c in _INT:
        df[c] = df[c].astype("int64")
    df["redundant"] = df["redundant"].astype("boolean")
    df.attrs["run_flags"] = list(run_flags)
    return df


def _closure(step: Mapping[str, Iterable[str]]) -> Callable[[str], frozenset[str]]:
    """Transitive closure (excluding the start) of an acyclic successor map."""
    memo: dict[str, frozenset[str]] = {}

    def walk(t: str) -> frozenset[str]:
        if t not in memo:
            out: set[str] = set()
            for s in step.get(t, ()):
                out |= {s} | walk(s)
            memo[t] = frozenset(out)
        return memo[t]

    return walk


def _depth(parents: Mapping[str, Iterable[str]]) -> Callable[[str], int]:
    """Shortest number of ``parents`` steps up to a term without parents."""
    memo: dict[str, int] = {}

    def depth(t: str) -> int:
        if t not in memo:
            ps = list(parents.get(t, ()))
            memo[t] = 1 + min(depth(p) for p in ps) if ps else 0
        return memo[t]

    return depth


@dataclass(frozen=True, slots=True)
class _Main:
    coll: GeneSets
    prep: Prepared
    res: SimpleResult | None
    q: F64


_Ctx = tuple[F64, list[str], dict[str, str], GseaConfig, np.random.Generator]


def _retest(
    m: _Main,
    term: str,
    removed: frozenset[str],
    ctx: _Ctx,
) -> tuple[float, int]:
    """Re-test ``term`` without ``removed`` genes (min size 1): (p, size in TUs). An empty set
    is not tested (p NaN, size 0); an unestimated re-test has p NaN."""
    stats, labels, tu_of, cfg, rng = ctx
    tus = _members(tu_of, m.coll.sets[term] - removed)
    if not tus:
        return math.nan, 0
    prep = prepare(stats, labels, {term: sorted(tus)}, 1, len(labels))
    res = fgsea_simple(prep, cfg.n_perm_simple, cfg.sample_size, rng)
    return float(res.p[0]), len(tus)


def _go_hierarchy(
    m: _Main, alpha: float, min_size: int, ctx: _Ctx
) -> tuple[dict[str, Any], set[str], dict[str, str]]:
    """``redundant`` (True / False / NA), the ``size_filtered`` terms and ``parent_group`` of a
    GO collection's significant terms; re-tests run in term order."""
    assert m.res is not None
    coll, p, names = m.coll, m.res.p, m.prep.names
    sig = {t for t, q in zip(names, m.q, strict=True) if q <= alpha}
    redundant: dict[str, Any] = dict.fromkeys(names, False)
    size_filtered: set[str] = set()
    group: dict[str, str] = {}
    if coll.collection not in GO_NAMESPACES or not sig:
        return redundant, size_filtered, group
    children: dict[str, set[str]] = {}
    for t, ps in coll.parents.items():
        for parent in ps:
            children.setdefault(parent, set()).add(t)
    below, above, depth = _closure(children), _closure(coll.parents), _depth(coll.parents)
    for i, t in enumerate(names):
        if t not in sig:
            continue
        group[t] = min({t} | (above(t) & sig), key=lambda a: (depth(a), a))
        desc = below(t) & sig
        if not desc:
            continue
        p_re, k = _retest(m, t, frozenset().union(*(coll.sets[d] for d in desc)), ctx)
        if k < min_size:
            size_filtered.add(t)
        if k == 0:  # the parent is the union of its significant descendants
            redundant[t] = True
        elif math.isnan(p_re):  # lead decision 7.4: undecided
            redundant[t] = pd.NA
        else:
            sub = np.array(p, dtype=np.float64)
            sub[i] = p_re
            redundant[t] = bool(bh(sub)[i] > alpha)
    return redundant, size_filtered, group


def gsea_table(
    genes: pd.DataFrame,
    colls: Sequence[GeneSets],
    size: tuple[int, int],
    cfg: GseaConfig,
    alpha: float,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """The enrichment_gsea rows of every collection (§9.3); run flags in ``attrs["run_flags"]``."""
    if _no_abundance(genes):
        return _frame([], ["layer_missing:abundance"])
    rank = tu_ranking(genes)
    n = len(rank.tu_ids)
    if n < 2:
        return _frame([], ["too_few_ranked_tus"])
    perm = rng.permutation(n)  # seeded tie order (M8-02 lead decision 7)
    labels = [rank.tu_ids[int(i)] for i in perm]
    stats = np.asarray(rank.scores, dtype=np.float64)[perm]
    tu_of = _tus_of(rank)

    mains: list[_Main] = []
    for coll in colls:
        sets = {t: sorted(s) for t, s in tu_sets(rank, coll).items()}
        prep = prepare(stats, labels, sets, size[0], size[1])
        res = fgsea_simple(prep, cfg.n_perm_simple, cfg.sample_size, rng) if prep.names else None
        q = bh(res.p) if res is not None else np.empty(0)
        mains.append(_Main(coll, prep, res, q))

    ctx = (stats, labels, tu_of, cfg, rng)
    rows: list[dict[str, Any]] = []
    for m in mains:
        if m.res is None:
            continue
        coll, res, names = m.coll, m.res, m.prep.names
        redundant, size_filtered, group = _go_hierarchy(m, alpha, size[0], ctx)
        for i, t in enumerate(names):
            p = float(res.p[i])
            es = float(res.es[i])
            edge_tus = [labels[int(j)] for j in res.leading_edge[i]] if es != 0 else []
            gs = coll.sets[t]
            edge = [g for tu in edge_tus for g in rank.tu_genes[tu] if g in gs]
            flags = ["ci_unestimated"]
            if math.isnan(p):
                flags.append("unestimated")
            if t in size_filtered:
                flags.append("size_filtered")
            rows.append(
                {
                    "collection": coll.collection,
                    "term": t,
                    "name": coll.names.get(t, "NA"),
                    "size_tus": int(m.prep.sizes[i]),
                    "es": es,
                    "nes": float(res.nes[i]),
                    "log2err": float(res.log2err[i]),
                    "leading_edge": edge,
                    "redundant": redundant[t],
                    "parent_group": group.get(t),
                    "coverage": float(coll.coverage),
                    "stat": es,
                    "effect": float(res.nes[i]),
                    "ci_low": math.nan,
                    "ci_high": math.nan,
                    "p_value": p,
                    "q_value": float(m.q[i]),
                    "method": METHOD,
                    "null_model": NULL_MODEL,
                    "n_perm": int(cfg.n_perm_simple),
                    "flags": flags,
                }
            )
    return _frame(rows, [])
