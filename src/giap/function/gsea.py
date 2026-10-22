"""Abundance GSEA on a TU-collapsed ranking (spec §7.7.4, §8.1, §8.2; issue M8-05).

Pure. The ranked units are TUs: a TU's score is the mean ``abundance_log2`` of its detected,
eligible genes (imputed values included), median-centred (ADR-B), and a TU is in a set when any
of its scoring genes is. Sets are sized in TUs. Per collection, M8-02's ``fgsea_simple`` runs on
the TUs in a seeded random order (M8-02 lead decision 7: ties are never broken by TU ID), then
M8-11's ``fgsea_multilevel`` refines the sets it flags (``method`` fgsea_multilevel on every row;
a finite p with NaN ``log2err`` is flagged ``log2err_unestimated``); q is BH within the
collection. A significant GO parent is ``redundant`` when, re-tested without the genes
of its significant descendants, its p substituted into the collection's BH vector gives q > alpha
(lead decisions 3, 6 and 7: a re-test with 0 TUs is redundant, an unestimated re-test leaves
``redundant`` NA); re-tests run the same two stages (M8-11 lead decision 3). The writer joins
the list cells ``leading_edge`` and ``flags``.

Every main-run row also carries the position flags of lead decision 9 (ADR-0025, ADR-0026), both
pure and RNG-free: ``genomic_cluster`` when the set's AR(1) variance inflation delta along the
ranked TU order of each replicon is >= ``CLUSTER_TAU``, and ``detection_confounded`` when the
hypergeometric z of the set's detected gene count is >= ``DETECTION_Z0`` in size with the sign
that would produce the row's ES. Redundancy re-tests compute no flags.
"""

from __future__ import annotations

import math
import re
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, TypeAlias

import numpy as np
import numpy.typing as npt
import pandas as pd

from giap.config.model import GseaConfig
from giap.function.collections import GeneSets
from giap.stats.fgsea import Prepared, SimpleResult, fgsea_multilevel, fgsea_simple, prepare
from giap.stats.multitest import bh

__all__ = [
    "CLUSTER_TAU",
    "COLUMNS",
    "DETECTION_Z0",
    "TuRanking",
    "cluster_delta",
    "cluster_deltas",
    "detection_z",
    "gsea_table",
    "is_detection_confounded",
    "is_genomic_cluster",
    "tu_ranking",
    "tu_sets",
]

F64 = npt.NDArray[np.float64]
GO_NAMESPACES = frozenset({"go_bp", "go_mf", "go_cc"})
METHOD = "fgsea_multilevel"
NULL_MODEL = "tu_set_permutation"
CLUSTER_TAU = 0.15  # decision 9a: genomic_cluster when delta >= tau
DETECTION_Z0 = 1.0  # decision 9b: detection_confounded needs |z| >= z0
_PHI0 = 0.5  # AR(1) kernel w(d) = phi0**d for 1 <= d <= _DMAX
_DMAX = 10
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


def _eligible(genes: pd.DataFrame) -> pd.Series:
    """Genes the ranking could use if detected: non-partial, finite ``p_detect`` (the ORA
    universe of §7.7.3)."""
    ok = pd.Series(True, index=genes.index)
    if "partial" in genes:
        ok &= ~genes["partial"].fillna(False).astype(bool)
    if "p_detect" in genes:
        ok &= np.isfinite(pd.to_numeric(genes["p_detect"], errors="coerce").to_numpy(float))
    return ok


def _scoring(genes: pd.DataFrame) -> pd.Series:
    """Detected, eligible (non-partial, finite ``p_detect``) genes with a non-NA abundance."""
    detected = genes["detected"].fillna(False).astype(bool) & genes["abundance_log2"].notna()
    return detected & _eligible(genes)


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


def _weights(d: npt.NDArray[np.int64]) -> F64:
    return np.where((d >= 1) & (d <= _DMAX), _PHI0 ** np.minimum(d, _DMAX), 0.0)


def _pair_weight(pos: Sequence[int], n_r: int, circular: bool) -> float:
    """Sum of w(d) over the member pairs of one replicon. Sorted positions are distinct, so a
    pair more than ``_DMAX`` members apart is more than ``_DMAX`` TUs apart; on a ring the pairs
    across the seam are taken by their wrap-around lag."""
    p = np.sort(np.asarray(pos, dtype=np.int64))
    m = p.size
    total = 0.0
    for lag in range(1, min(_DMAX, m - 1) + 1):
        d = p[lag:] - p[:-lag]
        if circular:  # a straight gap > _DMAX is left to the seam pass below
            d = np.where(d <= _DMAX, np.minimum(d, n_r - d), 0)
        total += float(_weights(d).sum())
        if circular:  # seam pairs (p[m - lag + i], p[i]) whose straight gap is > _DMAX
            wrap = p[:lag] + n_r - p[m - lag :]
            total += float(_weights(np.where(n_r - wrap > _DMAX, wrap, 0)).sum())
    return total


def _all_pairs_weight(n_r: int, circular: bool) -> float:
    """P_w of one replicon in closed form: w(d) summed over every pair of its ranked TUs."""
    if circular:
        return sum((n_r if 2 * d < n_r else n_r / 2) * _PHI0**d
                   for d in range(1, min(_DMAX, n_r // 2) + 1))  # fmt: skip
    return sum((n_r - d) * _PHI0**d for d in range(1, min(_DMAX, n_r - 1) + 1))


def cluster_delta(
    positions: Mapping[str, Sequence[int]],
    n_ranked: Mapping[str, int],
    circular: Mapping[str, bool],
) -> float:
    """Decision 9a: delta = (C - E[C]) / (1 - k/n) of the members at ``positions`` (0-based
    indices in each replicon's ranked TU order). ``n_ranked`` holds every replicon with ranked
    TUs; a replicon missing from ``circular`` is linear. delta = 0 when k < 2 or k >= n."""
    k = sum(len(v) for v in positions.values())
    n = sum(n_ranked.values())
    if k < 2 or k >= n:
        return 0.0
    s = sum(_pair_weight(p, n_ranked[r], circular.get(r, False)) for r, p in positions.items())
    pw = sum(_all_pairs_weight(nr, circular.get(r, False)) for r, nr in n_ranked.items())
    c = 2.0 * s / k
    ec = 2.0 * (k - 1) * pw / (n * (n - 1))
    return float((c - ec) / (1.0 - k / n))


Replicons: TypeAlias = Mapping[str, str] | pd.DataFrame | None


def _topologies(replicons: Replicons) -> dict[str, str]:
    if replicons is None:
        return {}
    if isinstance(replicons, pd.DataFrame):
        return dict(zip(replicons["contig"].astype(str), replicons["topology"].astype(str),
                        strict=True))  # fmt: skip
    return {str(c): str(t) for c, t in replicons.items()}


def cluster_deltas(
    rank: TuRanking,
    genes: pd.DataFrame,
    replicons: Replicons,
    sets: Mapping[str, Iterable[str]],
) -> dict[str, float]:
    """delta per term from the member TUs in ``sets``. The ranked TUs are split by the
    ``contig`` of their scoring genes (one linear replicon without that column) and keep their
    natural TU-ID order; only ``"circular"`` topologies (C1 or a contig -> topology map) wrap."""
    contig_of: dict[str, str] = {}
    if "contig" in genes:
        contig_of = dict(zip(genes["locus_tag"].astype(str), genes["contig"].astype(str),
                             strict=True))  # fmt: skip
    where: dict[str, tuple[str, int]] = {}
    n_ranked: dict[str, int] = {}
    for tu in rank.tu_ids:
        r = contig_of.get(rank.tu_genes[tu][0], "")
        where[tu] = (r, n_ranked.get(r, 0))
        n_ranked[r] = n_ranked.get(r, 0) + 1
    topo = _topologies(replicons)
    circular = {r: topo.get(r) == "circular" for r in n_ranked}
    out: dict[str, float] = {}
    for term, tus in sets.items():
        positions: dict[str, list[int]] = {}
        for tu in tus:
            if tu in where:
                r, i = where[tu]
                positions.setdefault(r, []).append(i)
        out[term] = cluster_delta(positions, n_ranked, circular)
    return out


def _hypergeometric_z(M: int, N: int, n: int, x: int) -> float:
    if M <= 1:
        return 0.0
    f = N / M
    v = n * f * (1 - f) * (M - n) / (M - 1)
    return 0.0 if v <= 0 else (x - n * f) / math.sqrt(v)


@dataclass(frozen=True, slots=True)
class _Detection:
    eligible: frozenset[str]
    scoring: frozenset[str]

    @classmethod
    def of(cls, genes: pd.DataFrame) -> _Detection:
        tags = genes["locus_tag"].astype(str)
        return cls(frozenset(tags[_eligible(genes)]), frozenset(tags[_scoring(genes)]))

    def z(self, members: Iterable[str]) -> float:
        s = set(members)
        M, N = len(self.eligible), len(self.scoring)
        return _hypergeometric_z(M, N, len(s & self.eligible), len(s & self.scoring))


def detection_z(genes: pd.DataFrame, members: Iterable[str]) -> float:
    """Decision 9b: hypergeometric z of the set's scoring-gene count among its eligible genes,
    against the eligible universe; 0 when the variance is 0."""
    return _Detection.of(genes).z(members)


def is_genomic_cluster(delta: float) -> bool:
    return bool(delta >= CLUSTER_TAU)


def is_detection_confounded(z: float, es: float) -> bool:
    """|z| >= z0 with z * ES < 0; for an ES that is NaN or 0, |z| >= z0 alone."""
    if math.isnan(es) or es == 0:
        return bool(abs(z) >= DETECTION_Z0)
    return bool(z * es < 0 and abs(z) >= DETECTION_Z0)


@dataclass(frozen=True, slots=True)
class _Main:
    coll: GeneSets
    tus: dict[str, frozenset[str]]
    prep: Prepared
    res: SimpleResult | None
    q: F64


_Ctx = tuple[F64, list[str], dict[str, str], GseaConfig, np.random.Generator]


def _engine(prep: Prepared, cfg: GseaConfig, rng: np.random.Generator) -> SimpleResult:
    """``fgsea_simple``, then ``fgsea_multilevel`` on the same ``prep`` (spec M8-11 item 6)."""
    simple = fgsea_simple(prep, cfg.n_perm_simple, cfg.sample_size, rng)
    return fgsea_multilevel(prep, simple, cfg.sample_size, cfg.eps, rng)


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
    return float(_engine(prep, cfg, rng).p[0]), len(tus)


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
        anc = above(t) & sig  # the term itself only when no ancestor is significant
        group[t] = min(anc, key=lambda a: (depth(a), a)) if anc else t
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
    replicons: Replicons = None,
) -> pd.DataFrame:
    """The enrichment_gsea rows of every collection (§9.3); run flags in ``attrs["run_flags"]``.
    ``replicons`` (C1 frame or contig -> topology) marks circular contigs; others are linear."""
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
        tus = tu_sets(rank, coll)
        prep = prepare(stats, labels, {t: sorted(s) for t, s in tus.items()}, size[0], size[1])
        res = _engine(prep, cfg, rng) if prep.names else None
        q = bh(res.p) if res is not None else np.empty(0)
        mains.append(_Main(coll, tus, prep, res, q))

    ctx = (stats, labels, tu_of, cfg, rng)
    detection = _Detection.of(genes)
    rows: list[dict[str, Any]] = []
    for m in mains:
        if m.res is None:
            continue
        coll, res, names = m.coll, m.res, m.prep.names
        redundant, size_filtered, group = _go_hierarchy(m, alpha, size[0], ctx)
        deltas = cluster_deltas(rank, genes, replicons, {t: m.tus[t] for t in names})
        for i, t in enumerate(names):
            p = float(res.p[i])
            es = float(res.es[i])
            edge_tus = [labels[int(j)] for j in res.leading_edge[i]] if es != 0 else []
            gs = coll.sets[t]
            edge = [g for tu in edge_tus for g in rank.tu_genes[tu] if g in gs]
            flags = ["ci_unestimated"]
            if math.isnan(p):
                flags.append("unestimated")
            elif math.isnan(float(res.log2err[i])):
                flags.append("log2err_unestimated")
            if t in size_filtered:
                flags.append("size_filtered")
            if is_genomic_cluster(deltas[t]):
                flags.append("genomic_cluster")
            if is_detection_confounded(detection.z(gs), es):
                flags.append("detection_confounded")
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
