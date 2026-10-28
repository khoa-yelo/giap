"""Abundance GSEA on a TU-collapsed ranking (spec §7.7.4, §8.1, §8.2; issue M8-05).

Pure. The ranked units are TUs: a TU's score is the mean ``abundance_log2`` of its detected,
eligible genes (imputed values included), median-centred (ADR-B), and a TU is in a set when any
of its scoring genes is. Sets are sized in TUs. The statistics come from an injected engine
(M8-17 §7, ADR-0038): the step passes R ``fgseaMultilevel`` (``gippa.io.fgsea_r``), the default
is the M8-02 ``fgsea_simple`` port. One call holds a ``main:<collection>`` job per collection
with in-range sets, the TUs in a seeded random order (M8-02 lead decision 7: ties are never
broken by TU ID) and one seed drawn after that order; a set covering every ranked TU gets an NA
row (M8-17 lead decision 3). q is BH within the collection. A significant GO parent is
``redundant`` when, re-tested without the genes of its significant descendants, its p
substituted into the collection's BH vector gives q > alpha (lead decisions 3, 6 and 7: a
re-test with 0 TUs is redundant, an unestimated re-test leaves ``redundant`` NA); every re-test
goes into one second call. The writer joins the list cells ``leading_edge`` and ``flags``.

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
from typing import Any, Protocol, TypeAlias

import numpy as np
import numpy.typing as npt
import pandas as pd

from gippa.config.model import GseaConfig
from gippa.function.collections import GeneSets
from gippa.stats.fgsea import fgsea_simple, prepare
from gippa.stats.multitest import bh

__all__ = [
    "CLUSTER_TAU",
    "COLUMNS",
    "DETECTION_Z0",
    "SIMPLE_ENGINE",
    "EngineCall",
    "EngineProvenance",
    "EngineResult",
    "EngineRow",
    "GseaEngine",
    "GseaJob",
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
METHOD = "fgsea_simple"
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


# ---------------------------------------------------------------- the engine interface (M8-17 §7)


@dataclass(frozen=True, slots=True)
class GseaJob:
    """One engine job: ``main:<collection>`` or ``retest:<collection>:<term>`` (M8-17 §7)."""

    job_id: str
    sets: Mapping[str, tuple[str, ...]]  # term -> member TU IDs, natural order
    min_size: int
    max_size: int


@dataclass(frozen=True, slots=True)
class EngineCall:
    """One engine call: a ranking in the seeded tie order (§5), its jobs, parameters and seed.
    Every score is finite and every TU ID unique, otherwise ``ValueError`` (§8 edge cases)."""

    tu_ids: tuple[str, ...]  # perm order (§5)
    scores: F64  # same order, all finite
    jobs: tuple[GseaJob, ...]
    params: GseaConfig
    seed: int

    def __post_init__(self) -> None:
        scores = np.asarray(self.scores, dtype=np.float64)
        if scores.ndim != 1 or scores.size != len(self.tu_ids):
            raise ValueError(f"{scores.size} scores for {len(self.tu_ids)} TUs")
        if not np.isfinite(scores).all():
            bad = [t for t, ok in zip(self.tu_ids, np.isfinite(scores), strict=True) if not ok]
            raise ValueError(f"non-finite GSEA scores for {len(bad)} TUs, e.g. {bad[0]}")
        if len(set(self.tu_ids)) != len(self.tu_ids):
            raise ValueError("duplicate TU IDs in the GSEA ranking")


@dataclass(frozen=True, slots=True)
class EngineRow:
    """One tested set of one job; NaN where the engine gives NA."""

    term: str
    size: int
    es: float
    nes: float
    p: float
    log2err: float
    leading_edge: tuple[str, ...]  # TU IDs, fgsea's order


EngineResult: TypeAlias = Mapping[str, tuple[EngineRow, ...]]  # job_id -> rows


@dataclass(frozen=True, slots=True)
class EngineProvenance:
    """What the engine reported about itself for one call (M8-17a lead decision 2)."""

    r: str
    fgsea: str
    warnings: tuple[tuple[str, str], ...]  # (job_id, message), distinct, in R's order
    seed: int = 0  # the seed R ran with (#204 review nit 2)
    jobs: int = 0  # the number of jobs R ran


class GseaEngine(Protocol):
    """A callable GSEA engine; ``method`` is written to the table's ``method`` column."""

    method: str

    def __call__(self, call: EngineCall) -> EngineResult: ...


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


class _SimpleEngine:
    """The M8-02 ``fgsea_simple`` port as a ``GseaEngine`` (ADR-0038 item 6: kept until the
    hardening pass). Each job draws from a fresh ``np.random.default_rng(call.seed)``, as R's
    ``set.seed`` runs before every job."""

    method = METHOD

    def __call__(self, call: EngineCall) -> EngineResult:
        stats = np.asarray(call.scores, dtype=np.float64)
        labels = list(call.tu_ids)
        p = call.params
        out: dict[str, tuple[EngineRow, ...]] = {}
        for job in call.jobs:
            prep = prepare(stats, labels, job.sets, job.min_size, job.max_size)
            if not prep.names:
                out[job.job_id] = ()
                continue
            rng = np.random.default_rng(call.seed)
            res = fgsea_simple(prep, p.n_perm_simple, p.sample_size, rng)
            out[job.job_id] = tuple(
                EngineRow(
                    term=t,
                    size=int(prep.sizes[i]),
                    es=float(res.es[i]),
                    nes=float(res.nes[i]),
                    p=float(res.p[i]),
                    log2err=float(res.log2err[i]),
                    leading_edge=tuple(labels[int(j)] for j in res.leading_edge[i]),
                )
                for i, t in enumerate(prep.names)
            )
        return out


SIMPLE_ENGINE: GseaEngine = _SimpleEngine()


@dataclass(frozen=True, slots=True)
class _Main:
    """One collection's main run: the engine rows of its tested terms, their q, and the terms
    covering every ranked TU (no job, NA row; lead decision 3)."""

    coll: GeneSets
    tus: dict[str, frozenset[str]]
    rows: tuple[EngineRow, ...]
    q: F64
    whole: frozenset[str]


@dataclass(frozen=True, slots=True)
class _Hierarchy:
    redundant: dict[str, Any]
    size_filtered: set[str]
    group: dict[str, str]
    retests: dict[str, tuple[int, int]]  # job_id -> (row index, k) awaiting the re-test call


_Natural = Callable[[Iterable[str]], tuple[str, ...]]
_NA_ROW = EngineRow("", 0, math.nan, math.nan, math.nan, math.nan, ())


def _checked(job: GseaJob, rows: Sequence[EngineRow]) -> tuple[EngineRow, ...]:
    """The engine's rows of ``job``; a term not sent, a size not Python's or a sent term without
    its row is a ``ValueError`` (a dropped row would shrink BH's count)."""
    for r in rows:
        if r.term not in job.sets or r.size != len(set(job.sets[r.term])):
            raise ValueError(f"GSEA engine row {job.job_id} {r.term} (size {r.size}) does not "
                             "match the sent set")  # fmt: skip
    missing = [t for t in job.sets if t not in {r.term for r in rows}]
    if missing:
        raise ValueError(f"GSEA engine returned no row for {job.job_id} {', '.join(missing)}")
    return tuple(rows)


def _go_hierarchy(
    m: _Main, alpha: float, min_size: int, tu_of: Mapping[str, str], natural: _Natural, n: int
) -> tuple[_Hierarchy, list[GseaJob]]:
    """``redundant`` (True / False, decided later for re-tested terms), ``size_filtered`` and
    ``parent_group`` of a GO collection's significant terms, plus their re-test jobs."""
    coll, names = m.coll, [r.term for r in m.rows]
    sig = {t for t, q in zip(names, m.q, strict=True) if q <= alpha}
    h = _Hierarchy(dict.fromkeys(names, False), set(), {}, {})
    jobs: list[GseaJob] = []
    if coll.collection not in GO_NAMESPACES or not sig:
        return h, jobs
    children: dict[str, set[str]] = {}
    for t, ps in coll.parents.items():
        for parent in ps:
            children.setdefault(parent, set()).add(t)
    below, above, depth = _closure(children), _closure(coll.parents), _depth(coll.parents)
    for i, t in enumerate(names):
        if t not in sig:
            continue
        anc = above(t) & sig  # the term itself only when no ancestor is significant
        h.group[t] = min(anc, key=lambda a: (depth(a), a)) if anc else t
        desc = below(t) & sig
        if not desc:
            continue
        removed = frozenset().union(*(coll.sets[d] for d in desc))
        tus = _members(tu_of, coll.sets[t] - removed)
        if len(tus) < min_size:
            h.size_filtered.add(t)
        if not tus:  # the parent is the union of its significant descendants
            h.redundant[t] = True
            continue
        job = GseaJob(f"retest:{coll.collection}:{t}", {t: natural(tus)}, 1, n)
        jobs.append(job)
        h.retests[job.job_id] = (i, len(tus))
    return h, jobs


def _resolve(m: _Main, h: _Hierarchy, alpha: float, result: EngineResult) -> None:
    """``redundant`` of each re-tested term from its re-test p (NaN or no row: NA)."""
    p = np.array([r.p for r in m.rows], dtype=np.float64)
    for job_id, (i, _) in h.retests.items():
        got = result.get(job_id, ())
        p_re = float(got[0].p) if got else math.nan
        t = m.rows[i].term
        if math.isnan(p_re):  # lead decision 7.4: undecided
            h.redundant[t] = pd.NA
        else:
            sub = p.copy()
            sub[i] = p_re
            h.redundant[t] = bool(bh(sub)[i] > alpha)


def gsea_table(
    genes: pd.DataFrame,
    colls: Sequence[GeneSets],
    size: tuple[int, int],
    cfg: GseaConfig,
    alpha: float,
    rng: np.random.Generator,
    replicons: Replicons = None,
    *,
    engine: GseaEngine = SIMPLE_ENGINE,
) -> pd.DataFrame:
    """The enrichment_gsea rows of every collection (§9.3); run flags in ``attrs["run_flags"]``.
    ``replicons`` (C1 frame or contig -> topology) marks circular contigs; others are linear.

    ``engine`` runs at most two calls (M8-17 §8): one ``main:<collection>`` job per collection
    with in-range sets, then one call with every GO re-test job. ``rng`` gives the seeded tie
    order and then the one engine seed (§5). Raises ``ValueError`` when an engine row does not
    match its sent set."""
    if _no_abundance(genes):
        return _frame([], ["layer_missing:abundance"])
    rank = tu_ranking(genes)
    n = len(rank.tu_ids)
    if n < 2:
        return _frame([], ["too_few_ranked_tus"])
    perm = rng.permutation(n)  # seeded tie order (M8-02 lead decision 7)
    seed = int(rng.integers(1, 2**31 - 1))  # one R seed for both calls (§5)
    labels = tuple(rank.tu_ids[int(i)] for i in perm)
    stats = np.asarray(rank.scores, dtype=np.float64)[perm]
    tu_of = _tus_of(rank)
    pos = {t: i for i, t in enumerate(rank.tu_ids)}

    def natural(tus: Iterable[str]) -> tuple[str, ...]:
        return tuple(sorted(tus, key=pos.__getitem__))

    lo, hi = size
    tested: list[tuple[GeneSets, dict[str, frozenset[str]], GseaJob | None, frozenset[str]]] = []
    for coll in colls:
        tus = tu_sets(rank, coll)
        sets = {t: natural(s) for t, s in tus.items() if lo <= len(s) <= hi and len(s) < n}
        job = GseaJob(f"main:{coll.collection}", sets, lo, hi) if sets else None
        # a collection with no job adds no row, whole-ranking sets included (lead decision 6a)
        whole = frozenset(t for t, s in tus.items() if len(s) == n and lo <= n <= hi and sets)
        tested.append((coll, tus, job, whole))
    jobs = tuple(j for _, _, j, _ in tested if j is not None)
    if not jobs:
        return _frame([], [])
    result = engine(EngineCall(labels, stats, jobs, cfg, seed))

    mains: list[_Main] = []
    for coll, tus, job, whole in tested:
        rows = _checked(job, result.get(job.job_id, ())) if job is not None else ()
        q = bh(np.array([r.p for r in rows], dtype=np.float64)) if rows else np.empty(0)
        mains.append(_Main(coll, tus, rows, q, whole))
    hier: list[_Hierarchy] = []
    retests: list[GseaJob] = []
    for m in mains:
        h, more = _go_hierarchy(m, alpha, lo, tu_of, natural, n)
        hier.append(h)
        retests += more
    if retests:
        jobs_re = tuple(retests)
        got = engine(EngineCall(labels, stats, jobs_re, cfg, seed))
        checked = {j.job_id: _checked(j, got.get(j.job_id, ())) for j in jobs_re}
        for m, h in zip(mains, hier, strict=True):
            _resolve(m, h, alpha, checked)

    detection = _Detection.of(genes)
    out: list[dict[str, Any]] = []
    for m, h in zip(mains, hier, strict=True):
        by_term = {r.term: (r, float(q)) for r, q in zip(m.rows, m.q, strict=True)}
        terms = [t for t in m.coll.sets if t in by_term or t in m.whole]
        deltas = cluster_deltas(rank, genes, replicons, {t: m.tus[t] for t in terms})
        for t in terms:
            r, qv = by_term.get(t, (_NA_ROW, math.nan))
            out.append(_out_row(m, t, r, qv, h, deltas[t], detection, rank, engine.method, cfg))
    table = _frame(out, [])
    table.attrs["rank_labels"] = labels  # the engine's TU order (M9-17 running sums)
    return table


def _out_row(
    m: _Main,
    t: str,
    r: EngineRow,
    q: float,
    h: _Hierarchy,
    delta: float,
    detection: _Detection,
    rank: TuRanking,
    method: str,
    cfg: GseaConfig,
) -> dict[str, Any]:
    """One table row from an engine row (§8 row mapping, flags, leading edge in genes)."""
    gs = m.coll.sets[t]
    edge_tus = r.leading_edge if r.es != 0 else ()
    edge = [g for tu in edge_tus for g in rank.tu_genes[tu] if g in gs]
    flags = ["ci_unestimated"]
    if math.isnan(r.p):
        flags.append("unestimated")
    elif math.isnan(r.log2err):
        flags.append("log2err_unestimated")
    if t in h.size_filtered:
        flags.append("size_filtered")
    if is_genomic_cluster(delta):
        flags.append("genomic_cluster")
    if is_detection_confounded(detection.z(gs), r.es):
        flags.append("detection_confounded")
    return {
        "collection": m.coll.collection,
        "term": t,
        "name": m.coll.names.get(t, "NA"),
        "size_tus": len(m.tus[t]),
        "es": float(r.es),
        "nes": float(r.nes),
        "log2err": float(r.log2err),
        "leading_edge": edge,
        "redundant": h.redundant.get(t, False),
        "parent_group": h.group.get(t),
        "coverage": float(m.coll.coverage),
        "stat": float(r.es),
        "effect": float(r.nes),
        "ci_low": math.nan,
        "ci_high": math.nan,
        "p_value": float(r.p),
        "q_value": q,
        "method": method,
        "null_model": NULL_MODEL,
        "n_perm": int(cfg.n_perm_simple),
        "flags": flags,
    }
