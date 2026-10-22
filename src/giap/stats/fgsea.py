"""fgsea port I and II: stat preparation, ES and leading edge, the simple permutation stage, NES,
the simple p and ``log2err``, the multilevel decision (issue M8-02) and the multilevel splitting
stage (issue M8-11; GIAP_SPEC.md §7.7.4, §2.4 G15, §8.2, §8.3 row ``fgsea_multilevel``, §6.6).

Ported from fgsea 1.36.2 (Bioconductor 3.22, git b7862b7): ``R/fgsea.R``
(``preparePathwaysAndStats``, ``calcGseaStat``, ``fgseaSimpleImpl``), ``R/util.R``
(``preparePathways``) and ``R/fgseaMultilevel.R`` (simple stage, ``multilevelImpl``,
``multilevelError``, ``denomProb``, ``eps``, ``isCpGeHalf``); scoreType std and gseaParam 1
only. fgsea is MIT-licensed:

    Copyright (c) 2016-2019 Alexey Sergushichev

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software
    and associated documentation files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions: The above copyright
    notice and this permission notice shall be included in all copies or substantial portions of
    the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND (full text in
    NOTICE.md).

fgsea 1.36.2 computes ES on the float |stats| (no integer scaling; M8-02 lead decision 5). All
randomness comes from the caller's generator (stream ``function.gsea``): permutations run in
chunks of 100, each drawn from a child of ``rng.spawn``; counts and sums are accumulated chunk by
chunk in chunk order, so results do not depend on the numba thread count. The multilevel stage
spawns one child per size group (ascending size) and, in each group, one child per sign.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from numpy.typing import ArrayLike
from scipy.special import betaincinv, digamma, polygamma

from giap.stats._fgsea_kernels import MAX_LEVELS, es_walk, multilevel_rulers, random_es_batch

__all__ = [
    "MultilevelResult",
    "Prepared",
    "SimpleResult",
    "fgsea_multilevel",
    "fgsea_simple",
    "gsea_es",
    "prepare",
]

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]
CHUNK = 100
_LN2 = math.log(2.0)


@dataclass(frozen=True, slots=True)
class Prepared:
    """Sorted |stats| (decreasing, stable), set members as sorted 1-based ranks."""

    stats: F64  # |stats| in rank order (float64; 1.36.2 has no integer scaling)
    order: I64  # input index of each sorted position
    names: tuple[str, ...]
    members: tuple[I64, ...]
    sizes: I64
    dropped: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SimpleResult:
    es: F64
    nes: F64
    p: F64
    log2err: F64
    n_more_extreme: I64
    mode_fraction: I64
    needs_multilevel: npt.NDArray[np.bool_]
    leading_edge: tuple[I64, ...]  # input indices, fgsea's order


def prepare(
    stats: ArrayLike,
    labels: Sequence[str],
    sets: Mapping[str, Sequence[str]],
    min_size: int,
    max_size: int,
) -> Prepared:
    """``preparePathwaysAndStats`` + ``preparePathways`` (scoreType std, gseaParam 1)."""
    x = np.asarray(stats, dtype=np.float64)
    if x.ndim != 1 or x.shape[0] != len(labels):
        raise ValueError("stats must be 1-D with one value per label")
    if not np.all(np.isfinite(x)):
        raise ValueError("Not all stats values are finite numbers")
    if len(set(labels)) != len(labels):
        raise ValueError("duplicate labels are not allowed")
    order = np.argsort(-x, kind="stable").astype(np.int64)
    rank = {labels[int(i)]: r + 1 for r, i in enumerate(order)}
    lo, hi = max(min_size, 1), min(max_size, len(labels) - 1)
    names: list[str] = []
    members: list[I64] = []
    dropped: list[str] = []
    for name, genes in sets.items():
        ranks = np.unique(np.array([rank[g] for g in genes if g in rank], dtype=np.int64))
        if lo <= ranks.shape[0] <= hi:
            names.append(name)
            members.append(ranks)
        else:
            dropped.append(name)
    return Prepared(
        stats=np.abs(x[order]),
        order=order,
        names=tuple(names),
        members=tuple(members),
        sizes=np.array([m.shape[0] for m in members], dtype=np.int64),
        dropped=tuple(dropped),
    )


def gsea_es(stats: F64, members: I64) -> tuple[float, I64]:
    """``calcGseaStat`` (std): ES and leading edge (1-based ranks, fgsea's order) of ``members``.

    ES is 0 with an empty edge on the ``maxP == -minP`` tie; a negative ES's edge runs from the
    bottom hit up to the first minimum.
    """
    r = np.asarray(stats, dtype=np.float64)
    S = np.sort(np.asarray(members, dtype=np.int64))
    n, m = r.shape[0], S.shape[0]
    if m == n:
        raise ValueError("GSEA statistic is not defined when all genes are selected")
    if m == 0 or S[0] < 1 or S[-1] > n or np.any(S[1:] == S[:-1]):
        raise ValueError("members must be distinct 1-based ranks")
    es, side, idx = es_walk(S, np.abs(r[S - 1]), n)
    if side > 0:
        return float(es), S[: idx + 1].copy()
    if side < 0:
        return float(es), S[idx:][::-1].copy()
    return 0.0, np.empty(0, dtype=np.int64)


def _trigamma(x: ArrayLike) -> F64:
    return np.asarray(polygamma(1, x), dtype=np.float64)


def _needs_multilevel(nme: I64, n_perm: int, sample_size: int) -> npt.NDArray[np.bool_]:
    """fgseaMultilevel: ``multilevelError(p0, sampleSize) < simpleError`` (Beta borders)."""
    k = nme.astype(np.float64)
    with np.errstate(divide="ignore"):  # qbeta with a zero shape is a point mass (0 or 1)
        q_lo = betaincinv(np.maximum(k, 1), n_perm - k + 1, 0.025)
        q_hi = betaincinv(k + 1, np.maximum(n_perm - k, 1), 0.975)
        left = np.where(nme == 0, -np.inf, np.log2(q_lo))
        right = np.where(nme == n_perm, 0.0, np.log2(q_hi))
    crude = np.log2((k + 1) / (n_perm + 1))
    simple_err = 0.5 * np.maximum(crude - left, right - crude)
    level = np.floor(-crude + 1)
    single = _trigamma((sample_size + 1) / 2) - _trigamma(sample_size + 1)
    mult_err = np.sqrt(level * single) / _LN2
    return np.asarray(mult_err < simple_err, dtype=np.bool_)


def _random_samples(rng: np.random.Generator, n: int, k: int, n_perm: int) -> list[I64]:
    """One (chunk, K) array of 1-based random ranks in draw order per chunk of 100."""
    chunks = [min(CHUNK, n_perm - c) for c in range(0, n_perm, CHUNK)]
    out = []
    for child, size in zip(rng.spawn(len(chunks)), chunks, strict=True):
        rows = [child.choice(n, size=k, replace=False) for _ in range(size)]
        out.append(np.stack(rows).astype(np.int64) + 1)
    return out


def fgsea_simple(
    prep: Prepared, n_perm: int, sample_size: int, rng: np.random.Generator
) -> SimpleResult:
    """``fgseaSimpleImpl`` via the cumulative batch, then fgseaMultilevel's simple-stage rules.

    NES = ES / mean(max(ES_rand, 0) | ES_rand >= 0) if ES > 0, else ES / |mean(min(ES_rand, 0) |
    ES_rand <= 0)|; p = min((1 + nLeEs)/(1 + nLeZero), (1 + nGeEs)/(1 + nGeZero)); sets with
    modeFraction < 10 get p, NES and log2err NaN. ``needs_multilevel`` marks the sets M8-11 refines.
    """
    if n_perm < 1:
        raise ValueError("n_perm must be >= 1")
    if sample_size < 3:
        raise ValueError("sample_size must be >= 3")
    if sample_size % 2 == 0:  # fgseaMultilevel rounds an even sampleSize up
        sample_size += 1
    k = len(prep.names)
    stats = np.asarray(prep.stats, dtype=np.float64)
    obs = [gsea_es(stats, m) for m in prep.members]
    es = np.array([e for e, _ in obs], dtype=np.float64)
    edges = tuple(np.asarray(prep.order)[le - 1].astype(np.int64) for _, le in obs)
    le_es, ge_es, le_zero, ge_zero = (np.zeros(k, dtype=np.int64) for _ in range(4))
    le_sum, ge_sum = np.zeros(k), np.zeros(k)
    if k:
        cols = np.asarray(prep.sizes, dtype=np.int64) - 1
        for samples in _random_samples(rng, stats.shape[0], int(cols.max()) + 1, n_perm):
            rand = random_es_batch(stats, samples)[:, cols]
            le_es += (rand <= es).sum(axis=0)
            ge_es += (rand >= es).sum(axis=0)
            le_zero += (rand <= 0).sum(axis=0)
            ge_zero += (rand >= 0).sum(axis=0)
            le_sum += np.minimum(rand, 0.0).sum(axis=0)
            ge_sum += np.maximum(rand, 0.0).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        denom = np.where(es > 0, ge_sum / ge_zero, np.abs(le_sum / le_zero))
        ok = np.isfinite(denom) & (denom != 0)
        nes = np.where(ok, es / np.where(ok, denom, 1.0), np.nan)
        p = np.minimum((1 + le_es) / (1 + le_zero), (1 + ge_es) / (1 + ge_zero))
    nme = np.where(es > 0, ge_es, le_es).astype(np.int64)
    mode = np.where(es >= 0, ge_zero, le_zero).astype(np.int64)
    estimated = ok & (mode >= 10)
    nes = np.where(estimated, nes, np.nan)
    p = np.where(estimated, p, np.nan)
    log2err = np.sqrt(_trigamma(nme + 1.0) - _trigamma(n_perm + 1.0)) / _LN2
    log2err = np.where(estimated, log2err, np.nan)
    needs = _needs_multilevel(nme, n_perm, sample_size) & estimated
    return SimpleResult(es, nes, p, log2err, nme, mode, needs, edges)


@dataclass(frozen=True, slots=True)
class MultilevelResult:
    """Multilevel p, ``log2err`` and fgsea's ``isCpGeHalf`` of the refined sets."""

    p: F64
    log2err: F64
    is_cp_ge_half: npt.NDArray[np.bool_]


def _multilevel_error(p: F64, sample_size: int) -> F64:
    """fgsea ``multilevelError``: sqrt(floor(-log2 p + 1) * single-level variance) / ln 2."""
    single = _trigamma((sample_size + 1) / 2) - _trigamma(sample_size + 1)
    return np.asarray(np.sqrt(np.floor(-np.log2(p) + 1) * single) / _LN2, dtype=np.float64)


def _beta_mean_log(a: ArrayLike, b: int) -> F64:
    """fgsea ``betaMeanLog``: E log of a Beta(a, b - a + 1) variable."""
    return np.asarray(digamma(a) - digamma(b + 1), dtype=np.float64)


def _n_perm(simple: SimpleResult, need: npt.NDArray[np.bool_]) -> int:
    """``nPermSimple`` recovered from the simple ``log2err`` (trigamma(nPerm + 1) =
    trigamma(nMoreExtreme + 1) - (ln 2 * log2err)^2), inverted by Newton steps."""
    err = np.asarray(simple.log2err, dtype=np.float64)[need]
    if not np.all(np.isfinite(err)):
        raise ValueError("a set needing multilevel has no finite simple log2err")
    tg = _trigamma(np.asarray(simple.n_more_extreme)[need] + 1.0) - (err * _LN2) ** 2
    y = 1.0 / tg + 0.5
    for _ in range(8):
        y = y - (_trigamma(y) - tg) / np.asarray(polygamma(2, y), dtype=np.float64)
    return int(np.rint(np.median(y))) - 1


def _rulers(
    prep: Prepared, simple: SimpleResult, need: npt.NDArray[np.bool_], denom: F64,
    sample_size: int, eps: float, rng: np.random.Generator,
) -> tuple[I64, I64]:  # fmt: skip
    """Run one ruler per (size, sign) of the ``need`` sets; (set index, kernel row) per target.

    Size groups take ``rng.spawn`` children in ascending size order; each group's child spawns
    one stream per sign (positive ES on the ranks, negative on the reversed ranks)."""
    sizes, es = np.asarray(prep.sizes), np.asarray(simple.es)
    groups = sorted({int(k) for k in sizes[need]})
    ids: list[I64] = []
    meta: list[tuple[int, int, int, np.uint64]] = []
    for k, child in zip(groups, rng.spawn(len(groups)), strict=True):
        in_group = need & (sizes == k)
        rounds = MAX_LEVELS
        if eps > 0:  # fgsea stops extending once 2^-levels < eps * min(denomProb) / 2
            eps_group = eps * float(np.min(denom[in_group]))
            rounds = min(MAX_LEVELS, math.ceil(-math.log2(0.5 * eps_group)))
        for side, (sign, stream) in enumerate(zip((es > 0, es < 0), child.spawn(2), strict=True)):
            sel = np.flatnonzero(in_group & sign)
            if sel.size:
                ids.append(sel)
                seed = np.uint64(stream.integers(0, 2**63, dtype=np.int64))
                meta.append((side, k, rounds, seed))
    stats = np.asarray(prep.stats, dtype=np.float64)
    flat = np.concatenate(ids)
    out = multilevel_rulers(
        np.stack([stats, stats[::-1].copy()]),
        np.array([m[0] for m in meta], dtype=np.int64),
        np.array([m[1] for m in meta], dtype=np.int64),
        sample_size,
        np.abs(es[flat]).astype(np.float64),
        np.concatenate([[0], np.cumsum([i.size for i in ids])]).astype(np.int64),
        np.array([m[2] for m in meta], dtype=np.int64),
        np.array([m[3] for m in meta], dtype=np.uint64),
    )
    return flat, out


def fgsea_multilevel(
    prep: Prepared,
    simple: SimpleResult,
    sample_size: int,
    eps: float,
    rng: np.random.Generator,
) -> SimpleResult:
    """fgseaMultilevel's multilevel stage for the sets in ``simple.needs_multilevel``.

    Per target ES, after r levels: log cppMPval = r * betaMeanLog(h, s) + betaMeanLog(c + 1, s)
    + digamma(c_pos + 1) - digamma(c + 1), with h = (s + 1)/2, c the samples of the resolving
    level above ES and c_pos those with a positive ES (the sign correction; isCpGeHalf = 2 c_pos
    >= c). p = min(1, cppMPval / denomProb); p < eps > 0 gives p = eps and log2err NaN; otherwise
    log2err = multilevelError(p, s), NaN when isCpGeHalf is false. Other fields are unchanged.
    """
    if sample_size < 3:
        raise ValueError("sample_size must be >= 3")
    if sample_size % 2 == 0:  # fgseaMultilevel rounds an even sampleSize up
        sample_size += 1
    p = np.array(simple.p, dtype=np.float64)
    log2err = np.array(simple.log2err, dtype=np.float64)
    need = np.asarray(simple.needs_multilevel, dtype=np.bool_) & np.isfinite(p)
    if not need.any():
        return dataclasses.replace(simple, p=p, log2err=log2err)
    mode = np.asarray(simple.mode_fraction, dtype=np.float64)
    denom = (mode + 1.0) / (_n_perm(simple, need) + 1.0)
    flat, out = _rulers(prep, simple, need, denom, sample_size, eps, rng)
    rounds, c, c_pos = out[:, 0], out[:, 1], out[:, 2]
    half = (sample_size + 1) // 2
    log_cpp = (
        rounds * _beta_mean_log(half, sample_size)
        + _beta_mean_log(c + 1, sample_size)
        + digamma(c_pos + 1)
        - digamma(c + 1)
    )
    ml = MultilevelResult(
        p=np.minimum(1.0, np.exp(log_cpp) / denom[flat]),
        log2err=np.empty(flat.size),
        is_cp_ge_half=2 * c_pos >= c,
    )
    floored = ml.p < eps
    ml.log2err[:] = np.where(ml.is_cp_ge_half & ~floored, _multilevel_error(ml.p, sample_size),
                             np.nan)  # fmt: skip
    p[flat] = np.where(floored, eps, ml.p)
    log2err[flat] = ml.log2err
    return dataclasses.replace(simple, p=p, log2err=log2err)
