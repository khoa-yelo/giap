"""Numba kernels of the fgsea port (issues M8-02, M8-11; GIAP_SPEC.md §7.7.4, §6.6).

Ported from fgsea 1.36.2 (Bioconductor 3.22, git b7862b7): ``R/fgsea.R`` ``calcGseaStat``,
``src/fastGSEA.cpp`` ``calcGseaStatCumulative`` / ``calcGseaStatCumulativeBatch`` (scoreType std,
gseaParam 1) and the ``EsRuler`` of ``src/fgseaMultilevelSupplement.cpp`` (adaptive multilevel
splitting, Korotkevich et al. 2021). fgsea is MIT-licensed:

    Copyright (c) 2016-2019 Alexey Sergushichev

    Permission is hereby granted, free of charge, to any person obtaining a copy of this software
    and associated documentation files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions: The above copyright
    notice and this permission notice shall be included in all copies or substantial portions of
    the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND (full text in
    NOTICE.md).

fgsea's C++ finds each prefix's running-sum extremes with a block convex hull; this port walks the
prefix directly (O(K^2) per sample) with the exact arithmetic of ``calcGseaStat``, so the random
ES of size k equals ``gsea_es`` on the same ranks. The batch kernel parallelises over samples
only; every sample is computed by one thread, so results do not depend on the thread count. The
multilevel kernel parallelises over rulers the same way (one ruler, one thread, one seeded
stream).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from numba import njit, prange

F64 = npt.NDArray[np.float64]
I64 = npt.NDArray[np.int64]


@njit(cache=True, nogil=True)
def es_walk(
    S: I64, w: F64, n: int
) -> tuple[float, int, int]:  # pragma: no cover - compiled by numba; via gsea_es
    """calcGseaStat (std) on sorted 1-based hits ``S`` with weights ``w``: (ES, side, index).

    side is +1 (leading edge ``S[:index + 1]``), -1 (``S[index:]`` reversed) or 0 (the
    ``maxP == -minP`` tie). Sums are accumulated sequentially in rank order.
    """
    m = S.shape[0]
    nm = n - m
    nr = 0.0
    for j in range(m):
        nr += w[j]
    cum = 0.0
    max_p, min_p = -np.inf, np.inf
    i_max, i_min = 0, 0
    for j in range(m):
        cum += w[j]
        rc = (j + 1) / m if nr == 0 else cum / nr
        top = rc - (S[j] - (j + 1)) / nm
        bottom = top - 1 / m if nr == 0 else top - w[j] / nr
        if top > max_p:
            max_p, i_max = top, j
        if bottom < min_p:
            min_p, i_min = bottom, j
    if max_p == -min_p:
        return 0.0, 0, 0
    if max_p > -min_p:
        return max_p, 1, i_max
    return min_p, -1, i_min


@njit(cache=True, nogil=True)
def random_es_prefixes(
    stats: F64, sample: I64
) -> F64:  # pragma: no cover - compiled by numba; via random_es_batch
    """ES of ``sort(sample[:k])`` for k = 1..K (calcGseaStatCumulative, std).

    As in fgsea 1.36.2, zero stats weigh ``statEps = min(1e-5, min non-zero |s| of the
    sample) / 1024`` in this kernel only.
    """
    n = stats.shape[0]
    k_max = sample.shape[0]
    eps = 1e-5
    for i in range(k_max):
        x = abs(stats[sample[i] - 1])
        if x > 0:
            eps = min(x, eps)
    eps /= 1024
    S = np.empty(k_max, dtype=np.int64)
    w = np.empty(k_max, dtype=np.float64)
    out = np.empty(k_max, dtype=np.float64)
    for k in range(k_max):
        s = sample[k]
        x = max(abs(stats[s - 1]), eps)
        j = k
        while j > 0 and S[j - 1] > s:  # insertion keeps the hits sorted
            S[j] = S[j - 1]
            w[j] = w[j - 1]
            j -= 1
        S[j] = s
        w[j] = x
        out[k] = es_walk(S[: k + 1], w[: k + 1], n)[0]
    return out


@njit(cache=True, parallel=True)
def random_es_batch(
    stats: F64, samples: I64
) -> F64:  # pragma: no cover - compiled by numba; via fgsea_simple
    """:func:`random_es_prefixes` for every row of ``samples`` (one sample per thread)."""
    n_samples, k_max = samples.shape
    out = np.empty((n_samples, k_max), dtype=np.float64)
    for p in prange(n_samples):
        out[p] = random_es_prefixes(stats, samples[p])
    return out


# --------------------------------------------------------------------------- multilevel (M8-11)
# The MCMC moves draw from splitmix64 streams seeded from the caller's numpy generator (one seed
# per ruler); numba's own np.random state is per thread and is never used.
_GAMMA = np.uint64(0x9E3779B97F4A7C15)
_MIX1 = np.uint64(0xBF58476D1CE4E5B9)
_MIX2 = np.uint64(0x94D049BB133111EB)
MAX_LEVELS = 1000  # 0.4976 ** 1000 ~ 1e-303: the splitting never reaches p = 0
_ATTEMPT_CAP = 10  # a level stalls after 10 * sampleSize * K attempts (acceptance < 10 %)


@njit(cache=True, nogil=True)
def _below(state: npt.NDArray[np.uint64], n: int) -> int:  # pragma: no cover - numba
    """Uniform integer in [0, n) from the splitmix64 state ``state[0]``."""
    state[0] += _GAMMA
    z = state[0]
    z = (z ^ (z >> np.uint64(30))) * _MIX1
    z = (z ^ (z >> np.uint64(27))) * _MIX2
    z ^= z >> np.uint64(31)
    return min(int(float(z >> np.uint64(11)) * (1.0 / 9007199254740992.0) * n), n - 1)


@njit(cache=True, nogil=True)
def _extremes(S: I64, w: F64) -> tuple[float, float]:  # pragma: no cover - numba
    """(maxP, minP) of the running sum of sorted 1-based hits ``S`` with weights ``w[S - 1]``,
    in the exact arithmetic of :func:`es_walk`; the sample's ES is positive iff maxP > -minP."""
    m = S.shape[0]
    nm = w.shape[0] - m
    nr = 0.0
    for j in range(m):
        nr += w[S[j] - 1]
    cum = 0.0
    max_p, min_p = -np.inf, np.inf
    for j in range(m):
        x = w[S[j] - 1]
        cum += x
        rc = (j + 1) / m if nr == 0 else cum / nr
        top = rc - (S[j] - (j + 1)) / nm
        bottom = top - 1 / m if nr == 0 else top - x / nr
        max_p = max(max_p, top)
        min_p = min(min_p, bottom)
    return max_p, min_p


@njit(cache=True, nogil=True)
def _reaches(S: I64, w: F64, level: float) -> bool:  # pragma: no cover - numba
    """maxP >= ``level``, walked from the last hit back: top_j <= cum_j / nr, and cum_j only
    falls going back, so the walk stops at the first cum_j / nr < level (cum_j as nr minus the
    later weights; it can differ from the forward sum by rounding)."""
    m = S.shape[0]
    nm = w.shape[0] - m
    nr = 0.0
    for j in range(m):
        nr += w[S[j] - 1]
    if nr == 0:
        return _extremes(S, w)[0] >= level
    later = 0.0
    for j in range(m - 1, -1, -1):
        rc = (nr - later) / nr
        if rc < level:
            return False
        if rc - (S[j] - (j + 1)) / nm >= level:
            return True
        later += w[S[j] - 1]
    return False


@njit(cache=True, nogil=True)
def _perturb(
    S: I64, w: F64, level: float, state: npt.NDArray[np.uint64], backup: I64
) -> int:  # pragma: no cover - numba
    """One MCMC move: swap a random member for a random non-member and keep it when the sample's
    maxP stays >= ``level``; returns 1 when accepted."""
    n, k = w.shape[0], S.shape[0]
    x = _below(state, n) + 1
    i = np.searchsorted(S, x)
    while i < k and S[i] == x:
        x = _below(state, n) + 1
        i = np.searchsorted(S, x)
    backup[:] = S
    for j in range(_below(state, k), k - 1):  # drop one member, keep the rest sorted
        S[j] = S[j + 1]
    j = k - 1
    while j > 0 and S[j - 1] > x:
        S[j] = S[j - 1]
        j -= 1
    S[j] = x
    if _reaches(S, w, level):
        return 1
    S[:] = backup
    return 0


@njit(cache=True, nogil=True)
def _ruler(
    w: F64, k: int, s: int, targets: F64, max_rounds: int, seed: np.uint64, out: I64
) -> None:  # pragma: no cover - compiled by numba; via multilevel_rulers
    """EsRuler of one (size, sign): ``s`` random K-sets, then levels at the sample median maxP
    (the top half kept twice, the median once, then ``s * K`` accepted MCMC moves above the
    level) until every target is below the median, ``max_rounds`` levels or a stalled level.
    ``out[t]`` = (levels below target t, #samples with maxP > t, #those with a positive ES)."""
    n, h = w.shape[0], (s + 1) // 2
    state = np.array([seed], dtype=np.uint64)
    pop = np.empty((s, k), dtype=np.int64)
    for i in range(s):
        for j in range(k):  # k distinct ranks by rejection
            x = _below(state, n) + 1
            while np.any(pop[i, :j] == x):
                x = _below(state, n) + 1
            pop[i, j] = x
        pop[i] = np.sort(pop[i])
    vals, pos = np.empty(s), np.empty(s, dtype=np.bool_)
    backup = np.empty(k, dtype=np.int64)
    done = np.zeros(targets.shape[0], dtype=np.bool_)
    rounds = 0
    while True:
        for i in range(s):
            mx, mn = _extremes(pop[i], w)
            vals[i], pos[i] = mx, mx > -mn
        order = np.argsort(vals, kind="mergesort")
        med = vals[order[h - 1]]
        for t in range(targets.shape[0]):
            if not done[t] and (med > targets[t] or rounds >= max_rounds):
                above = vals > targets[t]
                out[t, 0], out[t, 1] = rounds, above.sum()
                out[t, 2] = (above & pos).sum()
                done[t] = True
        if done.all():
            return
        new = np.empty_like(pop)
        new[0] = pop[order[h - 1]]
        for i in range(h, s):
            new[2 * (i - h) + 1] = pop[order[i]]
            new[2 * (i - h) + 2] = pop[order[i]]
        pop = new
        rounds += 1
        moves, attempts, iters = 0, 0, max(1, k // 10)
        while moves < s * k and attempts < _ATTEMPT_CAP * s * k:
            for i in range(s):
                for _ in range(iters):
                    moves += _perturb(pop[i], w, med, state, backup)
                attempts += iters
        if moves < s * k:  # the level cannot be mixed: resolve the rest at this level
            max_rounds = rounds


@njit(cache=True, parallel=True)
def multilevel_rulers(
    w2: F64,
    sides: I64,
    sizes: I64,
    s: int,
    targets: F64,
    offsets: I64,
    max_rounds: I64,
    seeds: npt.NDArray[np.uint64],
) -> I64:  # pragma: no cover - compiled by numba; via fgsea_multilevel
    """:func:`_ruler` for every ruler g (weights ``w2[sides[g]]``, targets
    ``targets[offsets[g]:offsets[g + 1]]``), one ruler per thread: thread-count independent."""
    out = np.zeros((targets.shape[0], 3), dtype=np.int64)
    for g in prange(sizes.shape[0]):
        a, b = offsets[g], offsets[g + 1]
        _ruler(w2[sides[g]], sizes[g], s, targets[a:b], max_rounds[g], seeds[g], out[a:b])
    return out
