"""Numba kernels of the fgsea port (issue M8-02; GIAP_SPEC.md §7.7.4, §6.6).

Ported from fgsea 1.36.2 (Bioconductor 3.22, git b7862b7): ``R/fgsea.R`` ``calcGseaStat`` and
``src/fastGSEA.cpp`` ``calcGseaStatCumulative`` / ``calcGseaStatCumulativeBatch`` (scoreType std,
gseaParam 1). fgsea is MIT-licensed:

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
only; every sample is computed by one thread, so results do not depend on the thread count.
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
