"""Slow oracles for the fgsea port (issue M8-02), transcribed from fgsea 1.36.2 (MIT;
Copyright 2016-2019 Alexey Sergushichev; Bioconductor 3.22 tarball ``fgsea_1.36.2.tar.gz``):
``R/fgsea.R`` ``calcGseaStat`` (scoreType std, gseaParam 1) and the simple-stage error terms of
``R/fgseaMultilevel.R``.

``reference_gsea_stat`` uses float64 with every sum accumulated sequentially in rank order (spec
item 2, "accumulated in rank order"), so it reproduces the implementation's tie decisions
``maxP == -minP`` bit for bit. ``exact_es`` repeats the computation in rational arithmetic; it is
used to recognise inputs whose ES is an exact tie, which float rounding may decide either way.
"""

from __future__ import annotations

import math
from fractions import Fraction

import numpy as np
from scipy.special import polygamma
from scipy.stats import beta


def reference_gsea_stat(r: np.ndarray, selected: np.ndarray) -> tuple[float, list[int]]:
    """(ES, leading edge as 1-based ranks in fgsea's order) of the hits ``selected``."""
    S = np.sort(np.asarray(selected, dtype=np.int64))
    m, N = len(S), len(r)
    if m == N:
        raise ValueError("all selected")
    r_adj = np.abs(np.asarray(r, dtype=np.float64)[S - 1])
    cum = np.cumsum(r_adj)  # sequential, rank order
    NR = cum[-1]
    seq = np.arange(1, m + 1)
    r_cum = seq / m if NR == 0 else cum / NR
    tops = r_cum - (S - seq) / (N - m)
    bottoms = tops - 1 / m if NR == 0 else tops - r_adj / NR
    max_p, min_p = tops.max(), bottoms.min()
    if max_p == -min_p:
        return 0.0, []
    if max_p > -min_p:
        return float(max_p), [int(x) for x in S[: int(np.argmax(tops)) + 1]]
    return float(min_p), [int(x) for x in S[int(np.argmin(bottoms)) :][::-1]]


def exact_es(r: np.ndarray, selected: np.ndarray) -> tuple[Fraction, Fraction, Fraction]:
    """(ES, maxP, minP) in exact rational arithmetic (floats converted exactly)."""
    S = sorted(int(x) for x in selected)
    m, N = len(S), len(r)
    w = [Fraction(abs(float(r[s - 1]))) for s in S]
    NR = sum(w, Fraction(0))
    tops, bottoms = [], []
    cum = Fraction(0)
    for j, s in enumerate(S, start=1):
        cum += w[j - 1]
        rc = Fraction(j, m) if NR == 0 else cum / NR
        top = rc - Fraction(s - j, N - m)
        tops.append(top)
        bottoms.append(top - (Fraction(1, m) if NR == 0 else w[j - 1] / NR))
    max_p, min_p = max(tops), min(bottoms)
    if max_p == -min_p:
        return Fraction(0), max_p, min_p
    return (max_p if max_p > -min_p else min_p), max_p, min_p


def is_exact_tie(r: np.ndarray, selected: np.ndarray) -> bool:
    """True when maxP == -minP in exact arithmetic (ES 0 by fgsea's rule, fragile in floats)."""
    _, max_p, min_p = exact_es(r, selected)
    return max_p == -min_p


def trigamma(x: float) -> float:
    return float(polygamma(1, x))


def simple_log2err(n_more_extreme: int, n_perm: int) -> float:
    """fgseaMultilevel: 1/log(2) * sqrt(trigamma(nMoreExtreme + 1) - trigamma(nPermSimple + 1))."""
    return math.sqrt(trigamma(n_more_extreme + 1) - trigamma(n_perm + 1)) / math.log(2)


def simple_and_multilevel_errors(
    n_more_extreme: int, n_perm: int, sample_size: int
) -> tuple[float, float]:
    """(simpleError, multError) of fgseaMultilevel; R's qbeta with a zero shape is a point mass."""
    k = n_more_extreme
    left = -math.inf if k == 0 else math.log2(beta.ppf(0.025, k, n_perm - k + 1))
    right = 0.0 if k == n_perm else math.log2(beta.ppf(0.975, k + 1, n_perm - k))
    crude = math.log2((k + 1) / (n_perm + 1))
    simple_err = 0.5 * max(crude - left, right - crude)
    p0 = (k + 1) / (n_perm + 1)
    level = math.floor(-math.log2(p0) + 1)
    single = trigamma((sample_size + 1) / 2) - trigamma(sample_size + 1)
    return simple_err, math.sqrt(level * single) / math.log(2)
