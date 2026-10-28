"""Hypothesis properties of the fgsea port I (issue M8-02; GIPPA_SPEC.md §7.7.4, §8.3).

- ES in [-1, 1]; the leading edge is a non-empty subset of the set iff ES != 0;
- negating every stat negates the ES and keeps the leading-edge genes in the same order (top of
  the walk first). Stats are continuous (no ties, so the stable tie order plays no part), and
  inputs whose ES is an exact rational tie (maxP == -minP) are skipped: IEEE rounding decides
  those, and it need not round the mirrored walk the same way;
- the batch random ES of each prefix size k equals ``gsea_es`` on the first k drawn ranks.

The batch kernel is pinned as ``gippa.stats._fgsea_kernels.random_es_prefixes(stats, sample)``:
``stats`` the prepared float64 |stats| (length N), ``sample`` int64 1-based ranks in draw order
(length K < N, distinct); returns float64[K] with element k - 1 the ES of ``sort(sample[:k])``
(fgsea 1.36.2 ``calcGseaStatCumulative``, scoreType std, ties to 0). 1.36.2 gives zero stats a
weight of statEps in this kernel only, so the property uses non-zero stats; exact ties are skipped
as above.
"""

from __future__ import annotations

import numpy as np
from helpers.fgsea_oracle import is_exact_tie
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from gippa.stats._fgsea_kernels import random_es_prefixes
from gippa.stats.fgsea import gsea_es, prepare

FAST = settings(max_examples=200, deadline=None)
TOL = 1e-12

STAT = st.one_of(
    st.integers(-6, 6).map(float),  # equal magnitudes make tops ties likely
    st.floats(-1e3, 1e3, allow_nan=False, allow_infinity=False),
)


@st.composite
def stats_and_set(draw: st.DrawFn) -> tuple[list[float], list[int]]:
    n = draw(st.integers(2, 40))
    stats = draw(st.lists(STAT, min_size=n, max_size=n))
    idx = draw(st.lists(st.integers(0, n - 1), min_size=1, max_size=n - 1, unique=True))
    return stats, idx


@st.composite
def continuous_stats_and_set(draw: st.DrawFn) -> tuple[np.ndarray, list[int]]:
    n = draw(st.integers(2, 60))
    rng = np.random.default_rng(draw(st.integers(0, 2**32 - 1)))
    scale = 10.0 ** rng.integers(-4, 4, size=n) if draw(st.booleans()) else 1.0
    stats = rng.normal(size=n) * scale
    idx = draw(st.lists(st.integers(0, n - 1), min_size=1, max_size=n - 1, unique=True))
    return stats, idx


def _labels(n: int) -> list[str]:
    return [f"g{i}" for i in range(n)]


@FAST
@given(stats_and_set())
def test_es_bounded_and_edge_iff_nonzero(case: tuple[list[float], list[int]]) -> None:
    stats, idx = case
    lab = _labels(len(stats))
    prep = prepare(stats, lab, {"s": [lab[i] for i in idx]}, 1, len(stats) - 1)
    es, le = gsea_es(prep.stats, prep.members[0])
    assert -1 - TOL <= es <= 1 + TOL
    le_list = [int(x) for x in le]
    members = {int(x) for x in prep.members[0]}
    assert set(le_list) <= members
    assert len(set(le_list)) == len(le_list)
    assert (len(le_list) > 0) == (es != 0)
    if es > 0:
        assert le_list == sorted(le_list)
        assert le_list[0] == min(members)
    elif es < 0:
        assert le_list == sorted(le_list, reverse=True)
        assert le_list[0] == max(members)


@FAST
@given(continuous_stats_and_set())
def test_negating_stats_negates_es_and_keeps_edge(case: tuple[np.ndarray, list[int]]) -> None:
    stats, idx = case
    n = len(stats)
    lab = _labels(n)
    sets = {"s": [lab[i] for i in idx]}
    out = []
    for s in (stats, -stats):
        prep = prepare(s, lab, sets, 1, n - 1)
        assume(not is_exact_tie(np.asarray(prep.stats, dtype=float), prep.members[0]))
        es, le = gsea_es(prep.stats, prep.members[0])
        out.append((es, [int(prep.order[r - 1]) for r in le]))
    assert abs(out[0][0] + out[1][0]) <= TOL
    assert out[0][1] == out[1][1]


@st.composite
def prefix_case(draw: st.DrawFn) -> tuple[np.ndarray, np.ndarray]:
    n = draw(st.integers(3, 60))
    seed = draw(st.integers(0, 2**32 - 1))
    rng = np.random.default_rng(seed)
    x = rng.normal(size=n)
    r = np.abs(np.sort(x)[::-1])
    k = draw(st.integers(1, n - 1))
    sample = (rng.permutation(n)[:k] + 1).astype(np.int64)
    return r, sample


@FAST
@given(prefix_case())
def test_batch_prefix_es_equals_gsea_es(case: tuple[np.ndarray, np.ndarray]) -> None:
    r, sample = case
    r0, s0 = r.copy(), sample.copy()
    got = np.asarray(random_es_prefixes(r, sample), dtype=float)
    assert got.shape == (len(sample),)
    for k in range(1, len(sample) + 1):
        if is_exact_tie(r, sample[:k]):
            continue
        es, _ = gsea_es(r, np.sort(sample[:k]))
        assert abs(got[k - 1] - es) <= TOL, (k, got[k - 1], es)
    np.testing.assert_array_equal(r, r0)
    np.testing.assert_array_equal(sample, s0)
