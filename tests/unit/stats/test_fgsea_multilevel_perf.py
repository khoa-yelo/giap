"""Performance of the fgsea multilevel stage (issue M8-11 acceptance criterion 3): 2,500 TUs and
1,500 sets of which 100 need multilevel, in < 60 s on 4 CPUs (``slow``: run on Slurm).

The 100 refined sets are planted (members drawn from the top or bottom of the ranking, sizes
10-200, so several size groups and both signs); ``needs_multilevel`` is restricted to them so that
exactly 100 sets are refined whatever the null sets draw. The timed call is
``fgsea_multilevel(prep, simple, 101, 0.0, rng)`` after a small warm-up call that compiles the
numba kernels.
"""

from __future__ import annotations

import dataclasses
import time

import numpy as np
import pytest

from giap.rng import for_step
from giap.stats.fgsea import fgsea_simple, prepare

pytestmark = pytest.mark.slow

N_TUS = 2500
N_SETS = 1500
N_PLANTED = 100
BUDGET_S = 60.0


def test_100_multilevel_sets_of_1500_in_under_a_minute_on_4_cpus() -> None:
    numba = pytest.importorskip("numba")
    from giap.stats.fgsea import fgsea_multilevel

    data = np.random.default_rng(2500)  # test data, not a GIAP stream
    x = data.standard_normal(N_TUS)
    order = np.argsort(-x, kind="stable")
    lab = [f"TU_{i}" for i in range(N_TUS)]
    sets: dict[str, list[str]] = {}
    for j in range(N_PLANTED):
        k = int(data.integers(10, 201))
        pool = order[: 3 * k] if j % 2 == 0 else order[::-1][: 3 * k]
        sets[f"planted{j:03d}"] = [lab[int(i)] for i in data.choice(pool, size=k, replace=False)]
    for j in range(N_SETS - N_PLANTED):
        k = int(data.integers(10, 201))
        sets[f"null{j:04d}"] = [lab[int(i)] for i in data.choice(N_TUS, size=k, replace=False)]
    prep = prepare(x, lab, sets, 10, 500)
    assert len(prep.names) == N_SETS
    simple = fgsea_simple(prep, 1000, 101, for_step("function.gsea", 1))
    planted = np.array([n.startswith("planted") for n in prep.names])
    needs = np.asarray(simple.needs_multilevel) & planted
    assert int(needs.sum()) == N_PLANTED
    simple = dataclasses.replace(simple, needs_multilevel=needs)

    before = numba.get_num_threads()
    try:
        numba.set_num_threads(min(4, int(numba.config.NUMBA_NUM_THREADS)))
        top15 = [lab[int(i)] for i in np.argsort(-x[:200], kind="stable")[:15]]
        small = prepare(x[:200], lab[:200], {"w": top15}, 10, 500)
        warm = fgsea_simple(small, 100, 101, for_step("function.gsea", 2))
        fgsea_multilevel(small, warm, 101, 0.0, for_step("function.gsea", 3))
        t0 = time.perf_counter()
        out = fgsea_multilevel(prep, simple, 101, 0.0, for_step("function.gsea", 4))
        elapsed = time.perf_counter() - t0
    finally:
        numba.set_num_threads(before)
    assert np.all(np.isfinite(out.p[needs]))
    assert elapsed < BUDGET_S, elapsed
