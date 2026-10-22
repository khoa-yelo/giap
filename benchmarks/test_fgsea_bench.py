"""Performance budget for the fgsea simple stage (docs/specs/M8-02 acceptance).

``fgsea_simple`` on 2,500 ranked TUs x 1,500 sets (sizes uniform in ``function.gsea_size``
[10, 500], one of size 500 so K = 500) x 1,000 permutations: < 20 s on 4 CPUs (4 numba threads),
observed ES and the random draws included. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

from typing import Any

import numba
import numpy as np

from giap.rng import for_step
from giap.stats.fgsea import fgsea_simple, prepare

from ._unit_helpers import timed

N_TU, N_SETS, N_PERM, SAMPLE_SIZE = 2_500, 1_500, 1_000, 101
THREADS = 4
BUDGET_S = 20.0


def test_fgsea_simple_2500_tus_1500_sets_1000_perms_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(82)
    stats = rng.standard_normal(N_TU)
    labels = [f"tu{i:05d}" for i in range(N_TU)]
    sizes = rng.integers(10, 501, size=N_SETS)
    sizes[0] = 500
    sets = {
        f"set{j:04d}": [labels[i] for i in rng.choice(N_TU, size=int(s), replace=False)]
        for j, s in enumerate(sizes)
    }
    prep = prepare(stats, labels, sets, 10, 500)
    assert len(prep.names) == N_SETS
    fgsea_simple(prep, 10, SAMPLE_SIZE, for_step("function.gsea", 0))  # compile outside timing
    before = numba.get_num_threads()
    numba.set_num_threads(THREADS)  # the job runs with 4 CPUs

    def run() -> float:
        secs, res = timed(
            lambda: fgsea_simple(prep, N_PERM, SAMPLE_SIZE, for_step("function.gsea", 1))
        )
        assert res.es.shape == (N_SETS,) and np.isfinite(res.p).all()
        return secs

    try:
        secs = benchmark.pedantic(run, rounds=1, iterations=1)
    finally:
        numba.set_num_threads(before)
    benchmark.extra_info["secs"] = secs
    assert secs < BUDGET_S, f"fgsea_simple took {secs:.2f} s"
