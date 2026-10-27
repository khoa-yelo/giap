"""Performance budget for the conditional Poisson-binomial primitive (issue M8-04a; docs/specs/
M8-04 §1 step 1 and the acceptance benchmark).

The ORA acceptance benchmark (M8-04b) is 3,000 categories x both tails plus one elim re-test each
at T = 2,500 TUs and N = 1,250 in < 10 s on 1 CPU. The primitive's share of it is measured here:
one ``tilt`` per run, then per category one ``conditional_poisson_binomial`` call on the full set
(both tails come from one call) and one on an elim-reduced set (half its TUs), category sizes
uniform in ``function.ora_size`` [5, 500]. Budget: the same 10 s on 1 CPU (numba warmed up
first). The spec-writer measured about 0.4 ms per call, i.e. about 2.4 s here.

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numba
import numpy as np

from gippa.stats.poisson_binomial import conditional_poisson_binomial, tilt

from ._unit_helpers import timed

N_TU, N_DETECTED, N_CATEGORIES = 2_500, 1_250, 3_000
SIZE_RANGE = (5, 500)
BUDGET_S = 10.0


def test_cpb_3000_categories_full_and_elim_2500_tus_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(84)
    q = rng.beta(2.0, 2.0, N_TU)
    d = np.zeros(N_TU, dtype=bool)
    d[rng.choice(N_TU, size=N_DETECTED, replace=False)] = True
    sizes = rng.integers(SIZE_RANGE[0], SIZE_RANGE[1] + 1, size=N_CATEGORIES)
    sizes[0] = SIZE_RANGE[1]
    cats = [rng.choice(N_TU, size=int(s), replace=False) for s in sizes]
    reduced = [c[: max(1, c.size // 2)] for c in cats]
    qt0 = tilt(q, N_DETECTED)
    conditional_poisson_binomial(qt0[:5], qt0[5:], int(d[:5].sum()), N_DETECTED)  # compile

    def work() -> int:
        qt = tilt(q, N_DETECTED)
        n_calls = 0
        for c in (*cats, *reduced):
            m = np.zeros(N_TU, dtype=bool)
            m[c] = True
            t = conditional_poisson_binomial(qt[m], qt[~m], int(d[m].sum()), N_DETECTED)
            assert np.isfinite(t.p_over) and np.isfinite(t.p_under)
            n_calls += 1
        return n_calls

    before = numba.get_num_threads()
    numba.set_num_threads(1)

    def run() -> float:
        secs, n_calls = timed(work)
        assert n_calls == 2 * N_CATEGORIES
        return secs

    try:
        secs = benchmark.pedantic(run, rounds=1, iterations=1)
    finally:
        numba.set_num_threads(before)
    benchmark.extra_info["secs"] = secs
    benchmark.extra_info["ms_per_call"] = 1e3 * secs / (2 * N_CATEGORIES)
    assert secs < BUDGET_S, f"conditional_poisson_binomial x {2 * N_CATEGORIES} took {secs:.2f} s"
