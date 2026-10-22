"""Performance budget for M8-01 Wallenius ORA (docs/specs/M8-01 acceptance, lead decision 3).

3,000 categories x ``wallenius_tails`` on a universe of M = 4,500 genes: < 5 s on 1 CPU.
Category sizes are log-uniform on [5, 500] (GO-like), N = 1,500 detected, odds in [0.5, 2] and k
drawn near its central expectation. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from giap.stats.wallenius import wallenius_tails

from ._unit_helpers import timed

CATEGORIES, M, N = 3_000, 4_500, 1_500
BUDGET_S = 5.0


def test_wallenius_tails_3000_categories_m4500_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(81)
    sizes = np.exp(rng.uniform(np.log(5), np.log(500), size=CATEGORIES)).astype(int)
    odds = rng.uniform(0.5, 2.0, size=CATEGORIES)
    ks = np.minimum(rng.binomial(sizes, N / M), np.minimum(sizes, N))
    cases = [(int(k), int(n), float(w)) for k, n, w in zip(ks, sizes, odds, strict=True)]

    def run() -> float:
        secs, out = timed(lambda: [wallenius_tails(k, M, n, N, w) for k, n, w in cases])
        assert len(out) == CATEGORIES and all(0.0 <= o <= 1.0 for o, _ in out)
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < BUDGET_S, f"wallenius_tails x {CATEGORIES} took {secs:.2f} s"
