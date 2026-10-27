"""Performance budget for the M5-02 autocorrelation engine (docs/specs/M5-02 acceptance, lead
decision 12).

``morans_i_matrix`` on B = 1,000 gene-label permutations of a 4,500-gene circular series at
lags 1 ... 20: < 5 s on 1 CPU (the permutation draw is excluded). Not part of the per-commit
suite; run with ``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from gippa.stats.autocorr import label_permutations, morans_i_matrix

from ._unit_helpers import timed

B, N, LAGS = 1_000, 4_500, list(range(1, 21))
BUDGET_S = 5.0


def test_morans_i_matrix_b1000_n4500_20_lags_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(52)
    x = rng.normal(size=N)
    X = x[label_permutations(N, B, rng)]

    def run() -> float:
        secs, out = timed(lambda: morans_i_matrix(X, LAGS, True))
        assert out.shape == (B, len(LAGS)) and np.isfinite(out).all()
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < BUDGET_S, f"morans_i_matrix took {secs:.2f} s"
