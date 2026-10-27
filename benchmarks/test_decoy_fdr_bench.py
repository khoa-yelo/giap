"""Performance budget for the M6-01 decoy FDR (docs/specs/M6-01-stats-tdist-decoy-fdr.md).

The acceptance budget is 50,000 targets x 2 decoy sets x 5 distance bins < 2 s; lead decision 11
keeps it here rather than in ``tests/``. One ``max_decoy_fdr`` call does everything the budget
covers: the two reweighted step functions, the element-wise maximum, the monotone q and the two
local-fdr KDEs. Decoy sets are the same size as the target list, which is the worst case the
pair table can produce. Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from gippa.stats.decoy_fdr import bin_index, max_decoy_fdr

from ._unit_helpers import timed

N_TARGET = 50_000
N_DECOY = 50_000
#: ``operons.distance_bins_bp``: five half-open bins
EDGES = (-np.inf, 0.0, 50.0, 150.0, 300.0, np.inf)


def _pairs(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Scores ``-|Delta|`` and distance bins of ``n`` pairs, bins drawn from real-ish distances."""
    scores = -np.abs(rng.normal(size=n))
    distance = rng.normal(loc=60.0, scale=200.0, size=n)
    return scores, bin_index(distance, EDGES)


def test_decoy_fdr_50k_targets_2_decoy_sets_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(601)
    target_s, target_bin = _pairs(rng, N_TARGET)
    decoys = {name: _pairs(rng, N_DECOY) for name in ("convergent", "sandwich")}

    def run() -> float:
        secs, res = timed(lambda: max_decoy_fdr(target_s, target_bin, decoys, n_bins=5))
        assert res.q.shape == (N_TARGET,)
        assert np.all((res.q >= 0.0) & (res.q <= 1.0))
        assert set(res.which.tolist()) <= set(decoys)
        assert np.all(np.isnan(res.lfdr) | ((res.lfdr >= 0.0) & (res.lfdr <= 1.0)))
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    benchmark.extra_info["n_target"] = N_TARGET
    benchmark.extra_info["n_decoy"] = N_DECOY
    assert secs < 2.0, f"max_decoy_fdr took {secs:.2f} s"
