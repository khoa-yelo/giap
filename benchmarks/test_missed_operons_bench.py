"""Performance budget for M6-06 missed operons (docs/specs/M6-06-missed-operons.md "Acceptance").

20,000 NO targets + 20,000 decoys (10,000 CV and 10,000 SW) over the five default distance bins:
one ``missed_operons`` call < 2 s - target and decoy selection, binning, the M6-01 engine (both
reweighted FDRs, the maximum, monotone q, both lfdr KDEs) and the row frame. The frame is built
with the unit tests' ``missed_data`` helpers, so the inputs cannot drift from the tests. Run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
from helpers import missed_data as md

from gippa.operons.missed import missed_operons

N_TARGET = 20_000
N_DECOY = 10_000  # per decoy set, two sets
BUDGET_S = 2.0


def _frame(rng: np.random.Generator) -> Any:
    parts = []
    for cls, n in (("NO", N_TARGET), ("CV", N_DECOY), ("SW", N_DECOY)):
        b = rng.integers(0, 5, size=n)
        d = rng.normal(size=n)
        parts += [md.block(cls, d[b == k], int(k)) for k in range(5)]
    return md.stack(*parts)


def test_missed_operons_20k_targets_20k_decoys_budget(benchmark: Any) -> None:
    pairs = _frame(np.random.default_rng(606))
    cfg = md.cfg()

    def run() -> float:
        secs, res = timed(lambda: missed_operons(pairs, cfg))
        assert res.n_targets == N_TARGET
        assert sum(res.n_decoys.values()) == 2 * N_DECOY
        q = res.rows["q_value"].to_numpy(dtype=float)
        assert np.all((q >= 0.0) & (q <= 1.0))
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    benchmark.extra_info["n_target"] = N_TARGET
    benchmark.extra_info["n_decoy"] = 2 * N_DECOY
    assert secs < BUDGET_S, f"missed_operons took {secs:.2f} s"
