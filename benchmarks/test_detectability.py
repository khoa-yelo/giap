"""Runtime of one D1 replicate of the M4-13 detectability calibration (lead decision 8,
docs/specs/M4-13-detectability.md).

One replicate = ``detectability_data.sim_geneset`` (M4-07 layout, 2,000 chromosome genes, the D1
fast size) + ``fit_detectability`` (5 TU folds, inner C search). Only the fit is timed. The accepted
numbers (per-replicate and the whole fast D1 run, R = 20) live in ``benchmarks/baseline.json``
under ``detectability``; the budget below is a loose guard. Run with
``pytest benchmarks --benchmark-only`` (needs ``tests/sim`` on PYTHONPATH for giapsim).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from giap.config.model import DetectabilityConfig
from giap.rng import for_step

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
import detectability_data as dd

N_GENES = 2000
BUDGET_S = 30.0


def test_detectability_d1_replicate_budget(benchmark: Any) -> None:
    from giap.genetable.detectability import fit_detectability

    gs = dd.sim_geneset(np.random.default_rng(20261015), N_GENES)

    def run() -> float:
        secs, res = timed(
            lambda: fit_detectability(
                gs.detected,
                gs.cov,
                gs.tu_of_gene,
                DetectabilityConfig(),
                for_step("detectability", 1),
            )
        )
        assert res.constant is False
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= BUDGET_S, f"one D1 replicate took {secs:.2f} s (budget {BUDGET_S} s)"
