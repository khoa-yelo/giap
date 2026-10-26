"""Runtime of the determinants permutation null on worker processes (GIAP_SPEC.md §7.8.2 step 4;
ADR-0063).

The spec example 1 world of ``tests/helpers/determinants_models_data.py`` (400 genes, TUs of
1-5 genes, the six planted columns, signal 2), built once and not timed; ``null_p`` with
B = 400 permutations (>= ``_PROCESS_MIN_B``, so the worker-process path) of both models
(shipped ``DeterminantsConfig``) on 8 threads. The measured time belongs under
``determinants_null_b400_n400_t8`` in ``benchmarks/baseline.json``; the budget is about 3 x it.

Run with ``pytest benchmarks/test_determinants_null_processes_bench.py --benchmark-only``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
import determinants_models_data as d

SEED = 20261023
B = 400
THREADS = 8
BUDGET_S = 140.0


def test_determinants_null_b400_n400_t8_budget(benchmark: Any) -> None:
    from giap.config.model import DeterminantsConfig
    from giap.determinants.null import null_p

    X, y, tu, _ = d.planted(np.random.default_rng(SEED))
    cfg = DeterminantsConfig()
    folds = d.folds_for(tu, cfg, SEED)
    observed = {"lightgbm": 0.5, "elastic_net": 0.5}

    def run() -> float:
        secs, res = timed(lambda: null_p(X, y, tu, folds[0], observed, cfg, B, SEED, THREADS))
        assert all(r.null_r2.shape == (B,) for r in res.values()), res
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= BUDGET_S, f"null_p took {secs:.1f} s (budget {BUDGET_S})"
