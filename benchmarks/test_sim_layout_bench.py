"""Performance budget for the giapsim sequence-free layout (docs/specs/M4-07-sim-layout-tables.md).

``simulate_layout(LayoutConfig(), default_rng(1))`` (1,000 genes): best of 3 < 0.2 s (M4-07
acceptance "Speed: 1,000 genes < 0.2 s ... (``benchmarks/``)"). Moved from
``tests/unit/test_sim_layout.py::test_speed_1000_genes``; the wall-clock assert was flaky on
loaded Slurm nodes. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

from typing import Any

from ._unit_helpers import load_test_module, timed

_L = load_test_module("unit/test_sim_layout.py")


def test_sim_layout_1000_genes_budget(benchmark: Any) -> None:
    def run() -> float:
        return min(timed(lambda: _L.simulate(1))[0] for _ in range(3))

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["best_secs"] = best
    assert best < 0.2, f"best of 3: {best:.3f} s for 1,000 genes"
