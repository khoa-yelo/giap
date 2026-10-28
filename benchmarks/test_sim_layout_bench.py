"""Performance budgets for the gippasim sequence-free layout (docs/specs/M4-07-sim-layout-tables.md
acceptance "Speed: 1,000 genes < 0.2 s, 4,500 genes < 2 s (``benchmarks/``)"), seed 1.

- ``LayoutConfig()`` (1,000 genes): best of 3 < 0.2 s. From
  ``tests/unit/test_sim_layout.py::test_speed_1000_genes``.
- ``LayoutConfig(n_genes_chromosome=4500)`` (example 7): best of 3 < 2 s. From
  ``tests/unit/test_sim_layout.py::test_example7_speed_4500_genes``.

The wall-clock asserts were flaky on loaded Slurm nodes. Not part of the per-commit suite; run with
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


def test_sim_layout_example7_4500_genes_budget(benchmark: Any) -> None:
    def run() -> float:
        return min(timed(lambda: _L.simulate(1, n_genes_chromosome=4500))[0] for _ in range(3))

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["best_secs"] = best
    assert best < 2.0, f"best of 3: {best:.3f} s for 4,500 genes"
