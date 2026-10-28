"""Performance budgets for the gippasim simulator (tests/sim; spec example 10).

- 100-kb genome, defaults: <= 2 s (from tests/unit/test_sim_basic.py).
- 100-kb genome, all options: best of 3 <= 1 s (from tests/unit/test_sim_features.py).
- 5-Mb genome, all options (example 10): <= 30 s (from tests/unit/test_sim_features.py).

The unit tests keep loose 10x guards. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import gippasim

from ._unit_helpers import load_test_module, timed

_F = load_test_module("unit/test_sim_features.py")


def test_sim_100kb_default_budget(benchmark: Any, tmp_path: Path) -> None:
    cfg = gippasim.SimConfig(genome_length=100_000)

    def run() -> float:
        return timed(lambda: gippasim.simulate(cfg, seed=5, outdir=tmp_path))[0]

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 2.0, secs


def test_sim_100kb_all_options_budget(benchmark: Any, tmp_path: Path) -> None:
    cfg = gippasim.SimConfig(**_F.all_options(100_000))

    def once(i: int) -> float:
        return timed(lambda: gippasim.simulate(cfg, seed=9, outdir=tmp_path / str(i)))[0]

    def run() -> float:
        return min(once(i) for i in range(3))

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["best_secs"] = best
    assert best < 1.0, best


def test_sim_example10_5mb_all_options_budget(benchmark: Any, tmp_path: Path) -> None:
    opts = _F.all_options(5_000_000) | {"n_rrna_operons": 7, "n_trna": 86, "n_sorfs": 40}
    cfg = gippasim.SimConfig(**opts)

    def run() -> float:
        return timed(lambda: gippasim.simulate(cfg, seed=10, outdir=tmp_path))[0]

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 30.0, secs
