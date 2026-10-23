"""Runtime budgets of the M7-12 giapsim scenarios (docs/specs/M7-12-sim-variant-novel-scenarios.md
rule 4 and lead decision 5): each scenario best of 3 under 2 s at 100 kb and under 20 s at 1 Mb
(``variant_f1`` and ``novel_orf_power``, the two with a ``size``), seed 20261022.

Not part of the per-commit suite; run with ``pytest benchmarks --benchmark-only``. Recorded in
``benchmarks/baseline.json`` (``sim_m7_scenarios``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import giapsim
import giapsim.scenarios as sc

from ._unit_helpers import timed

SEED = 20261022
_CASES = [
    ("variant_f1", "100k", 2.0),
    ("novel_orf_power", "100k", 2.0),
    ("same_strain_insilico", None, 2.0),
    ("isobaric_attack", None, 2.0),
    ("variant_f1", "1M", 20.0),
    ("novel_orf_power", "1M", 20.0),
]


@pytest.mark.parametrize(("name", "size", "budget"), _CASES)
def test_sim_m7_scenario_budget(
    benchmark: Any, tmp_path: Path, name: str, size: str | None, budget: float
) -> None:
    fn = getattr(sc, name)
    cfg = fn(SEED) if size is None else fn(SEED, size)

    def once(i: int) -> float:
        return timed(lambda: giapsim.simulate(cfg, SEED, tmp_path / str(i)))[0]

    def run() -> float:
        return min(once(i) for i in range(3))

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["best_secs"] = best
    assert best < budget, f"best of 3: {best:.3f} s for {name} {size or '100k'}"
