"""Performance budget for the giapsim operon scenarios (docs/specs/M6-03-sim-operon-scenarios.md
acceptance "Runtime: a 1,800-gene scenario in < 5 s (the giapsim budget)"), seed 1.

One ``simulate_operons`` replicate of each of the five scenarios at 1,800 genes (layout,
expression, pairs and truth frames): best of 3 < 5 s. The scenario arguments are those of
tests/unit/sim/test_operon_scenarios.py examples 1-6. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import pytest

from ._unit_helpers import timed

N = 1800
BUDGET_S = 5.0
CASES = [
    ("operon_null", (N,)),
    ("missed_operon_power", (N, 50, 0.1)),
    ("discordance_power", (N, 40, 2.0)),
    ("operon_attack_ar1", (N, 0.9)),
    ("operon_attack_detection", (N, "prophage")),
]


@pytest.mark.parametrize("name,args", CASES, ids=[c[0] for c in CASES])
def test_sim_operon_1800_genes_budget(benchmark: Any, name: str, args: tuple[Any, ...]) -> None:
    import giapsim.scenarios as sc

    cfgs = getattr(sc, name)(*args)

    def run() -> float:
        return min(timed(lambda: sc.simulate_operons(cfgs, 1))[0] for _ in range(3))

    best = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["best_secs"] = best
    assert best < BUDGET_S, f"best of 3: {best:.3f} s for {name} at 1,800 genes"
