"""Runtime of the M4 chain on ``sim4500`` (M4-20 acceptance "Benchmarks recorded
(``benchmarks/baseline.json``): gene table, detectability, context stats on ``sim4500``").

One ``build_m4_fixture(size="sim4500")`` work directory (built once, not timed), then each step
timed on its own: ``run_gene_table``, ``run_detectability`` and ``run_context_stats`` at the
shipped ``stats.permutations`` / ``stats.bootstrap`` = 1,000. The measured medians belong under
``m4_chain`` in ``benchmarks/baseline.json`` (recorded by whoever runs the benchmark, like the
other budgets); the budgets below are deliberately generous first numbers for the M4-20 PR and
should be tightened to about 3 x the measured time once recorded.

Run with ``pytest benchmarks/test_m4_chain_bench.py --benchmark-only`` (needs ``tests/sim`` on
PYTHONPATH for gippasim).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pytest

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
import m4_chain as mc

SEED = 20261021
BUDGETS_S = {"gene_table": 120.0, "detectability": 300.0, "context_stats": 600.0}


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Any]:
    root = tmp_path_factory.mktemp("m4_bench")
    w = mc.build(root / "work", seed=SEED, size="sim4500")
    return w, mc.test_db(root / "db"), mc.config(SEED)


@pytest.mark.parametrize("step", sorted(BUDGETS_S))
def test_m4_step_budget_sim4500(benchmark: Any, work: tuple[Path, Path, Any], step: str) -> None:
    from gippa.steps.context_stats import run_context_stats
    from gippa.steps.detectability import run_detectability
    from gippa.steps.gene_table import run_gene_table

    w, db, cfg = work

    def run() -> float:
        with mc.db_env(db):
            if step == "gene_table":
                secs, _ = timed(lambda: run_gene_table(w, cfg))
            elif step == "detectability":
                run_gene_table(w, cfg)
                secs, _ = timed(lambda: run_detectability(w, cfg))
            else:
                run_gene_table(w, cfg)
                run_detectability(w, cfg)
                secs, _ = timed(lambda: run_context_stats(w, cfg, None))
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    budget = BUDGETS_S[step]
    assert secs <= budget, f"{step} on sim4500 took {secs:.1f} s (budget {budget} s)"
