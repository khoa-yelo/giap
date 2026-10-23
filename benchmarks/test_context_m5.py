"""Runtime and peak memory of ``context_stats`` with its M5 part on the 4,500-gene fixture
(M5-11 acceptance: "``context_stats`` M5 part on the 4,500-gene fixture (B = 1,000): wall time and
peak memory recorded in ``benchmarks/baseline.json``; budget < 2 min on 4 CPUs").

One ``build_m4_fixture(size="sim4500", oric_frac=0.8)`` work directory with the island of
``m5_chain.ISLAND`` planted (so the scan, the effect CI and the ori-ter bootstrap all have work),
``run_gene_table`` and ``run_detectability`` once (not timed), then ``run_context_stats`` at the
shipped ``stats.permutations`` / ``stats.bootstrap`` = 1,000: timed once (the M4 part of the step
is included, so the M5 budget holds a fortiori; the M4 part alone is ``m4_chain.context_stats``
in the baseline, about 10 s), and once more under ``tracemalloc`` for the peak of Python-tracked
allocations (NumPy reports to tracemalloc) plus the process high-water mark ``ru_maxrss``.

The numbers belong under ``context_m5`` in ``benchmarks/baseline.json`` (``secs``, ``peak_mb``,
``maxrss_mb``, ``budget_s`` = 120, ``cpus`` = 4, node, Slurm job, commit, date) --
tests/integration/test_m5_gate.py checks that entry. Run on Slurm with 4 CPUs:
``pytest benchmarks/test_context_m5.py --benchmark-only``.
"""

from __future__ import annotations

import resource
import sys
import tracemalloc
from pathlib import Path
from typing import Any

import pytest

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
import m4_chain as mc
import m5_chain as m5

SEED = 20261021
BUDGET_S = 120.0
CASE = m5.IslandCase("sim4500_bench", "sim4500", 4500, 20, 3.0, m5.ISLAND.detection)


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Any]:
    from giap.steps.detectability import run_detectability
    from giap.steps.gene_table import run_gene_table

    root = tmp_path_factory.mktemp("m5_bench")
    w = m5.build_island(root / "work", SEED, CASE)
    db, cfg = mc.test_db(root / "db"), mc.config(SEED)
    with mc.db_env(db):
        run_gene_table(w, cfg)
        run_detectability(w, cfg)
    return w, db, cfg


def test_context_stats_m5_budget_sim4500(benchmark: Any, work: tuple[Path, Path, Any]) -> None:
    from giap.steps.context_stats import run_context_stats

    w, db, cfg = work

    def run() -> float:
        with mc.db_env(db):
            secs, _ = timed(lambda: run_context_stats(w, cfg, None))
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    tracemalloc.start()
    try:
        with mc.db_env(db):
            run_context_stats(w, cfg, None)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    maxrss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    benchmark.extra_info.update(
        secs=secs, peak_mb=peak / 2**20, maxrss_mb=maxrss_kb / 1024, n_genes=CASE.n_genes
    )
    assert (w / "tables" / "islands.tsv").is_file() and (
        w / "tables" / "autocorrelation.tsv"
    ).is_file()
    assert secs < BUDGET_S, (
        f"context_stats with M5 on sim4500 took {secs:.1f} s (budget {BUDGET_S} s)"
    )
