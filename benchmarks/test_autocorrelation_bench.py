"""Performance budget for the M5-05 autocorrelation table (docs/specs/M5-05 acceptance
"Performance", lead decision 30).

``giap.context.autocorrelation.autocorrelation_rows`` on a 4,500-gene circular replicon with
B = 1,000 permutations per null, both nulls, both statistics and lags 1 ... 20: **< 30 s on
4 CPUs** (the giapsim simulation of the gene table is excluded from the timing). The dominant
cost is the ADR-0020 re-detrend of every replicate, one (B x n) (n x n) smoother product per
null, so the number scales with the detected-gene count and with the BLAS thread count: run it
with ``OMP_NUM_THREADS=4`` / ``srun -c 4``.

Not part of the per-commit suite; run with ``pytest benchmarks --benchmark-only``. Recorded in
``benchmarks/baseline.json`` (entry ``autocorrelation_rows_b1000_n4500_20_lags``).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "sim"))

B, N, LAGS = 1_000, 4_500, 20
BUDGET_S = 30.0


def _inputs() -> tuple[Any, Any, Any]:
    """The gene table, C1 replicons and config of one 4,500-gene TU-null replicate."""
    import autocorr_calib as acal
    import autocorr_data as ad
    import region_tests_data as rd

    sc = acal.tu_null(N)
    rng = np.random.default_rng(7)
    tables = rd.sim_tables(rng, acal.layout_config(N), sc.expr_kw, p_source="true")
    return tables.genes, tables.replicons, ad.config(B, lags=(1, LAGS))


def test_autocorrelation_rows_b1000_n4500_20_lags_budget(benchmark: Any) -> None:
    from giap.context.autocorrelation import autocorrelation_rows
    from giap.rng import for_step

    genes, replicons, cfg = _inputs()

    def run() -> float:
        secs, out = timed(
            lambda: autocorrelation_rows(
                genes, replicons, cfg, for_step("context_stats.autocorrelation", 11)
            )
        )
        assert len(out) == 2 * 2 * LAGS
        assert int(out["n_perm"].max()) == B
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    benchmark.extra_info["threads"] = os.environ.get("OMP_NUM_THREADS", "unset")
    assert secs < BUDGET_S, f"autocorrelation_rows took {secs:.2f} s (budget {BUDGET_S} s)"
