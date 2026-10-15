"""Performance budget for the protein covariate table (issue M4-10 acceptance).

``covariate_table`` on a 1,000-gene C2 table (M4-07 layout size, median 280 aa): < 2 s. Inputs are
built by tests/unit/genetable/test_covariates.py (``build_genes``), which keeps a loose 10x guard.
Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

from ._unit_helpers import load_test_module

_T = load_test_module("unit/genetable/test_covariates.py")
BUDGET_S = 2.0


def test_covariate_table_1000_genes_budget(benchmark: Any) -> None:
    genes = _T.build_genes(1000)
    secs = benchmark.pedantic(lambda: _T.time_thousand(genes), rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= BUDGET_S, f"1,000 genes took {secs:.2f} s (budget {BUDGET_S} s)"
