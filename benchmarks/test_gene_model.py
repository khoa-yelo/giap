"""Performance budgets for ``GeneModel.overlapping`` on a 5-Mb, 5,000-gene model.

100k random queries: <= 1.5 s; with a 1-Mb gap feature: <= 3 s. Moved from
tests/unit/test_gene_model_review.py (which keeps loose 10x guards). Not part of the per-commit
suite; run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import pytest

from ._unit_helpers import load_test_module

_T = load_test_module("unit/test_gene_model_review.py")


@pytest.mark.parametrize(("with_gap", "budget"), [(False, 1.5), (True, 3.0)])
def test_overlapping_100k_queries_budget(benchmark: Any, with_gap: bool, budget: float) -> None:
    model, queries = _T._big_model(with_gap=with_gap), _T._queries()
    secs = benchmark.pedantic(lambda: _T._time_queries(model, queries), rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < budget, f"100k overlapping queries took {secs:.2f} s (budget {budget} s)"
