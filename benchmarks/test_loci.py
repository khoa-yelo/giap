"""Performance budget for locus classification (issue M3-14 acceptance; GIAP_SPEC.md §7.2.6).

500,000 accepted hits on a 5-Mb, ~4,500-gene model through ``annotate_peptide_hits``: <= 20 s
(interval lookup by sorted arrays). Inputs are built by tests/unit/test_loci.py (``big_case``),
which keeps a loose 10x guard. Not part of the per-commit suite; run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

from ._unit_helpers import load_test_module

_T = load_test_module("unit/test_loci.py")
BUDGET_S = 20.0


def test_annotate_500k_hits_budget(benchmark: Any) -> None:
    case = _T.big_case()
    secs = benchmark.pedantic(lambda: _T.time_big_case(case), rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= BUDGET_S, f"500k hits took {secs:.2f} s (budget {BUDGET_S} s)"
