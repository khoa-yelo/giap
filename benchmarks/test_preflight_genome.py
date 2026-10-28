"""Performance budget for genome preflight on a pathological 5-Mb alternating R/Y contig.

``summarize_genome`` + ``genome_checks``: <= 2 s (moved from
tests/unit/preflight/test_genome_rules_review.py, which keeps a loose 10x guard).

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

from gippa.preflight.genome import LIST_CAP, genome_checks, summarize_genome

from ._unit_helpers import load_test_module, timed

_T = load_test_module("unit/preflight/test_genome_rules_review.py")


def test_alternating_ry_5mb_preflight_budget(benchmark: Any) -> None:
    n = 5_000_000
    config = _T.cfg(*_T.SMALL)
    records = [_T.rec("a", "RY" * (n // 2))]

    def run() -> float:
        def body() -> Any:
            s = summarize_genome(records, config.preflight.n_run_min_bp)
            return s, {c.rule_id: c for c in genome_checks(s, config)}

        secs, (s, out) = timed(body)
        assert len(s.ambiguous) <= LIST_CAP
        assert out["genome.gaps"].level == "warn"
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 2.0, f"took {secs:.2f} s"
