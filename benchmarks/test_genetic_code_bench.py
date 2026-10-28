"""Performance budget for the M11-14w genetic-code check: a 10-Mb synthetic genome in <= 5 s on
one core (docs/specs/M11-14w-wrong-genetic-code-warning.md AC 8; reported only, not a CI gate).

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

from helpers.genetic_code_fixture import synthetic_genome

from gippa.preflight.genetic_code import genetic_code_check

from ._unit_helpers import timed


def test_genetic_code_10mb_budget(benchmark: Any) -> None:
    # ~10 Mb code-4 fixture (8,000 genes x 400 codons + spacers): thousands of long runs, so the
    # coverage path is exercised (random DNA has no 300-codon runs and returns early; #423 N2)
    seq = synthetic_genome(7, 0.30, 4, n_genes=8_000)
    records = [seq[:9_000_000], seq[9_000_000:9_800_000], seq[9_800_000:]]

    def run() -> float:
        secs, check = timed(lambda: genetic_code_check(records, 11))
        assert check.level == "warn"
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 5.0, f"took {secs:.2f} s"
