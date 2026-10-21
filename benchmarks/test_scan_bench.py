"""Performance budget for the M5-03 window scan (docs/specs/M5-03-stats-scan.md "Acceptance").

``scan_null_max`` on B = 1,000 bootstrap rows of n = 2,500 TUs, windows of 2-20 TUs, circular:
< 2 s on 1 CPU (plain numpy, nothing to compile). Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from giap.stats.scan import scan_null_max

from ._unit_helpers import timed


def test_scan_null_max_b1000_n2500_budget(benchmark: Any) -> None:
    Z = np.random.default_rng(503).normal(size=(1000, 2500))

    def run() -> float:
        secs, m = timed(lambda: scan_null_max(Z, 2, 20, circular=True))
        assert m.shape == (1000,) and np.isfinite(m).all()
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 2.0, f"scan_null_max took {secs:.2f} s"
