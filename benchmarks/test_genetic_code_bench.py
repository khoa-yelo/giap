"""Performance budget for the M11-14w genetic-code check: a 10-Mb synthetic genome in <= 5 s on
one core (docs/specs/M11-14w-wrong-genetic-code-warning.md AC 8; reported only, not a CI gate).

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from giap.preflight.genetic_code import genetic_code_check

from ._unit_helpers import timed


def _genome(n: int, gc: float, seed: int) -> str:
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    codes = np.random.default_rng(seed).choice(4, size=n, p=p)
    return np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")


def test_genetic_code_10mb_budget(benchmark: Any) -> None:
    # three records (a chromosome and two plasmids) at low GC, where S4 runs are longest
    seq = _genome(10_000_000, 0.30, 7)
    records = [seq[:9_000_000], seq[9_000_000:9_800_000], seq[9_800_000:]]

    def run() -> float:
        secs, check = timed(lambda: genetic_code_check(records, 11))
        assert check.level in {"ok", "info", "warn"}
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 5.0, f"took {secs:.2f} s"
