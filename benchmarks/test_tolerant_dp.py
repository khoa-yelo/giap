"""Performance budget for the banded verifier (issue M3-06 acceptance; GIAP_SPEC.md §7.2.4).

10**6 ``verify`` calls on 12-aa queries (k = 2, a mix of exact, planted 1-2 difference and
random candidates on a 100,000-aa random frame) <= 5 s on one core, numba warmed up first. Run
in a fresh single-threaded subprocess.

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from typing import Any

SCRIPT = textwrap.dedent(
    """
    import json, time
    import numpy as np
    from giap.mapping.diffops import block_tables
    from giap.mapping.tolerant_dp import verify

    tables = block_tables()
    rng = np.random.default_rng(0)
    frame = rng.integers(0, 19, size=100_000).astype(np.uint8)
    frame[rng.random(frame.shape[0]) < 0.02] = 19  # some stops
    cases = []
    for i in range(1000):
        s = int(rng.integers(10, frame.shape[0] - 30))
        q = frame[s : s + 12].copy()
        q[q >= 19] = 0
        mode = i % 4
        if mode == 1:
            q[int(rng.integers(12))] = rng.integers(0, 19)
        elif mode == 2:
            j = int(rng.integers(11))
            q[j], q[j + 1] = q[j + 1], q[j]
            q[int(rng.integers(12))] = rng.integers(0, 19)
        elif mode == 3:
            q = rng.integers(0, 19, size=12).astype(np.uint8)
        cases.append((np.ascontiguousarray(q), s + int(rng.integers(-1, 2))))
    for q, d0 in cases[:50]:
        verify(q, frame, d0, 2, tables)  # warm-up (compile or load the numba cache)
    found = 0
    t0 = time.perf_counter()
    for _ in range(1000):
        for q, d0 in cases:
            found += len(verify(q, frame, d0, 2, tables))
    secs = time.perf_counter() - t0
    assert found >= 400_000, found
    print(json.dumps({"secs": secs, "calls": 1_000_000, "found": found}))
    """
)


def _run_budget_script() -> dict[str, float]:
    env = dict(os.environ, NUMBA_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(
        [sys.executable, "-c", SCRIPT], capture_output=True, text=True, env=env, timeout=900
    )
    assert proc.returncode == 0, proc.stderr
    result: dict[str, float] = json.loads(proc.stdout.strip().splitlines()[-1])
    return result


def test_verify_1e6_calls_12aa_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(_run_budget_script, rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] <= 5.0, result
