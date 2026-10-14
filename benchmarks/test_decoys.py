"""Performance budget for decoy translations (issue M3-03 acceptance; §6.6).

Both decoy sets (``reverse`` and ``shuffle``, window 50) for a 5-Mb random genome: <= 3 s on one
core, excluding ``build_sixframe`` and numba compilation (warmed up on a tiny genome first).
Measured in a fresh single-threaded subprocess.

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
    from giap.mapping.decoys import build_decoys
    from giap.mapping.sixframe import build_sixframe
    from giap.model import Replicon
    from giap.rng import for_step

    tiny = build_sixframe({"w": "ATGAAATAG" * 50}, {"w": Replicon("w", 450, "circular")})
    build_decoys(tiny, "reverse", 50, for_step("sixframe_index", 0))
    build_decoys(tiny, "shuffle", 50, for_step("sixframe_index", 0))
    L = 5_000_000
    codes = np.random.default_rng(0).integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    idx = build_sixframe({"chr": seq}, {"chr": Replicon("chr", L, "circular")})
    rng = for_step("sixframe_index", 1)
    t0 = time.perf_counter()
    rev = build_decoys(idx, "reverse", 50, rng)
    shu = build_decoys(idx, "shuffle", 50, rng)
    secs = time.perf_counter() - t0
    assert len(rev) == len(shu) == 6
    print(json.dumps({"secs": secs}))
    """
)


def _run_budget_script() -> dict[str, float]:
    env = dict(os.environ, NUMBA_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(
        [sys.executable, "-c", SCRIPT], capture_output=True, text=True, env=env, timeout=600
    )
    assert proc.returncode == 0, proc.stderr
    result: dict[str, float] = json.loads(proc.stdout.strip().splitlines()[-1])
    return result


def test_build_decoys_5mb_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(_run_budget_script, rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] <= 3.0, result
