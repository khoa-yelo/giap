"""Performance budget for the six-frame index (issue M3-02 acceptance; decision 20).

``build_sixframe`` on a 5-Mb random genome: <= 3 s and <= 100 MB peak on one core. Measured in a
fresh single-threaded subprocess; numba compilation is warmed up on a tiny genome first and the
peak RSS increase is taken relative to the state just before the build. ``ru_maxrss`` is in KiB on
Linux but in bytes on macOS, so it is normalised to bytes before converting to MB.

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
    import json, resource, sys, time
    import numpy as np
    from gippa.mapping.sixframe import build_sixframe
    from gippa.model import Replicon

    def peak_bytes():
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss if sys.platform == "darwin" else rss * 1024

    build_sixframe({"w": "ATGAAATAG" * 50}, {"w": Replicon("w", 450, "circular")})
    L = 5_000_000
    codes = np.random.default_rng(0).integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    genome = {"chr": seq}
    reps = {"chr": Replicon("chr", L, "circular")}
    before = peak_bytes()
    t0 = time.perf_counter()
    idx = build_sixframe(genome, reps)
    secs = time.perf_counter() - t0
    after = peak_bytes()
    assert sum(len(f.aa) for f in idx.frames) > 9_000_000
    print(json.dumps({"secs": secs, "peak_mb": (after - before) / 2**20}))
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


def test_build_sixframe_5mb_budget(benchmark: Any) -> None:
    # One round: the subprocess is the unit; the budget is checked on the in-process timings.
    result = benchmark.pedantic(_run_budget_script, rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] <= 3.0, result
    assert result["peak_mb"] <= 100.0, result
