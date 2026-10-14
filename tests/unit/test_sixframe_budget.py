"""Performance budget for the six-frame index (issue M3-02 acceptance; decision 20).

``build_sixframe`` on a 5-Mb random genome: <= 3 s and <= 100 MB peak on one core. Measured in a
fresh single-threaded subprocess; numba compilation is warmed up on a tiny genome first and the
peak RSS increase is taken relative to the state just before the build.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap

SCRIPT = textwrap.dedent(
    """
    import json, resource, time
    import numpy as np
    from giap.mapping.sixframe import build_sixframe
    from giap.model import Replicon

    build_sixframe({"w": "ATGAAATAG" * 50}, {"w": Replicon("w", 450, "circular")})
    L = 5_000_000
    codes = np.random.default_rng(0).integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    genome = {"chr": seq}
    reps = {"chr": Replicon("chr", L, "circular")}
    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()
    idx = build_sixframe(genome, reps)
    secs = time.perf_counter() - t0
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert sum(len(f.aa) for f in idx.frames) > 9_000_000
    print(json.dumps({"secs": secs, "peak_mb": (after - before) / 1024}))
    """
)


def test_build_sixframe_5mb_budget() -> None:
    env = dict(os.environ, NUMBA_NUM_THREADS="1", OMP_NUM_THREADS="1")
    proc = subprocess.run(
        [sys.executable, "-c", SCRIPT], capture_output=True, text=True, env=env, timeout=600
    )
    assert proc.returncode == 0, proc.stderr
    result = json.loads(proc.stdout.strip().splitlines()[-1])
    assert result["secs"] <= 3.0, result
    assert result["peak_mb"] <= 100.0, result
