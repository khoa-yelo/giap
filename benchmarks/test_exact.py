"""Performance budget for exact matching (issue M3-04 acceptance; GIPPA_SPEC.md §7.2.2).

100,000 random tryptic peptides (7-40 aa, ending in K/R; half cut from the genome's frames so
that hits exist, half random) on a 5-Mb random circular genome, target + reverse decoy:
<= 30 s and <= 2 GB peak memory on one core. Index and decoy construction are excluded from the
time (numba warmed up first) but included in the peak. Measured in a fresh single-threaded
subprocess; peak from ``ru_maxrss`` (KiB on Linux, bytes on macOS).

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
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.mapping.sixframe import build_sixframe, decode
    from gippa.model import Replicon
    from gippa.rng import for_step

    tiny = build_sixframe({"w": "ATGAAATAG" * 50}, {"w": Replicon("w", 450, "circular")})
    build_decoys(tiny, "reverse", 50, for_step("sixframe_index", 0))
    L = 5_000_000
    rng = np.random.default_rng(0)
    codes = rng.integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    idx = build_sixframe({"chr": seq}, {"chr": Replicon("chr", L, "circular")})
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 1))}
    frames = [decode(f.aa) for f in idx.frames]
    residues = "ACDEFGHIKLMNPQRSTVWY"
    peps = set()
    while len(peps) < 50_000:
        fr = frames[int(rng.integers(6))]
        k = int(rng.integers(7, 41))
        a = int(rng.integers(0, len(fr) - k))
        win = fr[a : a + k]
        if "*" not in win and "X" not in win:
            peps.add(win)
    while len(peps) < 100_000:
        k = int(rng.integers(7, 41))
        body = "".join(residues[i] for i in rng.integers(0, 20, size=k - 1))
        peps.add(body + "KR"[int(rng.integers(2))])
    del frames
    queries = sorted(peps)
    decoy_of = dict.fromkeys(queries, "reverse")
    t0 = time.perf_counter()
    df = match_exact(queries, idx, decoys, decoy_of)
    secs = time.perf_counter() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = rss if sys.platform == "darwin" else rss * 1024
    assert (df["search"] == "target").sum() >= 50_000
    print(json.dumps({"secs": secs, "peak_bytes": peak, "rows": len(df)}))
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


def test_match_exact_100k_on_5mb_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(_run_budget_script, rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] <= 30.0, result
    assert result["peak_bytes"] <= 2 * 1024**3, result
