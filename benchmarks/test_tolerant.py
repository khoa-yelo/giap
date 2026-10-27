"""Performance budget for tolerant matching (issue M3-07 acceptance; GIPPA_SPEC.md §7.2.4).

- 10,000 leftover peptides (10-30 aa) on a 5-Mb random circular genome, target + one decoy type
  (reverse): <= 5 min on 8 CPUs and <= 2 GB peak memory;
- 1,000 peptides on a 100-kb genome: <= 5 s (the calibration suite depends on it).

Peptides: 60 % cut from the genome's six frames with 0-2 substitutions (hits exist), 39 % random,
1 % low-complexity (poly-A, AG/GA repeats; many candidate diagonals). Six-frame index and decoy
construction are excluded from the time (numba warmed up first on a tiny genome) but included in
the peak; the q-gram indexes are built inside ``match_tolerant`` (lazily) and are timed. Measured
in a fresh subprocess with ``NUMBA_NUM_THREADS=8``; peak from ``ru_maxrss`` (KiB on Linux, bytes on
macOS).

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
    from gippa.config.model import TolerantConfig
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.sixframe import build_sixframe, decode
    from gippa.mapping.tolerant import match_tolerant
    from gippa.model import Replicon
    from gippa.rng import for_step

    L, N = int(sys.argv[1]), int(sys.argv[2])
    cfg = TolerantConfig()
    tiny = build_sixframe({"w": "ATGAAAGCTGGTTGGTAG" * 40}, {"w": Replicon("w", 720, "circular")})
    tdec = {"reverse": build_decoys(tiny, "reverse", 50, for_step("sixframe_index", 0))}
    warm = ["MKAGWMKAGWMK", "KAGWMKAGWM", "AAAAAAAAAAAA"]
    match_tolerant(warm, tiny, tdec, dict.fromkeys(warm, "reverse"), set(warm), cfg)
    rng = np.random.default_rng(0)
    codes = rng.integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    idx = build_sixframe({"chr": seq}, {"chr": Replicon("chr", L, "circular")})
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 1))}
    frames = [decode(f.aa) for f in idx.frames]
    residues = "ACDEFGHKLMNPQRSTVWY"
    peps = set()
    n_low = max(N // 100, 3)
    for i in range(n_low):
        n = 10 + i % 21
        peps.add(["A" * n, ("AG" * n)[:n], ("GA" * n)[:n]][i % 3] + "K" * (i // 3 % 2))
    n_cut = n_low + int(0.6 * N)
    while len(peps) < n_cut:
        fr = frames[int(rng.integers(6))]
        k = int(rng.integers(10, 31))
        a = int(rng.integers(0, len(fr) - k))
        win = list(fr[a : a + k])
        if "*" in win or "X" in win:
            continue
        for _ in range(int(rng.integers(0, 3 if k >= 12 else 2))):
            win[int(rng.integers(k))] = residues[int(rng.integers(19))]
        peps.add("".join(win))
    while len(peps) < N:
        k = int(rng.integers(10, 31))
        peps.add("".join(residues[i] for i in rng.integers(0, 19, size=k)))
    del frames
    queries = sorted(peps)
    decoy_of = dict.fromkeys(queries, "reverse")
    t0 = time.perf_counter()
    df = match_tolerant(queries, idx, decoys, decoy_of, set(queries), cfg)
    secs = time.perf_counter() - t0
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = rss if sys.platform == "darwin" else rss * 1024
    n_tgt = df.loc[df["search"] == "target", "peptide"].nunique()
    assert n_tgt >= int(0.6 * N), n_tgt
    print(json.dumps({"secs": secs, "peak_bytes": peak, "rows": len(df), "peptides": N}))
    """
)


def _run_budget_script(genome_bp: int, n_peptides: int) -> dict[str, float]:
    env = dict(os.environ, NUMBA_NUM_THREADS="8", OMP_NUM_THREADS="8")
    proc = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(genome_bp), str(n_peptides)],
        capture_output=True,
        text=True,
        env=env,
        timeout=3600,
    )
    assert proc.returncode == 0, proc.stderr
    result: dict[str, float] = json.loads(proc.stdout.strip().splitlines()[-1])
    return result


def test_match_tolerant_10k_on_5mb_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(
        _run_budget_script, args=(5_000_000, 10_000), rounds=1, iterations=1
    )
    benchmark.extra_info.update(result)
    assert result["secs"] <= 300.0, result
    assert result["peak_bytes"] <= 2 * 1024**3, result


def test_match_tolerant_1k_on_100kb_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(_run_budget_script, args=(100_000, 1_000), rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] <= 5.0, result
