"""Performance budget for frameshift splits (issue M3-08 acceptance; GIPPA_SPEC.md §7.2.4).

"Benchmark: adds <= 50 % to the M3-07 benchmark runtime." The M3-07 inputs (benchmarks/
test_tolerant.py: 10,000 leftover peptides of 10-30 aa on a 5-Mb random circular genome and 1,000
on 100 kb, target + one reverse decoy type) with 5 % of the peptides replaced by planted ±1-nt
splits (both strands, δ = ±1). In one fresh subprocess with ``NUMBA_NUM_THREADS=8`` both matchers
are warmed up on a tiny genome, then ``match_tolerant`` and ``match_frameshift`` run on the same
arguments; the split matcher's wall time must be <= 0.5 x the tolerant matcher's.

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path
from typing import Any

SCRIPT = textwrap.dedent(
    """
    import json, resource, sys, time
    sys.path.insert(0, sys.argv[3])
    import numpy as np
    from fs_genomes import codons_for, join_split, rc
    from gippa.config.model import TolerantConfig
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.frameshift import match_frameshift
    from gippa.mapping.sixframe import build_sixframe, decode
    from gippa.mapping.tolerant import match_tolerant
    from gippa.model import Replicon
    from gippa.rng import for_step

    L, N = int(sys.argv[1]), int(sys.argv[2])
    cfg = TolerantConfig()
    tiny = build_sixframe({"w": "ATGAAAGCTGGTTGGTAG" * 40}, {"w": Replicon("w", 720, "circular")})
    tdec = {"reverse": build_decoys(tiny, "reverse", 50, for_step("sixframe_index", 0))}
    warm = ["MKAGWMKAGWMK", "KAGWMKAGWM", "AAAAAAAAAAAA"]
    for fn in (match_tolerant, match_frameshift):
        fn(warm, tiny, tdec, dict.fromkeys(warm, "reverse"), set(warm), cfg)
    rng = np.random.default_rng(0)
    codes = rng.integers(0, 4, size=L, dtype=np.uint8)
    seq = np.frombuffer(b"ACGT", dtype=np.uint8)[codes].tobytes().decode("ascii")
    del codes
    residues = "ACDEFGHKLMNPQRSTVWY"
    planted = []
    n_fs = max(N // 20, 5)
    step = L // (n_fs + 1)
    for j in range(n_fs):
        while True:
            gp = "".join(rng.choice(list(residues), size=int(rng.integers(4, 16))))
            gs = "".join(rng.choice(list(residues), size=int(rng.integers(4, 16))))
            nt = join_split(codons_for(gp, rng), codons_for(gs, rng), (1, -1)[j % 2])
            if nt is not None and 10 <= len(gp) + len(gs) <= 30:
                break
        nt = nt if j % 3 else rc(nt)
        pos = step * (j + 1)
        seq = seq[:pos] + nt + seq[pos + len(nt):]
        planted.append(gp + gs)
    idx = build_sixframe({"chr": seq}, {"chr": Replicon("chr", L, "circular")})
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 1))}
    frames = [decode(f.aa) for f in idx.frames]
    peps = set(planted)
    n_low = max(N // 100, 3)
    for i in range(n_low):
        n = 10 + i % 21
        peps.add(["A" * n, ("AG" * n)[:n], ("GA" * n)[:n]][i % 3] + "K" * (i // 3 % 2))
    n_cut = len(peps) + int(0.6 * N)
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
    tol = match_tolerant(queries, idx, decoys, decoy_of, set(queries), cfg)
    t1 = time.perf_counter()
    fs = match_frameshift(queries, idx, decoys, decoy_of, set(queries), cfg)
    t2 = time.perf_counter()
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak = rss if sys.platform == "darwin" else rss * 1024
    found = set(fs.loc[fs["search"] == "target", "peptide"])
    assert len(found & set(planted)) >= 0.95 * len(planted), (len(found), len(planted))
    print(json.dumps({"tolerant_secs": t1 - t0, "frameshift_secs": t2 - t1, "peak_bytes": peak,
                      "rows": len(fs), "peptides": N}))
    """
)


def _run(genome_bp: int, n_peptides: int) -> dict[str, float]:
    env = dict(os.environ, NUMBA_NUM_THREADS="8", OMP_NUM_THREADS="8")
    helpers = str(Path(__file__).resolve().parents[1] / "tests" / "helpers")
    sim = str(Path(__file__).resolve().parents[1] / "tests" / "sim")
    env["PYTHONPATH"] = os.pathsep.join(p for p in (sim, env.get("PYTHONPATH", "")) if p)
    proc = subprocess.run(
        [sys.executable, "-c", SCRIPT, str(genome_bp), str(n_peptides), helpers],
        capture_output=True,
        text=True,
        env=env,
        timeout=3600,
    )
    assert proc.returncode == 0, proc.stderr
    result: dict[str, float] = json.loads(proc.stdout.strip().splitlines()[-1])
    return result


def test_frameshift_adds_at_most_half_on_5mb(benchmark: Any) -> None:
    result = benchmark.pedantic(_run, args=(5_000_000, 10_000), rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["frameshift_secs"] <= 0.5 * result["tolerant_secs"], result


def test_frameshift_adds_at_most_half_on_100kb(benchmark: Any) -> None:
    result = benchmark.pedantic(_run, args=(100_000, 1_000), rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["frameshift_secs"] <= 0.5 * result["tolerant_secs"], result
