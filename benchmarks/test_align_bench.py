"""Performance budget for the residue aligner (issue M7-01 acceptance; GIAP_SPEC.md §7.5.2).

100,000 tolerant hits (7-30 aa, both strands, one substitution, two separated substitutions or
one adjacent swap planted in a random 200-kb circular genome) aligned by ``align_hits`` in < 30 s
on one core, numba warmed up first. Run in a fresh single-threaded subprocess.

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
    import pandas as pd
    from giap.model import Replicon
    from giap.variants.align import align_hits, translate_locus

    rng = np.random.default_rng(0)
    n = 200_000
    seq = "".join(rng.choice(list("ACGT"), size=n))
    rep = Replicon("c1", n, "circular")
    aa = "ACDEFGHKLMNPQRSTVWY"
    rows = []
    while len(rows) < 100_000:
        k = int(rng.integers(7, 31))
        s = int(rng.integers(0, n))
        strand = "+" if rng.random() < 0.5 else "-"
        g = translate_locus(seq, rep, s, s + 3 * k, strand)
        if "*" in g or "X" in g:
            continue
        p = list(g)
        mode = len(rows) % 3
        if mode == 2:
            j = int(rng.integers(0, k - 1))
            if p[j].replace("I", "L") == p[j + 1].replace("I", "L"):
                continue
            p[j], p[j + 1] = p[j + 1], p[j]
            nd = 1
        else:
            js = [int(rng.integers(0, k))] if mode == 0 else [1, k - 2]
            for j in js:
                p[j] = rng.choice([c for c in aa if c != p[j].replace("I", "L")])
            nd = len(js)
        rows.append({"peptide": "".join(p), "contig": "c1", "start": s, "end": s + 3 * k,
                     "strand": strand, "route": "tolerant", "n_diff": nd,
                     "diff_types": "substitution"})
    ph = pd.DataFrame(rows)
    align_hits(ph.head(50), None, {"c1": seq}, {"c1": rep})  # warm-up (numba cache)
    t0 = time.perf_counter()
    got = align_hits(ph, None, {"c1": seq}, {"c1": rep})
    secs = time.perf_counter() - t0
    assert len(got) == len(ph)
    print(json.dumps({"secs": secs, "hits": len(got)}))
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


def test_align_hits_100k_tolerant_budget(benchmark: Any) -> None:
    result = benchmark.pedantic(_run_budget_script, rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["secs"] < 30.0, result
