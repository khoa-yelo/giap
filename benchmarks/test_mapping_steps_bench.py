"""Performance budget for the mapping steps (issue M3-17 acceptance; GIAP_SPEC.md §6.4, §7.2).

The whole mapping chain (``parse_bakta`` → ``sixframe_index`` → ``map_exact`` → ``map_tolerant``
→ ``placement`` → ``loci``) through ``giap.steps.registry.run_step`` on the 100-kb giapsim
genome with 5,000 query peptides (the simulator's M3-19 peptides, topped up with windows of the
genome's six-frame translation and random tryptic peptides, all ``q = 0.001``): <= 20 s wall
time (local executor, no external tools). Simulation and preflight are excluded; numba is warmed
up first on the 50-kb genome.

Not part of the per-commit suite (``testpaths = ["tests"]``); run with
``pytest benchmarks --benchmark-only``. Recorded in ``benchmarks/baseline.json``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))

from helpers.mapping_steps import (
    SimRun,
    drive,
    fresh_steps,
    roots_of,
    simulate_with_proteomics,
)

from giap.config import resolve_config
from giap.io.fasta import read_genome_fasta
from giap.steps.mapping import select_queries

N_PEPTIDES = 5_000
BUDGET_S = 20.0
_REV = str.maketrans("ACGT", "TGCA")
_CODE = dict(
    zip(
        [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        strict=True,
    )
)


def _frames(seq: str) -> list[str]:
    rc = seq.translate(_REV)[::-1]
    return [
        "".join(_CODE.get(s[i : i + 3], "X") for i in range(f, len(s) - 2, 3))
        for s in (seq, rc)
        for f in range(3)
    ]


def _five_thousand(run: SimRun, out: Path, seed: int) -> SimRun:
    """The simulator's peptides plus genomic windows and random peptides, 5,000 queries."""
    df = pd.read_parquet(run.peptides)
    rng = np.random.default_rng(seed)
    frames = [f for r in read_genome_fasta(run.genome) for f in _frames(r.sequence)]
    have = set(df["peptide"])
    extra: list[str] = []
    n_query = len(select_queries(df, resolve_config(run.config_yaml)))
    while n_query + len(extra) < N_PEPTIDES:
        k = int(rng.integers(7, 31))
        if len(extra) % 2 == 0:
            fr = frames[int(rng.integers(len(frames)))]
            a = int(rng.integers(0, len(fr) - k))
            pep = fr[a : a + k]
            if "*" in pep or "X" in pep:
                continue
        else:
            body = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=k - 1))
            pep = body + "KR"[int(rng.integers(2))]
        if pep not in have:
            have.add(pep)
            extra.append(pep)
    base = df[~df["is_decoy"]].iloc[0].to_dict()
    rows = [{**base, "peptide": p, "q_value": 0.001, "pep": 0.001} for p in extra]
    path = out / "peptides.5k.parquet"
    pd.concat([df, pd.DataFrame(rows)], ignore_index=True).to_parquet(path, index=False)
    return SimRun(run.res, run.genome, path, run.proteins_long, run.config_yaml)


def _chain(tmp: Path, run: SimRun, seed: int) -> float:
    cfg = resolve_config(run.config_yaml, [f"seed={seed}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        t0 = time.perf_counter()
        drive(reg, cfg, tmp, roots_of(run))
        return time.perf_counter() - t0


def _run(tmp: Path) -> dict[str, float]:
    warm = simulate_with_proteomics(tmp / "warm", 1)
    _chain(tmp / "warm_work", warm, 1)
    run = simulate_with_proteomics(tmp / "in", 5, length=100_000)
    run = _five_thousand(run, tmp, 5)
    secs = _chain(tmp / "work", run, 5)
    n = len(select_queries(pd.read_parquet(run.peptides), resolve_config(run.config_yaml)))
    return {"secs": secs, "queries": n}


def test_mapping_chain_100kb_5000_peptides_budget(benchmark: Any, tmp_path: Path) -> None:
    result = benchmark.pedantic(lambda: _run(tmp_path), rounds=1, iterations=1)
    benchmark.extra_info.update(result)
    assert result["queries"] >= N_PEPTIDES
    assert result["secs"] <= BUDGET_S, result
