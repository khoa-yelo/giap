"""Performance budget for the M5-06 island step (docs/specs/M5-06-islands.md "Acceptance").

4,500-gene circular chromosome (about 2,500 scored TUs), B = 1,000 bootstrap replicates and
1,000 block-length permutations: ``tu_series`` -> ``scan_replicons`` -> ``islands_frame`` < 20 s
(4 CPUs). The input is the hand-built null of the calibration suite (TU-level AR(1) abundance,
geometric TU sizes, detection independent of abundance), so the benchmark needs no annotation
tool. Run with ``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests" / "helpers"))
import islands_data as idat


def test_islands_4500_genes_b1000_budget(benchmark: Any) -> None:
    rng = np.random.default_rng(506)
    contigs = [idat.sim_contig(rng, "chr1", 4500, rho=0.5, tu_p=1 / 1.8)]
    cfg = idat.config(bootstrap=1000, permutations=1000, alpha=0.05)
    genes, tus, replicons = idat.frames(contigs)
    n_tus: list[int] = []

    def one() -> Any:
        series = idat.build_series(genes, tus, replicons, cfg)  # Blom scores, TU means, LOESS
        n_tus.append(len(series[0].tu_ids))
        return idat.step(genes, series, cfg, np.random.default_rng(1))

    def run() -> float:
        secs, out = timed(one)
        assert set(out[0].columns) == set(idat.ISLAND_COLUMNS)
        assert 2000 <= n_tus[-1] <= 3000
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    benchmark.extra_info["n_tus"] = n_tus[-1]
    assert secs < 20.0, f"island step took {secs:.1f} s"
