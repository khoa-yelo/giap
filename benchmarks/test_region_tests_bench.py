"""Performance budget for the M4-17 region tests (docs/specs/M4-17-region-tests.md "Performance").

4,500-gene circular chromosome (gippasim M4-07 layout + M4-08 truth), 30 regions (8 prophages,
8 BGCs, 8 CRISPR arrays, 6 AMR sets), exhaustive shifts (K = 4,499 per region and statistic),
B = 1,000 TU-block bootstrap replicates: ``test_all_regions`` < 60 s (4 CPUs). Run with
``pytest benchmarks --benchmark-only``.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests" / "helpers"))
import region_tests_data as rd


def _regions() -> tuple[Any, ...]:
    import gippasim.layout as la

    spec = la.RegionSpec
    return (
        *(spec("prophage", 40) for _ in range(8)),
        *(spec("bgc", 25) for _ in range(8)),
        *(spec("crispr", 6) for _ in range(8)),
        *(spec("amr", 10, scattered=True) for _ in range(6)),
    )


def test_region_tests_4500_genes_30_regions_budget(benchmark: Any) -> None:
    import gippasim.layout as la

    lc = la.LayoutConfig(n_genes_chromosome=4500, regions=_regions())
    tables = rd.sim_tables(np.random.default_rng(45), lc, {"ar1_rho": 0.5})
    assert int((tables.regions["region_type"] != "plasmid").sum()) == 30
    cfg = rd.config(bootstrap=1000, permutations=1000)

    def run() -> float:
        secs, res = timed(lambda: rd.run_region_tests(tables, cfg, np.random.default_rng(1)))
        assert len(res) == 60
        assert (res["n_perm"] == 4499).sum() >= 30  # every O - E row is exhaustive
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs < 60.0, f"region tests took {secs:.1f} s"
