"""Fast coverage calibration of the TU-block percentile bootstrap CI (issue M4-04; §7.3.5(a),
§8.2, §8.4; decision 14: B = 1000 by default).

Gated: under a null with TU-structured noise (each TU shares a random effect across its 1-4
genes, mean zero, so the true expected sum is 0), the 95 % TU-block percentile CI for the sum of
gene values over 100 TUs covers 0 in a share of replicates within the binomial band
0.95 +- k * sqrt(0.95 * 0.05 / R), with R = ``replicates_fast.regions`` (the region-effect CIs of
§7.3.5 use this bootstrap; lead decision 6: reuse existing counts), k = ``se_multiplier``, seed
``seeds.calibration_base``. Gene-level resampling would under-cover badly here (positive
within-TU correlation), so the band also checks that whole TUs are resampled.

Recorded, not gated (spec M4-04 acceptance): coverage for a sum over 20 TUs of iid N(0, 1) TU
effects, R = 200; a warning flags it when outside 0.95 +- 0.04 (percentile CIs under-cover with
few TUs; see docs/methods/bootstrap.md).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import pytest

from giap.rng import for_step
from giap.stats.bootstrap import tu_block_bootstrap

pytestmark = pytest.mark.calibration

LEVEL = 0.95
B = 1000
N_TU = 100


def _sum(d: pd.DataFrame) -> float:
    return float(d["v"].to_numpy().sum())


def _tu_null(rng: np.random.Generator, n_tu: int) -> pd.DataFrame:
    sizes = rng.integers(1, 5, n_tu)
    tu = np.repeat(np.arange(n_tu), sizes)
    effect = rng.normal(0.0, 1.0, n_tu)[tu]  # shared within TU
    noise = rng.normal(0.0, 0.5, tu.size)
    return pd.DataFrame({"tu": tu, "v": effect + noise})


def test_tu_block_ci_coverage_within_binomial_band(tolerances: dict[str, object]) -> None:
    reps = tolerances["replicates_fast"]
    seeds = tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    r = int(reps["regions"])
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    gen = for_step("calibration.tu_block_bootstrap", int(seeds["calibration_base"]))
    covered = 0
    for _ in range(r):
        df = _tu_null(gen, N_TU)
        ci = tu_block_bootstrap(df, "tu", _sum, B=B, rng=gen, level=LEVEL)
        assert ci.n_valid == B
        covered += ci.ci_low <= 0.0 <= ci.ci_high
    cov = covered / r
    half = k * math.sqrt(LEVEL * (1 - LEVEL) / r)
    print(f"TU-block percentile CI coverage (100 TUs, R={r}, B={B}): {cov:.3f}")
    assert LEVEL - half <= cov <= LEVEL + half, (
        f"coverage {cov:.3f} outside {LEVEL} +- {half:.3f} over R={r}"
    )


def test_recorded_coverage_twenty_tus(tolerances: dict[str, object]) -> None:
    seeds = tolerances["seeds"]
    assert isinstance(seeds, dict)
    r = 200
    gen = for_step("calibration.tu_block_bootstrap.20tu", int(seeds["calibration_base"]))
    covered = 0
    for _ in range(r):
        df = pd.DataFrame({"tu": np.arange(20), "v": gen.normal(0.0, 1.0, 20)})
        ci = tu_block_bootstrap(df, "tu", _sum, B=B, rng=gen, level=LEVEL)
        covered += ci.ci_low <= 0.0 <= ci.ci_high
    cov = covered / r
    print(f"recorded: TU-block percentile CI coverage (20 iid TUs, R={r}, B={B}): {cov:.3f}")
    if abs(cov - LEVEL) > 0.04:
        warnings.warn(
            f"M4-04 recorded coverage {cov:.3f} outside 0.95 +- 0.04 (20 TUs)", stacklevel=1
        )
