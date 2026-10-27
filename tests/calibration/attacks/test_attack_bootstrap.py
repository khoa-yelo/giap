"""Attack nulls for the M4-04 TU-block bootstrap (GIPPA_SPEC.md §8.2, §8.5; statistician-owned).

Each generator stresses the TU-block percentile CI (docs/methods/bootstrap.md) while staying inside
its assumptions (TUs independent, exchangeable within a stratum):

- ``few_large_tus``: 30 TUs of 5-15 genes with intra-TU correlation 0.9; statistic = mean over
  genes (a ratio estimator, since TU sizes vary). Effective sample size is ~30, not ~300.
- ``hl_two_group``: the §7.3.5(b) use: Hodges-Lehmann shift of a region's genes (12 TUs of 2-6
  genes) against the outside genes (80 TUs), intra-TU correlation 0.8, both groups resampled
  within strata. True shift 0.

Criterion. Percentile CIs from n independent units have a known small-sample deficit: for a mean
of n normal units the bootstrap variance is (n-1)/n of the plug-in one and the pivot is t, not z,
so the expected coverage is c(n) = P(|t_{n-1}| <= z_{0.975} sqrt((n-1)/n)) (0.936 at n = 30,
0.913 at n = 12). Gate: coverage >= c(n) - k sqrt(c(n)(1 - c(n))/R), with n the smallest number
of independently resampled TUs, R = ``replicates_*.regions`` and k = ``se_multiplier``. The gate
is the level check on CI misses at alpha = 1 - c(n) (the same inequality) through
``_criteria.calibrate`` (ADR-0015 two-block confirmation, sharding), stream
``regions.attack_bootstrap.<attack>``.

``test_gene_level_bootstrap_undercovers`` is the teeth check: resampling genes instead of TUs on
the ``few_large_tus`` data must fall far below the band: a confirmed ``fail`` through ``calibrate``
(stream ``regions.attack_bootstrap.naive_gene``, PR #69 decision 3). ``test_recorded_few_tus``
records (does not gate) coverage at 5 and 10 TUs against c(n), the documented regime where
callers flag CIs.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import NOT_FAILED
from scipy.stats import t as student_t

from gippa.rng import for_step
from gippa.stats.bootstrap import BootCI, tu_block_bootstrap
from gippa.stats.robust import hodges_lehmann

LEVEL = 0.95
B = 1000
Z = 1.959963984540054

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


def expected_percentile_coverage(n: int, level: float = LEVEL) -> float:
    """c(n): coverage of the percentile CI for a mean of n iid normal units."""
    return float(2.0 * student_t.cdf(Z * math.sqrt((n - 1) / n), n - 1) - 1.0)


def _clustered(rng: np.random.Generator, n_tu: int, lo: int, hi: int, icc: float) -> pd.DataFrame:
    sizes = rng.integers(lo, hi + 1, n_tu)
    tu = np.repeat(np.arange(n_tu), sizes)
    v = math.sqrt(icc) * rng.normal(size=n_tu)[tu] + math.sqrt(1 - icc) * rng.normal(size=tu.size)
    return pd.DataFrame({"tu": tu, "gene": np.arange(tu.size), "v": v})


def _mean(d: pd.DataFrame) -> float:
    return float(d["v"].to_numpy().mean())


def _hl(d: pd.DataFrame) -> float:
    inside = d["inside"].to_numpy()
    v = d["v"].to_numpy()
    return hodges_lehmann(v[inside], v[~inside])


@dataclass(frozen=True)
class Attack:
    n_min: int  # smallest number of independently resampled TUs (sets c(n))
    run: Callable[[np.random.Generator], BootCI]


def _few_large_tus(rng: np.random.Generator) -> BootCI:
    return tu_block_bootstrap(_clustered(rng, 30, 5, 15, 0.9), "tu", _mean, B=B, rng=rng)


def _hl_two_group(rng: np.random.Generator) -> BootCI:
    n_in, n_out = 12, 80
    d = _clustered(rng, n_in + n_out, 2, 6, 0.8)
    d["inside"] = d["tu"].to_numpy() < n_in
    return tu_block_bootstrap(d, "tu", _hl, B=B, rng=rng, strata_col="inside")


ATTACKS: dict[str, Attack] = {
    "few_large_tus": Attack(30, _few_large_tus),
    "hl_two_group": Attack(12, _hl_two_group),
}


def _settings(tolerances: dict[str, object], suite: str) -> tuple[int, int, float]:
    reps, seeds, se_k = (tolerances[key] for key in (suite, "seeds", "se_multiplier"))
    assert isinstance(reps, dict) and isinstance(seeds, dict) and isinstance(se_k, (int, float))
    return int(reps["regions"]), int(seeds["calibration_base"]), float(se_k)


def _lower_bound(c: float, r: int, k: float) -> float:
    return c - k * math.sqrt(c * (1 - c) / r)


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_tu_block_ci_coverage_under_attack(
    tolerances: dict[str, object], attack: str, suite: str
) -> None:
    r, _, _ = _settings(tolerances, suite)
    a = ATTACKS[attack]
    c = expected_percentile_coverage(a.n_min)

    def missed(rng: np.random.Generator) -> float:
        """0 (a 'rejection') when the CI misses the true value 0, else 1."""
        ci = a.run(rng)
        assert ci.n_valid == B and math.isfinite(ci.ci_low) and ci.ci_low <= ci.ci_high
        return 1.0 if ci.ci_low <= 0.0 <= ci.ci_high else 0.0

    # coverage >= c - k sqrt(c(1-c)/R)  <=>  miss rate <= (1-c) + k sqrt(c(1-c)/R) = level bound
    rec = C.calibrate(
        "regions", "level", 1.0 - c, missed, r, stream=f"regions.attack_bootstrap.{attack}"
    )
    print(f"attack {attack}: c(n) = {c:.3f}, R = {r}: {rec}")
    assert rec["status"] in NOT_FAILED, rec


@pytest.mark.parametrize("suite", SUITES)
def test_gene_level_bootstrap_undercovers(
    tolerances: dict[str, object], suite: str, teeth_log: Path
) -> None:
    """Teeth: resampling genes (each gene its own 'TU') ignores intra-TU correlation; through
    ``calibrate`` it is a confirmed ``fail`` (both blocks, PR #69 decision 3)."""
    r, _, _ = _settings(tolerances, suite)
    c = expected_percentile_coverage(30)

    def missed(rng: np.random.Generator) -> float:
        ci = tu_block_bootstrap(_clustered(rng, 30, 5, 15, 0.9), "gene", _mean, B=B, rng=rng)
        return 1.0 if ci.ci_low <= 0.0 <= ci.ci_high else 0.0

    rec = C.calibrate(
        "regions", "level", 1.0 - c, missed, r, stream="regions.attack_bootstrap.naive_gene",
        log_dir=teeth_log,
    )  # fmt: skip
    print(f"gene-level bootstrap on few_large_tus: {rec}")
    assert rec["status"] == "fail", f"gene-level bootstrap not a confirmed fail: {rec}"


@pytest.mark.calibration
@pytest.mark.parametrize("n_tu", [5, 10])
def test_recorded_few_tus(tolerances: dict[str, object], n_tu: int) -> None:
    """Recorded, not gated: percentile under-coverage with very few TUs tracks c(n)."""
    _, seed, _ = _settings(tolerances, "replicates_fast")
    rng = for_step(f"attack.bootstrap.few_tus.{n_tu}", seed)
    r = 200
    covered = 0
    for _ in range(r):
        ci = tu_block_bootstrap(_clustered(rng, n_tu, 1, 4, 0.9), "tu", _mean, B=B, rng=rng)
        covered += ci.ci_low <= 0.0 <= ci.ci_high
    cov = covered / r
    c = expected_percentile_coverage(n_tu)
    print(f"recorded: {n_tu} TUs, coverage {cov:.3f}, c(n) = {c:.3f}, R = {r}")
    if cov < _lower_bound(c, r, 2.0):
        warnings.warn(f"{n_tu} TUs: coverage {cov:.3f} below c(n) = {c:.3f}", stacklevel=1)
