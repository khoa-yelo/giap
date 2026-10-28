"""Fast calibration of the M6-12 cut-off statistics (spec §8.4; issue M6-12 part A; methods
docs/methods/operons.md "K-12 cut-off calibration (M6-12)").

Claims checked:
1. **Directon AUROC CI** (gated). Pairs cluster in directons: a shared directon effect on the
   scores and directon-homogeneous labels (an operon gives several positive pairs). The true
   AUROC is the between-directon value Φ(μ/√(2 + 2 s_u^2)); within-directon pos/neg pairs are a
   negligible share. The 95 % directon CI covers it in a share within the binomial band
   0.95 ± k·√(0.95·0.05/R), R = ``replicates_fast.regions``, k = ``se_multiplier``.
   Teeth: the single-pair CI (each pair its own cluster) under-covers below that band.
2. **Refit consistency** (gated). On Bernoulli data from (d₀, s) = (40, 15), the mean refitted d₀
   and s over R replicates lie within k standard errors of the truth.
   Recorded, not gated: the percentile coverage of the d₀ and s CIs (warning outside
   0.95 ± 0.04; a scratch run at R = 300 gave 0.93 / 0.92, the known small-sample deficit of
   percentile CIs for a ratio estimator).
3. **Wilson 95 %** (gated, exact). The exact binomial coverage of ``wilson_ci`` averaged over
   p ∈ [0.05, 0.95] is 0.95 ± 0.005 at n ∈ {50, 100, 727, 3217} (727 and 3,217 are the
   design-check TerminatorNet sizes).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from scipy.stats import binom, norm

from gippa.rng import for_step
from gippa.stats.cutoffs import auroc_ci, distance_fit_ci, wilson_ci

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("regions")]

LEVEL = 0.95
MU, SIGMA_U = 1.5, 1.0
TRUE_AUROC = float(norm.cdf(MU / math.sqrt(2.0 + 2.0 * SIGMA_U**2)))


def _ints(tolerances: dict[str, object]) -> tuple[int, float, int]:
    reps, seeds = tolerances["replicates_fast"], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    return int(reps["regions"]), float(tolerances["se_multiplier"]), int(seeds["calibration_base"])  # type: ignore[arg-type]


def clustered_pairs(
    g: np.random.Generator, n_dir: int = 150
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sizes = g.integers(1, 7, n_dir)
    cl = np.repeat(np.arange(n_dir), sizes)
    cls = (g.random(n_dir) < 0.6)[cl]
    flip = g.random(cl.size) < 0.1
    y = np.where(flip, ~cls, cls).astype(int)
    u = g.normal(0.0, SIGMA_U, n_dir)[cl]
    s = MU * y + u + g.normal(0.0, 1.0, cl.size)
    return s, y, cl


def test_directon_auroc_ci_coverage_and_iid_teeth(tolerances: dict[str, object]) -> None:
    r, k, seed = _ints(tolerances)
    gen = for_step("calibration.m6_12.auroc_ci", seed)
    cov_dir = cov_iid = 0
    for _ in range(r):
        s, y, cl = clustered_pairs(gen)
        a = auroc_ci(s, y, cl, rng=gen, B=400)
        b = auroc_ci(s, y, np.arange(y.size), rng=gen, B=400)
        cov_dir += a.ci_low <= TRUE_AUROC <= a.ci_high
        cov_iid += b.ci_low <= TRUE_AUROC <= b.ci_high
    half = k * math.sqrt(LEVEL * (1 - LEVEL) / r)
    cd, ci = cov_dir / r, cov_iid / r
    print(f"M6-12 AUROC CI coverage (R={r}): directon {cd:.3f}, single-pair {ci:.3f}")
    assert LEVEL - half <= cd <= LEVEL + half, f"directon coverage {cd:.3f} outside band"
    assert ci < LEVEL - half, f"single-pair coverage {ci:.3f} should under-cover (teeth)"


def test_refit_consistent_and_ci_coverage_recorded(tolerances: dict[str, object]) -> None:
    r, k, seed = _ints(tolerances)
    gen = for_step("calibration.m6_12.refit", seed)
    d0s: list[float] = []
    ss: list[float] = []
    cov0 = cov1 = 0
    for _ in range(r):
        n_dir = 400
        cl = np.repeat(np.arange(n_dir), gen.integers(1, 5, n_dir))
        d = gen.integers(-20, 201, cl.size)
        y = (gen.random(cl.size) < 1.0 / (1.0 + np.exp((d - 40.0) / 15.0))).astype(int)
        a, b = distance_fit_ci(d, y, cl, rng=gen, B=100)
        assert a.n_valid == b.n_valid == 100
        d0s.append(a.estimate)
        ss.append(b.estimate)
        cov0 += a.ci_low <= 40.0 <= a.ci_high
        cov1 += b.ci_low <= 15.0 <= b.ci_high
    for name, vals, truth in (("d0", d0s, 40.0), ("s", ss, 15.0)):
        v = np.asarray(vals)
        se = float(v.std(ddof=1)) / math.sqrt(r)
        print(f"M6-12 refit {name}: mean {v.mean():.3f} (truth {truth}), SE {se:.3f}")
        assert abs(float(v.mean()) - truth) <= k * se + 1e-9, f"{name} biased"
    c0, c1 = cov0 / r, cov1 / r
    print(f"recorded: M6-12 refit CI coverage (R={r}, B=100): d0 {c0:.3f}, s {c1:.3f}")
    for name, c in (("d0", c0), ("s", c1)):
        if abs(c - LEVEL) > 0.04:
            warnings.warn(
                f"M6-12 refit {name} CI coverage {c:.3f} outside 0.95 +- 0.04", stacklevel=1
            )


@pytest.mark.parametrize("n", [50, 100, 727, 3217])
def test_wilson_average_exact_coverage(n: int) -> None:
    k = np.arange(n + 1)
    bounds = [wilson_ci(int(j), n) for j in k]
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    ps = np.linspace(0.05, 0.95, 181)
    cov = np.array([binom.pmf(k, n, p)[(lo <= p) & (p <= hi)].sum() for p in ps])
    print(f"Wilson n={n}: mean exact coverage {cov.mean():.4f}, min {cov.min():.4f}")
    assert abs(float(cov.mean()) - LEVEL) <= 0.005
