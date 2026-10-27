"""Calibration of the coincidence-null Monte Carlo SE (spec §8.4; PR #264 review item 1; lead
decision of 2026-10-23: c = mean of one call's R per-placement rates, SE = sd(ddof=1)/sqrt(R)).

Claim (gated): the normal 95 % interval c ± 1.96 SE covers the exact null expectation c₀ in a
share inside the binomial band 0.95 ± k·√(0.95·0.05/M), where M = ``replicates_fast.regions``
and k = ``se_multiplier``.

The toy uses the five null windows of ``tests/unit/stats/test_cutoffs_coincidence_rates.py``
(+: gA [950, 1500), gB [2350, 2950); -: gC [2400, 3000), h1 [4000, 5050), h2 [6000, 6550)).
It has three calls, each placed independently among the windows long enough to hold it:
- [1000, 1030) +, 50 bp past gA's lo, fits gA and gB. Only the gA placement [1000, 1030) is
  within 50 bp of truth, so p = 1/2.
- [5000, 5020) -, 30 bp before h1's hi. Of [2950, 2970), [5000, 5020) and [6500, 6520), only
  h1's matches, so p = 1/3.
- [1130, 1160) +, 180 bp past lo. Truth sits only at gB's placement [2530, 2560), so p = 1/2.

So c₀ = (1/2 + 1/3 + 1/2) / 3 = 4/9 exactly.
Teeth: an SE that divides by R instead of √R under-covers far below the band.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene
from gippa.rng import for_step
from gippa.stats.cutoffs import coincidence_rate, downstream_windows

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("regions")]

C = "chr1"
LIN = Replicon(C, 10_000, "linear")
LEVEL = 0.95
R_PLACEMENTS = 40
C0 = 4 / 9


def _gene(tag: str, start: int, end: int, strand: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(C, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def _iv(start: int, end: int, strand: str) -> Interval:
    return Interval(C, start, end, strand)  # type: ignore[arg-type]


GENES = [
    _gene("gA", 100, 1000, "+"),
    _gene("gB", 1500, 2400, "+"),
    _gene("gC", 2950, 4000, "-"),
    _gene("h1", 5000, 6000, "-"),
    _gene("h2", 6500, 7000, "-"),
]
CALLS = [_iv(1000, 1030, "+"), _iv(5000, 5020, "-"), _iv(1130, 1160, "+")]
CONF = [0.9, 0.8, 0.7]
TRUTH = [_iv(1000, 1030, "+"), _iv(5000, 5020, "-"), _iv(2530, 2560, "+")]


def test_coincidence_mc_se_coverage_and_teeth(tolerances: dict[str, object]) -> None:
    reps, seeds = tolerances["replicates_fast"], tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    m, k = int(reps["regions"]), float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    gen = for_step("calibration.m6_12.coincidence_se", int(seeds["calibration_base"]))
    windows = downstream_windows(GENES, LIN)
    z = 1.959963984540054
    cov = cov_teeth = 0
    cs, ses = [], []
    for _ in range(m):
        rates = coincidence_rate(
            CALLS, CONF, windows, TRUTH, LIN, 0.0, rng=gen, R=R_PLACEMENTS, return_rates=True
        )
        c = float(np.mean(rates))
        sd = float(np.std(rates, ddof=1))
        se, se_bad = sd / math.sqrt(R_PLACEMENTS), sd / R_PLACEMENTS
        cs.append(c)
        ses.append(se)
        cov += abs(c - C0) <= z * se
        cov_teeth += abs(c - C0) <= z * se_bad
    half = k * math.sqrt(LEVEL * (1 - LEVEL) / m)
    share, share_bad = cov / m, cov_teeth / m
    ratio = float(np.mean(ses) / np.std(cs, ddof=1))
    print(
        f"M6-12 coincidence MC SE (M={m}, R={R_PLACEMENTS}): coverage {share:.3f}, "
        f"sd/R teeth {share_bad:.3f}, mean SE / sd(c) {ratio:.3f}"
    )
    assert LEVEL - half <= share <= LEVEL + half, f"coverage {share:.3f} outside band"
    assert share_bad < LEVEL - half, f"sd/R coverage {share_bad:.3f} should under-cover (teeth)"
