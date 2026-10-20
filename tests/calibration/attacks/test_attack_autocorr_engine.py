"""Attack nulls for the M5-02 autocorrelation engine (GIAP_SPEC.md §8.5; statistician-owned).

Proposed by the statistician review of PR #116 for tests/calibration/attacks/. Each generator stays
inside the null each engine claims and tries to push the lag test above level:

- ``jc_low_p_islands_label`` / ``jc_low_p_islands_tu``: detection d_i ~ Bernoulli(p_i) independently
  given p (the join-count null), with p ~ U(0.3, 0.9) except 30-gene islands covering 20 % of a
  1,000-gene circular replicon where p ~ U(0.01, 0.05) (§8.5 attack 2: detectability confounded with
  region; pathogenic E. coli carry 15-25 % of genes in low-expression islands). ADR-0020's residual
  Σ(d - p)(d - p) has variance Σ v_i v_j with v = p(1 - p); the permutation null averages v over
  random pairs, so clustered v makes the null too narrow (0.072-0.074 at α 0.05, R = 2,000,
  Slurm 1861073). The Pearson residual (d - p)/sqrt(p(1 - p)) is at level (0.051-0.054).
- ``moran_tu_large_heavy``: TU sizes 1 + Geometric(0.15) (mean ≈ 6.7), t3 TU effects and t3 gene
  noise, the series rotated so a TU crosses the origin; TU-block null, lags 1 and 20.
- ``moran_linear_missing``: linear contig, 600 genes, 45 % missing; TU null (gaps move with TUs)
  and label null with ``movable`` = genes with values; lags 1 and 20.
- ``moran_small_plasmid``: 40 genes in ≈ 13 TUs, circular, TU null (coarse permutation support).

Criterion: §8.4 level through ``_gate.calibrate_levels`` (ADR-0015 two-block confirmation),
R = replicates_fast/full.autocorrelation; B = 199 permutations.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.stats import autocorr as ac

B = 199
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
pytestmark = pytest.mark.calsuite("autocorrelation")


def _tus(rng: np.random.Generator, n: int, pg: float) -> tuple[np.ndarray, int]:
    sizes: list[int] = []
    while sum(sizes) < n:
        sizes.append(int(rng.geometric(pg)))
    sizes[-1] -= sum(sizes) - n
    if sizes[-1] == 0:
        sizes.pop()
    return np.repeat(np.arange(len(sizes)), sizes), len(sizes)


def _island_p(rng: np.random.Generator, n: int, frac: float = 0.2, width: int = 30) -> np.ndarray:
    p = rng.uniform(0.3, 0.9, size=n)
    starts = rng.choice(n // width, size=int(round(frac * n / width)), replace=False) * width
    for s in starts:
        p[s : s + width] = rng.uniform(0.01, 0.05, size=width)
    return p


def _jc(tu_null: bool) -> Callable[[np.random.Generator], float]:
    def rep(rng: np.random.Generator) -> float:
        n = 1_000
        p = _island_p(rng, n)
        d = rng.uniform(size=n) < p
        if tu_null:
            blk, _ = _tus(rng, n, 0.55)
            orders = ac.tu_block_permutations(blk, B, rng)
        else:
            orders = ac.label_permutations(n, B, rng)
        null = ac.join_count_matrix(d[orders], [1], True, P=p[orders])[:, 0]
        return ac.lag_test(ac.join_count(d, 1, True, p=p), null).p_value

    return rep


def _moran_tu_large_heavy(lag: int) -> Callable[[np.random.Generator], float]:
    def rep(rng: np.random.Generator) -> float:
        n = 1_000
        blk, nt = _tus(rng, n, 0.15)
        x = rng.standard_t(3, size=nt)[blk] + 0.5 * rng.standard_t(3, size=n)
        shift = int(rng.integers(0, n))
        x, blk = np.roll(x, shift), np.roll(blk, shift)
        orders = ac.tu_block_permutations(ac.tu_blocks([str(b) for b in blk], True), B, rng)
        null = ac.morans_i_matrix(x[orders], [lag], True)[:, 0]
        return ac.lag_test(ac.morans_i(x, lag, True), null).p_value

    return rep


def _moran_linear_missing(lag: int, tu_null: bool) -> Callable[[np.random.Generator], float]:
    def rep(rng: np.random.Generator) -> float:
        n = 600
        blk, nt = _tus(rng, n, 0.55)
        miss = rng.uniform(size=n) < 0.45
        if tu_null:
            x = rng.normal(size=nt)[blk] + 0.5 * rng.normal(size=n)
            x[miss] = np.nan
            orders = ac.tu_block_permutations(ac.tu_blocks([str(b) for b in blk], False), B, rng)
        else:
            x = rng.normal(size=n)
            x[miss] = np.nan
            orders = ac.label_permutations(n, B, rng, movable=~miss)
        null = ac.morans_i_matrix(x[orders], [lag], False)[:, 0]
        return ac.lag_test(ac.morans_i(x, lag, False), null).p_value

    return rep


def _moran_small_plasmid(rng: np.random.Generator) -> float:
    n = 40
    blk, nt = _tus(rng, n, 0.3)
    x = rng.normal(size=nt)[blk] + 0.5 * rng.normal(size=n)
    orders = ac.tu_block_permutations(ac.tu_blocks([str(b) for b in blk], True), B, rng)
    null = ac.morans_i_matrix(x[orders], [1], True)[:, 0]
    return ac.lag_test(ac.morans_i(x, 1, True), null).p_value


ATTACKS: dict[str, Callable[[np.random.Generator], float]] = {
    "jc_low_p_islands_label": _jc(False),
    "jc_low_p_islands_tu": _jc(True),
    "moran_tu_large_heavy_lag1": _moran_tu_large_heavy(1),
    "moran_tu_large_heavy_lag20": _moran_tu_large_heavy(20),
    "moran_linear_missing_tu_lag1": _moran_linear_missing(1, True),
    "moran_linear_missing_tu_lag20": _moran_linear_missing(20, True),
    "moran_linear_missing_label_lag1": _moran_linear_missing(1, False),
    "moran_linear_missing_label_lag20": _moran_linear_missing(20, False),
    "moran_small_plasmid": _moran_small_plasmid,
}


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_attack_autocorr_engine_level(attack: str, full: bool) -> None:
    r = C.replicates("autocorrelation", full=full)
    fn = ATTACKS[attack]

    def checked(rng: np.random.Generator) -> float:
        p = fn(rng)
        assert 0 < p <= 1, p  # every attack null is non-degenerate
        return p

    gate = calibrate_levels(
        "autocorrelation",
        "level",
        C.load_tolerances().levels,
        checked,
        r,
        stream=f"autocorrelation.attack_engine.{attack}",
    )
    gate.assert_ok()
