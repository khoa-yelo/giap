"""Fast level calibration of the M4-03 permutation engine (GIPPA_SPEC.md §8.4 level criterion).

For alpha in ``levels``, the null rejection rate is <= alpha + k * sqrt(alpha(1-alpha)/R), with
k = ``se_multiplier``, R = ``replicates_fast.regions`` (tests/calibration/tolerances.yaml).
Discrete permutation p-values use the level check only (no KS). Nulls are exchangeable:
two groups of 30 iid N(0, 1) with label permutations, and iid N(0, 1) on a circular replicon with
exhaustive circular shifts. Gated by ``_criteria.calibrate`` (ADR-0015 two-block confirmation,
sharding): replicate i of stream ``regions.stats.<test>`` draws from
``replicate_rng(stream, i)`` (seed ``SeedSequence([calibration_base, crc32(stream), i])``).
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from gippa.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

pytestmark = pytest.mark.calibration

SUITE = "regions"
N_GROUP = 30


def _check_level(step: str, replicate: Callable[[np.random.Generator], float]) -> None:
    def checked(rng: np.random.Generator) -> float:
        p = replicate(rng)
        assert np.isfinite(p) and 0 < p <= 1, p
        return p

    r = C.replicates(SUITE, full=False)
    levels = C.load_tolerances().levels
    calibrate_levels(SUITE, "level", levels, checked, r, stream=f"{SUITE}.stats.{step}").assert_ok()


def _diff(x: np.ndarray, labels: np.ndarray) -> float:
    return float(x[labels].mean() - x[~labels].mean())


def test_perm_p_level_label_permutation() -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = _diff(x, labels)
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(obs, null)

    _check_level("perm_p", rep)


@pytest.mark.parametrize("alternative", ["greater", "less"])
def test_perm_p_one_sided_level(alternative: str) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = _diff(x, labels)
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(obs, null, alternative)  # type: ignore[arg-type]

    _check_level(f"perm_p.{alternative}", rep)


def test_circular_shift_level_iid() -> None:
    n, width = 1000, 20

    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=n)
        start = int(rng.integers(0, n))
        idx = (start + np.arange(width)) % n
        valid = circular_shifts(n, idx, circular=True)
        shifts, exhaustive = choose_shifts(valid, 20_000, n, 1_000, rng)
        assert bool(exhaustive) and len(shifts) == n - 1
        obs = float(x[idx].sum())
        null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts])
        return perm_p(obs, null)

    _check_level("circular_shifts", rep)


def test_linear_shift_level_iid() -> None:
    n, width = 400, 10

    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=n)
        start = int(rng.integers(0, n - width + 1))
        idx = start + np.arange(width)
        shifts, _ = choose_shifts(circular_shifts(n, idx, circular=False), 20_000, n, 1_000, rng)
        obs = float(x[idx].sum())
        null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts])
        return perm_p(obs, null)

    _check_level("linear_shifts", rep)


def test_besag_clifford_level() -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = abs(_diff(x, labels))  # two-sided: caller passes centred absolute statistics

        def stat(g: np.random.Generator) -> float:
            return abs(_diff(x, g.permutation(labels)))

        p, b_used = besag_clifford(stat, obs, 10, 1_000, rng)
        assert 10 <= b_used <= 1_000
        return p

    _check_level("besag_clifford", rep)
