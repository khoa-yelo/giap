"""Fast level calibration of the M4-03 permutation engine (GIAP_SPEC.md §8.4 level criterion).

For alpha in ``levels``, the null rejection rate is <= alpha + k * sqrt(alpha(1-alpha)/R), with
k = ``se_multiplier``, R = ``replicates_fast.regions``, seed ``seeds.calibration_base``
(tests/calibration/tolerances.yaml).
Discrete permutation p-values use the level check only (no KS). Nulls are exchangeable:
two groups of 30 iid N(0, 1) with label permutations, and iid N(0, 1) on a circular replicon with
exhaustive circular shifts. All randomness via ``giap.rng.for_step``.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

from giap.rng import for_step
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

pytestmark = pytest.mark.calibration

N_GROUP = 30


def _settings(tolerances: dict[str, object]) -> tuple[list[float], int, int, float]:
    levels = [float(x) for x in tolerances["levels"]]  # type: ignore[union-attr]
    reps = tolerances["replicates_fast"]
    seeds = tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    return (
        levels,
        int(reps["regions"]),
        int(seeds["calibration_base"]),
        float(tolerances["se_multiplier"]),  # type: ignore[arg-type]
    )


def _check_level(
    tolerances: dict[str, object], step: str, replicate: Callable[[np.random.Generator], float]
) -> None:
    levels, r, seed, k = _settings(tolerances)
    rng = for_step(step, seed)
    p = np.array([replicate(rng) for _ in range(r)])
    assert np.isfinite(p).all()
    assert ((p > 0) & (p <= 1)).all()
    for a in levels:
        rate = float((p <= a).mean())
        bound = a + k * np.sqrt(a * (1 - a) / r)
        assert rate <= bound, f"alpha={a}: rejection rate {rate:.4f} > bound {bound:.4f}, R={r}"


def _diff(x: np.ndarray, labels: np.ndarray) -> float:
    return float(x[labels].mean() - x[~labels].mean())


def test_perm_p_level_label_permutation(tolerances: dict[str, object]) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = _diff(x, labels)
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(obs, null)

    _check_level(tolerances, "calibration.perm_p", rep)


@pytest.mark.parametrize("alternative", ["greater", "less"])
def test_perm_p_one_sided_level(tolerances: dict[str, object], alternative: str) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = _diff(x, labels)
        null = np.array([_diff(x, rng.permutation(labels)) for _ in range(199)])
        return perm_p(obs, null, alternative)  # type: ignore[arg-type]

    _check_level(tolerances, f"calibration.perm_p.{alternative}", rep)


def test_circular_shift_level_iid(tolerances: dict[str, object]) -> None:
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

    _check_level(tolerances, "calibration.circular_shifts", rep)


def test_linear_shift_level_iid(tolerances: dict[str, object]) -> None:
    n, width = 400, 10

    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=n)
        start = int(rng.integers(0, n - width + 1))
        idx = start + np.arange(width)
        shifts, _ = choose_shifts(circular_shifts(n, idx, circular=False), 20_000, n, 1_000, rng)
        obs = float(x[idx].sum())
        null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts])
        return perm_p(obs, null)

    _check_level(tolerances, "calibration.linear_shifts", rep)


def test_besag_clifford_level(tolerances: dict[str, object]) -> None:
    def rep(rng: np.random.Generator) -> float:
        x = rng.normal(size=2 * N_GROUP)
        labels = np.arange(2 * N_GROUP) < N_GROUP
        obs = abs(_diff(x, labels))  # two-sided: caller passes centred absolute statistics

        def stat(g: np.random.Generator) -> float:
            return abs(_diff(x, g.permutation(labels)))

        p, b_used = besag_clifford(stat, obs, 10, 1_000, rng)
        assert 10 <= b_used <= 1_000
        return p

    _check_level(tolerances, "calibration.besag_clifford", rep)
