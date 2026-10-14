"""Fast FDR calibration of the M4-01 engine (GIAP_SPEC.md §8.4 FDR criterion).

Mean FDP over R replicates <= q + 2 SE(FDP) for q in ``levels``, R = ``replicates_fast.gene_fdr``,
seed ``seeds.calibration_base`` (tests/calibration/tolerances.yaml). Three generators, all with
known truth: BH on independent p-values (null U(0, 1), alternatives Beta), separate-search
target-decoy with as many null decoys as targets, and picked gene competition (Savitski 2015).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

from giap.stats.multitest import bh
from giap.stats.tda import picked_fdr, tda_qvalues

pytestmark = pytest.mark.calibration

N = 2000
PI0 = 0.8


def _settings(tolerances: dict[str, object]) -> tuple[list[float], int, int, float]:
    levels = [float(x) for x in tolerances["levels"]]  # type: ignore[union-attr]
    reps = tolerances["replicates_fast"]
    seeds = tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    return (
        levels,
        int(reps["gene_fdr"]),
        int(seeds["calibration_base"]),
        float(tolerances["se_multiplier"]),  # type: ignore[arg-type]
    )


def _check_fdr(
    tolerances: dict[str, object],
    replicate: Callable[[np.random.Generator], tuple[np.ndarray, np.ndarray]],
) -> None:
    """``replicate`` returns (q per unit, is_false per unit); NaN q = not called."""
    levels, r, seed, k = _settings(tolerances)
    rng = np.random.default_rng(seed)
    fdp = np.zeros((r, len(levels)))
    for i in range(r):
        q, false = replicate(rng)
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            fdp[i, j] = (false & called).sum() / n if n else 0.0
    mean = fdp.mean(axis=0)
    se = fdp.std(axis=0, ddof=1) / np.sqrt(r)
    for j, a in enumerate(levels):
        assert mean[j] <= a + k * se[j], (
            f"q={a}: mean FDP {mean[j]:.4f} > {a} + {k}*{se[j]:.4f} over R={r}"
        )


def test_bh_controls_fdr_on_independent_p(tolerances: dict[str, object]) -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        null = rng.random(N) < PI0
        p = np.where(null, rng.random(N), rng.beta(0.1, 8.0, N))
        return bh(p), null

    _check_fdr(tolerances, rep)


@pytest.mark.parametrize("use_true_pi", [False, True])
def test_tda_separate_search_controls_fdr(tolerances: dict[str, object], use_true_pi: bool) -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        null = rng.random(N) < PI0
        targets = np.where(null, rng.normal(0.0, 1.0, N), rng.normal(3.5, 1.0, N))
        decoys = rng.normal(0.0, 1.0, N)  # one decoy per target, from the null distribution
        pi = float(null.mean()) if use_true_pi else 1.0
        return tda_qvalues(targets, decoys, pi=pi), null

    _check_fdr(tolerances, rep)


def test_picked_gene_fdr_is_controlled(tolerances: dict[str, object]) -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        absent = rng.random(N) < PI0
        target = np.where(absent, rng.normal(0.0, 1.0, N), rng.normal(3.5, 1.0, N))
        decoy = rng.normal(0.0, 1.0, N)
        # some genes have no target or no decoy evidence at all
        target[rng.random(N) < 0.05] = np.nan
        decoy[rng.random(N) < 0.05] = np.nan
        return picked_fdr(target, decoy).q, absent

    _check_fdr(tolerances, rep)
