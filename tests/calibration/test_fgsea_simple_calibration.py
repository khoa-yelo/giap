"""Level of the fgsea simple-stage p (issue M8-02; GIPPA_SPEC.md §7.7.4, §8.4 level criterion).

Claim checked: the simple-stage p, min((1 + nLeEs)/(1 + nLeZero), (1 + nGeEs)/(1 + nGeZero)), is a
valid (+1, sign-conditional) permutation p. Under the null the tested set is a uniform random set
of the same size as the random sets, so its ES is exchangeable with theirs and, given its sign,
P(p <= alpha) <= alpha. The level check follows §8.4: rejection rate <= alpha + k sqrt(alpha
(1 - alpha)/R) with R = ``replicates_fast.gsea`` from tests/calibration/tolerances.yaml.

Two nulls: iid N(0, 1) stats (balanced signs) and the median-centred skewed stats M8-05 feeds fgsea
(exp(N(0, 1)) minus its median), each with one random set of 15 among 500.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from gippa.stats.fgsea import fgsea_simple, prepare

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("gsea")]

SUITE = "gsea"
N_STATS = 500
SET_SIZE = 15
N_PERM = 999
SAMPLE_SIZE = 101


def _replicate(
    draw_stats: Callable[[np.random.Generator], np.ndarray],
) -> Callable[[np.random.Generator], float]:
    labels = [f"g{i}" for i in range(N_STATS)]

    def rep(rng: np.random.Generator) -> float:
        stats = draw_stats(rng)
        members = [labels[i] for i in rng.choice(N_STATS, size=SET_SIZE, replace=False)]
        prep = prepare(stats, labels, {"s": members}, SET_SIZE, SET_SIZE)
        res = fgsea_simple(prep, N_PERM, SAMPLE_SIZE, rng)
        p = float(res.p[0])
        assert np.isfinite(p) and 1 / (N_PERM + 1) <= p <= 1, p
        return p

    return rep


def _check(step: str, rep: Callable[[np.random.Generator], float]) -> None:
    r = C.replicates(SUITE, full=False)
    levels = C.load_tolerances().levels
    calibrate_levels(SUITE, "level", levels, rep, r, stream=f"{SUITE}.simple.{step}").assert_ok()


def test_simple_p_level_normal_stats() -> None:
    _check("normal", _replicate(lambda rng: rng.standard_normal(N_STATS)))


def test_simple_p_level_centred_skewed_stats() -> None:
    def draw(rng: np.random.Generator) -> np.ndarray:
        x = np.exp(rng.standard_normal(N_STATS))
        return x - np.median(x)

    _check("centred_lognormal", _replicate(draw))
