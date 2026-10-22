"""Fast level calibration of Wallenius ORA under its own model (M8-01; GIAP_SPEC.md §7.7.3, §8.4).

Generator with known truth: a universe of M = 400 genes, a category C of n = 40, p_detect = 0.6
inside C and 0.3 outside (so the true odds are ω = 2 and the mean-odds approximation is exact).
N = 120 genes are "detected" by sequential weighted sampling without replacement (each remaining
gene drawn with probability ∝ p_detect), drawn as the N smallest of Exp(1) / p_detect keys
(exponential races: the same law as sequential draws), so k ~ Wallenius(M, n, N, ω) exactly.
Claim: ``wallenius_ora``'s ``p_over`` and ``p_under`` are each valid p-values, i.e. the null
rejection rate at alpha in ``levels`` is <= alpha + k*sqrt(alpha(1 - alpha)/R) (§8.4 level
criterion; discrete, so typically below alpha). R = ``replicates_fast.ora``; gated by
``_gate.calibrate_levels`` (ADR-0015).
The suite-level ``ora`` gate on the synthetic ontology and attack 2 are M8-15's.
"""

from __future__ import annotations

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.stats.wallenius import wallenius_ora

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("ora")]

SUITE = "ora"
M, N_CAT, N_DET = 400, 40, 120
P_IN, P_OUT = 0.6, 0.3


def _draw(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    in_cat = np.arange(M) < N_CAT
    p = np.where(in_cat, P_IN, P_OUT)
    keys = rng.exponential(size=M) / p
    detected = np.zeros(M, dtype=bool)
    detected[np.argsort(keys, kind="stable")[:N_DET]] = True
    return detected, p, in_cat


@pytest.mark.parametrize("tail", ["p_over", "p_under"])
def test_wallenius_ora_level_under_urn_model(tail: str) -> None:
    def rep(rng: np.random.Generator) -> float:
        detected, p, in_cat = _draw(rng)
        r = wallenius_ora(detected, p, in_cat)
        assert r.odds == pytest.approx(P_IN / P_OUT, rel=1e-12)
        value = float(getattr(r, tail))
        assert 0.0 < value <= 1.0, value
        return value

    r = C.replicates(SUITE, full=False)
    levels = C.load_tolerances().levels
    calibrate_levels(
        SUITE, "level", levels, rep, r, stream=f"{SUITE}.stats.wallenius.{tail}"
    ).assert_ok()
