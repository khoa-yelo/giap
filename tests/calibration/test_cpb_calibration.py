"""Level of the conditional Poisson-binomial tails (issue M8-04a; docs/specs/M8-04 §1, ADR-0029;
GIAP_SPEC.md §8.4 level criterion).

Claim checked: under the model the primitive assumes (independent units, unit t detected with
probability q_t, a category fixed independently of the detections) the tails ``p_over`` and
``p_under`` of ``conditional_poisson_binomial(q[C], q[not C], k, N)`` are exact conditional
p-values given S = N, so P(p <= alpha) <= alpha in each direction. Checked with the pooled level
criterion (``check_level_pooled``): R = ``replicates_fast.ora`` replicates, m = 20 disjoint
categories per replicate, levels from tolerances.yaml.

Generator: T = 300 units with heterogeneous q ~ Beta(1.2, 2.5) (as TU sizes make q unequal); the
categories are the 20 blocks of 15 consecutive q-ranks, so some categories hold only high-q units
and some only low-q units (the size bias of spec arm 6). The tilt (``tilt(q, N)``, computed once
per replicate as in ORA) is applied before the call; the law is invariant to it.

Teeth: the same data tested with the unweighted null (all q equal, the TU hypergeometric) over-
rejects: the pooled ``p_over`` rate at alpha = 0.05 exceeds 2 alpha, since high-q categories are
detected more often than the hypergeometric allows.

M8-15 runs the gated ORA arms on real gene frames; this file calibrates the primitive alone.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.stats.poisson_binomial import conditional_poisson_binomial, tilt

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("ora")]

SUITE = "ora"
T = 300
M = 20
SIZE = T // M


def _tails(rng: np.random.Generator, weighted: bool) -> tuple[np.ndarray, np.ndarray]:
    q = rng.beta(1.2, 2.5, T)
    d = rng.random(T) < q
    N = int(d.sum())
    order = np.argsort(q, kind="stable")
    q_used = tilt(q, N) if weighted else np.full(T, 0.5)
    over, under = np.empty(M), np.empty(M)
    for i in range(M):
        m = np.zeros(T, dtype=bool)
        m[order[i * SIZE : (i + 1) * SIZE]] = True
        t = conditional_poisson_binomial(q_used[m], q_used[~m], int(d[m].sum()), N)
        over[i], under[i] = t.p_over, t.p_under
    return over, under


def _shares(direction: int, levels: list[float]) -> Callable[[np.random.Generator], np.ndarray]:
    def rep(rng: np.random.Generator) -> np.ndarray:
        p = _tails(rng, weighted=True)[direction]
        assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1 + 1e-12))
        return np.array([np.count_nonzero(p <= a) / M for a in levels])

    return rep


@pytest.mark.parametrize(("direction", "name"), [(0, "over"), (1, "under")])
def test_cpb_tails_level_under_the_conditional_null(direction: int, name: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    gate = calibrate_levels(
        SUITE, "level", levels, _shares(direction, levels), r, stream=f"{SUITE}.cpb.{name}", m=M
    )
    gate.assert_ok()


def test_teeth_unweighted_null_over_rejects_high_q_categories() -> None:
    r = C.replicates(SUITE, full=False)
    alpha = 0.05
    rejected = 0
    for i in range(r):
        over, _ = _tails(C.replicate_rng(f"{SUITE}.cpb.teeth", i), weighted=False)
        rejected += int(np.count_nonzero(over <= alpha))
    rate = rejected / (r * M)
    assert rate > 2 * alpha, f"unweighted null rejected {rate:.3f} <= 2 alpha"
