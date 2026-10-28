"""Level of the conditional Poisson-binomial tails on ``tilt_pair`` with complements (issue M8-04b;
docs/specs/M8-04 lead decision 12, ADR-0035 decision 3; GIPPA_SPEC.md §8.4 level criterion).

Replaces ``test_cpb_calibration.py`` (M8-04a), which tilts with ``tilt`` and drops the
complements. Same claim, generator and teeth, plus the decision-12 tilt-stress arm:

Claim: under the model the primitive assumes (independent units, unit t detected with
probability q_t, categories fixed independently of the detections) the tails ``p_over`` and
``p_under`` of ``conditional_poisson_binomial`` on ``tilt_pair(q, N, comp=1 - q)`` (complements
passed through ``comp_in``/``comp_out``) are exact conditional p-values given S = N, so
P(p <= alpha) <= alpha in each direction. Pooled level criterion (``check_level_pooled``):
R = ``replicates_fast.ora``, m = 20 disjoint categories per replicate (negatively dependent given
N, which makes the pooled bound conservative), levels from tolerances.yaml.

- **Base arm.** T = 300 units, q ~ Beta(1.2, 2.5); the 20 categories are the blocks of 15
  consecutive q-ranks (some all-high-q, some all-low-q: spec arm 6's size bias).
- **Tilt-stress arm (decision 12).** The tested q is badly scaled and partly saturated: every
  log-odds is shifted by -8 (odds x e^-8) after 5 % of the units (15) are put at log-odds
  U(45, 52), so their scaled q is 1.0 in float and only the complement carries them. Detection
  is drawn from the unscaled q; the conditional law is invariant to a common odds factor, so the
  scaled q is still the exact null. The pairs (expit(x), expit(-x)) go to
  ``tilt_pair(q, N, comp=)``. Categories: 20 random disjoint blocks of 15, so near-1 units are
  spread over them.
- **Teeth.** The base data under the unweighted null (all q equal, the TU hypergeometric)
  over-rejects ``p_over`` at alpha = 0.05 beyond 2 alpha.
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from scipy.special import expit

from gippa.stats.poisson_binomial import conditional_poisson_binomial, tilt_pair

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("ora")]

SUITE = "ora"
T = 300
M = 20
SIZE = T // M
N_NEAR_ONE = 15


def _test(q: np.ndarray, comp: np.ndarray, d: np.ndarray, blocks: list[np.ndarray]
          ) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    N = int(d.sum())
    qt, ct = tilt_pair(q, N, comp=comp)
    over, under = np.empty(len(blocks)), np.empty(len(blocks))
    for i, b in enumerate(blocks):
        m = np.zeros(q.size, dtype=bool)
        m[b] = True
        t = conditional_poisson_binomial(
            qt[m], qt[~m], int(d[m].sum()), N, comp_in=ct[m], comp_out=ct[~m]
        )
        over[i], under[i] = t.p_over, t.p_under
    return over, under


def _base(rng: np.random.Generator, weighted: bool) -> tuple[np.ndarray, np.ndarray]:
    q = rng.beta(1.2, 2.5, T)
    d = rng.random(T) < q
    order = np.argsort(q, kind="stable")
    blocks = [order[i * SIZE : (i + 1) * SIZE] for i in range(M)]
    q_used = q if weighted else np.full(T, 0.5)
    return _test(q_used, 1.0 - q_used, d, blocks)


def _stress(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    b = rng.beta(1.2, 2.5, T)
    x0 = np.log(b) - np.log1p(-b)
    near = rng.choice(T, N_NEAR_ONE, replace=False)
    x0[near] = rng.uniform(45.0, 52.0, N_NEAR_ONE)
    d = rng.random(T) < expit(x0)
    x = x0 - 8.0
    q, comp = expit(x), expit(-x)
    assert np.all(q[near] == 1.0) and np.all(comp[near] > 0.0)  # the saturated case is exercised
    perm = rng.permutation(T)
    blocks = [perm[i * SIZE : (i + 1) * SIZE] for i in range(M)]
    return _test(q, comp, d, blocks)


def _shares(
    gen: Callable[[np.random.Generator], tuple[np.ndarray, np.ndarray]],
    direction: int,
    levels: list[float],
) -> Callable[[np.random.Generator], np.ndarray]:
    def rep(rng: np.random.Generator) -> np.ndarray:
        p = gen(rng)[direction]
        assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1 + 1e-12))
        return np.array([np.count_nonzero(p <= a) / M for a in levels])

    return rep


@pytest.mark.parametrize(("direction", "name"), [(0, "over"), (1, "under")])
def test_cpb_pair_tails_level_under_the_conditional_null(direction: int, name: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    gate = calibrate_levels(
        SUITE,
        "level",
        levels,
        _shares(lambda g: _base(g, weighted=True), direction, levels),
        r,
        stream=f"{SUITE}.cpb_pair.{name}",
        m=M,
    )
    gate.assert_ok()


@pytest.mark.parametrize(("direction", "name"), [(0, "over"), (1, "under")])
def test_cpb_pair_tails_level_under_tilt_stress(direction: int, name: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    gate = calibrate_levels(
        SUITE,
        "level",
        levels,
        _shares(_stress, direction, levels),
        r,
        stream=f"{SUITE}.cpb_pair.stress.{name}",
        m=M,
    )
    gate.assert_ok()


def test_teeth_unweighted_null_over_rejects_high_q_categories() -> None:
    r = C.replicates(SUITE, full=False)
    alpha = 0.05
    rejected = 0
    for i in range(r):
        over, _ = _base(C.replicate_rng(f"{SUITE}.cpb_pair.teeth", i), weighted=False)
        rejected += int(np.count_nonzero(over <= alpha))
    rate = rejected / (r * M)
    assert rate > 2 * alpha, f"unweighted null rejected {rate:.3f} <= 2 alpha"
