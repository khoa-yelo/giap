"""Every ORA and CPB tail is a probability in [0, 1] (statistician review of PR #179, B1).

``conditional_poisson_binomial`` normalises by the support sum: p_over = w[i:].sum() / w.sum().
numpy's pairwise summation groups the two sums differently, so when k sits far below the mode
(an under-detected category) p_over can be 1 + 2^-52. ``detection_enrichment`` passes it to
``bh``, which rejects p > 1 and aborts the whole run. In the review's search, 663 of 4,000
strongly under-detected categories (T = 1,300 units, 60-300 per category) hit it, and so did
1 of 500 4b c = 2.5 replicates with ~114-TU categories.

Proposed repo path: tests/unit/function/test_ora_pvalue_domain.py
"""

from __future__ import annotations

import numpy as np
import pytest
from helpers.ora_calsim import collection, simulate, whole_tu_sets
from helpers.ora_data import ora

from giap.stats.poisson_binomial import CpbTails, conditional_poisson_binomial, tilt_pair

ONE_PLUS_ULP = float(np.nextafter(1.0, 2.0))


def _under_detected(seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    rng = np.random.default_rng(seed)
    T, n = 1300, int(rng.integers(60, 300))
    q = rng.beta(1.5, 1.0, T)
    d = rng.random(T) < q
    d[:n] = rng.random(n) < 0.5 * q[:n]
    N = int(d.sum())
    qt, ct = tilt_pair(q, N, comp=1 - q)
    return qt, ct, d, np.arange(T) < n, int(d[:n].sum()), N


@pytest.mark.parametrize("seed", [1, 2, 6])
def test_cpb_tails_never_exceed_one_far_below_the_mode(seed: int) -> None:
    qt, ct, _, m, k, N = _under_detected(seed)
    t = conditional_poisson_binomial(qt[m], qt[~m], k, N, comp_in=ct[m], comp_out=ct[~m])
    assert 0.0 <= t.p_over <= 1.0 and 0.0 <= t.p_under <= 1.0, (t.p_over, t.p_under)


def test_detection_enrichment_survives_a_tail_of_one_plus_ulp(monkeypatch) -> None:  # type: ignore[no-untyped-def]
    """The ORA layer never hands BH (or the table) a p above 1, whatever the primitive returns."""
    m = ora()
    real = m.conditional_poisson_binomial

    def plus_ulp(*a, **kw):  # type: ignore[no-untyped-def]
        t = real(*a, **kw)
        return CpbTails(t.k, t.N, ONE_PLUS_ULP, t.p_under, t.mean, t.lo, t.hi)

    monkeypatch.setattr(m, "conditional_poisson_binomial", plus_ulp)
    rng = np.random.default_rng(3)
    sim = simulate(rng)
    sets = whole_tu_sets(sim, rng, "random", 10)
    out = m.detection_enrichment(sim.genes, [collection(sets)], (5, 500), 0.05, folds=sim.folds)
    for c in ("p_value", "p_noelim", "q_value"):
        v = out[c].to_numpy(float)
        assert np.all((v >= 0.0) & (v <= 1.0)), c
