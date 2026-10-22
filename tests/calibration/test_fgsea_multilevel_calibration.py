"""Accuracy and level of the fgsea multilevel p (issue M8-11; GIAP_SPEC.md §7.7.4, §8.2, §8.4;
docs/specs/M8-11-stats-fgsea-multilevel.md acceptance criterion 2).

Claim 1 (accuracy, spec acceptance criterion): for a set with a planted shift whose exact p is
known by brute-force permutation at 10^6 draws (p ~ 1e-4), |log2 p_hat - log2 p| <= 2 log2err in
>= 90 % of R = ``replicates_fast.gsea`` (50) seeds. The brute-force p is fgsea's conditional p,
#{ES_rand >= ES} / #{ES_rand >= 0}, which is what ``cppMPval / denomProb`` estimates. Its own
Monte-Carlo error (about 70 exceedances, sd(log2) ~ 0.17) is small against 2 log2err ~ 1.1.
A teeth check: the same estimates against a reference 8x off (3 log2 units) must fail the
criterion, so the test can tell a right p from a wrong one.

Claim 2 (level, §8.2/§8.4): under the null (a uniform random set) the upgraded p is a valid p,
P(p <= alpha) <= alpha within the §8.4 bound, at the R of ``replicates_fast.gsea``. The tail
(alpha = 0.01) inspection at full R is M8-15's job; this is the fast subset.

Interface: ``fgsea_multilevel(prep, simple, sample_size, eps, rng) -> SimpleResult``.
"""

from __future__ import annotations

import math
from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from giap.rng import for_step
from giap.stats._fgsea_kernels import random_es_batch
from giap.stats.fgsea import Prepared, fgsea_simple, gsea_es, prepare

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("gsea")]

SUITE = "gsea"
STREAM = "function.gsea"
N_STATS = 500
SET_SIZE = 15
SHIFT = 1.2  # tuned so that the exact p is ~1.4e-4
N_PERM = 1000
SAMPLE_SIZE = 101
N_REFERENCE = 10**6
COVER = 0.90


def _planted() -> Prepared:
    rng = np.random.default_rng(20261022)  # test data, not a GIAP stream
    x = rng.standard_normal(N_STATS)
    x[:SET_SIZE] += SHIFT
    lab = [f"g{i}" for i in range(N_STATS)]
    return prepare(x, lab, {"planted": lab[:SET_SIZE]}, SET_SIZE, SET_SIZE)


def _brute_force_p(prep: Prepared) -> float:
    es, _ = gsea_es(prep.stats, prep.members[0])
    assert es > 0
    rng = for_step("calibration.gsea.multilevel.reference", 1)
    ge = ge0 = 0
    chunk = 20_000
    for _ in range(N_REFERENCE // chunk):
        # the SET_SIZE smallest of iid uniforms: a uniform random subset of ranks
        idx = np.argpartition(rng.random((chunk, N_STATS)), SET_SIZE, axis=1)[:, :SET_SIZE]
        rand = random_es_batch(prep.stats, idx.astype(np.int64) + 1)[:, -1]
        ge += int((rand >= es).sum())
        ge0 += int((rand >= 0).sum())
    assert ge >= 30, ge  # enough exceedances for a reference
    return ge / ge0


@pytest.fixture(scope="module")
def estimates() -> tuple[float, list[tuple[float, float]]]:
    from giap.stats.fgsea import fgsea_multilevel

    prep = _planted()
    p_exact = _brute_force_p(prep)
    assert 3e-5 < p_exact < 5e-4, p_exact
    out = []
    for rep in range(C.replicates(SUITE, full=False)):
        rng = for_step(STREAM, C.suite_seed(f"{SUITE}.multilevel.accuracy", rep))
        simple = fgsea_simple(prep, N_PERM, SAMPLE_SIZE, rng)
        assert bool(simple.needs_multilevel[0])
        res = fgsea_multilevel(prep, simple, SAMPLE_SIZE, 0.0, rng)
        p, err = float(res.p[0]), float(res.log2err[0])
        assert 0 < p <= 1
        out.append((p, err))
    return p_exact, out


def _coverage(p_ref: float, est: list[tuple[float, float]]) -> float:
    """Share of seeds with |log2 p_hat - log2 p_ref| <= 2 log2err (NaN log2err counts as a miss)."""
    hits = [
        math.isfinite(err) and abs(math.log2(p) - math.log2(p_ref)) <= 2 * err for p, err in est
    ]
    return sum(hits) / len(hits)


def test_multilevel_p_within_two_log2err_of_the_brute_force_p(
    estimates: tuple[float, list[tuple[float, float]]],
) -> None:
    p_exact, est = estimates
    cover = _coverage(p_exact, est)
    assert cover >= COVER, (p_exact, cover, est[:5])


def test_teeth_a_reference_8x_off_fails_the_criterion(
    estimates: tuple[float, list[tuple[float, float]]],
) -> None:
    p_exact, est = estimates
    assert _coverage(p_exact * 8, est) < COVER
    assert _coverage(p_exact / 8, est) < COVER


def _null_replicate() -> Callable[[np.random.Generator], float]:
    from giap.stats.fgsea import fgsea_multilevel

    lab = [f"g{i}" for i in range(N_STATS)]

    def rep(rng: np.random.Generator) -> float:
        stats = rng.standard_normal(N_STATS)
        members = [lab[i] for i in rng.choice(N_STATS, size=SET_SIZE, replace=False)]
        prep = prepare(stats, lab, {"s": members}, SET_SIZE, SET_SIZE)
        simple = fgsea_simple(prep, N_PERM, SAMPLE_SIZE, rng)
        p = float(fgsea_multilevel(prep, simple, SAMPLE_SIZE, 0.0, rng).p[0])
        assert np.isfinite(p) and 0 < p <= 1, p
        return p

    return rep


def test_multilevel_p_level_under_the_null() -> None:
    r = C.replicates(SUITE, full=False)
    levels = C.load_tolerances().levels
    calibrate_levels(
        SUITE, "level", levels, _null_replicate(), r, stream=f"{SUITE}.multilevel.null"
    ).assert_ok()
