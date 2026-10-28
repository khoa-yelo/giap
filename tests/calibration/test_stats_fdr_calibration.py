"""Fast FDR calibration of the M4-01 engine (GIPPA_SPEC.md §8.4 FDR criterion).

Mean FDP over R replicates <= q + 2 SE(FDP) for q in ``levels``, R = ``replicates_fast.gene_fdr``
(tests/calibration/tolerances.yaml). Three generators, all with known truth: BH on independent
p-values (null U(0, 1), alternatives Beta), separate-search target-decoy with as many null decoys
as targets, and picked gene competition (Savitski 2015). Gated by ``_criteria.calibrate``
(ADR-0015 two-block confirmation, sharding): replicate i of stream ``gene_fdr.stats.<test>`` draws
from ``replicate_rng(stream, i)`` (seed ``SeedSequence([calibration_base, crc32(stream), i])``).
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

from gippa.stats.multitest import bh
from gippa.stats.tda import picked_fdr, tda_qvalues

pytestmark = pytest.mark.calibration

SUITE = "gene_fdr"
N = 2000
PI0 = 0.8


def _check_fdr(
    name: str, replicate: Callable[[np.random.Generator], tuple[np.ndarray, np.ndarray]]
) -> None:
    """``replicate`` returns (q per unit, is_false per unit); NaN q = not called."""
    levels = C.load_tolerances().levels

    def fdp(rng: np.random.Generator) -> np.ndarray:
        q, false = replicate(rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (false & called).sum() / n if n else 0.0
        return out

    r = C.replicates(SUITE, full=False)
    calibrate_levels(SUITE, "fdr", levels, fdp, r, stream=f"{SUITE}.stats.{name}").assert_ok()


def test_bh_controls_fdr_on_independent_p() -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        null = rng.random(N) < PI0
        p = np.where(null, rng.random(N), rng.beta(0.1, 8.0, N))
        return bh(p), null

    _check_fdr("bh", rep)


@pytest.mark.parametrize("use_true_pi", [False, True])
def test_tda_separate_search_controls_fdr(use_true_pi: bool) -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        null = rng.random(N) < PI0
        targets = np.where(null, rng.normal(0.0, 1.0, N), rng.normal(3.5, 1.0, N))
        decoys = rng.normal(0.0, 1.0, N)  # one decoy per target, from the null distribution
        pi = float(null.mean()) if use_true_pi else 1.0
        return tda_qvalues(targets, decoys, pi=pi), null

    _check_fdr(f"tda_separate.true_pi_{use_true_pi}", rep)


def test_picked_gene_fdr_is_controlled() -> None:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        absent = rng.random(N) < PI0
        target = np.where(absent, rng.normal(0.0, 1.0, N), rng.normal(3.5, 1.0, N))
        decoy = rng.normal(0.0, 1.0, N)
        # some genes have no target or no decoy evidence at all
        target[rng.random(N) < 0.05] = np.nan
        decoy[rng.random(N) < 0.05] = np.nan
        return picked_fdr(target, decoy).q, absent

    _check_fdr("picked", rep)
