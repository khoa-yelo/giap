"""Level of the TU-counted ORA on gene frames (issue M8-04b; docs/specs/M8-04 calibration arms 1
and 6, lead decision 12; ADR-0029, ADR-0035; GIPPA_SPEC.md §8.4).

Claim: with gene detection ~ Bernoulli(p_true), p̂ = p_true (calibrated) and categories fixed
independently of the detections, the ``p_value`` of ``detection_enrichment(genes, [GeneSets],
size, alpha, folds=)`` (raw q recalibrated on M4-13-style TU folds, decision 12) is a valid
p-value in each direction. Pooled level criterion per (category type, direction), m = 10
categories per replicate, R = ``replicates_fast.ora``, levels from tolerances.yaml.

Generator (``helpers.ora_calsim.simulate``, 3,000 genes, TU sizes min(Geometric(0.45), 8),
c = 0, ``oracle=True``): the spec's ``ora_null`` (M8-14 will own the scenario; this is the
runnable stand-in on a gene frame plus ``GeneSets``). Category types:

- ``whole_tu``: 10 disjoint categories of 15 random whole TUs (spec arm 4 style);
- ``genes``: 10 disjoint categories of 40 random genes (arm 1/6: TU membership by the
  any-eligible-gene rule, so q varies with the covered TUs' sizes).

p̂ here is the truth, not a fitted M4-13 p_detect (arm 9 style; M8-15 runs the fitted version):
the recalibration should then be near the identity and must not break level.

Teeth (spec arm 6): on the ``genes`` categories, the unweighted TU null (all q equal, the TU
hypergeometric) over-rejects ``over`` at alpha = 0.05 beyond 2 alpha (spec pre-check 0.40),
because categories of genes cover larger-than-average TUs.

The replicate simulation is shared by the four gates (cached by replicate index); each gate has
its own stream name for the log.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import cache

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from helpers.ora_calsim import collection, random_gene_sets, simulate, whole_tu_sets
from scipy.stats import hypergeom

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calibration_nightly,
    pytest.mark.calsuite("ora"),
]

SUITE = "ora"
SIM_STREAM = f"{SUITE}.tu_null.sim"
M = 10
SIZE = (5, 500)
KINDS = ("whole_tu", "genes")


@cache
def _replicate(i: int) -> dict[tuple[str, str], np.ndarray]:
    from gippa.function import ora

    rng = C.replicate_rng(SIM_STREAM, i)
    sim = simulate(rng, oracle=True)
    sets = {
        "whole_tu": whole_tu_sets(sim, rng, "random", M),
        "genes": random_gene_sets(sim, rng, M),
    }
    out: dict[tuple[str, str], np.ndarray] = {}
    for kind, s in sets.items():
        assert len(s) == M
        tab = ora.detection_enrichment(sim.genes, [collection(s)], SIZE, 0.05, folds=sim.folds)
        assert "tu_recal_skipped" not in tab.attrs["run_flags"]
        for direction in ("over", "under"):
            sub = tab[tab["direction"] == direction].sort_values("term")
            assert len(sub) == M, (kind, direction, len(sub))
            out[(kind, direction)] = sub["p_value"].to_numpy(float)
        # teeth: the TU hypergeometric (equal q) on the same TU counts
        sub = tab[tab["direction"] == "over"].sort_values("term")
        n_u, n_d = int(sub["n_universe"].iloc[0]), int(sub["n_detected"].iloc[0])
        out[(kind, "hyper_over")] = np.array(
            [
                hypergeom(n_u, int(n), n_d).sf(int(k) - 1)
                for n, k in zip(sub["n_in_set"], sub["k_detected_in_set"], strict=True)
            ]
        )
    return out


def _fn(stream: str, R: int, cell: tuple[str, str], levels: list[float]
        ) -> Callable[[np.random.Generator], np.ndarray]:  # fmt: skip
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def rep(rng: np.random.Generator) -> np.ndarray:
        i = index[int(rng.bit_generator.seed_seq.entropy)]
        p = _replicate(i)[cell]
        assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1 + 1e-12))
        return np.array([np.count_nonzero(p <= a) / M for a in levels])

    return rep


@pytest.mark.parametrize("direction", ["over", "under"])
@pytest.mark.parametrize("kind", KINDS)
def test_tu_ora_level_under_the_null(kind: str, direction: str) -> None:
    r = C.replicates(SUITE, full=False)
    levels = list(C.load_tolerances().levels)
    stream = f"{SUITE}.tu_null.{kind}.{direction}"
    gate = calibrate_levels(
        SUITE, "level", levels, _fn(stream, r, (kind, direction), levels), r, stream=stream, m=M
    )
    gate.assert_ok()


def test_teeth_tu_hypergeometric_over_rejects_gene_categories() -> None:
    r = C.replicates(SUITE, full=False) // 2
    alpha = 0.05
    p = np.concatenate([_replicate(i)[("genes", "hyper_over")] for i in range(r)])
    rate = float(np.mean(p <= alpha))
    assert rate > 2 * alpha, f"TU hypergeometric rejected {rate:.3f} <= 2 alpha"
