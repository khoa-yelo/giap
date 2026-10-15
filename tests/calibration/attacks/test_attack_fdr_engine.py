"""Attack nulls for the M4-01 FDR engine (GIAP_SPEC.md §8.5; statistician-owned).

Each generator is built to push ``bh``, ``tda_qvalues`` or ``picked_fdr`` above level while staying
inside the method's stated assumptions (docs/methods/target-decoy.md, multiple-testing.md):

- picked, ``sizes_corr_ties``: null genes of very different sizes (best of k_g null peptides,
  k_g log-normal), a gene effect shared by a gene's target and decoy, scores rounded to 0.1 (many
  target-decoy ties) and some genes with only target or only decoy evidence;
- picked, ``coarse_ties``: scores rounded to integers, so a large share of genes tie;
- separate search, ``discrete_batch``: exponential scores on a 0.5 grid plus a replicate-wide batch
  shift shared by targets and decoys;
- global nulls, picked and separate search, large and small (B.6-sized) strata: every call is false;
- BH, ``equicorrelated_ties``: one-sided z with pairwise correlation 0.5 (PRDS) and p rounded to
  three decimals.

Criterion (§8.4): mean FDP <= q + 2 SE(FDP) over R replicates, q in ``levels``. R is
``replicates_fast.gene_fdr`` (marker ``calibration``) or ``replicates_full.gene_fdr``
(marker ``calibration_full``). Gated by ``_criteria.calibrate`` (ADR-0015 two-block confirmation,
sharding); replicate i of stream ``gene_fdr.attack_fdr_engine.<attack>`` draws from
``replicate_rng(stream, i)`` (seeds from ``seeds.calibration_base``).
"""

from __future__ import annotations

from collections.abc import Callable

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels
from scipy.stats import norm

from giap.stats.multitest import bh
from giap.stats.tda import picked_fdr, tda_qvalues

Replicate = Callable[[np.random.Generator], tuple[np.ndarray, np.ndarray]]
N = 2000


def _best_of(rng: np.random.Generator, k: np.ndarray) -> np.ndarray:
    """Maximum of k_i standard normals, drawn exactly via U^(1/k)."""
    return np.asarray(norm.ppf(rng.random(k.size) ** (1.0 / k)))


def _sizes_corr_ties(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    k = np.maximum(1, np.round(rng.lognormal(1.5, 1.0, N))).astype(int)
    absent = rng.random(N) < 0.8
    gene = 0.7 * rng.normal(0.0, 1.0, N)
    target = _best_of(rng, k) + gene + np.where(absent, 0.0, 3.0)
    decoy = _best_of(rng, k) + gene
    target, decoy = np.round(target, 1), np.round(decoy, 1)
    target[rng.random(N) < 0.05] = np.nan
    decoy[rng.random(N) < 0.05] = -np.inf
    return picked_fdr(target, decoy).q, absent


def _coarse_ties(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    absent = rng.random(N) < 0.8
    target = np.round(np.where(absent, rng.normal(0.0, 1.0, N), rng.normal(3.5, 1.0, N)))
    decoy = np.round(rng.normal(0.0, 1.0, N))
    return picked_fdr(target, decoy).q, absent


def _discrete_batch(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n = 500
    null = rng.random(n) < 0.8
    batch = rng.normal(0.0, 1.0)
    target = rng.exponential(1.0, n) + np.where(null, 0.0, 4.0) + batch
    decoy = rng.exponential(1.0, n) + batch
    return tda_qvalues(np.floor(target * 2) / 2, np.floor(decoy * 2) / 2), null


def _global_null(n: int, picked: bool) -> Replicate:
    def rep(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
        target, decoy = rng.normal(0.0, 1.0, n), rng.normal(0.0, 1.0, n)
        q = picked_fdr(target, decoy).q if picked else tda_qvalues(target, decoy)
        return q, np.ones(n, dtype=bool)

    return rep


def _bh_equicorrelated_ties(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n = 1000
    null = rng.random(n) < 0.8
    z = np.sqrt(0.5) * rng.normal() + np.sqrt(0.5) * rng.normal(0.0, 1.0, n)
    p = np.round(norm.sf(z + np.where(null, 0.0, 3.0)), 3)
    return bh(p), null


ATTACKS: dict[str, Replicate] = {
    "picked_sizes_corr_ties": _sizes_corr_ties,
    "picked_coarse_ties": _coarse_ties,
    "tda_discrete_batch": _discrete_batch,
    "picked_global_null_2000": _global_null(2000, picked=True),
    "picked_global_null_25": _global_null(25, picked=True),
    "tda_global_null_2000": _global_null(2000, picked=False),
    "tda_global_null_25": _global_null(25, picked=False),
    "bh_equicorrelated_ties": _bh_equicorrelated_ties,
}

SUITES = [
    pytest.param("replicates_fast", marks=pytest.mark.calibration, id="fast"),
    pytest.param("replicates_full", marks=pytest.mark.calibration_full, id="full"),
]


@pytest.mark.parametrize("suite", SUITES)
@pytest.mark.parametrize("attack", sorted(ATTACKS))
def test_attack_stays_at_level(attack: str, suite: str) -> None:
    levels = C.load_tolerances().levels
    r = C.replicates("gene_fdr", full=suite == "replicates_full")

    def fdp(rng: np.random.Generator) -> np.ndarray:
        q, false = ATTACKS[attack](rng)
        out = np.zeros(len(levels))
        for j, a in enumerate(levels):
            called = ~np.isnan(q) & (q <= a)
            n = int(called.sum())
            out[j] = (false & called).sum() / n if n else 0.0
        return out

    stream = f"gene_fdr.attack_fdr_engine.{attack}"
    calibrate_levels("gene_fdr", "fdr", levels, fdp, r, stream=stream).assert_ok()
