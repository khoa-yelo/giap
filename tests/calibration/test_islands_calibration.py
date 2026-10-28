"""FWER calibration of the M5-06 island step (GIPPA_SPEC.md §7.3.5(d), §8.1, §8.4; ADR-0021 with
Amendment 1).

Validity domain (ADR-0021 Amendment 1, accepted 2026-10-20): the claim holds for the *dependence
regime* the calibration covers, not for a replicon size, and no replicon is flagged
``uncalibrated``; the run report carries l-hat, l' and the lag-1 Moran's I per replicon instead
(decision 48). The gated cells below are rho <= 0.5, where every size measured is calibrated; the
breaches at rho 0.9 are recorded in attacks/test_attack_islands.py (arm D).

Claim (§8.1 "Islands | maximum statistic | FWER | 0.05"): under the complete null (no island) the
share of runs that report at least one island is at most alpha. The family is genome-wide - M_b is
the maximum over every scanned replicon (lead decision 21) - so the claim covers the whole run and
not 1 - (1 - alpha)^k over k replicons; ``test_genome_wide_family_over_two_replicons`` checks a
chromosome and a plasmid together.

Each replicate draws a replicon under the null with ``islands_data.sim_contig`` (TU-level AR(1)
log2 abundance, gene noise, detection ``Bernoulli(0.55)`` independent of abundance, so the
left-censored Blom scores carry no island) and runs the real step ``tu_series`` ->
``scan_replicons`` -> ``islands_frame`` with B = 199 bootstrap replicates and 199 block-length
permutations. The replicate value is the indicator "at least one island" at every alpha in
``tolerances.yaml``; the criterion is ``fwer`` with R = ``replicates_fast/full.islands``.

Measured with the spec oracle before the implementation existed (R = 50, bound 0.112 at
alpha 0.05): rho 0 -> 0.00, rho 0.5 -> 0.02, chromosome + plasmid -> 0.02; the statistician's
gippasim re-measurement at R = 1,000 gives 0.020 / 0.031 / 0.037 for ``ar1_background`` at rho 0.5
on 557 / 1,001 / 2,501 scored TUs (bound 0.0638). The attack nulls (the M5-04 ori-ter scenario,
autocorrelation at rho 0.9, large TUs, draft assemblies, many replicons and the structural
low-direction ceiling) live in attacks/test_attack_islands.py.

Teeth: the spec's block length l-hat (instead of ADR-0021's l' = min(20, 2 l-hat + 1, T)) must be
a confirmed ``fail`` of the same criterion (ADR-0021: FWER 0.205-0.27 at rho >= 0.5; 0.62 here at
rho 0.9, R = 50).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import islands_data as idat

SUITE = "islands"
B, PERM = 199, 199
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def null_genome(rng: np.random.Generator, n_genes: int = 1000, rho: float = 0.5) -> list[Any]:
    return [idat.sim_contig(rng, "chr1", n_genes, rho=rho)]


@pytest.mark.parametrize("full", MODES)
@pytest.mark.parametrize("rho", [0.0, 0.5])
def test_island_fwer_at_level(rho: float, full: bool) -> None:
    """No island is planted, so every island reported is a false discovery."""
    levels = C.load_tolerances().levels

    def rep(rng: np.random.Generator) -> np.ndarray:
        return idat.any_island(
            rng, null_genome(rng, 1000, rho), levels, bootstrap=B, permutations=PERM
        )

    calibrate_levels(
        SUITE, "fwer", levels, rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.step.null.n1000.rho{rho:g}",
    ).assert_ok()  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_genome_wide_family_over_two_replicons(full: bool) -> None:
    """Decision 21: with a chromosome and a plasmid in one family the run-level rate stays at
    alpha; a per-replicon family would allow up to 1 - (1 - alpha)^2."""
    levels = C.load_tolerances().levels

    def rep(rng: np.random.Generator) -> np.ndarray:
        contigs = [
            idat.sim_contig(rng, "chr1", 1000, rho=0.5),
            idat.sim_contig(rng, "p1", 60, rho=0.5),
        ]
        return idat.any_island(rng, contigs, levels, bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", levels, rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.step.null.chr_plus_plasmid",
    ).assert_ok()  # fmt: skip


@pytest.mark.calibration
def test_teeth_spec_block_length_is_anticonservative(teeth_log: Path) -> None:
    """Teeth: the spec's l-hat as the bootstrap block length (without ADR-0021's widening) must be
    a confirmed ``fail`` of the same criterion at rho 0.9."""
    from gippa.stats.autocorr import choose_block_len

    def rep(rng: np.random.Generator) -> float:
        cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.05)
        genes, tus, replicons = idat.frames(null_genome(rng, 1000, 0.9))
        series = idat.build_series(genes, tus, replicons, cfg)
        l_hat = choose_block_len(series[0].detrended, series[0].circular, PERM, rng, 0.05, 20)
        spec = idat.config(bootstrap=B, permutations=PERM, alpha=0.05, block_len=int(l_hat))
        _, scan = idat.step(genes, series, spec, rng)
        assert scan.block_len == (int(l_hat),)
        return float(idat.min_p(scan) <= 0.05)

    rec = C.calibrate(
        SUITE, "fwer", 0.05, rep, C.replicates(SUITE, full=False),
        stream=f"{SUITE}.step.teeth.spec_block_len", log_dir=teeth_log,
    )  # fmt: skip
    print(f"spec block length l-hat: {rec}")
    assert rec["status"] == "fail", f"the spec's l-hat is not a confirmed fail: {rec}"
