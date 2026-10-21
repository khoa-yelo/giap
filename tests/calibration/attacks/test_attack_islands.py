"""Attack nulls for the M5-06 island block bootstrap (GIAP_SPEC.md §8.5; statistician-owned).

Each generator stays inside the null the step claims (no island planted) and tries to push the
genome-wide FWER above level. They are the statistician's list for this issue: autocorrelated
abundance, large TUs, the ori-ter gradient, draft assemblies and multiple replicons. Criterion:
``fwer`` through ``_gate.calibrate_levels`` (ADR-0015 two-block confirmation), R =
``replicates_fast/full.islands``, B = 199 bootstrap replicates, 199 block-length permutations.

Rates measured with the spec oracle before the implementation existed (R = 50, bounds
[0.038, 0.112, 0.185] at levels [0.01, 0.05, 0.10]):

* ``large TUs`` (sizes 1 + Geometric(0.15), mean ~6.7 genes, 2,000 genes, rho 0.5):
  [0.00, 0.00, 0.02] - gated. TU means of many genes have a smaller variance, so the windows of
  single-gene TUs dominate M_b (docs/methods/islands.md "Why S/sqrt(n_w) is not standardised").
* ``draft assembly`` (500-gene circular chromosome + ten 30-gene linear contigs):
  [0.00, 0.00, 0.04] - gated; the small contigs are never scanned (decision 40, ``few_values``).
* ``many replicons`` (chromosome + three plasmids): [0.00, 0.02, 0.04] - gated; one genome-wide
  family over four replicons.
* ``autocorrelation rho 0.9`` on 1,000 single-gene TUs: [0.02, 0.16, 0.30] - **recorded, not
  gated**: this is ADR-0021 decision 2's uncalibrated corner (the ADR measured 0.13 at
  1,000 giapsim genes and 0.023 on a 4,500-gene chromosome). The islands carry ``uncalibrated``.
* ``ori-ter gradient`` (cosine, 2 log2 peak to trough, rho 0.5): [0.00, 0.08, 0.12] - **recorded,
  not gated**: the gradient raises the rate from 0.02 to 0.08 at alpha 0.05 even though LOESS
  detrending removes the trend itself. ADR-0021 measured 0.04 for the giapsim gradient; the
  gated version of this attack waits for the M5-04 oriC/ori-ter scenario (PR #119) and a
  statistician decision.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import islands_data as idat

B, PERM = 199, 199

SUITE = "islands"
MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]


def _levels() -> tuple[float, ...]:
    return C.load_tolerances().levels


def _record(name: str, make: Any, r: int) -> np.ndarray:
    """Run ``r`` replicates of a non-gated attack and return the rate per level."""
    levels = _levels()
    hits = np.zeros(len(levels))
    for i in range(r):
        rng = C.replicate_rng(f"{SUITE}.attack.{name}", i)
        hits += idat.any_island(rng, make(rng), levels, bootstrap=B, permutations=PERM)
    rate = hits / r
    print(f"recorded attack {name}: R = {r}, FWER {np.round(rate, 3).tolist()} at {list(levels)}"
          f" (bounds {[round(C.level_bound(a, r), 3) for a in levels]})")  # fmt: skip
    return rate


# --------------------------------------------------------------------------- gated attacks


@pytest.mark.parametrize("full", MODES)
def test_attack_large_tus(full: bool) -> None:
    """TU sizes 1 + Geometric(0.15) (mean ~6.7 genes): heteroscedastic TU scores."""

    def rep(rng: np.random.Generator) -> np.ndarray:
        return idat.any_island(
            rng,
            [idat.sim_contig(rng, "chr1", 2000, rho=0.5, tu_p=0.15)],
            _levels(),
            bootstrap=B,
            permutations=PERM,
        )

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.large_tus",
    ).assert_ok()  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_attack_draft_assembly(full: bool) -> None:
    """A draft assembly: one 500-gene circular chromosome and ten 30-gene linear contigs."""

    def rep(rng: np.random.Generator) -> np.ndarray:
        return idat.any_island(rng, _draft(rng), _levels(), bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.draft_assembly",
    ).assert_ok()  # fmt: skip


def _draft(rng: np.random.Generator) -> list[Any]:
    return [idat.sim_contig(rng, "chr1", 500, rho=0.5)] + [
        idat.sim_contig(rng, f"ctg{i}", 30, rho=0.5, circular=False) for i in range(10)
    ]


@pytest.mark.calibration
def test_draft_contigs_are_never_scanned() -> None:
    """Decision 40/41: the 30-gene linear contigs give a span of k < 4, so they are flagged
    ``few_values``, contribute nothing to M_b and never carry an island."""
    rng = C.replicate_rng(f"{SUITE}.attack.draft_assembly.structure", 0)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.10)
    genes, tus, replicons = idat.frames(_draft(rng))
    out, scan = idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg, rng)
    assert scan.scanned[0] and not any(scan.scanned[1:])
    assert all(f == (idat.FEW_VALUES,) for f in scan.flags[1:])
    assert all(p.size == 0 for p in scan.p_value[1:])
    assert set(out["contig"]) <= {"chr1"}


@pytest.mark.parametrize("full", MODES)
def test_attack_many_replicons(full: bool) -> None:
    """One chromosome and three plasmids in a single genome-wide family."""

    def rep(rng: np.random.Generator) -> np.ndarray:
        contigs = [idat.sim_contig(rng, "chr1", 1000, rho=0.5)] + [
            idat.sim_contig(rng, f"p{i}", n, rho=0.5) for i, n in enumerate((60, 100, 200))
        ]
        return idat.any_island(rng, contigs, _levels(), bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.many_replicons",
    ).assert_ok()  # fmt: skip


# --------------------------------------------------------------------------- recorded attacks


@pytest.mark.calibration
def test_recorded_autocorrelation_rho09() -> None:
    """Recorded (ADR-0021 decision 2, "Open limitation"): strong TU autocorrelation on a replicon
    at the edge of the calibrated domain inflates the FWER; those islands are ``uncalibrated``."""
    r = C.replicates(SUITE, full=False)
    rate = _record("autocorr_rho09", lambda rng: [idat.sim_contig(rng, "chr1", 1000, rho=0.9)], r)
    alpha = 0.05
    i = _levels().index(alpha)
    if rate[i] > C.level_bound(alpha, r):
        warnings.warn(
            f"rho 0.9 on 1,000 TUs: FWER {rate[i]:.3f} > {C.level_bound(alpha, r):.3f}"
            " (ADR-0021 decision 2: recorded, flagged uncalibrated, not gated)",
            stacklevel=1,
        )


@pytest.mark.calibration
def test_recorded_oriter_gradient() -> None:
    """Recorded: a smooth ori-ter gradient of 2 log2 peak to trough, which the LOESS detrending
    is meant to remove. Gating waits for the M5-04 ori-ter scenario (PR #119)."""
    r = C.replicates(SUITE, full=False)
    rate = _record(
        "oriter_gradient",
        lambda rng: [idat.sim_contig(rng, "chr1", 1000, rho=0.5, gradient=2.0)],
        r,
    )
    alpha = 0.05
    i = _levels().index(alpha)
    if rate[i] > C.level_bound(alpha, r):
        warnings.warn(
            f"ori-ter gradient: FWER {rate[i]:.3f} > {C.level_bound(alpha, r):.3f}"
            " (statistician decision pending; ADR-0021 measured 0.04 on giapsim)",
            stacklevel=1,
        )


@pytest.mark.calibration
def test_recorded_detrending_removes_the_gradient() -> None:
    """The detrended scores of a gradient replicon keep no trend: the LOESS fit explains the
    gradient (correlation of the fit with the planted shape > 0.9) and the residual mean over the
    first and second half of the replicon differ by less than 0.1 Blom units."""
    rng = C.replicate_rng(f"{SUITE}.attack.oriter_gradient.shape", 0)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.05)
    genes, tus, replicons = idat.frames([idat.sim_contig(rng, "chr1", 1000, rho=0.0,
                                                        gradient=4.0)])  # fmt: skip
    series = idat.build_series(genes, tus, replicons, cfg)
    ts = series[0]
    planted = idat.oriter_gradient(len(ts.tu_ids), 4.0)
    fit = np.asarray(ts.score) - np.asarray(ts.detrended)
    corr = float(np.corrcoef(fit, planted)[0, 1])
    half = len(ts.tu_ids) // 2
    drift = abs(float(np.mean(ts.detrended[:half]) - np.mean(ts.detrended[half:])))
    print(f"gradient: corr(fit, planted) = {corr:.3f}, residual half-drift = {drift:.3f}")
    assert corr > 0.9
    assert drift < 0.1
