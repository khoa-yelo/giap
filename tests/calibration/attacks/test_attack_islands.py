"""Attack nulls for the M5-06 island block bootstrap (GIPPA_SPEC.md §8.5; statistician-owned).

Each generator stays inside the null the step claims (no island planted) and tries to push the
genome-wide FWER above level. Arms A-E are the statistician's proposal
(local/scratch/m5-06-stat/work/proposed_attacks.py.txt, lead decision 51); the large-TU, draft
and many-replicon arms are the rest of the §8.5 list. Criterion: ``fwer`` through
``_gate.calibrate_levels`` (ADR-0015 two-block confirmation), R = ``replicates_fast/full.islands``,
B = 199 bootstrap replicates and 199 block-length permutations.

Amplitudes are quoted in **background SD** as well as log2 (decision 49): gippasim's AR(1)
innovation SD is 2 log2, the hand-built ``sim_contig`` background is 1, so a 2-log2 gippasim
gradient is 0.87 background SD peak to trough and the hand-built ``gradient=2.0`` is the same
point as a 4-log2 gippasim gradient (1.73 SD).

Rates measured by the statistician with this branch's oracle at R = 1,000, B = PERM = 199,
``level_bound(0.05, 1000) = 0.0638`` (levels [0.01, 0.05, 0.10]):

* A ``oriter_gradient(1800, 0.5, 2.0)`` (2 log2 = 0.87 SD, a 4-fold ori:ter ratio):
  [0.003, 0.042, 0.094] - **gated**; the same scenario without the gradient is
  [0.001, 0.029, 0.065].
* B the same with ``detection = "abundance"``: [0.000, 0.017, 0.057] - **gated**; detection then
  carries the TU structure, Moran's I sees it (l-hat 3) and the block grows to 7.
* C 4 log2 (1.73 SD, a 16-fold ratio, outside any growth regime): [0.009, 0.062, 0.136], and the
  hand-built equivalent [0.009, 0.066, 0.156] - **recorded**, not gated (decision 49).
* D ``ar1_background`` at rho 0.9: [0.020, 0.113, 0.218] at 557 TUs, [0.018, 0.095, 0.206] at
  1,001 TUs, [0.002, 0.024, 0.085] at 2,501 TUs - **recorded** (ADR-0021 Amendment 1: the
  dependence regime, not the TU count, decides; at rho 0.5 every size is calibrated).
* E the structural low-direction ceiling (decision 50): silencing whole TUs cannot produce a low
  island, because sqrt(w)*|Phi^-1(q/2)| ~ 3.4 at q = 0.45, w = 20 against a critical value ~ 4.3.
* large TUs (sizes 1 + Geometric(0.15), mean ~6.7 genes, 2,000 genes, rho 0.5): [0, 0, 0.02] at
  R = 50 - gated; draft assembly ``draft_ar1(0.5)`` (500-gene chromosome + 133 linear 30-gene
  contigs): gated, and the contigs are never scanned (decision 40, ``few_values``); chromosome
  plus three plasmids: [0, 0.02, 0.04] at R = 50 - gated, one genome-wide family.
"""

from __future__ import annotations

import dataclasses
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
    print(
        f"recorded attack {name}: R = {r}, FWER {np.round(rate, 3).tolist()} at {list(levels)}"
        f" (bounds {[round(C.level_bound(a, r), 3) for a in levels]})"
    )
    return rate


def _warn_above_bound(name: str, rate: np.ndarray, r: int, note: str) -> None:
    i = _levels().index(0.05)
    if rate[i] > C.level_bound(0.05, r):
        warnings.warn(f"{name}: FWER {rate[i]:.3f} > {C.level_bound(0.05, r):.3f} ({note})",
                      stacklevel=2)  # fmt: skip


def _abundance(pair: Any) -> Any:
    """The scenario with detection driven by abundance (``scenarios.ABUNDANCE_BETA``)."""
    from gippasim import scenarios as sc

    lcfg, ecfg = pair
    return lcfg, dataclasses.replace(ecfg, detection="abundance", beta=sc.ABUNDANCE_BETA)


# ----------------------------------------------------- A. gated: the M5-04 ori-ter scenario


@pytest.mark.parametrize("full", MODES)
def test_attack_oriter_gradient_scenario(full: bool) -> None:
    """M5-04's ori-ter scenario at a realistic amplitude: 2 log2 (0.87 background SD) ori minus
    ter, a 4-fold copy-number ratio, oriC at 0.8, 1,800 genes (~1,000 scored TUs), rho 0.5.

    Left censoring compresses the Blom scores of the low-abundance half against the censored atom
    while the high half spreads into the upper tail, so the detrended series is heteroscedastic
    along the replicon and the stationary block bootstrap is slightly too narrow (decision 49).
    """
    from gippasim import scenarios as sc

    def rep(rng: np.random.Generator) -> np.ndarray:
        contigs = idat.scenario_contigs(rng, sc.oriter_gradient(1800, 0.5, 2.0))
        return idat.any_island(rng, contigs, _levels(), bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.oriter_gradient_scenario",
    ).assert_ok()  # fmt: skip


# ------------------------------------- B. gated: the gradient with detection from abundance


@pytest.mark.parametrize("full", MODES)
def test_attack_oriter_gradient_detection_abundance(full: bool) -> None:
    """The same gradient with ``detection = "abundance"``: censoring then clusters at the
    terminus, the realistic worst case for a rank score. Measured more conservative than
    covariates-only detection, because detection carries the TU structure, Moran's I sees it
    (l-hat 3) and the block grows to 7."""
    from gippasim import scenarios as sc

    def rep(rng: np.random.Generator) -> np.ndarray:
        contigs = idat.scenario_contigs(rng, _abundance(sc.oriter_gradient(1800, 0.5, 2.0)))
        return idat.any_island(rng, contigs, _levels(), bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.oriter_gradient_abundance",
    ).assert_ok()  # fmt: skip


# ----------------------------------------- C. recorded: the gradient beyond the real range


@pytest.mark.calibration
def test_recorded_oriter_gradient_extreme() -> None:
    """4 log2 (1.73 background SD) ori minus ter, a 16-fold ratio outside any growth regime, and
    the hand-built equivalent (``sim_contig(gradient=2.0)``, whose background SD is 1): the same
    point, both within the bound at alpha 0.05 and above it at alpha 0.10 (decision 49).
    Recorded with the mechanism, not gated."""
    from gippasim import scenarios as sc

    r = C.replicates(SUITE, full=False)
    rate = _record(
        "oriter_gradient_4log2",
        lambda rng: idat.scenario_contigs(rng, sc.oriter_gradient(1800, 0.5, 4.0)),
        r,
    )
    _warn_above_bound("ori-ter gradient 4 log2", rate, r, "decision 49: recorded, not gated")
    hand = _record(
        "oriter_gradient_handbuilt",
        lambda rng: [idat.sim_contig(rng, "chr1", 1000, rho=0.5, gradient=2.0)],
        r,
    )
    _warn_above_bound("ori-ter gradient 1.73 SD (hand-built)", hand, r, "decision 49")


@pytest.mark.calibration
def test_recorded_detrending_removes_the_gradient() -> None:
    """The mechanism is not a leftover trend: the LOESS fit explains the planted gradient
    (correlation > 0.9) and the residual means of the two halves differ by < 0.1 Blom units, yet
    the detrended series is heteroscedastic (decision 49)."""
    rng = C.replicate_rng(f"{SUITE}.attack.oriter_gradient.shape", 0)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.05)
    genes, tus, replicons = idat.frames([idat.sim_contig(rng, "chr1", 1000, rho=0.0, gradient=4.0)])
    series = idat.build_series(genes, tus, replicons, cfg)
    ts = series[0]
    planted = idat.oriter_gradient(len(ts.tu_ids), 4.0)
    fit = np.asarray(ts.score) - np.asarray(ts.detrended)
    corr = float(np.corrcoef(fit, planted)[0, 1])
    half = len(ts.tu_ids) // 2
    drift = abs(float(np.mean(ts.detrended[:half]) - np.mean(ts.detrended[half:])))
    scale = np.array([float(np.std(ts.detrended[i : i + 50])) for i in range(0, half * 2 - 50, 25)])
    print(
        f"gradient: corr(fit, planted) = {corr:.3f}, residual half-drift = {drift:.3f},"
        f" local-SD q90/q10 = {np.quantile(scale, 0.9) / np.quantile(scale, 0.1):.2f}"
    )
    assert corr > 0.9
    assert drift < 0.1


# ------------------------------- D. recorded: strong autocorrelation at three replicon sizes


@pytest.mark.calibration
@pytest.mark.parametrize("n_genes", [1000, 1800, 4500])
def test_recorded_autocorrelation_rho09_sizes(n_genes: int) -> None:
    """ADR-0021 Amendment 1: the failure is the dependence, not the TU count. gippasim
    ``ar1_background`` at rho 0.9 breaches at 557 and 1,001 scored TUs but not at 2,501, where
    l-hat ~ 8 makes l' = 17 span the 20-TU windows; at rho 0.5 every size is calibrated."""
    from gippasim import scenarios as sc

    r = C.replicates(SUITE, full=False)
    rate = _record(
        f"ar1_rho09.n{n_genes}",
        lambda rng: idat.scenario_contigs(rng, sc.ar1_background(n_genes, 0.9)),
        r,
    )
    _warn_above_bound(
        f"ar1_background rho 0.9, {n_genes} genes", rate, r,
        "ADR-0021 Amendment 1: outside the claimed dependence regime, recorded, not gated",
    )  # fmt: skip


# ----------------------------------- E. structural: the low direction cannot be reported


@pytest.mark.calibration
@pytest.mark.parametrize("silenced", [10, 20, 40])
def test_low_direction_ceiling(silenced: int) -> None:
    """Decision 50, quantified. Left censoring caps a low window of w TUs at
    ``sqrt(w) * |Phi^-1(q / 2)|`` with q the undetected share: 3.4 at q = 0.45, w = 20, against a
    genome-wide critical value of about 4.3. Silence *every* gene of ``silenced`` consecutive TUs
    of an rho = 0.5 null: no low island is reported even at infinite effect. The ceiling and the
    95th percentile of M_b are printed so M5-07's low-power row can be read as structural."""
    from gippa.stats.scan import scan_max

    rng = C.replicate_rng(f"{SUITE}.attack.low_ceiling.{silenced}", 0)
    c = idat.sim_contig(rng, "chr1", 1800, rho=0.5)  # one gene per TU
    d = list(c.d)
    for i in range(300, 300 + silenced):  # silence every gene of `silenced` consecutive TUs
        d[i] = 0
    contigs = [idat.Contig(c.name, d=d, p=c.p, a=c.a, circular=c.circular, tu=c.tu)]
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.05)
    genes, tus, replicons = idat.frames(contigs)
    out, scan = idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg, rng)
    ts = scan.series[0]
    q = idat.censored_share(genes, "chr1")
    w = idat.max_window(ts, cfg)
    ceiling = idat.ceiling_from_q(q, w)
    crit = float(np.quantile(scan.null_max, 0.95))
    raw = scan_max(ts.score, cfg.context.islands.min_tus, cfg.context.islands.max_tus, ts.circular)
    worst_raw, worst = float(np.min(raw.s)), float(np.min(scan.windows[0].s))
    neg = scan.windows[0].s < 0
    p_low = float(np.min(scan.p_value[0][neg])) if neg.any() else float("nan")
    print(
        f"low ceiling: {silenced} TUs silenced, q = {q:.2f}, w = {w}, ceiling {ceiling:.2f},"
        f" M_b q95 {crit:.2f}, most negative S {worst_raw:.2f} raw / {worst:.2f} detrended,"
        f" best low p {p_low:.3f}"
    )
    assert worst_raw >= -ceiling  # exact: every TU score is at least the censored score
    assert ceiling < crit, "the low direction is structurally unreachable at this censored share"
    assert "low" not in set(out["direction"])
    assert p_low > cfg.stats.alpha.context


# --------------------------------------------------------------------------- the §8.5 list


@pytest.mark.parametrize("full", MODES)
def test_attack_large_tus(full: bool) -> None:
    """TU sizes 1 + Geometric(0.15) (mean ~6.7 genes): heteroscedastic TU scores, because a TU
    mean of many genes has a smaller variance than a single-gene TU
    (docs/methods/islands.md "Why S/sqrt(n_w) is not standardised")."""

    def rep(rng: np.random.Generator) -> np.ndarray:
        contigs = [idat.sim_contig(rng, "chr1", 2000, rho=0.5, tu_p=0.15)]
        return idat.any_island(rng, contigs, _levels(), bootstrap=B, permutations=PERM)

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.large_tus",
    ).assert_ok()  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_attack_draft_assembly(full: bool) -> None:
    """``scenarios.draft_ar1(0.5)``: a 500-gene circular chromosome and 133 linear 30-gene
    contigs (ADR-0021 measured 0.003)."""
    from gippasim import scenarios as sc

    def rep(rng: np.random.Generator) -> np.ndarray:
        return idat.any_island(
            rng, idat.scenario_contigs(rng, sc.draft_ar1(0.5)), _levels(),
            bootstrap=B, permutations=PERM,
        )  # fmt: skip

    calibrate_levels(
        SUITE, "fwer", _levels(), rep, C.replicates(SUITE, full=full),
        stream=f"{SUITE}.attack.draft_assembly",
    ).assert_ok()  # fmt: skip


@pytest.mark.calibration
def test_draft_contigs_are_never_scanned() -> None:
    """Decisions 40/41: the 30-gene linear contigs give a span of k < 4, so they are flagged
    ``few_values``, contribute nothing to M_b and never carry an island."""
    from gippasim import scenarios as sc

    rng = C.replicate_rng(f"{SUITE}.attack.draft_assembly.structure", 0)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=0.10)
    genes, tus, replicons = idat.frames(idat.scenario_contigs(rng, sc.draft_ar1(0.5)))
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
