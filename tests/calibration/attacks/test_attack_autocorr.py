"""Attack nulls for the M5-05 autocorrelation table (GIAP_SPEC.md §8.5 attack 1; ADR-0020
decision 19, M5-05 lead decisions 19 and 30; statistician-owned).

All four arms of ADR-0020 decision 19 are here (lead decision 33: M5-04's
``giapsim.scenarios`` is merged, so (a) and (b) have their ori-ter gradient). Each stays inside
the null the table claims and tries to push the lag tests above level:

* **(a) TU null plus a 2 log2 ori-ter cosine** (gated): ``scenarios.tu_null(oriter_log2=2.0)``
  -- within-TU structure, an oriC at ``scenarios.ORIC_FRAC`` and a cosine of 2 log2 peak to
  trough over the replichores, with ``detection = "abundance"`` so the gradient reaches
  ``p_detect`` too. The 10 % LOESS must remove the gradient: both statistics stay at level under
  the TU-block null at lags 1 and 20. Without re-detrending the same arm is anti-conservative
  (ADR-0020: pooled 0.34).
* **(b) label null plus the gradient** (gated): ``scenarios.label_null(oriter_log2=2.0)`` -- the
  gradient is then the *only* structure, so the gene-label null still holds after detrending
  (ADR-0020 evidence: pooled 0.010 / 0.051 / 0.104 at R = 400). Gated at lags 1 and 20 for both
  statistics under the gene-label null.
* **(c) short-gene prophage under the TU null** (gated): the TU-null scenario
  (``ar1_rho = 0``, ``tu_sd = 2``, decay 0.2, ``detection = "abundance"``) on a layout whose
  prophage genes are drawn short (M4-07 ``short_genes_in``), so low ``p_detect`` clusters along
  the replicon while abundance has no structure beyond operons (§8.5 attack 2's confounder inside
  attack 1's null). The Pearson residual join count (ADR-0020 Amendment 1) must stay at level
  under the TU-block null at lags 1 and 20; the raw count rejected 0.082-0.092 there.
* **(d) no re-detrending: teeth** (must fail, recorded as evidence that the gate bites): the
  superseded §7.3.5(c) procedure permutes the LOESS residuals **without** recomputing the trend
  (M5-05's original "The detrended x is permuted as is"). The smoother writes a negative
  autocorrelation into the observed residuals that the permuted residuals keep only partly, so
  the TU-block test over-rejects: on a 1,000-gene circular TU-null series (R = 200, B = 199)
  0.030 / 0.105 / 0.190 at lag 1 and 0.076 / 0.230 / 0.352 pooled over 20 lags at
  alpha 0.01 / 0.05 / 0.10, against 0.015 / 0.065 / 0.105 and 0.009 / 0.046 / 0.098 for the
  ADR-0020 procedure. The two-block gate must confirm the failure at alpha 0.05.
  Arm (d) builds its series from the engines alone, so it is independent of
  ``giap.context.autocorrelation`` -- it is the control for the arms that are not.

Criterion: §8.4 level through ``_gate.calibrate_levels`` / ``_criteria.calibrate`` (ADR-0015
two-block confirmation), R = ``replicates_fast``/``replicates_full``.``autocorrelation``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from _gate import calibrate_levels

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import autocorr_calib as acal
import autocorr_data as ad

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
pytestmark = pytest.mark.calsuite("autocorrelation")
FAST_GENES, FULL_GENES = 1000, 4500
B = acal.PERMUTATIONS
LAGS = list(acal.LAGS)
#: lags gated in every arm: the first and the last of ``context.moran_lags`` (lead decision 18)
GATED_LAGS = (1, 20)
#: ori-ter difference of arms (a) and (b), in log2 (ADR-0020 decision 19)
ORITER_LOG2 = 2.0


def _n(full: bool) -> int:
    return FULL_GENES if full else FAST_GENES


def _level_gates(sc: acal.Scenario, statistics: tuple[str, ...], null_model: str,
                 R: int) -> list[str]:  # fmt: skip
    """``check_level`` at every tolerance alpha for each statistic at the gated lags."""
    levels = list(C.load_tolerances().levels)
    failed: list[str] = []
    for statistic in statistics:
        for lag in GATED_LAGS:
            stream = f"{sc.stream}.{statistic}.{null_model}.lag{lag}"
            fn = acal.shared_fn(
                sc, stream, R, lambda r, s=statistic, m=null_model, k=lag: acal.lag_p(r, s, m, k)
            )
            gate = calibrate_levels(acal.SUITE, "level", levels, fn, R, stream=stream)
            failed += gate.detail() if gate.failed else []
    return failed


# --------------------------------------------------------------------------- arms (a) and (b)


def oriter_tu_null(n_genes: int) -> acal.Scenario:
    """Arm (a): the TU null with the 2 log2 ori-ter cosine on top."""
    return acal.tu_null(n_genes, oriter_log2=ORITER_LOG2,
                        name=f"attack_oriter_tu_null_{n_genes}")  # fmt: skip


def oriter_label_null(n_genes: int) -> acal.Scenario:
    """Arm (b): the label null with the 2 log2 ori-ter cosine as its only structure."""
    return acal.label_null(n_genes, oriter_log2=ORITER_LOG2,
                           name=f"attack_oriter_label_null_{n_genes}")  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_attack_oriter_gradient_tu_block_at_level(full: bool) -> None:
    """Arm (a): the LOESS must remove the ori-ter gradient, so the TU-block null holds."""
    R = C.replicates(acal.SUITE, full)
    failed = _level_gates(oriter_tu_null(_n(full)), ad.STATISTICS, "tu_block", R)
    assert not failed, failed


@pytest.mark.parametrize("full", MODES)
def test_attack_oriter_gradient_gene_label_at_level(full: bool) -> None:
    """Arm (b): with the gradient as the only structure the gene-label null holds too."""
    R = C.replicates(acal.SUITE, full)
    failed = _level_gates(oriter_label_null(_n(full)), ad.STATISTICS, "gene_label", R)
    assert not failed, failed


def test_attack_ab_plant_a_two_log2_oriter_cosine() -> None:
    """The arms really plant the gradient: an oriC, and a 2 log2 ori-minus-ter difference of the
    TU-level truth along the replichore position."""
    import giapsim.expression as ex
    import giapsim.layout as la
    import giapsim.scenarios as scen

    for sc in (oriter_tu_null(FAST_GENES), oriter_label_null(FAST_GENES)):
        assert sc.layout_cfg.oric_frac == scen.ORIC_FRAC
        assert sc.expr_cfg.oriter_log2 == ORITER_LOG2 and sc.expr_cfg.oriter_shape == "cosine"
        rng = np.random.default_rng(3)
        lay = la.simulate_layout(sc.layout_cfg, rng)
        assert (lay.features["type"] == "oric").sum() == 1
        genes = ex.simulate_expression(lay, sc.expr_cfg, rng).genes
        pos = genes["replichore_position_true"].to_numpy(float)
        a = genes["abundance_log2_true"].to_numpy(float)
        near_ori, near_ter = pos < 0.1, pos > 0.9
        drop = float(np.nanmean(a[near_ori]) - np.nanmean(a[near_ter]))
        assert 1.0 < drop < 3.0, (sc.name, drop)


# --------------------------------------------------------------------------- arm (c)


def short_prophage_tu_null(n_genes: int) -> acal.Scenario:
    """Arm (c): the TU-null scenario on a layout whose prophage genes are short (M4-07 attack
    null 2, ``short_genes_in``)."""
    return acal.tu_null(n_genes, short_genes_in="prophage",
                        name=f"attack_short_prophage_tu_null_{n_genes}")  # fmt: skip


@pytest.mark.parametrize("full", MODES)
def test_attack_short_prophage_join_count_at_level(full: bool) -> None:
    """Arm (c): clustered detectability must not make the join count reject under the TU null."""
    R = C.replicates(acal.SUITE, full)
    failed = _level_gates(short_prophage_tu_null(_n(full)), ("join_count_bb",), "tu_block", R)
    assert not failed, failed


def test_attack_c_layout_draws_short_prophage_genes() -> None:
    """The arm really plants the confounder: prophage genes are shorter than the rest."""
    import giapsim.layout as la

    sc = short_prophage_tu_null(FAST_GENES)
    lay = la.simulate_layout(sc.layout_cfg, np.random.default_rng(0))
    ids = lay.region_genes["region_id"].astype(str).str.startswith("prophage")
    tags = set(lay.region_genes.loc[ids, "locus_tag"])
    length = lay.genes.set_index("locus_tag")["aa"].str.len().astype(float)
    inside, outside = length[length.index.isin(tags)], length[~length.index.isin(tags)]
    assert len(inside) > 10
    assert inside.median() < 0.6 * outside.median(), (inside.median(), outside.median())


# --------------------------------------------------------------------------- arm (d)


def _tu_series(rng: np.random.Generator, n: int, p_geom: float = 0.55, tu_sd: float = 2.0,
               gene_sd: float = 0.5) -> tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """A TU-null series: iid TU effects over geometric TU sizes plus iid gene noise."""
    sizes: list[int] = []
    while sum(sizes) < n:
        sizes.append(int(rng.geometric(p_geom)))
    sizes[-1] -= sum(sizes) - n
    if sizes[-1] == 0:
        sizes.pop()
    blocks = np.repeat(np.arange(len(sizes)), sizes)
    y = tu_sd * rng.normal(size=len(sizes))[blocks] + gene_sd * rng.normal(size=n)
    return y, blocks


def _no_redetrend_p(rng: np.random.Generator, lag: int = 1, n: int = FAST_GENES) -> float:
    """Two-sided p of Moran's I at ``lag`` under the superseded procedure: detrend once, then
    permute the residuals by TU block without recomputing the trend."""
    from giap.stats import autocorr as ac

    y, blocks = _tu_series(rng, n)
    x = ad.oracle_values(y, True)
    orders = ac.tu_block_permutations(ac.tu_blocks([str(b) for b in blocks], True), B, rng)
    null = ac.morans_i_matrix(x[orders], [lag], True)[:, 0]
    return float(ac.lag_test(ac.morans_i(x, lag, True), null).p_value)


@pytest.mark.calibration
def test_teeth_permuting_residuals_without_redetrending_fails(teeth_log: Path) -> None:
    """Arm (d): the spec's "detrend then permute" must be a confirmed fail at alpha 0.05, so the
    level gates the ADR-0020 procedure passes have teeth (M5-05 lead decision 19)."""
    R = C.replicates(acal.SUITE, full=True)
    rec = C.calibrate(
        acal.SUITE,
        "level",
        0.05,
        _no_redetrend_p,
        R,
        stream="autocorrelation.attack.teeth_no_redetrend_tu_block_lag1",
        log_dir=teeth_log,
    )
    assert rec["status"] == "fail", f"no re-detrending is not a confirmed fail: {rec}"


def test_teeth_arm_d_uses_a_true_tu_block_null() -> None:
    """The arm's series has no structure beyond its TUs, so a rejection is the procedure's
    fault: the ADR-0020 values (re-detrended per replicate) are at level on the same draws."""
    from giap.stats import autocorr as ac
    from giap.stats.smooth import loess_detrend, loess_matrix

    n = 300
    rng = np.random.default_rng(4)
    y, blocks = _tu_series(rng, n)
    pos = np.arange(n, dtype=float)
    wide = loess_detrend(y - y.mean(), pos, ad.WIDE_SPAN, circular=True, period=float(n)).resid
    h = loess_matrix(pos, ad.SPAN, circular=True, period=float(n))
    orders = ac.tu_block_permutations(ac.tu_blocks([str(b) for b in blocks], True), 49, rng)
    moved = wide[orders]
    fixed = ac.morans_i_matrix(moved - moved @ h.T, [1], True)[:, 0]
    broken = ac.morans_i_matrix((wide - h @ wide)[orders], [1], True)[:, 0]
    assert np.isfinite(fixed).all() and np.isfinite(broken).all()
    # the re-detrended replicates carry the smoother's negative lag-1 bias, the permuted ones less
    assert float(np.mean(fixed)) < float(np.mean(broken))
