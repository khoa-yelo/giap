"""Attack nulls for the M5-05 autocorrelation table (GIAP_SPEC.md §8.5 attack 1; ADR-0020
decision 19, M5-05 lead decisions 19 and 30; statistician-owned).

Two of the four arms are here; each stays inside the null the table claims and tries to push the
lag tests above level:

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

**Pending (M5-04, PR #119):** ADR-0020 arms **(a)** TU null plus a 2 log2 ori-ter cosine and
**(b)** label null plus the same gradient need ``giapsim.scenarios.oriter_gradient`` (an oriC
feature and ``ExpressionConfig.oriter_log2``), which is not merged. They are the arms that show
the LOESS removes a smooth gradient, so this file is incomplete until #119 lands; add them then.

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
#: lags gated in arm (c): the first and the last of ``context.moran_lags`` (lead decision 18)
GATED_LAGS = (1, 20)


# --------------------------------------------------------------------------- arm (c)


def short_prophage_tu_null(n_genes: int) -> acal.Scenario:
    """The TU-null scenario on a layout whose prophage genes are short (M4-07 attack null 2)."""
    import giapsim.layout as la

    base = acal.tu_null(n_genes)
    return acal.Scenario(
        f"attack_short_prophage_{base.name}",
        n_genes,
        base.expr_kw,
        {"regions": la.DEFAULT_REGIONS, "short_genes_in": "prophage"},
        nulls=("tu_block",),
    )


@pytest.mark.parametrize("full", MODES)
def test_attack_short_prophage_join_count_at_level(full: bool) -> None:
    """Arm (c): clustered detectability must not make the join count reject under the TU null."""
    R = C.replicates(acal.SUITE, full)
    sc = short_prophage_tu_null(FULL_GENES if full else FAST_GENES)
    levels = list(C.load_tolerances().levels)
    failed: list[str] = []
    for lag in GATED_LAGS:
        stream = f"{sc.stream}.join_count_bb.tu_block.lag{lag}"
        fn = acal.shared_fn(
            sc, stream, R, lambda r, k=lag: acal.lag_p(r, "join_count_bb", "tu_block", k)
        )
        gate = calibrate_levels(acal.SUITE, "level", levels, fn, R, stream=stream)
        failed += gate.detail() if gate.failed else []
    assert not failed, failed


def test_attack_c_layout_draws_short_prophage_genes() -> None:
    """The arm really plants the confounder: prophage genes are shorter than the rest."""
    import giapsim.layout as la

    sc = short_prophage_tu_null(FAST_GENES)
    lay = la.simulate_layout(
        acal.layout_config(sc.n_genes, **dict(sc.layout_kw)), np.random.default_rng(0)
    )
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
