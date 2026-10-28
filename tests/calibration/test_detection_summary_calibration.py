"""Calibration of the M4-16 basic O - E confidence interval (§7.3.5(a), §8.2, §8.4).

Claim (gated through ``_criteria.calibrate``, ADR-0015 two-block confirmation): **the 95 %
TU-block bootstrap CI of ``o_minus_e`` reported by ``basic_rows`` for the chromosome replicon
group covers the true expectation 0 in at least 95 % of null replicates**, i.e. its non-coverage
rate is <= 0.05 + k sqrt(0.05 * 0.95 / R) (the ``fwer`` criterion applied to 0/1 miss indicators;
k = ``se_multiplier``, R = ``replicates_fast/full.regions``, lead decision 6: reuse existing R).

Null generator, per replicate: a gippasim M4-07 layout (one circular chromosome of 400 genes, real
TU structure and feature classes) and the M4-08 ``p_detect_true``; detection is drawn
comonotonically within each TU (one uniform u per TU, d_g = 1{u < p_g}), so each gene keeps
P(d_g = 1) = p_g exactly (E[O - E] = 0 with ``p_detect = p_detect_true``) while detections inside a
TU are strongly positively correlated. A gene-level bootstrap under-covers here; only resampling
whole TUs (§8.2) passes. B = 1000 (``stats.bootstrap`` default).
"""

from __future__ import annotations

import sys
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_summary_fixture as fx

MODES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
SUITE = "regions"
STREAM = "detection_summary.o_minus_e_ci.chromosome"
B = 1000
N_GENES = 400


def null_genes(rng: np.random.Generator) -> pd.DataFrame:
    from gippasim.expression import ExpressionConfig, simulate_expression
    from gippasim.layout import LayoutConfig, simulate_layout

    layout = simulate_layout(LayoutConfig(n_genes_chromosome=N_GENES), rng)
    expr = simulate_expression(layout, ExpressionConfig(), rng)
    g = layout.genes
    tags = g["locus_tag"].astype(str)
    tu = tags.map(layout.truth.set_index("locus_tag")["tu_id"].astype(str)).to_numpy()
    p = tags.map(expr.genes.set_index("locus_tag")["p_detect_true"]).astype(float).to_numpy()
    tus = sorted(set(tu))
    u = dict(zip(tus, rng.random(len(tus)), strict=True))
    d = np.array([u[t] < pi for t, pi in zip(tu, p, strict=True)], dtype=bool)
    partial = g["partial"].astype(bool).to_numpy()
    status = np.where(partial, "partial_excluded", np.where(d, "detected", "not_detected"))
    return pd.DataFrame(
        {
            "locus_tag": tags.to_numpy(),
            "contig": g["contig"].astype(str).to_numpy(),
            "feature_type": g["feature_type"].astype(str).to_numpy(),
            "hypothetical": g["hypothetical"].astype(bool).to_numpy(),
            "partial": partial,
            "tu_id": tu,
            "detected": d & ~partial,
            "detection_status": status,
            "p_detect": np.where(partial, np.nan, p),
            "cog_category": [None] * len(g),
        }
    )


def miss(rng: np.random.Generator) -> float:
    from gippa.context.summary import basic_rows

    genes = null_genes(rng)
    out = basic_rows(genes, B, np.random.default_rng(rng.integers(2**63)))
    r = fx.metric(out, "replicon", "chr1", "o_minus_e")
    lo, hi = float(r["ci_low"]), float(r["ci_high"])
    assert np.isfinite(lo) and np.isfinite(hi), r
    assert int(fx.metric(out, "replicon", "chr1", "n_tu")["value"]) >= 100
    return 0.0 if lo <= 0.0 <= hi else 1.0


@pytest.mark.parametrize("full", MODES)
def test_o_minus_e_ci_covers_zero_under_tu_correlated_null(full: bool) -> None:
    R = C.replicates(SUITE, full)
    rec = C.calibrate(SUITE, "fwer", 0.05, miss, R, stream=STREAM)
    assert rec["status"] in ("pass", "unconfirmed", C.SHARD_STATUS, C.SKIPPED_STATUS), rec


def test_null_generator_is_tu_correlated_and_unbiased() -> None:
    """Generator sanity (not a gate): E[d] = p per gene and strong within-TU correlation."""
    rng = np.random.default_rng(5)
    resid, pair = [], []
    for _ in range(30):
        g = null_genes(rng)
        g = g[g["detection_status"] != "partial_excluded"]
        r = g["detected"].astype(float) - g["p_detect"]
        resid.append(float(r.sum()) / np.sqrt(float((g["p_detect"] * (1 - g["p_detect"])).sum())))
        same = g.groupby("tu_id")["detected"].agg(["sum", "size"])
        multi = same[same["size"] > 1]
        pair.append(float(((multi["sum"] == 0) | (multi["sum"] == multi["size"])).mean()))
    assert abs(np.mean(resid)) < 3 * np.std(resid) / np.sqrt(len(resid)) + 0.5
    assert np.mean(pair) > 0.6  # most multi-gene TUs are all-or-nothing
