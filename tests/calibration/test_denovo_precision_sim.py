"""De novo precision on gippasim, through the mapping steps (issue M3-21; GIPPA_SPEC.md §7.2.5
"De novo precision, estimated from the genome", §8.4 "De novo precision | known precision per bin
| |bias| <= 0.05 | 50 | 300"; docs/methods/denovo-precision.md, forwarded to M3-21).

Replicate i (R = ``replicates_*.denovo_precision``, genome seed from ``replicate_rng`` of stream
``denovo_precision.sim``): 100-kb genome, M3-19 database table (pi_wrong 0.05, 3 runs), its M3-20 de
novo reads as input and the M3-20 reference peptides; the ``denovo_precision`` rows of
``mapping_fdr`` are compared with the truth.

**Estimand.** GIPPA bins the *distinct* de novo query peptides (preflight collapses rows by
peptide), so the precision of bin b is the share of correct peptides among the bin's queries
(score deciles, right-closed, recomputed here with gippasim's ``decile_bins`` on the queries' input
scores). ``truth/denovo_precision.tsv`` is the per-*read* precision of gippasim's deciles over all
reads: collapsing repeated correct reads of one peptide while errors stay distinct puts the
per-peptide precision 0.01-0.085 below it (bins 1-6 measured at 100 kb, 2026-10-20), so a
comparison with the table would fail by construction although the estimator is unbiased; the
table bias is recorded in the diagnostics, not gated.

Claims, k = ``se_multiplier``:

1. bias: per reported bin |mean(prec_b - realised_b)| <= ``denovo_precision_max_abs_bias``;
2. conservative π̂ (ADR-0010): per bin mean(π̂_b) >= mean(1 - realised_b) - k SE, π̂_b = 1 - ci_low.

Recorded (diagnostics jsonl; the two directions that make π̂ anti-conservative, forwarded by the
M3-11 review): c_ref - c_true (c_true = share of correct de novo queries with an exact target hit)
and the wrong reads' exact target / decoy match ratio.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
import yaml
from _placement_suite import write_diagnostics
from helpers.mapping_steps import read_tsv_df
from helpers.placement_sims import hits_of, run_chain

pytestmark = pytest.mark.calibration

SUITE = "denovo_precision"


def _replicate(seed: int, d: Path) -> dict[str, np.ndarray | float]:
    from gippa.steps.mapping import length_label, select_queries
    from gippasim.denovo import DenovoConfig, decile_bins
    from gippasim.proteomics import ProteomicsConfig

    ch = run_chain(
        "denovo", seed, d, proteomics=ProteomicsConfig(pi_wrong=0.05, search_db="other"),
        denovo=DenovoConfig(),
    )  # fmt: skip
    f = pd.read_parquet(ch.avail["mapping_fdr"])
    f = f[f["section"] == "denovo_precision"].set_index("score_bin")
    est = f["precision"].astype(float).reindex(range(10)).to_numpy()
    pi_hat = 1.0 - f["ci_low"].astype(float).reindex(range(10)).to_numpy()
    q = select_queries(pd.read_parquet(ch.avail["peptides"]), ch.cfg)
    classes = ch.cfg.mapping.length_classes
    q = q[[length_label(len(p), classes) is not None for p in q["peptide"]]]
    score = pd.to_numeric(q["score"], errors="coerce").to_numpy(dtype=float)
    keep = ~np.isnan(score)
    peps = q["peptide"].to_numpy()[keep]
    b = decile_bins(score[keep])
    ok = ~pd.Series(peps).isin(ch.truth.wrong).to_numpy()
    real = pd.Series(ok).groupby(b).mean().reindex(range(10)).to_numpy()
    table = read_tsv_df(ch.res.truth["denovo_precision"]).set_index("score_bin")
    true_tab = table["precision_true"].astype(float).reindex([str(i) for i in range(10)])
    h = hits_of(ch.avail)
    exact = h[h["n_diff"] == 0]
    tx = set(exact.loc[exact["search"] == "target", "peptide"])
    dx = set(exact.loc[exact["search"] == "decoy", "peptide"])
    p = pd.Series(peps)
    in_t, in_d = p.isin(tx).to_numpy(), p.isin(dx).to_numpy()
    c_true = float(in_t[ok].mean()) if ok.any() else np.nan
    c_ref = float(np.nanmean(f["c_value"].astype(float))) if len(f) else np.nan
    shutil.rmtree(d, ignore_errors=True)
    return {
        "est": est,
        "pi_hat": pi_hat,
        "real": real,
        "table": true_tab.to_numpy(),
        "c_gap": c_ref - c_true,
        "wrong_t": float(in_t[~ok].mean()) if (~ok).any() else np.nan,
        "wrong_d": float(in_d[~ok].mean()) if (~ok).any() else np.nan,
    }


def _run(full: bool, tmp: Path) -> None:
    tol = C.load_tolerances()
    R = C.replicates(SUITE, full)
    reps = [
        _replicate(int(C.replicate_rng(f"{SUITE}.sim", i).integers(0, 2**31 - 1)), tmp / f"r{i}")
        for i in C.replicate_range(R)
    ]
    est = np.vstack([r["est"] for r in reps])
    real = np.vstack([r["real"] for r in reps])
    pih = np.vstack([r["pi_hat"] for r in reps])
    tab = np.vstack([r["table"] for r in reps])
    bias = np.nanmean(est - real, axis=0)
    write_diagnostics(
        {
            "suite": f"{SUITE}.sim",
            "R_run": len(reps),
            "bias_vs_realised": bias.tolist(),
            "bias_vs_truth_table": np.nanmean(est - tab, axis=0).tolist(),
            "reported_share": np.mean(~np.isnan(est), axis=0).tolist(),
            "c_ref_minus_c_true": float(np.nanmean([r["c_gap"] for r in reps])),
            "wrong_target_exact": float(np.nanmean([r["wrong_t"] for r in reps])),
            "wrong_decoy_exact": float(np.nanmean([r["wrong_d"] for r in reps])),
        }
    )
    limit = tol.se_multiplier
    raw = yaml.safe_load(C.TOLERANCES_PATH.read_text())
    max_bias = float(raw["denovo_precision_max_abs_bias"])
    reported = np.sum(~np.isnan(est), axis=0) >= max(2, len(reps) // 2)
    assert reported.sum() >= 8, f"too few reported bins: {np.sum(~np.isnan(est), axis=0)}"
    bad = [(b, round(float(bias[b]), 4)) for b in range(10) if reported[b]
           and abs(bias[b]) > max_bias]  # fmt: skip
    assert not bad, f"|bias| > {max_bias} in bins {bad}; per-bin bias {np.round(bias, 4)}"
    d = pih - (1.0 - real)
    n = np.sum(~np.isnan(d), axis=0)
    mean_d = np.nanmean(d, axis=0)
    se = np.nanstd(d, axis=0, ddof=1) / np.sqrt(np.maximum(n, 1))
    low = [(b, round(float(mean_d[b]), 4)) for b in range(10) if reported[b]
           and mean_d[b] < -limit * se[b]]  # fmt: skip
    assert not low, f"π̂ = 1 - ci_low below the realised wrong share in bins {low}"


def test_denovo_precision_on_gippasim(tmp_path: Path) -> None:
    _run(False, tmp_path)


@pytest.mark.calibration_full
def test_denovo_precision_on_gippasim_full(tmp_path: Path) -> None:
    _run(True, tmp_path)
