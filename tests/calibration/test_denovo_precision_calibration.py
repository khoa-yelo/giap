"""Statistics-only calibration of the de novo precision estimator (M3-11; GIPPA_SPEC.md §7.2.5, §8.4
"De novo precision | known precision per bin | |bias| ≤ 0.05 | 50 | 300").

Generator under the estimator's own model: 10 score bins with known precision p_b; two length
classes ("8": c = 0.95, r = 0.05; "12-14": c = 0.99, r = 0.01) mixed 60/40 within each bin. A
correct peptide matches the target exactly with probability c_l, a wrong one with probability
r_l; every peptide matches the shuffled decoy with probability r_l. Scores are drawn inside
disjoint per-bin intervals with equal bin sizes, so the deciles recover the true bins exactly.

Claims over R replicates (``replicates_*.denovo_precision``), seed ``seeds.calibration_base``,
k = ``se_multiplier``:

1. bias: per bin |mean(prec_b) - p_b| ≤ ``denovo_precision_max_abs_bias`` (§8.4);
2. CI coverage: pooled share of percentile 95 % CIs containing p_b ≥ 0.95 - k·sqrt(.95·.05/N)
   (binomial band, N = R·10 intervals);
3. conservative π̂ (ADR-0010: π̂ handed to M3-10 must not be anti-conservative): pooled share
   of peptides' bins with π̂ = 1 - ci_low ≥ 1 - p_b is ≥ 0.975 - k·sqrt(.975·.025/N), and the
   mean π̂ per bin is ≥ 1 - p_b.

The genome-level version with gippasim de novo errors is M3-21.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gippa.mapping.denovo_precision import denovo_precision, pi_wrong_denovo, score_bins
from gippa.rng import for_step

P_BIN = np.linspace(0.2, 0.85, 10)
N_PER_BIN = 400
CLASSES = {"8": (0.95, 0.05, 0.6), "12-14": (0.99, 0.01, 0.4)}  # c, r, share


def _settings(tolerances: dict[str, object], key: str) -> tuple[int, int, float, float]:
    reps = tolerances[key]
    seeds = tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    return (
        int(reps["denovo_precision"]),
        int(seeds["calibration_base"]),
        float(tolerances["se_multiplier"]),  # type: ignore[arg-type]
        float(tolerances["denovo_precision_max_abs_bias"]),  # type: ignore[arg-type]
    )


def _simulate(rng: np.random.Generator) -> pd.DataFrame:
    parts = []
    labels = list(CLASSES)
    shares = np.array([CLASSES[k][2] for k in labels])
    for b, p in enumerate(P_BIN):
        n = N_PER_BIN
        cls_idx = rng.choice(len(labels), size=n, p=shares)
        c = np.array([CLASSES[labels[i]][0] for i in cls_idx])
        r = np.array([CLASSES[labels[i]][1] for i in cls_idx])
        correct = rng.random(n) < p
        target = np.where(correct, rng.random(n) < c, rng.random(n) < r)
        decoy = rng.random(n) < r
        parts.append(
            pd.DataFrame(
                {
                    "peptide": ["PEPTIDEK"] * n,
                    "score": b + rng.uniform(0.01, 0.99, size=n),
                    "length_class": [labels[i] for i in cls_idx],
                    "target_exact": target,
                    "decoy_exact": decoy,
                }
            )
        )
    df = pd.concat(parts, ignore_index=True)
    return df.sample(frac=1.0, random_state=int(rng.integers(1 << 31))).reset_index(drop=True)


def _run(r: int, seed: int, k: float, max_bias: float) -> None:
    sim_rng = for_step("denovo_precision_calibration", seed)
    c = {lab: v[0] for lab, v in CLASSES.items()}
    est = np.full((r, 10), np.nan)
    covered = np.zeros((r, 10), dtype=bool)
    conservative = np.zeros((r, 10), dtype=bool)
    pi_mean = np.full((r, 10), np.nan)
    for i in range(r):
        df = _simulate(sim_rng)
        true_bin = np.floor(df["score"].to_numpy()).astype(int)
        bins = score_bins(df["score"].to_numpy())
        np.testing.assert_array_equal(bins, true_bin)  # deciles recover the true bins
        rows = denovo_precision(df, c, False, for_step("placement", seed + i))
        assert [row.score_bin for row in rows] == list(range(10))
        pi = pi_wrong_denovo(rows, bins)
        for row in rows:
            b = row.score_bin
            assert row.n == N_PER_BIN
            assert row.precision is not None and row.ci_low is not None
            assert row.ci_high is not None
            est[i, b] = row.precision
            covered[i, b] = row.ci_low <= P_BIN[b] <= row.ci_high
            conservative[i, b] = 1.0 - row.ci_low >= 1.0 - P_BIN[b]
            pi_mean[i, b] = pi[bins == b].mean()

    # 1. bias per bin (§8.4)
    bias = est.mean(axis=0) - P_BIN
    assert np.all(np.abs(bias) <= max_bias), f"per-bin bias {np.round(bias, 4)}"

    n_int = r * 10
    # 2. two-sided 95 % coverage, binomial band
    cov = covered.mean()
    band = k * np.sqrt(0.95 * 0.05 / n_int)
    assert cov >= 0.95 - band, f"CI coverage {cov:.4f} < 0.95 - {band:.4f}"

    # 3. conservative π̂: one-sided 97.5 % lower bound
    cons = conservative.mean()
    band1 = k * np.sqrt(0.975 * 0.025 / n_int)
    assert cons >= 0.975 - band1, f"P(π̂ ≥ π) {cons:.4f} < 0.975 - {band1:.4f}"
    assert np.all(pi_mean.mean(axis=0) >= 1.0 - P_BIN), "mean π̂ below true π in some bin"


@pytest.mark.calibration
def test_denovo_precision_calibrated_fast(tolerances: dict[str, object]) -> None:
    _run(*_settings(tolerances, "replicates_fast"))


@pytest.mark.calibration_full
def test_denovo_precision_calibrated_full(tolerances: dict[str, object]) -> None:
    _run(*_settings(tolerances, "replicates_full"))
