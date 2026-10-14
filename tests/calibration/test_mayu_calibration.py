"""Fast calibration of the MAYU closed form under its own model (M4-02; GIAP_SPEC.md §7.3.2, §8.4).

Generator with known truth: G genes in two size groups (bins aligned with the groups), a fixed set
of truly present genes, and F false peptides landing on genes with probability ∝ a_g. A gene hit
only by false peptides is a false gene. Claims checked over R = ``replicates_fast.gene_fdr``
replicates, seed ``seeds.calibration_base``, k = ``se_multiplier``:

1. binomial band: the pooled share of truly absent genes that are hit lies within
   p̄ ± k·sqrt(p̄(1 - p̄)/(R·M)), p̄ = Poisson-binomial mean / M of the exact per-gene hit
   probabilities (equal to the binned p_b, checked); hits within a replicate are negatively
   correlated, so this band is conservative;
2. FP unbiased: mean(FP_est - FP_true) within ±k·SE;
3. FDR criterion (§8.4): mean FDP <= mean estimated FDR + k·SE of the paired difference.

Gene-FDR calibration of the full tier logic (picked, MAYU fallback, decoy loss per ADR-0009) is
M4-12's ``gene_fdr`` suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from giap.stats.mayu import mayu_gene_fdr
from giap.stats.poisson_binomial import poisson_binomial_pmf

pytestmark = pytest.mark.calibration

N_SMALL, N_LARGE = 1200, 800
A_SMALL, A_LARGE = 4.0, 16.0
N_BINS = 10  # 200 genes per bin: aligned with the two size groups
F_FALSE = 120
PRESENT = 0.4


def _settings(tolerances: dict[str, object]) -> tuple[int, int, float]:
    reps = tolerances["replicates_fast"]
    seeds = tolerances["seeds"]
    assert isinstance(reps, dict) and isinstance(seeds, dict)
    return int(reps["gene_fdr"]), int(seeds["calibration_base"]), float(tolerances["se_multiplier"])  # type: ignore[arg-type]


def test_mayu_closed_form_calibrated_under_its_model(tolerances: dict[str, object]) -> None:
    r, seed, k = _settings(tolerances)
    rng = np.random.default_rng(seed)
    sizes = np.r_[np.full(N_SMALL, A_SMALL), np.full(N_LARGE, A_LARGE)]
    g = sizes.size
    share = sizes / sizes.sum()
    present = rng.random(g) < PRESENT
    absent = ~present
    m = int(absent.sum())

    p_exact = 1.0 - (1.0 - share) ** F_FALSE  # per-gene hit probability under the model
    pmf = poisson_binomial_pmf(p_exact[absent])
    expected_false = float((np.arange(pmf.size) * pmf).sum())

    fp_est = np.empty(r)
    fp_true = np.empty(r)
    fdr_est = np.empty(r)
    fdp = np.empty(r)
    p_bar_fn = np.empty(r)
    for i in range(r):
        hit = np.zeros(g, dtype=bool)
        hit[rng.choice(g, size=F_FALSE, replace=True, p=share)] = True
        detected = present | hit
        false_gene = hit & absent
        best_q = np.where(detected, 0.01, np.nan)
        n_true = int(present.sum())
        peptide_q = np.full(n_true + F_FALSE, 0.01)
        pep = np.r_[np.zeros(n_true), np.ones(F_FALSE)]  # F(t) = number of false peptides
        res = mayu_gene_fdr(best_q, peptide_q, pep, sizes, N_BINS)
        row = res.curve.iloc[0]
        assert row["F"] == pytest.approx(F_FALSE)
        fp_est[i] = row["fp"]
        fdr_est[i] = row["fdr"]
        fp_true[i] = false_gene.sum()
        fdp[i] = false_gene.sum() / detected.sum()
        # binned p_b = 1 - (1 - w_b)^F from the function's bins (one size per bin here)
        p_gene = np.empty(g)
        a = np.maximum(sizes, 1.0)
        for b in range(N_BINS):
            in_b = res.bin_of_gene == b
            w_b = a[in_b].mean() / a.sum()
            p_gene[in_b] = 1.0 - (1.0 - w_b) ** F_FALSE
        p_bar_fn[i] = p_gene[absent].mean()

    # the function's bins reproduce each gene's exact hit probability
    np.testing.assert_allclose(p_bar_fn, p_exact[absent].mean(), rtol=1e-9)

    # 1. binomial band on the pooled hit share of truly absent genes
    p_bar = expected_false / m
    pooled = fp_true.sum() / (r * m)
    band = k * np.sqrt(p_bar * (1 - p_bar) / (r * m))
    assert abs(pooled - p_bar) <= band, f"pooled hit share {pooled:.5f} vs {p_bar:.5f} ± {band:.5f}"

    # 2. FP estimate unbiased for the false-gene count
    d = fp_est - fp_true
    se = d.std(ddof=1) / np.sqrt(r)
    assert abs(d.mean()) <= k * se + 1e-9, f"mean FP bias {d.mean():.3f} > {k}*{se:.3f}"

    # 3. FDR criterion: mean FDP <= mean estimated FDR + k SE
    e = fdp - fdr_est
    se_e = e.std(ddof=1) / np.sqrt(r)
    assert fdp.mean() <= fdr_est.mean() + k * se_e, (
        f"mean FDP {fdp.mean():.4f} > mean FDR {fdr_est.mean():.4f} + {k}*{se_e:.4f}"
    )
