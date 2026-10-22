# ADR-0033: Elastic-net ΔR² is the primary block attribution; LightGBM early stopping in hardening
- Status: accepted-by-khoa (2026-10-21, in session: "Take your recommendations")
- Date: 2026-10-21
- Spec sections: §7.8, §5.2 (LightGBM parameters), §9 `tables/determinants_model.tsv`
- Needs Khoa: yes (what the report presents as the attribution). Answered: option (c), plus (a) in hardening.

## Context
With the §5.2 LightGBM parameters at n = 400–600 genes:
- CV R² is −0.27 under the null.
- In a planted world, ΔR²(codon) is 0.97 against a true 0.80.
The elastic net is unbiased in both worlds, and the M8-09 permutation p stays valid for both models.

## Decision
- **Pass 1:** the elastic-net ΔR² is the primary block attribution in the report and V15.
- LightGBM ΔR² and SHAP shares are reported as secondary, nonlinear checks. Their rows are not gated against truth.
- `determinants_delta_r2_max_abs_bias` (ADR-0031) gates the elastic net.
- **Hardening:** LightGBM early stopping inside the inner CV, with the number of trees chosen per fold. The LightGBM bias is then re-measured, and if it passes the bias gate the LightGBM ΔR² becomes primary again.
- The §5.2 defaults stay unchanged until that evidence exists.

## Evidence (before → after: metrics, tests, simulations)
Statistician probe p3_cv.py (Slurm 1865752, 64 replicates) and the planted calibration world.

## Consequences
- M8-10 and M9 present the elastic net first.
- SHAP (M8-09) stays a LightGBM output, with a caveat.
