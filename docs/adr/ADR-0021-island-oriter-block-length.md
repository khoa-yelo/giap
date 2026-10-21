# ADR-0021: Block length for the island scan and the ori–ter CI

- Status: accepted-by-khoa (2026-10-20); **decision 2 amended and accepted-by-khoa 2026-10-20 evening**
  (option B: the `uncalibrated` flag is dropped). See "Amendment 1".
- Date: 2026-10-20
- Issues: M5-02, M5-06, M5-07, M5-08; statistician attack check, Slurm 1860199, 1860281, 1860305, 1860318
  (local/scratch/stat-m5/sims/)
- Spec sections: §7.3.5(d) step 5, §7.3.5(e), docs/methods/bootstrap.md "Block length"

## Context
With the spec's ℓ̂ (the first lag at which TU-level Moran's I is not significant, bounded to [1, 20]),
island FWER (share of replicates with ≥ 1 island) on a 4,500-gene chromosome + 60-gene plasmid
(B = 500, R = 400) is 0.028 / 0.205 / 0.26 / 0.27 at ρ 0 / 0.5 / 0.9 / 0.95, and 0.15 at the fast-gate
setting (1,000 genes, ρ 0.5). ℓ̂ marks where correlation stops being detectable, not the block length a
window-sum maximum needs; the maximum is driven by extreme 2-TU windows of single-gene TUs. Fixed blocks
≥ 10 TUs copy the observed extreme window into ≈ 63 % of replicates, so FWER and power go to 0. The same
ℓ̂ under-covers the ori–ter slope CI (0.88–0.93 at ρ 0.5). Rejected alternatives (Gaussian null with the
empirical ACF, sieve AR bootstrap 0.32 / 0.80, ACF-studentised scan 0.31 / 0.70 with V ≤ 0 on the plasmid,
TU-size-standardised scores 0.18 at ρ 0) all failed.

## Decision
1. Block length ℓ' = min(20, 2ℓ̂ + 1, T) per replicon (T = number of scored TUs), with ℓ̂ as in the spec,
   for both the island circular block bootstrap and the ori–ter slope CI. An integer
   `context.islands.block_len` overrides it unchanged.
2. Validity domain: islands on replicons with fewer than 1,000 scored TUs carry the §9.4 flag
   `uncalibrated` (statistician's option b). The M5-07 gate is claimed for replicons with ≥ 1,000 scored TUs;
   smaller replicons are recorded, not gated.

## Evidence
- Island FWER with ℓ': 0.003 / 0.040 / 0.023 / 0.023 at ρ 0 / 0.5 / 0.9 / 0.95 (R = 300; bound 0.075);
  0.013 at 1,000 genes, ρ 0.5.
- Attacks at level (R = 300): ori–ter gradient 0.04, large TUs 0.013, abundance-driven detection 0.04,
  short-gene prophage 0.04, draft assembly (500 + 133 × 30-gene linear contigs) 0.003 (spec rule 0.09).
- Ori–ter pooled-slope CI coverage 0.96 at ρ 0.5 (spec rule ≈ 0.905). **Corrected 2026-10-21:** that
  0.96 was an R = 50 fluctuation. Re-measured with the merged implementation at R = 300 (Slurm 1863166)
  the same cell gives **0.943** (SE 0.013); ρ 0 gives 0.925, ρ 0.9 0.900, no gradient 0.917, 1,000 genes
  0.905, and the arms 0.875–0.945. The estimator is unbiased (−2.4239 against the closed form −2.4317);
  the interval is 1–5 points narrow throughout.
- Power (4,500 genes, ρ 0.5, R = 100): 20 TUs × 2 SD 0.18, 20 × 1.5 0.06, 10 × 2 0.04, 5 × 2 0.00 — low
  because ≈ 45 % of genes are censored independently of abundance in giapsim defaults; power stays
  recorded, not gated, and is also run with `detection = "abundance"`.
- Open limitation: 1,000 genes at ρ 0.9 gives FWER 0.13 (spec rule 0.55), because ℓ̂ ≈ 1.3 there (Moran's
  I cannot detect the diluted correlation on ≈ 555 TUs); hence decision 2.

## Consequences
- M5-07's planning target (power ≥ 0.8) will not be met under giapsim defaults; the report says so.
- The effect CI stays the spec's TU-bootstrap CI; the methods note says it is conditional on the selected
  window (winner's curse) and ignores within-window autocorrelation; M5-07 records its coverage.
- The `sim1000_island` golden (M5-11) is checked against ℓ' before it is frozen.

## Alternatives
- (a) Accept the small-replicon limitation and gate only the full-gate sizes, without a flag.
- (c) ℓ = 20 on replicons with < 1,000 scored TUs: valid but no power there.
- The spec's ℓ̂ unchanged: anti-conservative (above).

## Amendment 1 (accepted-by-khoa, 2026-10-20): the validity domain is a dependence regime, not a size

Supersedes decision 2. Evidence: M5-06 statistician re-measurement, R = 1,000 per cell, B = 199, through
the M5-06 oracle on the merged engines (Slurm 1861482, 1861494, 1861551, 1861552; scripts and outputs
under local/scratch/m5-06-stat/). The original numbers reproduce (4,500 genes: 0.037 at ρ 0.5 against the
ADR's 0.040; 0.024 at ρ 0.9 against 0.023); the **inference** from them was wrong.

2. **Validity domain.** The island FWER is claimed for the dependence regime the calibration covers, not
   for a replicon size. It is claimed for giapsim `ar1_background` at ρ ≤ 0.5 at every size measured
   (557 / 1,001 / 2,501 scored TUs: 0.020 / 0.031 / 0.037), for the 4,500-gene chromosome at ρ = 0.9
   (0.024, where ℓ̂ ≈ 8 makes ℓ' = 17 span the 20-TU windows), and for the §8.5 attacks that stay in that
   regime. It is **not** claimed where strong TU-level autocorrelation is resolved by Moran's I only to a
   short lag (ℓ' ≤ 7): giapsim ρ = 0.9 gives 0.113 at 557 TUs and 0.095 at 1,001 TUs, and 1,000
   single-gene TUs give 0.068 / 0.110 / 0.208 / 0.327 at ρ = 0.7 / 0.8 / 0.9 / 0.95 (bound 0.0638).
   The scored-TU count decides nothing on its own — 400 TUs at ρ 0.5 give 0.005, 1,000 TUs at ρ 0.9 give
   0.095–0.208, 2,500 TUs at ρ 0.9 give 0.024 — it enters only through the power of the ℓ̂ test, together
   with ρ, the TU size distribution and the censored share (the same ρ = 0.9 replicon gives ℓ̂ = 9 and
   0.007 uncensored, ℓ̂ = 2 and 0.208 at 55 % detected).

   Two mechanisms, both measured at 1,000 TUs, ρ = 0.9, with the block forced: autocovariance beyond the
   block (uncensored, Î₁ = 0.78: 0.576 at block 5 against 0.016 at block 19) and local-scale heterogeneity
   at matched variance (45 % censored, variance ratio ≈ 1.0 yet 0.158 at block 5; block SD q90/q10 1.73 at
   ρ 0.9 against 1.24 at ρ 0). A long block fixes both and costs all the power, which decision 1 rejected.

2a. **No flag.** `uncalibrated` is **dropped** (Khoa, option B): no per-replicon statistic certifies the
   domain, so a flag would give false reassurance. ℓ̂, ℓ', the lag-1 Moran's I of the detrended TU scores,
   a Bartlett variance ratio of the window sums and a local-scale ratio all fail as keep/flag rules, and
   selecting on them is **adverse** — giapsim ρ = 0.9 at 557 TUs is 0.113 overall but **0.220** among the
   replicates that look calmest (ℓ̂ ≤ 1, flat local scale), because a small ℓ̂ is exactly what makes the
   block too short; the same rule flags 94 % of the perfectly calibrated 2,501-TU ρ = 0.5 cell.
   Instead: docs/methods/islands.md states the domain above, and the §9.4 run report carries ℓ̂, ℓ'
   (`block_len`) and Î₁ per replicon so a reader can see where their data sit.

### Measured and rejected
Local standardisation (dividing the detrended scores by a LOESS fit of |residual| before the scan and the
bootstrap) reduces the inflation (ρ 0.8: 0.106 → 0.032 at a 10 % span; ρ 0.9: 0.188 → 0.084, still a
breach, because the autocovariance mechanism remains) but absorbs the signal: power on a planted
20-TU × +2 SD island at ρ 0.5 falls from 0.103 to 0.003.

### Consequences
- M5-07 keeps its gate cells (ρ ∈ {0, 0.5, 0.9} at 4,500 genes) and adds **recorded, not gated** cells at
  1,000 and 1,800 genes × ρ = 0.9.
- docs/methods/islands.md gains the domain statement and the two mechanisms.
- M5-06 and M5-09 drop `uncalibrated`; the §9.4 run report gains the per-replicon ℓ̂, ℓ' and Î₁.
- Backlog: a scale-aware or variance-corrected island null.
