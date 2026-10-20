# ADR-0021: Block length for the island scan and the ori–ter CI

- Status: accepted-by-khoa (2026-10-20; small replicons flagged `uncalibrated`, option b)
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
- Ori–ter pooled-slope CI coverage 0.96 at ρ 0.5 (spec rule ≈ 0.905).
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
