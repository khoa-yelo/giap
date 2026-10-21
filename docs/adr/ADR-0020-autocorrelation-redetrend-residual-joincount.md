# ADR-0020: Autocorrelation nulls recompute the detrend; residual join count; no CI for lag statistics

- Status: accepted-by-khoa (2026-10-20)
- Date: 2026-10-20
- Issue: M5-05 (also M5-01, M5-02); statistician attack check, Slurm 1860199, 1860281, 1860305, 1860318
  (scripts and outputs under local/scratch/stat-m5/sims/)
- Spec sections: §7.3.5(c), §8.1, §8.2, §8.3

## Context
Simulations on giapsim M4-07/08 truth (B = 499–1,000):
- Permuting the LOESS residuals without recomputing the trend (§7.3.5(c) "detrend then permute", M5-05's
  "no re-detrending per permutation") is anti-conservative. At 1,000 genes, R = 1,000: gene-label null
  pooled rejection over 20 lags 0.026 / 0.096 / 0.162 at α 0.01 / 0.05 / 0.10, BH-family FDR 0.136;
  TU-block null pooled 0.155 at 0.05, family 0.27. At 4,500 genes still 0.062 / 0.085 at 0.05. The smoother
  puts a negative autocorrelation (mean observed I₁ ≈ −0.042) into the observed residuals that permuted
  residuals lack (null mean ≈ −0.002).
- The raw join count rejects 0.082–0.092 at lag 1 (0.084–0.086 pooled) at α 0.05 under both nulls when
  short genes cluster in a prophage (`short_genes_in = prophage`, R = 500): detectability, not
  expression, clusters. §8.2 already requires detection tests to adjust for detectability.
- No interval for a lag-k Moran's I reaches 95 % coverage (test inversion 0.85–0.93 under TU null and
  AR 0.5/0.9 at lags 1/5/20, R = 150; a bootstrap CI of I_k is uncalibrated too).

## Decision
1. Per replicon: y = centred log2 abundance of genes with values. A wide trend is removed once:
   y_w = y − LOESS(y; span 0.30). The statistic is Moran's I of x = y_w − LOESS(y_w; span
   `context.autocorr.loess_span` = 0.10) along the ordinal of genes with values (circular padding as in
   M5-01).
2. Each null replicate permutes y_w (gene label: among genes with values; TU block: whole TUs, values and
   gaps together) and applies the same 10 % LOESS before computing I_k, so observed and null are the same
   function of the data.
3. The join count test uses the residual statistic Σ over valid lag-k pairs of (d_i − p_i)(d_j − p_j), with
   p = cross-fitted `p_detect`, permuted together with d under both nulls. `value` reports the raw BB_k;
   `effect` = residual statistic − its null mean (raw effect, adjusted test, as ADR-0018 decision 2).
4. Autocorrelation rows report `effect` = statistic − null mean and `ci_low` / `ci_high` = NA with flag
   `ci_unestimated` — a deviation from §8.2 for `autocorrelation.tsv` only. The null envelope (`null_lo`,
   `null_hi`) is reported for the V10 figure.
5. Calibration gates per (statistic, null) and α ∈ {0.01, 0.05, 0.10}: lag-1 level, lag-20 level, and the
   BH family's any-rejection rate under the complete null (FWER-style bound). The pooled 20-lag rate is
   recorded, not gated (lags are positively correlated; pooling hid a lag-1 failure).

## Evidence
- Decisions 1–2: TU null lag-1 rates 0.005 / 0.052 / 0.092, pooled 0.0115 / 0.0488 / 0.0973; gene-label
  null with a 2 log2 ori–ter cosine gradient pooled 0.010 / 0.051 / 0.104 (R = 400).
- Alternatives measured: permute raw values then detrend each replicate fails the gradient attack (TU-null
  lag-1 0.946); double detrend without the wide pre-trend is 0.072 at lag 1 (bound 0.064, R = 1,000).
- Decision 3: residual join count pooled at α 0.05: 0.051 (short-gene prophage), 0.050–0.052 (label and TU
  nulls, R = 1,000).

## Consequences
- Runtime: one (B × n)·(n × n) smoother product per null (≈ 6 × 10⁹ flops at 4,500 genes, inside M5-05's
  30 s budget). `giap.stats.smooth.loess_detrend` must accept positions 0…n−1 so the smoother matrix is
  shared by every replicate.
- p-values on real data drop for lag-k structure the detrend created. Interpretation note: at ρ 0.9 the
  10 % detrend creates a real negative I₂₀ (effect −0.195 at 1,000 genes); large-lag negative values are a
  span artefact under strong short-range correlation.
- New flag `ci_unestimated` (vocabulary through M5-09).
- New statistician attack file tests/calibration/attacks/test_attack_autocorr.py: (a) TU null + 2 log2
  cosine, (b) label null + gradient, (c) short-gene prophage under the TU null (join count gated), all gated
  at level; (d) the no-detrend version of (a) recorded as teeth (pooled 0.34).

## Alternatives
- Spec as written (detrend once, permute residuals): anti-conservative (above).
- Permute raw, detrend per replicate: fails under an ori–ter gradient with weak TU structure.
- Double detrend without the wide pre-trend: 0.072 at lag 1.
- Bootstrap CI for I_k: needs blocks ≫ k, biased toward 0 at joins, uncalibrated.

## Amendment 1 (accepted-by-khoa, 2026-10-20): Pearson residual join count
Supersedes decision 3's statistic. Statistician review of PR #116 (Slurm 1861073, 1861078; scripts and
outputs under local/scratch/review-116-stat-sims/) found Σ(d − p)(d − p) anti-conservative when low
`p_detect` clusters along the replicon: the permutation null cannot reproduce the clustered p(1 − p)
variance. Circular n = 1,000, 30-gene islands with p ~ U(0.01, 0.05) amid p ~ U(0.3, 0.9), α = 0.05:
rejection 0.059 / 0.074 / 0.086 / 0.105 at 10 / 20 / 30 / 50 % island share (both nulls).

- The join count test uses r_i = (d_i − p_i)/√(p_i(1 − p_i)) and the statistic Σ over valid lag-k pairs of
  r_i r_j, with p unclipped and permuted together with d under both nulls. Genes with p ∈ {0, 1} (or NaN)
  have no residual and are treated as missing (they join no pair).
- Evidence: level 0.045–0.056 in every scenario above and in the original attack set; about 10 % power
  cost; at very small p (0.001–0.005) α = 0.01 gives 0.019 (recorded, not gated).
- A floor on p in the denominator is rejected: it reintroduces the clustering failure.
- `value` still reports the raw BB_k; `effect` = Pearson statistic − its null mean.
- The engine attack file tests/calibration/attacks/test_attack_autocorr_engine.py carries the low-p-island
  arms (pass with the Pearson statistic) and a must-fail arm (raw residual, 50 % islands).
