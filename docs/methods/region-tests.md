# Region tests: O − E and Hodges–Lehmann shift under circular gene-index shifts

Implements `giap.context.regions` (GIAP_SPEC.md §7.3.5(b), §8.1, Appendix B.4; issue M4-17,
docs/specs/M4-17-region-tests.md, ADR-0018 and lead decisions 6–10). Reuses `giap.stats.permutation` (M4-03), `giap.stats.bootstrap`
and `giap.stats.robust` (M4-04) and `giap.stats.multitest.bh` (M4-01). Statistician sign-off: see
the end of this note.

## Index

A replicon's index is its **eligible** genes — not partial and with a `p_detect` — sorted by start
(the origin-crossing gene, stored unwrapped, comes last), numbered 0 … N − 1. Partial genes occur
only at linear contig ends and are never indexed, so every shifted set has the same number of
eligible genes as the region (lead decision 5). A region's index set I = the indices of its
eligible members on that replicon. Members on two replicons (an AMR gene set spanning chromosome
and plasmid) give one test per replicon, each on its own shifts, flagged `split_contigs`.

## Statistics

- **O − E** S(I) = Σ_{i∈I} (d_i − p_i); `observed` = Σ d_i, `expected` = Σ p_i, `effect` = S
  (and its CI) on the O − E scale. The **test statistic** is the studentised sum
  T(I) = S(I)/√Σ_{i∈I} p_i(1 − p_i) (ADR-0018 decision 2), recomputed at every shift, so a region
  whose p_detect sits nearer 0.5 than the replicon's is not anti-conservative (Σ p(1 − p) = 0 gives
  T = 0 for S = 0, ±∞ otherwise). With one p on the whole replicon T and S give the same p.
- **Residualised abundance** (ADR-0018 decision 1, lead decision 10): per replicon, OLS with an
  intercept of log2 abundance on log(length_aa), log1p(n_observable_peptides) and GRAVY (the §9.3
  genes.tsv covariates), fitted on the replicon's indexed genes that are detected with a finite
  abundance; a non-finite transformed covariate takes that column's median over the fitted genes
  (minimum-norm solution for rank-deficient designs, so constant covariates only centre a). The
  residuals ã are fitted once per replicon, so every shift sees the same values.
- **HL shift**: the median of all pairwise differences ã_i − ã_j between the region's detected genes
  with a finite residualised log2 abundance and the replicon's other such genes (`hodges_lehmann`; mean of the two
  middle values for an even count). NaN when either side is empty. Detected genes without abundance
  (policy `exclude`) drop out on both sides. Comparing within the replicon only keeps a high-copy
  plasmid's dosage out of the shift (attack 4). Effect, CI and p are on the residual scale: log2
  units relative to the covariate expectation, labelled "detectability-adjusted" in the report.

## Null: circular gene-index shifts

I ↦ (I + k) mod N for k = 1 … N − 1 (`circular_shifts`); on linear contigs only shifts that do not
wrap. Exhaustive while N ≤ `context.exhaustive_shift_max_genes`, otherwise `stats.permutations`
shifts drawn without replacement (`choose_shifts`). Both statistics are recomputed at each shift;
HL shifts with an empty side are dropped, so HL `n_perm` counts valid shifts.

What the null preserves: region size, the spacing of its members (a block stays a block, a scattered
AMR set keeps its gaps), and the replicon's autocorrelation of d − p and abundance; the region is
moved, the genome is not rearranged. Under H₀ "the region's position is exchangeable with every
other position on the replicon", the observed set is one member of the shift orbit, so the +1-rule
p is exact. The circular null is not perfect on simulated data (M4-08's AR(1) restarts per replicon
and does not wrap at the origin, one uncorrelated join) — negligible for level.

**p** = `perm_p(S_obs, null, "two-sided")`: (1 + #{|S_k − c| ≥ |S_obs − c|})/(1 + K), c the pooled
mean of the observed and null statistics (exact under exchangeability; ties within 1e-12 count).
K = 0 (a linear contig with no non-wrapping shift) → p NA. p is discrete with floor 1/(1 + K):
rows with K < 19 (no p can reach 0.05) are flagged `few_shifts`.

### Exact HL at every shift

Recomputing HL by enumeration (m·(M − m) differences at each of up to N − 1 shifts) is too slow at
4,500 genes. `hl_null` sorts the replicon's M observed values once; for each shift the k-th smallest
inside-minus-outside difference is found by bisection on its value, counting pairs x_i − z_j ≤ t by
binary search in the sorted values and subtracting the inside-inside pairs, for all shifts at once;
once ≤ 64 candidates remain they are enumerated. Results equal `hodges_lehmann` to ~1e-14 (the
property suite compares with a naive per-shift loop).

## CIs: TU-block bootstrap

`tu_block_bootstrap`, B = `stats.bootstrap`, 95 % percentile. O − E resamples the region's TUs,
counting only region members within each TU (a region inside one TU has a degenerate CI). HL
resamples inside and outside TUs as two strata (on the residualised abundances; the fit is not
repeated per resample).

## Multiple testing and flags

BH within each (region_type, stat) family (lead decision 5, ADR-0018 decision 3; detection and
abundance answer different questions; at least as strict as "within type" for each statistic). Flags: `small_region` (< 3 eligible genes; tests still run when K > 0; with 0 eligible
genes all statistics are NA and `n_perm` = 0), `few_shifts`, `split_contigs`.

## Calibration evidence

Suite `regions` (tests/calibration/test_region_tests.py, attacks in
tests/calibration/attacks/test_region_tests_attacks.py): M4-07 layout + 60-gene plasmid with its
own AMR set, M4-08 AR(1) ρ = 0.5, no region effects, detection from covariates only.

- Fast run (Slurm 1849058, R = 100 + confirmation block): C2 FDR, A4 (20x plasmid) and every C1
  check pass except `crispr`/`hl_shift` at α = 0.01 (4/100 then 3/100, bound 2.99/100). An
  independent 600-replicate run of the same scenario (Slurm 1849100) gives rejection rates
  0.012 / 0.040 / 0.072 at α = 0.01 / 0.05 / 0.10 for that check (all 8 (type, stat) checks within
  0.003-0.017 / 0.032-0.077 / 0.072-0.127, also without origin crossing and with ρ = 0):
  the gate failure looks like seed chance, for the statistician to rule on.
- A2 (short prophage genes) gates on M4-13 cross-fitted `p_detect` and waits for PR #83.
- ADR-0018 implementation (Slurm 1849426 fast, 1849427 full): A2-mid O − E (studentised) passes at
  α = 0.05 / 0.10 (0.05 / 0.07, R = 100; α = 0.01 unconfirmed 0.03 / 0.02). **A2-HL still fails**:
  linear residualisation (lead decision 10) gives 0.07 / 0.04 (fast, α = 0.01), 0.17 / 0.11 (0.05),
  0.22 / 0.15 (0.10), and 0.052 / 0.144 / 0.221 at α = 0.01 / 0.05 / 0.10 in full mode (R = 1,000 per
  block). Side runs on the same scenario (600 replicates at 1,000 genes / 200 at 4,500): raw HL
  0.107 / 0.243 / 0.328 and 0.105 / 0.25 / 0.38, linear 0.033 / 0.112 / 0.185 and
  0.055 / 0.14 / 0.195, adding (log length)² 0.025 / 0.09 / 0.14 and 0.05 / 0.12 / 0.18.
  Residualising on covariates halves the bias but does not remove detection selection; open for
  the statistician and the lead.
- C3 power (recorded, not gated; BGC O − E, α = 0.05, R = 100): +10 % 0.11, +20 % 0.21,
  +30 % 0.46, +50 % 0.88 (`local/logs/calibration/regions_power.jsonl`).

## Performance

4,500-gene chromosome, 30 regions, exhaustive shifts, B = 1,000: ~35 s on one core
(benchmarks/test_region_tests_bench.py; budget 60 s); the TU-block bootstrap of HL dominates.

## Limitations (statistician)

- **Exactness needs exchangeable positions.** The shift p is an exact randomization p only when the
  region's position is independent of the per-gene data. Two violations matter:
  (i) *detectability confounded with membership* (§8.5 attack 2): the unstandardised O − E sum has
  variance Σ p(1 − p), so a region whose p_detect sits nearer 0.5 than the replicon's is
  anti-conservative (synthetic, R = 4,000: replicon p ~ Beta(4, 1.5), 40-gene region p ~ Beta(2, 2)
  gives 0.016 / 0.074 / 0.131 at α = 0.01 / 0.05 / 0.10; replicon Beta(0.3, 0.3), region p = 0.5
  gives 0.074 / 0.212 / 0.303). The studentised shift statistic Σ(d − p)/√Σ p(1 − p), recomputed at
  each shift, restores level (0.010 / 0.048 / 0.103 and 0.010 / 0.046 / 0.097); adopted (ADR-0018). The giapsim attack-2
  layout (short prophage genes, low p) is the conservative direction (0.004 / 0.019 / 0.047 with
  p_detect_true) and cannot detect this.
  (ii) *detection-selection bias of HL*: HL uses detected genes only. When detection depends on
  abundance (giapsim `detection = "abundance"`, β_abundance 0.5–1) and the region's genes are short,
  only its abundant genes are detected: prophage HL rejects 0.089 / 0.21–0.25 / 0.30–0.35 (R = 2,000,
  mean HL +0.7 to +0.8 log2; control layout without short genes 0.004–0.009 / 0.036–0.048 /
  0.080–0.100). Residualising log2 abundance on detectability covariates among the replicon's
  detected genes before HL (prototype: quadratic in log length) gives 0.017 / 0.067 / 0.103 (R = 300);
  adopted as the linear fit above (ADR-0018, lead decision 10). Residualising removes the part of
  a real region shift that the covariates explain (e.g. a prophage that is both short and silent).
  **Correction (review round 2):** the R = 300 prototype figure was noise; residualisation removes
  the *location* bias but not the failure. The remaining failure is a third violation:
  (iii) *detected-count heterogeneity of HL*: only ~6–7 of the prophage's 40 genes are detected
  (short genes), against ~22 at a shifted position, so the observed HL has a wider sampling
  distribution than the shift null (same mechanism as (i) for O − E). It fails even when detection
  does not depend on abundance (A2-HL-MAR below) and even with the exact oracle E[a | x, detected].
- **Replication dosage.** Real replicons have an ori → ter gene-dosage gradient (up to ~1–2 log2 in
  fast growth) that giapsim does not model; circular shifts move the region along it, so HL partly
  measures position. Report HL as "relative to random positions on this replicon", not as a
  region-specific regulatory effect.
- **Plasmid copy number** cancels exactly (within-replicon comparisons): results are bitwise identical
  for `plasmid_copy_log2` 0 and log2 20 on the same seed (3 seeds checked), so attack 4 holds by
  construction and is best tested as an invariance, not a calibration stream.

## Statistician sign-off

Statistician (M4-17 review, PR #89, 2026-10-15): **changes requested** (not signed).

- C1 `crispr`/`hl_shift` fast confirmed fail: **chance**. 20,000 null replicates of the fast scenario
  (Slurm 1849136; same layout, ρ = 0.5, p = p_detect_true) give 0.0083 / 0.048 / 0.096 for that check
  and 0.0083–0.0103 / 0.045–0.052 / 0.096–0.102 for every chromosome check; 10,000 without origin
  crossing 0.0080–0.0112 / 0.046–0.055 / 0.087–0.108. Per stream the two-block false alarm is 0.0063
  at α = 0.01, R = 100 (fail iff ≥ 3 rejections in both blocks). Seeds are fixed, so the fail repeats
  until the claim is restructured; no seed change and no tolerance change.
- False-alarm cost as drafted: ≈ 0.13 per fast run (21 level streams: C1 10, A4 10, A2 1) and ≈ 0.09
  per nightly run, against ≈ 0.04 for the whole nightly (ADR-0015). Restructured: C1 one pooled level
  stream per statistic over the four chromosome types (`check_level_pooled`, m = 4; rejection
  indicators of different region types correlate ≤ 0.06 over 20,000 replicates), plasmid AMR at
  α ∈ {0.05, 0.10} only (K = 59: p ≤ 0.01 is unattainable), O − E at one ρ in full mode (O − E does not
  depend on ρ under covariate-only detection), A4 as a deterministic invariance test: ≈ 0.018 fast and
  ≈ 0.017 nightly.
- Blocking: limitations (i) and (ii) are attack-2 failures (§8.5: a failing attack blocks the merge);
  the A2 gate errors until M4-13 (#83) is merged (not marked xfail as lead decision 4 allows).

### Review round 2 (statistician, 2026-10-15): A2-HL root cause and recommendation

Scripts and outputs: `local/scratch/review-m4-17/v2/` (exp2–exp5.py, an.py; Slurm 1849806, 1849807,
1849809, 1849826, 1849914; attack test run 1849915). All rates are prophage `hl_shift` rejection
shares at α = 0.01 / 0.05 / 0.10, exhaustive circular shifts, ρ = 0.5, β_abundance = 0.5 unless stated.

1. **Decision 10 as implemented fixes the location bias.** Mean residualised HL over replicates:
   −0.07 to −0.09 log2 (raw +0.74). Fit-set, pooling and imputation variants do not matter:
   pooled chr1 + p1 fit equals per-replicon to ±0.003; adding (log length)² or a cubic in a fitted
   logistic detectability index (or in the true index) changes rates by ≤ 0.04. The oracle
   residual a − E[a | true covariate index, detected] (exact under the simulator, Gauss–Hermite)
   still rejects 0.053 / 0.131 / 0.199. So no covariate adjustment can pass this attack.
2. **Cause: detected-count heterogeneity (limitation (iii)).** A2-HL-MAR (short prophage,
   detection from covariates only, so no selection on abundance): plain HL 0.058 / 0.147 / 0.225
   (1,000 genes, R = 2,000) and 0.062 / 0.151 / 0.229 (4,500 genes, R = 1,000). Residualised:
   0.051 / 0.131 / 0.207 and 0.063 / 0.152 / 0.227. β_abundance = 1 detects more prophage genes
   (10 versus 7) and fails less (0.012 / 0.074 / 0.139).
3. **Fix (recommendation (c)): a TU-cluster-studentised HL test statistic**, which parallels
   ADR-0018 decision 2. Keep decision 10's residualisation, gene-index circular shifts, and HL as
   the reported effect with its CI. The p-value comes from T = HL / √V, recomputed at every shift.
   V = σ_u²(Σ_t m_t²/m² + Σ_t m′_t²/m′²) + σ_e²(1/m + 1/m′), where:
   - m_t and m′_t are the detected genes with finite ã in TU t, inside and outside the set;
   - m and m′ are their totals;
   - σ_u² (between-TU) and σ_e² (within-TU) are one-way ANOVA moment estimates on the replicon's
     residualised detected genes, computed once per replicon, with σ_u² truncated at 0.
   Level (prophage HL, 1,000 genes / 4,500 genes):
   - A2-HL β 0.5: 0.005 / 0.019 / 0.046 (R = 4,000) and 0.003 / 0.021 / 0.047 (R = 2,000).
   - A2-HL β 1.0: 0.001 / 0.014 / 0.045 (1,000 genes, R = 1,000).
   - A2-HL-MAR: 0.003 / 0.016 / 0.038 and 0.001 / 0.014 / 0.037.
   - A2-mid: 0.004 / 0.025 / 0.052 (1,000 genes, R = 2,000).
   - C1 pooled over the 4 chromosome types: 0.010 / 0.045 / 0.092 (R = 4,000) and
     0.009 / 0.048 / 0.099 (R = 1,000).
   - C1 plasmid AMR (p1, R = 2,000): 0.000 / 0.043 / 0.095.
   - Under the attack it is conservative because AR(1) between TUs is not modelled.
   Power without confounding (prophage +1 log2, R = 500) at α 0.01 / 0.05 / 0.10:
   - TU-studentised: 0.058 / 0.166 / 0.252.
   - plain HL: 0.048 / 0.156 / 0.246.
   With confounding (short prophage + abundance detection, +1 log2), TU-studentised gives
   0.010 / 0.070 / 0.134 against plain HL's invalid 0.056 / 0.182 / 0.258.

   Alternatives tested, all at level:
   - shifting the region's detected genes within the detected-only gene sequence:
     0.004 / 0.019 / 0.044. This changes the spec's null.
   - scaling by √(m m′/M): 0.002 / 0.010 / 0.028.
   - Wilcoxon z: 0.001 / 0.009 / 0.026.

   Both scaling and Wilcoxon z are over-conservative because they ignore TU clustering. Decision 10
   itself needs no amendment. Option (b) is not needed: the attack is claimable once the statistic
   is studentised.
4. **O − E re-check.** The studentised O − E statistic passes:
   - A2-mid (1,000 genes, p = truth, R = 6,000): 0.009 / 0.050 / 0.103, which confirms α = 0.01.
   - A2-HL data: 0.012 / 0.054 / 0.101 (R = 4,000).
   - A2-MAR data: 0.008 / 0.051 / 0.104.
5. **C1 pooled family.** BH within (type, stat) and pooled level per statistic over the 4 chromosome
   types:
   - O − E: 0.010 / 0.048 / 0.101 (1,000 genes, R = 4,000) and 0.009 / 0.047 / 0.099 (4,500 genes).
   - HL (decision 10): 0.008 / 0.045 / 0.094 and 0.010 / 0.049 / 0.099.

   All are within the bounds (R = 4,000: 0.0131 / 0.0569 / 0.1095). A single unpooled cell reached
   0.072 at α 0.05 (4,500 genes, AMR O − E, R = 1,000; bound 0.064), but an independent seed for the
   same cell gives 0.052. This is chance, so the pooled restructure stands.
6. **New attack null** `test_a2_short_prophage_hl_at_level_without_abundance_selection`
   (tests/calibration/attacks). It runs the A2-low data with HL and must stay at level. The current
   code is a confirmed fail in fast mode (Slurm 1849915), as is A2-HL.

**Sign-off state: changes requested (blocking).** Merge needs the studentised HL statistic
(item 3), with A2-HL and A2-HL-MAR passing. This changes the §7.3.5(b) HL test, so it needs an
ADR-0018 amendment (decision 4) accepted by Khoa.
