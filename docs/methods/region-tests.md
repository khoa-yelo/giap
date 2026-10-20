# Region tests: O − E and Hodges–Lehmann shift under circular gene-index shifts

Implements `giap.context.regions` (GIAP_SPEC.md §7.3.5(b), §8.1, Appendix B.4; issue M4-17,
docs/specs/M4-17-region-tests.md, ADR-0018 and lead decisions 6–10) and the `regions_stats.tsv`
rows of `giap.context.regions_table` (issue M4-18, docs/specs/M4-18-regions-stats-table.md,
sections "Replicon-level effects" and "Descriptive rows"). Reuses `giap.stats.permutation` (M4-03), `giap.stats.bootstrap`
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
- **HL test statistic** (ADR-0018 amendments 1–2, lead decisions 11–16): the shift p comes from
  the studentised T = HL/√V, recomputed at every shift (`studentised_hl_null`). V is label-free
  (amendment 2): V = Var(mean_in − mean_out) = S_in,in/m² + S_out,out/m′² − 2 S_in,out/(m m′), S_A,B
  the sum of γ(d(i, j)) over detected genes with finite ã in A × B, m, m′ their counts inside /
  outside the set (`hl_variance`). γ(h), h = 0 … 30 (`LAG_MAX`), is the replicon's autocovariance
  of ã at gene-index lag h (centred on the mean of ã, each lag divided by its number of detected
  pairs), estimated once per replicon and set to 0 from the first h ≥ 1 with γ(h) ≤ 0
  (`lag_autocovariance`). Lags and d are circular on circular replicons and plain |i − j| (pairs
  inside the sequence) on linear ones (decision 16). Implementation: c_i = Σ_j ok_j γ(d(i, j)) once
  per replicon (≤ 61 offsets), then per shift S_in,all = Σ_{i∈set} c_i, S_in,in by pairwise
  distances inside the set, S_in,out = S_in,all − S_in,in, S_out,out = S_all − 2 S_in,all + S_in,in.
  Guard (decision 12): when γ(0) ≤ (10⁻⁷ · max |raw log2 abundance|)², i.e. (numerically) constant
  residuals, V counts as 0 and T = 0 at every shift. Otherwise, a finite V ≤ 0 at the observed or
  any null shift (the truncated γ need not be positive definite) makes the region's HL p NA with
  flag `hl_variance_nonpositive` (fail closed, decision 16). TU labels do not enter V; a missing
  `tu_id` on a tested replicon raises ValueError (decision 13; C6 makes it non-null). The reported
  effect and CI stay HL. Why: a region's detected-gene count can differ strongly from that at
  shifted positions (≈ 7 vs ≈ 22 of 40 short prophage genes), so raw HL has the wrong spread under
  the shift null (limitation (iii)); label-free V is immune to TU-label errors (attack A-TU).

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
repeated per resample). A missing `tu_id` raises ValueError before any test (decision 13), so
every block is a real TU.

## Multiple testing and flags

BH within each (region_type, stat) family (lead decision 5, ADR-0018 decision 3; detection and
abundance answer different questions; at least as strict as "within type" for each statistic). Flags: `small_region` (< 3 eligible genes; tests still run when K > 0; with 0 eligible
genes all statistics are NA and `n_perm` = 0), `few_shifts`, `split_contigs`,
`hl_variance_nonpositive` (HL row only; p NA, see above).

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
  Residualising on covariates halves the bias but does not remove detection selection; resolved
  by ADR-0018 amendment 1 (next item).
- Studentised HL (dd663c0; suites run on that tree, Slurm 1850200 fast / 1850201 full; all pass).
  Full mode (R = 1,000, α = 0.01 / 0.05 / 0.10):
  - A2-HL 0.002 / 0.024 / 0.064, A2-HL-MAR 0.003 / 0.016 / 0.042;
  - A2-mid O − E 0.005 / 0.034 / 0.073, A2-low O − E (M4-13 p) 0.005 / 0.044 / 0.076;
  - C1 pooled (m = 4, ρ = 0.5): HL 0.007 / 0.049 / 0.096, O − E 0.009 / 0.050 / 0.103; pooled HL
    at ρ = 0 and 0.9 also pass, as does plasmid AMR at α 0.05 / 0.10.
  Fast: everything passes; A2-mid O − E at α = 0.01 is unconfirmed (0.03 / 0.02).
- Label-free lag V (ADR-0018 amendment 2, 04fc8f7; Slurm 1850777 fast, 1850778 full; all pass).
  Full mode (R = 1,000 per block, α = 0.01 / 0.05 / 0.10):
  - A-TU (region TUs fragmented, pooled m = 4) 0.0115 / 0.0525 / 0.0985 (failed with TU-label V);
  - A2-HL 0.009 / 0.050 / 0.098, A2-HL-MAR 0.012 / 0.041 / 0.090;
  - A2-mid O − E 0.005 / 0.034 / 0.073, A2-low O − E 0.005 / 0.044 / 0.076;
  - C1 pooled HL at ρ = 0 / 0.5 / 0.9: 0.0105 / 0.0505 / 0.102, 0.0085 / 0.047 / 0.096,
    0.009 / 0.047 / 0.096; O − E (ρ = 0.5) 0.009 / 0.050 / 0.103; plasmid AMR passes.
  Fast: everything passes (C1 pooled HL at α = 0.01 and A2-mid O − E at α = 0.01 unconfirmed:
  0.02 / 0.0175 and 0.03 / 0.02). `hl_variance_nonpositive` occurred in 3 of 59,000 HL rows in
  full mode (all in the ρ = 0.9 null) and 0 of 7,300 in fast mode.
- A-LIN (statistician, round 4; `test_alin_linear_replicon_hl_at_level`): no origin crossing, every
  replicon declared linear, ρ = 0.9, chromosome HL pooled over prophage / BGC / CRISPR (m = 3).
  Fast (Slurm 1851283, R = 100) 0.0133 / 0.070 / 0.097, full (1851308, R = 1,000)
  0.0083 / 0.0453 / 0.0993: pass.
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
- **TU labels enter only the CIs (amendment 2).** V is built from gene-index lags, so TU-label
  errors cannot move the HL p (attack A-TU at level). The TU-block bootstrap still uses the labels,
  so labels that fragment a region's operons narrow that region's CI; the p is unaffected.
- **V ≤ 0 (decision 16).** The truncated γ is not always positive definite. In 60,000+ chromosome
  HL rows (ρ up to 0.95, ori–ter cosine up to 2 log2) it never happened; on the 60-gene plasmid at
  ρ ≥ 0.9 it hits 0.15–0.35 % of rows (lags up to 30 cover half the replicon). Those rows lose
  their HL p, but level is unaffected: under exhaustive shifts the flag is a function of the whole
  shift orbit, so the same rows would be NA whichever orbit member was observed.

## Replicon-level effects (plasmids, M4-18)

Plasmids are not tested (§7.3.5(b), §8.5 attack 4). A plasmid region is the whole replicon, so no
circular shift moves it, and its contrast with the chromosome mixes gene dosage (copy number) with
regulation: nothing in proteomics data separates the two. Each `plasmid_*` region gets two rows with
`p_value = q_value = NA`, `null_model = none`, `n_perm = 0`, `method = replicon_contrast`,
`label = dosage_plus_regulation`:

- `o_minus_e`: the plasmid's own O − E over its eligible genes (E = Σ p_detect, which already
  carries the detectability model; lead decision 1). No chromosome-rate row.
- `hl_shift`: HL of the plasmid's **raw** log2 abundance against the reference replicons (lead
  decision 2). Not M4-17's per-replicon residuals: those centre every replicon and give HL ≈ 0
  whatever the dosage. Only detected, eligible genes with an observed value enter; **imputed genes
  are excluded** on both sides, because the genome-wide median fill (docs/methods/abundance.md)
  places imputed plasmid genes at the chromosome median and pulls the contrast towards 0.
- Reference (lead decision 3): the `kind = chromosome` replicons; without one (draft genomes) every
  replicon that is neither `kind = plasmid` nor a plasmid region's contig, flag
  `no_reference_replicon` (HL NA with the same flag if that set has no usable gene).
- CIs: 95 % TU-block bootstrap; for HL with strata plasmid / reference, so every replicate has both
  sides. Each plasmid draws from its own child stream keyed by its `region_id` (lead decision 7), so
  adding a plasmid changes no other row.
- Fewer than `SMALL_REGION` (3) eligible genes: every statistic NA, flag `small_region`.

**Depth offset.** With `--replicon-depth` (§5.3), an extra `hl_shift` row (`label =
dosage_adjusted`, flag `depth_adjusted`) reports HL − log2(plasmid depth / reference depth), the
reference depth being the length-weighted mean over the same reference replicons. **Depth is
treated as known**: the CI is the HL CI shifted by the same constant, so it omits the uncertainty
of the depth ratio (deferred, docs/backlog.md). The row is omitted, with a warning from the step,
when the plasmid or a reference replicon has no depth, or when HL is NA.

**Bias at the null (recorded, not gated).** The HL compares detected genes only. When detection
depends on abundance, the low-abundance genes of each replicon drop out, and the two replicons
lose different shares, so the raw contrast carries a selection bias at the true dosage; the lead
decision cites ≈ +0.3 log2 in simulation for abundance-driven detection. Under covariate-only
detection (`tests/calibration/attacks/test_regions_stats_attacks.py`, AR(1) ρ = 0.5, 1,000
chromosome genes, 60-gene plasmid at log2 20 copies, R = 100, B = 200; Slurm 1860163) the mean
HL − log2(20) is +0.013 log2 (sd 0.79), but the 95 % TU-block CI covers log2(20) in only 78 % of
replicates: AR(1) dependence between neighbouring TUs and the replicon-level abundance level of
the simulated plasmid are not captured by resampling independent TUs. Read the CI as a lower bound
on the uncertainty. The result is an effect estimate, never a test.

## Descriptive rows (prophage subsets, M4-18)

Per prophage and §7.1.7 function subset present in C8 (`integrase`, `repressor`, `structural`,
`lysis`, `cargo`): `subset` = the name, `stat = o_minus_e`, `observed` = detected eligible members,
`expected` = Σ p_detect, `o_minus_e` (= `effect`), `o_over_e` (NA when E = 0). `method = label =
descriptive`, flag `descriptive`, CI, p and q NA, `null_model = none`. Subsets are small (often one
to a few genes) and nested in the tested prophage; they are not tested in version 1. Fewer than
3 eligible genes add `small_region`; 0 eligible genes give NA statistics.

**Table rules.** M4-17 rows enter with `subset = all` and the detection columns of both stat rows
(`observed`, `expected`, `o_minus_e`, `o_over_e`, `hl_shift_log2`), q-values copied (BH is never
recomputed). Every row satisfies `null_model = none` ⇒ p, q NA and `n_perm = 0`; a p-value ⇒
`null_model ≠ none`; a q-value ⇒ a p-value (checked by the assembler).

## Statistician sign-off

Statistician (M4-17 review, PR #89): **signed 2026-10-16** (review round 4 below). Rounds 1–3
(2026-10-15, changes requested) are kept for the record.

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

### Review round 3 (statistician, 2026-10-15): studentised HL implementation (dd663c0)

Scripts and outputs: `local/scratch/review-m4-17/v3/` (cmp.py, tumis.py, tusplit.py, lagv.py, an.py;
Slurm 1850322, 1850439, 1850446, 1850492). Rates are rejection shares at α = 0.01 / 0.05 / 0.10,
1,000 genes, ρ = 0.5, exhaustive circular shifts. Bounds at R = 2,000: 0.0144 / 0.0597 / 0.1134.

1. **Formula: matches the reference.** `tu_variance_components`, `hl_variance` and
   `studentised_hl_null` were checked against the round-2 reference (exp5.py `tus_`) on the same
   seeds: 114 replicates, 455 regions (A2-HL at 1,000 and 4,500 genes, control, A2-MAR).
   max |ΔT| = 1.8e-15 over every shift, identical p, identical σ_u² and σ_e² (to 1e-12). The
   observed T (shift 0) is computed by the same path.
2. **Decision-12 guard: does not hide real regions.** It acts on the replicon's σ_u² + σ_e², not on
   a region's V. With log2 abundances of ~20–30 it triggers only below a residual SD of ~3e-6
   log2. Real proteomics quantification is coarser than that by orders of magnitude, so only
   (numerically) constant residuals reach it. Typical values: σ_u² ≈ 4.0, σ_e² ≈ 0.87, threshold
   7.5e-12. A small V caused by a large region never reaches the guard. σ_u² = σ_e² = 0 needs
   SSE = SSB = 0, i.e. constant data (all-singleton TUs give σ_e² = 0 but σ_u² > 0).
3. **Missing-TU rule: anti-conservative (change requested).** A gene without a TU label becomes its
   own TU. When the missing labels sit in the tested region, the region's V is too small, because
   its genes are still correlated. Region labels missing:
   - control: prophage 0.034 / 0.118 / 0.184, BGC 0.040 / 0.114 / 0.191, CRISPR 0.026 / 0.101 / 0.159;
   - A2-HL and A2-MAR: BGC 0.108–0.111 at α 0.05.

   All labels missing stays at level (0.009 / 0.048 / 0.102). C6 makes `tu_id` non-null, and M4-13
   and `tu_block_bootstrap` refuse NA TU IDs. Either refuse missing labels with a ValueError, or
   ignore labels on the whole replicon when any indexed gene lacks one. Attack test
   `test_atu_missing_region_tu_labels_do_not_change_the_hl_p` accepts either fix and fails on dd663c0.
4. **New weakness: region-specific TU fragmentation (attack A-TU).** The data are the control, with
   no region effect. Each region TU on chr1 is cut into 2-gene chunks; the TUs outside the regions
   are left whole. HL then rejects:
   - prophage 0.017 / 0.082 / 0.143, BGC 0.020 / 0.080 / 0.144, CRISPR 0.015 / 0.071 / 0.124;
   - for comparison, true labels give 0.010 / 0.051 / 0.097 and the same cut genome-wide gives
     0.008 / 0.049 / 0.104 (R = 2,000).

   Over-merging (adjacent TUs merged in pairs) is conservative. An operon predictor that splits
   prophage or BGC operons more than core operons (different intergenic spacing) would do this.
   Gate: `test_atu_region_specific_tu_fragmentation_hl_at_level`, pooled m = 4 as C1, fast and full.
   Slurm 1850492, fast (R = 100, pooled): 0.015 at α 0.01 (pass); 0.080 then 0.055 at 0.05 and
   0.133 then 0.120 at 0.10 (unconfirmed: R = 100 is too small to confirm a rate near 0.08).
   Full (4,500 genes, R = 1,000 per block, pooled m = 4): **confirmed fail at every α**,
   0.019 / 0.019, 0.071 / 0.072 and 0.130 / 0.126 against bounds 0.0131 / 0.0569 / 0.1095.
   The full stream costs ≈ 46 min per block on one core.
   **A label-free fix exists (prototype `lagv.py`).** It keeps T = HL/√V but builds V from
   positions instead of labels:
   - V = Var(mean_in − mean_out) = Σ_{in,in} γ/m² + Σ_{out,out} γ/m′² − 2 Σ_{in,out} γ/(m m′);
   - γ(h) is the replicon's autocovariance of residualised detected abundances at circular
     gene-index lag h, estimated once per replicon, h ≤ 30, truncated at its first non-positive
     value.
   - Results:

     | Scenario (R) | Lag V | Raw HL | TU-label V |
     |---|---|---|---|
     | Control (2,000) | 0.008–0.010 / 0.049–0.057 / 0.097–0.104 | | |
     | Control, ρ = 0.9 (1,000) | 0.006–0.012 / 0.043–0.052 / 0.087–0.110 | | |
     | A2-HL (2,000) | 0.009 / 0.032 / 0.074 | 0.034 / 0.109 / 0.178 | |
     | A2-MAR (2,000) | 0.010 / 0.045 / 0.087 | 0.061 / 0.136 / 0.216 | |
     | Power, prophage +1 log2, α 0.05 (500) | 0.172 | 0.168 | 0.178 |

   - It is immune to TU-label errors by construction. Its cost is the same order as `hl_variance`
     (pairwise distances inside the set).
   - It is a change to amendment 1, so it needs an ADR-0018 amendment 2 accepted by Khoa.
5. **Calibration log ("null_fast crispr hl_shift 0.03 at α 0.01: fail").** The record is stale and
   is not a gating error. `local/logs/calibration/regions.jsonl` is append-only and holds every run
   since 09:44. Lines 47–51 and 115–119 are the pre-pooling runs (Slurm 1849058 era, round 1).
   Jobs 1850200 and 1850201 wrote lines 495–726. They hold no per-type `null_fast` record, only
   `chr1_pooled`, and every check there is pass or unconfirmed. Records carry no run ID, commit or
   timestamp, so runs cannot be told apart. Harness follow-up for the backlog: add `git_head` and
   `slurm_job` / `ts` to every record, or write one file per run. The gating needs no fix.
   Separately, jobs 1850199–1850201 logged head 5a7da49 with the implementation uncommitted;
   dd663c0 was committed 1 min after 1850201 ended. The numbers stand only if the tree equals dd663c0.
6. **False-alarm budget, final gate set** (ADR-0015 two-block rule; per-check probability
   P(block fails)²):
   - Fast: C1 (2 pooled streams × 3 α, plus plasmid AMR 2 × 2 α) 0.010; A2 (4 streams × 3 α)
     0.006 at the measured rates (0.030 if every A2 stream sat at its nominal level); total
     ≈ 0.016.
   - Nightly: C1 (4 pooled + 4 AMR streams) 0.014; A2 0.004–0.009; total ≈ 0.018–0.023.
   - The two A-TU streams add ≤ 0.003 once they pass.

   Acceptable against ≈ 0.04 for the whole nightly. A2-mid O − E at α 0.01, R = 100, is the largest
   single term (0.006, near nominal). No change.
7. Notes: the bootstrap/V inconsistency for missing TU labels (CIs section) goes away with fix 3.
   Effect and CI are unchanged from round 2 (the HL point estimate on the residual scale).

**Sign-off state: changes requested (blocking).**
- (a) Missing-TU rule: item 3, a small code change, gated by the deterministic attack test.
- (b) A-TU: item 4. Either ADR-0018 amendment 2 (label-free lag V, the prototype passes every
  scenario above), or an explicit waiver by Khoa that keeps A-TU as a documented limitation and
  marks it xfail with that reason.

Everything else is signed: formula, guard, C1, A2-HL, A2-HL-MAR, A2-mid, A2-low, the A4 invariance,
FDR families and the false-alarm budget.

### Review round 4 (statistician, 2026-10-16): label-free lag V (04fc8f7, 651cf9c) — signed

Scripts and outputs: `local/scratch/review-m4-17/v4/` (rev.py, an.py; Slurm 1851233, 1851234 for the
side runs, 1851283 / 1851284 for the A-LIN gate). Rates at α = 0.01 / 0.05 / 0.10, 1,000 genes.

1. **Code equals the v3 reference on the same seeds.** `lag_autocovariance`, `hl_variance` and
   `studentised_hl_null` against `lagv.py` (round-3 seeds, chr1, every shift): 7,000 replicates,
   27,909 regions (control ρ 0.5 and 0.9, A2-HL, A2-MAR); max |Δγ| 1.1e-14, max relative |ΔT|
   6.6e-16, identical NaN pattern, identical p in every row (round-3 p reproduced). The prototype
   double-counts lag n/2 on an even replicon of ≤ 60 genes; the implementation counts it once
   (correct). Linear γ and V match a naive in-sequence-pair oracle to 7e-15 (400 random replicons,
   5–90 genes, every non-wrapping shift).
2. **Linear replicons (decision 16, not in v3): at level.** Every replicon set to linear (rates among
   rows with a p; regions crossing the origin get K = 0 by design): ρ 0.5 HL
   0.0084 / 0.0496 / 0.0985 and O − E 0.0084 / 0.046 / 0.098 (R = 4,000 × 4 types); ρ 0.9 HL
   0.0089 / 0.0488 / 0.0957; A2-HL 0.0108 / 0.0465 / 0.0968; A2-MAR 0.0076 / 0.0465 / 0.0978;
   4,500 genes 0.0057 / 0.0456 / 0.0937 (R = 400). The linear 60-gene plasmid is conservative
   (K ≈ 14). New gate **A-LIN** (`test_alin_linear_replicon_hl_at_level`, fast and full; no origin
   crossing, ρ = 0.9, HL pooled over prophage / BGC / CRISPR, m = 3; the scattered AMR set can span
   a linear chromosome, K = 0). Fast 0.0133 / 0.070 / 0.097 (bounds 0.0215 / 0.0752 / 0.1346), pass;
   full (4,500 genes, R = 1,000, Slurm 1851308) 0.0083 / 0.0453 / 0.0993 (bounds
   0.0136 / 0.0580 / 0.1110), pass; ≈ 57 min per block on one core (shard it in the nightly).
3. **V ≤ 0 → NA (decision 16): accepted.** Chromosome: 0 of ~60,000 HL rows (circular ρ 0.9 and
   0.95, R = 4,000 / 2,000; ori–ter cosine of amplitude 1 and 2 log2 at ρ 0.5 / 0.9). 60-gene
   plasmid: 6 / 4,000 at ρ 0.9, 7 / 2,000 at ρ 0.95, 2 / 2,000 linear ρ 0.9, 0–1 elsewhere — so the
   implementer's 3 / 59,000 are small-replicon rows. No bias: with exhaustive shifts the flag is
   decided on the observed set plus every shift, i.e. on the whole orbit, so it is the same whatever
   orbit member was observed and the p of the remaining rows stays exact (rates among non-NA rows are
   at level in every scenario above). The cost is lost HL p on a few small replicons at high
   autocorrelation; fail-closed is right (a silent fallback would need its own calibration).
4. **Dosage gradient (limitation "Replication dosage").** With a random region position, an ori–ter
   cosine leaves level intact (1 log2: 0.0097 / 0.049 / 0.102 at ρ 0.5, 0.0106 / 0.049 / 0.092 at
   ρ 0.9; 2 log2 at ρ 0.9: 0.0105 / 0.052 / 0.100; R = 2,000 / 1,000 × 4 types) and creates no
   V ≤ 0. The limitation (HL partly measures position) stands as written.
5. **A-TU 0.0525 at α 0.05 (bound 0.0569, n = 4,000 pooled): at level.** z = 0.7 against 0.05;
   the same statistic on the unfragmented control gives 0.0516 (R = 2,000 × 4) and 0.046–0.049 at
   ρ 0.9. A-TU and the control cannot differ for the lag V (labels do not enter it); the gap to
   round 3 (0.071) is the fix, not chance.
6. **Fast "unconfirmed" cells: chance, no action.** C1 pooled HL α 0.01: 8 / 400 then 7 / 400
   (bound 7.98 / 400) against full 0.0085 (n = 4,000) and 0.009–0.010 here (n ≈ 28,000). A2-mid O − E
   α 0.01: 3 / 100 then 2 / 100 against full 0.005 and round-2 0.009 (R = 6,000). Seeds are fixed,
   so both cells will keep reporting "unconfirmed" (not a failure) until the replicate streams
   change.
7. **False-alarm budget, final gate set** (two-block rule, P(block fails)², every stream at nominal
   level except A2-mid / A2-low O − E at their measured conservative rates): fast ≈ 0.034 (C1 0.008,
   plasmid AMR 0.002, A2 0.016, A-TU + A-LIN 0.008), nightly ≈ 0.023 (C1 0.008, AMR 0.006, A2 0.005,
   A-TU + A-LIN 0.004). Up from round 3 (0.016 / 0.018–0.023) because A2-HL and A2-MAR are no longer
   conservative under the lag V; nightly is within ≈ 0.04 (ADR-0015). The largest terms are single
   R = 100 streams at α 0.01 (0.0063 each); backlog suggestion, not a condition: pool A2-HL and
   A2-HL-MAR at α 0.01 in fast mode.
8. **Families and outputs (§8.1).** Unchanged: BH within (region_type, stat); every row with a p has
   a q, an effect and a TU-block 95 % CI on the same scale; NA p (K = 0, empty HL side, V ≤ 0) gets NA
   q and stays out of BH.

**Sign-off state: signed.** ADR-0018 amendment 2 and decisions 13–16 are implemented as specified;
every C1, A2, A4, A-TU and A-LIN gate passes. Rounds 1–3 conditions are all met.
