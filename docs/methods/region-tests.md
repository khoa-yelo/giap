# Region tests: O − E and Hodges–Lehmann shift under circular gene-index shifts

Implements `giap.context.regions` (GIAP_SPEC.md §7.3.5(b), §8.1, Appendix B.4; issue M4-17,
docs/specs/M4-17-region-tests.md). Reuses `giap.stats.permutation` (M4-03), `giap.stats.bootstrap`
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

- **O − E** S(I) = Σ_{i∈I} (d_i − p_i); `observed` = Σ d_i, `expected` = Σ p_i, `effect` = S.
- **HL shift**: the median of all pairwise differences a_i − a_j between the region's detected genes
  with a finite log2 abundance and the replicon's other such genes (`hodges_lehmann`; mean of the two
  middle values for an even count). NaN when either side is empty. Detected genes without abundance
  (policy `exclude`) drop out on both sides. Comparing within the replicon only keeps a high-copy
  plasmid's dosage out of the shift (attack 4).

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
resamples inside and outside TUs as two strata.

## Multiple testing and flags

BH within each (region_type, stat) family (lead decision 5; detection and abundance answer different
questions). Flags: `small_region` (< 3 eligible genes; tests still run when K > 0; with 0 eligible
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
- C3 power (recorded, not gated; BGC O − E, α = 0.05, R = 100): +10 % 0.11, +20 % 0.21,
  +30 % 0.46, +50 % 0.88 (`local/logs/calibration/regions_power.jsonl`).

## Performance

4,500-gene chromosome, 30 regions, exhaustive shifts, B = 1,000: ~35 s on one core
(benchmarks/test_region_tests_bench.py; budget 60 s); the TU-block bootstrap of HL dominates.

## Statistician sign-off

Pending (lead decision 5: claims C1–C3, A2, A4 reviewed before merge).
