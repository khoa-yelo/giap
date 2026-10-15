# ADR-0018: Region tests robust to detectability confounding

- Status: accepted-by-khoa (2026-10-15)
- Date: 2026-10-15
- Issue: M4-17 (PR #89); statistician review https://github.com/khoa-yelo/giap/pull/89#issuecomment-6198595819

## Context
Under the attack-2 scenario (detection depends on abundance; prophage genes are short), the region tests
of GIAP_SPEC.md §7.3.5 are anti-conservative (statistician, Slurm runs under local/scratch/review-m4-17):
- (b) HL abundance shift: prophage HL rejects 0.21–0.25 at α = 0.05 (bias ≈ +0.75 log2).
- O−E detection: when a region's p_detect sits nearer 0.5 than the replicon's, plain O−E rejects 0.074 at
  α = 0.05 (A2 only covered the conservative direction).
Also §8.1 says "BH within type"; the code (and the test-author) use BH within (type, stat).

## Decision
1. §7.3.5(b): abundance is residualised on the M4-13 detectability covariates (log length, n_obs, GRAVY)
   before the HL/circular-shift test (statistician: 0.067 at α = 0.05, R = 300).
2. The O−E test statistic is the studentised sum Σ(d − p)/√Σ p(1 − p) (0.048 at α = 0.05); the reported
   effect stays O−E.
3. BH families are (region type, statistic): two statistics per type are different hypotheses; this is
   at least as strict as "within type" for each statistic.

## Consequences
Region-test p-values change on real data where detectability differs between regions; methods note and
report text follow. Calibration gates (C1, A2 both directions) verify level under attack 2.
