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

## Amendment 1 (accepted-by-khoa 2026-10-16)
Statistician review round 2 (https://github.com/khoa-yelo/giap/pull/89#issuecomment-6205786444;
docs/methods/region-tests.md "Review round 2"; local/scratch/review-m4-17/v2): residualisation (decision 1)
removes the bias (prophage HL −0.08 log2 vs +0.74 raw) but HL stays anti-conservative because a region's
detected-gene count differs from that at shifted positions (≈ 7 vs ≈ 22 of 40), so the raw HL's spread under
the shift null is wrong; this holds even when detection ignores abundance (new attack null A2-HL-MAR: 0.147 /
0.151 at α 0.05), and no covariate fit fixes it (oracle correction 0.131).

Decision 1 is extended: the circular-shift p-value is computed from the studentised statistic
T = HL / √V at every shift, V a per-TU variance estimated once per replicon (formula in the methods note);
the reported HL effect and CI are unchanged (residual scale). Evidence at α 0.01/0.05/0.10: A2-HL
0.005/0.019/0.046 (1,000 genes, R 4000) and 0.003/0.021/0.047 (4,500, R 2000); A2-HL-MAR ≤ 0.038 at 0.10;
C1 pooled at level; power without confounding unchanged (0.166 vs 0.156).

## Amendment 2 (proposed, needs Khoa) — 2026-10-16
Statistician final review (https://github.com/khoa-yelo/giap/pull/89#issuecomment-6208970265;
local/scratch/review-m4-17/v3): amendment 1's per-TU variance V depends on the TU labels, so a region whose
TUs are fragmented more finely than the rest of the replicon (new attack A-TU) is anti-conservative:
0.019 / 0.071 / 0.128 at α 0.01 / 0.05 / 0.10 (full mode, Slurm 1850492).

Proposed: V is estimated label-free from the lag autocovariance of the residual abundance along the replicon
(formula in docs/methods/region-tests.md, statistician's v3 section); T = HL/√V as in amendment 1; reported
HL effect and CI unchanged. Statistician: passes every scenario (A2-HL, A2-HL-MAR, A2-mid, A-TU, C1) with
the same power. Alternative: waive A-TU (xfail, documented limitation).
