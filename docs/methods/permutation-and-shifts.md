# Permutation p-values, circular shifts and Besag–Clifford

Implements `gippa.stats.permutation` (`perm_p`, `circular_shifts`, `apply_shift`, `choose_shifts`,
`besag_clifford`; GIPPA_SPEC.md §7.3.5(b), §8.2, §8.3, decision 14, Appendix B.4; issue M4-03).
Statistician sign-off: see the end of this note.

## The +1 rule

For an observed statistic T₀ and B null statistics T₁ … T_B,

    p = (1 + b) / (1 + B),     b = #{ i ≥ 1 : T_i is at least as extreme as T₀ }.

If T₀, T₁, …, T_B are **exchangeable** under H₀, the rank of T₀ among the B + 1 values is uniform,
so P(p ≤ α) ≤ α for every α, exactly and for any B (Phipson & Smyth 2010). Consequences:

- p is never 0 and never exceeds 1: 1/(B + 1) ≤ p ≤ 1. The floor 1/(B + 1) is why many-test families
  use Besag–Clifford, Wallenius or fgsea multilevel instead (§8.2, decision 14).
- p is **discrete** on the grid {1, …, B + 1}/(B + 1). Under H₀ with a continuous statistic,
  P(p ≤ α) = ⌊α(B + 1)⌋/(B + 1) ≤ α; it is uniform on the grid, not on (0, 1). The calibration
  suite therefore uses the level check only, never KS (§8.4).
- Ties count as exceedances, so ties only make p larger (conservative).
- Dropping the "+1" (p = b/B) gives p = 0 and is anti-conservative; it is never used.

## Two-sided p: centring

Two-sided p counts null values whose distance from a centre c is at least the observed one:
b = #{ |T_i − c| ≥ |T₀ − c| }. The spec (§7.3.5(b), §8.3) centres on "the null mean". Exactness
requires the statistic |T − c| to treat T₀ and the T_i symmetrically, i.e. **c must be a symmetric
function of all B + 1 values**:

    c = (T₀ + T₁ + … + T_B) / (B + 1)      (mean of the randomization distribution, identity included).

For exhaustive shifts this is the mean over all n placements of the region (the shift k = 0 is a
member of the null orbit), which is what "the null mean" of a randomization distribution means.
Centring on the mean of T₁ … T_B **only** (the first implementation) makes T₀ look more extreme by
the factor (B + 1)/B and is anti-conservative when B is small:

| iid N(0, 1), 4 × 10⁵ replicates | α | null-only centre | pooled centre (exact) |
|---|---|---|---|
| B = 9 | 0.10 | 0.141 | 0.099 |
| B = 19 | 0.05 | 0.062 | 0.050 |
| B = 19 | 0.10 | 0.118 | 0.100 |
| B = 99 | 0.05 | 0.053 | 0.051 |
| B = 999 | 0.05 | 0.048 | 0.048 |

Small K is routine in GIPPA: plasmids with 10–20 genes (K = n − 1) and short linear contigs
(K = n − span). The attack nulls `tiny_circular_n10`, `tiny_circular_n20` and `tiny_linear_n12`
(tests/calibration/attacks/test_attack_permutation.py) fail with null-only centring (rates 0.14,
0.064 and 0.126 against bounds 0.109, 0.056, 0.109) and pass with pooled centring (0.104, 0.052,
0.091). **Implemented (5e80944):** `perm_p` two-sided centres on the pooled mean of the finite
values among T₀, T₁ … T_B; ±inf null values are kept in B and always count as extreme (an infinite
T₀ is left out of the centre). M4-03 example 1 is now
`perm_p(5.5, [1, 2, 3, 4, 5]) = 2/6` (pooled centre 3.417; |1 − 3.417| = 2.417 > |5.5 − 3.417|
= 2.083); examples 2 and 4 and B.4 are unchanged. Centring is translation-equivariant either way,
so the "invariant to adding a constant" property still holds.

Besag–Clifford has no pooled set to centre on; for two-sided use the caller must centre on a value
known **a priori** (e.g. 0 for a difference in means under label permutation, whose permutation
mean is exactly 0), never on a running estimate.

## Float ties

Exhaustive shifts sum equal multisets in different orders, so exact ties can differ by rounding.
`perm_p` counts T_i as an exceedance when it is within tol of the threshold:

    one-sided:  tol = 1e-12 · max(1, |T₀|)
    two-sided:  tol = 1e-12 · max(1, |T₀ − c|, |c|)

The implementer added |c| to the two-sided scale (a deviation from the spec text
`1e-12 · max(1, |x|)`). **Assessment: correct and needed.** |T_i − c| is computed after
subtracting c, so its rounding error scales with |c|, not with the deviation; without |c| a large
centre (log₂ intensities ≈ 20–30 summed over a region) would leave ties unresolved. The tolerance
only ever adds exceedances, so it can only make p larger (valid, at worst slightly less powerful).

Limitation (non-blocking): the rounding error of a sum scales with the **summands**, not with the
sum. With large summands that cancel (terms ±1000.1, width 50; or ±100000.1, width 200) a true tie can be missed:
in a probe, 1/300 and 7/300 replicates had float p below the exact rational p (one tie missed;
effect ≤ 1/(B + 1)). For GIPPA's region statistics (z-scores or log intensities, widths ≤ a few
hundred genes) the error stays far below 1e-12·scale; the attack `test_tie_tolerance_never_lowers_p`
(heavy ties on log₂-like decimals, width 50, exhaustive shifts) checks float p ≥ exact p. If a
caller ever sums large cancelling terms, it should pass rounded statistics or a scale.

## Circular shifts (§7.3.5(b), §8.1 "region effect")

A region is a set I of gene indices on a replicon of N genes. The null statistics are the same
statistic on I + k (mod N) for k = 1 … N − 1 (`circular_shifts(n, idx, circular=True)`;
`apply_shift`). B.4: N = 100 gives K = 99, and an observation more extreme than all 99 gives
p = 1/100.

**What the shift null preserves.** Region size, the gaps between its genes (non-contiguous AMR sets
move as a block), and the joint distribution of abundance along the replicon, in particular its
**autocorrelation**: each shifted set sees a stretch of the same autocorrelated series, so the
null variance of a region sum includes the positive covariances between neighbours, exactly as the
observed one does. **What it breaks:** only the region's position. H₀ is "the region's position
is unrelated to abundance", formalised as the abundance series being stationary on the circle
(rotation-invariant in distribution); then T₀ and the shifted T_k are exchangeable and the +1 rule
is exact.

**Why gene permutation is wrong here.** Permuting gene labels draws the null region from scattered
genes, whose sum has variance ≈ w·σ²; the observed contiguous region has variance
≈ w·σ²·(1 + ρ)/(1 − ρ) under AR(1). The test is then anti-conservative. Attack
`test_naive_permutation_is_anticonservative` (circular-stationary AR(1), ρ = 0.9, w = 20, N = 1000)
confirms the attack has teeth: label permutation rejects far above 0.05, while the shift null on
the same data stays at level (table below).

Caveats the callers own (M5): stationarity fails near replication-origin/terminus gradients
(gene dosage) — shifts then compare a region with regions at other dosage; GIPPA mitigates this by
detrending (§7.3.5) and by reporting plasmids as replicon-level effects (§8.5 attack 4).

**Linear replicons.** On a linear contig a wrapped set is not a region, so only translations that
keep I inside [0, N) are used: k with max(I) + k < N (right) or min(I) + k ≥ N (left by N − k).
With span = max(I) − min(I) + 1 this gives **K = N − span** shifts, and with the identity the
N − span + 1 placements of the block; the set of placements is the same from every placement, so
the exhaustive linear null is exact under translation-stationarity. K = 0 (the region spans the
contig) → no test; callers report `p_value = NA`.

**Exhaustive vs sampled.** Exhaustive (all K shifts) while N ≤ `context.exhaustive_shift_max_genes`
(default 20,000), otherwise B random shifts drawn **without replacement** from the valid set (all of
them when B ≥ K). Sampling a uniformly random subset of the other placements keeps T₀ exchangeable
with the chosen ones, so the +1 rule stays exact (Hemerik & Goeman 2018); without replacement
avoids duplicated null draws and is never worse. Determinism: the shifts depend only on the
`Generator` state (from `gippa.rng.for_step`) and are returned sorted.

## Besag–Clifford sequential p (§8.2, decision 14)

Draw null statistics one at a time; stop at the h-th draw ≥ T₀ (h = 10 by default). If that is
draw L, p = h/L; if `max_B` draws give b < h exceedances, p = (b + 1)/(max_B + 1). This is the
sequential Monte Carlo p-value of Besag & Clifford (1991, Biometrika 78:301): it is exactly valid,
P(p ≤ α) ≤ α under exchangeable draws, with ties counted as exceedances. It spends few draws on
clearly null tests (L ≈ h/p) and up to `max_B` on small p, so its floor is 1/(max_B + 1) rather than
1/(B + 1) for the fixed B, at a fraction of the cost. h < 1 or max_B < h is rejected.

## Simulation evidence

Level criterion (§8.4): rate ≤ α + 2√(α(1 − α)/M), M independent p-values.

Fast calibration (tests/calibration/test_stats_permutation_calibration.py, R = 100): 6/6 pass.
Larger R = 2,000 (same generators, scratch run), rates at α = 0.01 / 0.05 / 0.10 (bounds 0.0144 /
0.0597 / 0.1134):

| Null | rates |
|---|---|
| two-sided label permutation, B = 199 | 0.0110 / 0.0515 / 0.0955 |
| one-sided greater / less, B = 199 | 0.0080 / 0.048 / 0.0965; 0.0065 / 0.046 / 0.098 |
| circular exhaustive, N = 1000, w = 20 | 0.0095 / 0.0515 / 0.1105 |
| linear exhaustive, N = 400, w = 10 | 0.0060 / 0.051 / 0.1035 |
| Besag–Clifford h = 10, max_B = 1000 | 0.0100 / 0.0540 / 0.1070 |

Attack nulls (tests/calibration/attacks/test_attack_permutation.py), full suite (R = 1000; M in
brackets), α = 0.01 / 0.05 / 0.10, current code → with pooled centring:

| Attack | current | pooled |
|---|---|---|
| AR(1) ρ = 0.9, exhaustive shifts (10⁴) | 0.0117 / 0.050 / 0.099 | 0.0115 / 0.050 / 0.098 |
| AR(1), sampled B = 199 w/o replacement (10⁴) | 0.0108 / 0.051 / 0.102 | 0.0102 / 0.049 / 0.100 |
| heavy ties (10⁴) | 0.0072 / 0.041 / 0.080 | 0.0072 / 0.041 / 0.079 |
| Besag–Clifford, tied Bernoulli groups (10³) | 0.003 / 0.018 / 0.032 | same |
| tiny circular N = 10 (5 × 10⁴) | 0 / 0 / **0.140** | 0 / 0 / 0.098 |
| tiny circular N = 20 (5 × 10⁴) | 0 / **0.061** / **0.119** | 0 / 0.050 / 0.100 |
| tiny linear N = 12, w = 2 (5 × 10⁴) | 0 / 0 / **0.126** | 0 / 0 / 0.090 |

(Bold = above bound.) A 10⁵-replicate check of the sampled AR(1) attack at four seeds gave
0.0103 / 0.0513 / 0.1020 (null-only) and 0.0099 / 0.0496 / 0.1002 (pooled). Disclosure: the first
draft seeded fast and full suites from one stream; the full sampled-AR(1) case then exceeded the
α = 0.01 bound by 0.0004 under both centrings (2.5 SE, a chance excursion given the 10⁵ check).
The streams were split per suite once (`attack.permutation.<name>.<suite>`); no further reseeding.

## Missing and infinite values

NaN null values are dropped (B counts the rest); NaN T₀ or B = 0 gives NaN. ±inf null values are
kept in B and count as extreme in their direction (+inf for "greater", −inf for "less", both for
"two-sided"), so `perm_p(5, [1, 2, inf], "greater")` = 2/4. The first implementation dropped them,
which lowered p; fixed in 5e80944 (lead decision). `choose_shifts` reports `exhaustive=True`
whenever every valid shift is used, including n > n_max with B ≥ K.

## Sign-off

First review (2bee83c): BLOCK on null-only centring. Re-review of 5e80944: **APPROVE**. All 24
permutation calibration and attack tests pass (fast and full), and 144 unit/property tests pass.
An independent level check of two-sided `perm_p` (R = 20,000 per cell, iid normal and lognormal
B + 1 values) gives rates at α = 0.05 / 0.10 of 0.047–0.051 / 0.097–0.102 for B = 19, 0.099–0.102
at α = 0.10 for B = 9, and 0.0089–0.0098 / 0.049–0.051 / 0.096–0.102 at α = 0.01 / 0.05 / 0.10
for B = 99. All are within bounds. Approved: the +1 rule, one-sided and pooled two-sided p, tie
handling including the |c| scale, ±inf handling, circular and linear shift sets, sampling without
replacement, Besag–Clifford, and determinism. The tie-tolerance limitation above stays open and is
non-blocking.
