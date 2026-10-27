# TU-block bootstrap, circular block bootstrap and the Hodges–Lehmann shift

Implements `gippa.stats.bootstrap` (`tu_block_bootstrap`, `percentile_ci`,
`circular_block_bootstrap`) and `gippa.stats.robust.hodges_lehmann` (GIPPA_SPEC.md §7.3.5(a)–(b),
§8.2, §8.3, decision 14; issue M4-04). Statistician sign-off: see the end of this note.

## Why the resampling unit is the TU

Genes in one transcription unit share a promoter, a transcript and most of their regulation, so
their abundances, detection indicators and O − E contributions are positively correlated. With an
intra-TU correlation ρ and TUs of size m, the variance of a gene-level mean is inflated by the
design effect 1 + (m − 1)ρ relative to iid genes. A bootstrap that resamples genes treats them as
independent and estimates the iid variance, so its CI is too narrow by √(1 + (m − 1)ρ). TUs, not
genes, are the units that are (approximately) independent, so §8.2 rules "resample TUs, not genes".

`tu_block_bootstrap(df, tu_col, stat_fn, B, rng, level, strata_col)`:

- each replicate draws, uniformly with replacement, as many TU IDs as there are TUs, and
  concatenates the drawn TUs' rows **whole** (a TU drawn twice contributes all its rows twice; a TU
  is never split; row order inside a TU is kept). TU sizes vary, so the replicate gene count varies;
  `stat_fn` must therefore be a statistic of the rows (sum, mean, HL, O − E), not of a fixed n;
- **stratified** two-group statistics (`strata_col`): TUs are drawn within each stratum, as many as
  the stratum has, so every replicate keeps the inside/outside TU counts. This is what §7.3.5(b)
  asks for the HL shift ("resample the region's TUs and, for the HL shift, the outside TUs too"):
  both groups carry sampling error, and pooling them would let the region's size drift. A TU's
  stratum is the stratum of its first row; callers must not let a TU straddle a region boundary
  (M4-17 must assign each TU wholly to one side before calling);
- NaN replicate statistics are dropped; with fewer than B/2 valid replicates the CI is (NaN, NaN)
  and `n_valid` says why. `stat_fn` errors propagate.

Randomness: the caller passes the `Generator` from `gippa.rng.for_step` (§6.6); draws are made per
stratum in sorted stratum order, each a (B × n_s) block, so results are reproducible given the seed.

Assumption: TUs are independent (or weakly dependent). Long-range dependence along the replicon
(gene dosage, replichore gradients) is handled by the circular shift null for p-values (see
permutation-and-shifts.md) and by detrending; the bootstrap CI is for the effect size only.

## CI method: percentile

Spec §8.3 fixes the percentile CI: `np.quantile(boot, [(1 − level)/2, (1 + level)/2],
method="linear")`. Not BCa, not basic. Properties:

- transformation-respecting and always inside the range of the statistic (useful for bounded
  statistics such as fractions);
- first-order accurate only; it inherits the bootstrap's small-sample bias. For a mean of n
  independent units the bootstrap variance is (n − 1)/n of the plug-in variance and the true pivot
  is t_{n−1}, not z, so the expected coverage is

      c(n) = P(|t_{n−1}| ≤ z_{0.975} √((n − 1)/n)):
      n = 5: 0.846   10: 0.904   12: 0.913   20: 0.929   30: 0.936   100: 0.948.

  **Percentile CIs from fewer than ~20 TUs under-cover.** Downstream modules (M4-17 region tests)
  must report `n_tu` with every CI and flag CIs from < 10 TUs (`ci_few_units`); they are not to be
  read as 95 % intervals. BCa or a studentised bootstrap would reduce the deficit and is a later
  improvement (spec M4-04 lead question 1), not required for M4.
- Monte Carlo error: with B = 1000 each endpoint rests on ~25 tail replicates; the endpoint's MC
  SD is about 0.085 bootstrap SDs, enough for reporting. `--publication` uses B = 10000 (§5.1,
  `stats.bootstrap`), cutting it by √10. B = 1000 is the default (decision 14).

### Recorded coverage (fixed seeds; tests/calibration)

| Scenario | TUs | statistic | coverage | c(n) |
|---|---|---|---|---|
| test_stats_bootstrap_calibration (gated) | 100, 1–4 genes, TU effect + noise | sum | in band | 0.948 |
| attack `few_large_tus` (gated, R = 100) | 30, 5–15 genes, ρ = 0.9 | mean | ≥ bound | 0.936 |
| same data, **gene-level** resampling (teeth) | 300 genes | mean | 0.44 | — |
| attack `hl_two_group` (gated) | 12 in + 80 out, 2–6 genes, ρ = 0.8 | HL | ≥ bound | 0.913 |
| exploratory, R = 400 | 12 in + 80 out, ρ = 0.8 | HL | 0.925 | 0.913 |
| exploratory, R = 400 | 30 in + 80 out, ρ = 0.8 / ρ = 0 | HL | 0.945 / 0.945 | 0.936 |
| recorded, R = 200 | 5 / 10 TUs, ρ = 0.9 | mean | 0.86 / 0.89 | 0.846 / 0.904 |

The TU-block CI tracks c(n); the gene-level bootstrap covers 44 % on the same data, which is why
the TU is the resampling unit. The attacks gate coverage ≥ c(n) − 2√(c(n)(1 − c(n))/R)
(tests/calibration/attacks/test_attack_bootstrap.py).

## Circular block bootstrap (Politis–Romano)

`circular_block_bootstrap(x, block_len, B, rng)` returns a (B × n) array: per replicate,
⌈n/ℓ⌉ start positions uniform on [0, n), each block ℓ consecutive values wrapping modulo n,
concatenated and truncated to n. Wrapping makes every position equally likely to be sampled
(unlike the moving-block bootstrap, which under-samples the ends), which matches a circular
replicon. For a linear contig the wrap joins the two ends; with ℓ ≪ n the induced artefact is
O(ℓ/n) and is accepted (spec uses the same function for both).

## Block length

The caller chooses the block length; this function does not. The island circular block bootstrap
(§7.3.5(d) step 5) and the ori–ter slope CI (§7.3.5(e)) use the same rule, from ADR-0021:

1. **ℓ̂** = `gippa.stats.autocorr.choose_block_len(z, circular, B, rng)`: the smallest lag at which
   the TU-level Moran's I is not significant (raw two-sided p > 0.05 under the gene-label null,
   bounded to [1, 20]). See docs/methods/autocorrelation.md "Block length".
2. **ℓ'** = `block_len_prime(ℓ̂, T)` = min(20, 2ℓ̂ + 1, T), where T is the number of scored TUs.
   ℓ' is the block length actually used. An integer `context.islands.block_len` overrides it
   unchanged.

**Why ℓ' and not ℓ̂.** ℓ̂ marks where the correlation stops being detectable, not the block length that
a window-sum maximum needs (ADR-0021).
- With ℓ̂ as the block length, island FWER on a 4,500-gene chromosome is 0.205 / 0.26 / 0.27 at
  ρ 0.5 / 0.9 / 0.95, and the ori–ter CI covers only 0.88–0.93.
- With ℓ', FWER is 0.003–0.040 and the ori–ter CI covers **0.943** at ρ 0.5 (R = 300, Slurm 1863166;
  the 0.96 first reported was an R = 50 fluctuation). Coverage runs 0.875–0.945 across cells and arms,
  so the interval is slightly narrow everywhere — see docs/methods/oriter.md.

**Validity domain.** ADR-0021 **Amendment 1** supersedes the original "fewer than 1,000 scored TUs"
rule and **drops the `uncalibrated` flag**: re-measurement at R = 1,000 showed the scored-TU count
decides nothing on its own (400 TUs at ρ 0.5 give FWER 0.005; 1,000 TUs at ρ 0.9 give 0.095–0.208;
2,500 TUs at ρ 0.9 give 0.024). What matters is whether strong TU-level autocorrelation is resolved by
Moran's I to a long enough lag: every breach has ℓ' ≤ 7, and censoring is what hides the correlation
(the same ρ 0.9 replicon gives ℓ̂ = 9 and FWER 0.007 uncensored, ℓ̂ = 2 and 0.208 at 55 % detected).
No per-replicon statistic certifies the domain, and selecting on one is adverse, so the run report
carries ℓ̂, ℓ' and the lag-1 Moran's I per replicon instead of a flag. See docs/methods/islands.md,
"The island step (M5-06)".

**Joins between blocks.** ℓ = 1 is the iid bootstrap. Blocks join independent stretches, so the
bootstrap slightly under-represents dependence across block joins. With ℓ above the correlation
length this bias is small; the M5 island calibration checks FWER with this null.

Signed off by the statistician (PR #116 re-review).

## Hodges–Lehmann shift

HL(x, y) = median{x_i − y_j : all i, j}, the two-sample location estimator paired with the
Wilcoxon rank-sum test: it is the Δ that makes the rank-sum statistic of x − Δ vs y central.
Robust (breakdown ~29 %), invariant under common shifts, HL(x + c, y) = HL(x, y) + c and
HL(y, x) = −HL(x, y).

- **Ties**: all m·n differences count with multiplicity; tied values give tied differences and
  need no correction for the point estimate. Even m·n: mean of the two middle order statistics.
- **Missing**: NaN raises (callers drop undetected genes first, §7.3.5(b) uses detected genes only).
  Empty sample → NaN.
- **CI**: no analytic (Moses/Wilcoxon-inversion) CI is used, because genes are not independent;
  the CI comes from the stratified TU-block bootstrap above.
- **Exact computation**: direct enumeration (`np.median(np.subtract.outer(x, y))`) for
  m·n ≤ 5 × 10⁴; above, bisection on Δ counting pairs with x_i − y_j ≤ Δ by binary search on sorted
  y, then a partition of the few remaining candidates. Float subtraction is monotone, so the counts
  use the very same float differences as enumeration and the result is bit-identical to it. The
  spec's threshold is 10⁶; 5 × 10⁴ was verified by the statistician (below) and is ~6× faster at the
  M4-17 size (m = 200, n = 4500: 2.0 ms vs 12.4 ms).
- **Non-finite values** (NaN, ±inf from log2(0)) raise `ValueError` on both paths. Differences that
  overflow float64 fall back to enumeration up to m·n = 2 × 10⁷, else raise.

## Statistician sign-off (M4-04, PR #31)

Verified on branch m4/bootstrap (commit 75d7a92; re-verified at 5850e09). **APPROVED.**

1. Bisection HL vs `np.median(np.subtract.outer)` with the threshold forced to 0: 297 cases
   (sizes 1×1 … 1000×999; normal, integer ties, decimal ties, 70 % one-value ties, Cauchy × 10⁶,
   1e-300 scale, ±0.0) — all bit-identical. At the 5 × 10⁴ boundary (200×250, 250×201, 1×50001,
   50001×1, 2×25001, 224², 225×223) and 2000×3000, 200×4500: 66 cases, all identical.
   **The 5 × 10⁴ deviation is acceptable**; the spec text should read 5 × 10⁴.
2. **Defect (blocking):** non-finite inputs on the bisection path. x containing +inf (300×600):
   bisection 2.0 → inf vs enumeration 1.5; x containing −inf: 2.0 vs 0.5; differences that
   overflow (|values| ~ 1e308) or inf on both sides: infinite loop (mid = NaN, `lo`/`hi` never
   move). The direct path is fine. Fix: raise `ValueError` on non-finite input (like NaN), and fall
   back to enumeration (or raise) when `x.min() − y.max()` or `x.max() − y.min()` is not finite;
   add a regression test above the threshold. **Resolved in 5850e09**: all four reproducers now
   raise or match enumeration (overflow 400×400 → 0.0 = enumeration; 2.5 × 10⁷ overflowing pairs
   → ValueError). Exactness re-run: 264 forced-bisection cases + 10 × (200×4500) all equal;
   200×4500 takes 1.9 ms vs 12.2 ms by enumeration. Nit (not blocking): with only ±0.0 inputs the
   result can be −0.0 where `np.median` gives 0.0 (equal as numbers; "bit-identical" holds up
   to the sign of zero).
3. TU blocks: never split, internal order kept, stratum counts preserved (code + property tests);
   circular blocks wrap modulo n; RNG only from the caller's `Generator`. Correct.
4. Small-TU bias of the percentile CI: tracks c(n), documented above; callers flag < 10 TUs.
