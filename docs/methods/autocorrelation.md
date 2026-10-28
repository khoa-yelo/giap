# Autocorrelation of abundance along the replicon

Issue M5-01 adds the shared detrending helper and M5-02 the lag statistics, nulls and block length;
the per-replicon autocorrelation tables are M5-05 (GIPPA_SPEC.md §7.3.5(c), ADR-0020 with Amendment 1,
ADR-0021).

`gippa.stats.autocorr` provides `morans_i`, `join_count`, `label_permutations`, `tu_blocks`,
`tu_block_permutations`, `lag_test`, `choose_block_len` and `block_len_prime` (issue M5-02, lead
decisions 4–22; GIPPA_SPEC.md §7.3.5(d) step 5, §7.3.6, §8.3, Appendix B.4). Statistician sign-off for
the statistics is at the end of this note.

## Detrending

`gippa.stats.smooth.loess_detrend(y, pos, frac, *, circular, period=None, min_points=3)` and
`loess_matrix(...)` (the n × n smoother H with `loess_detrend(y, …).fit == H @ y`). M5-05 detrends
log2-abundance residuals along gene index and M5-06 detrends Blom scores along TU order, both with
span 10 % and circular padding on circular replicons, so both use this one helper.

- **Estimator.** Plain local-linear LOESS: at each point a weighted least-squares line over its k
  nearest neighbours by `pos`, tricube weights (1 − (d/radius)³)³ with radius the distance to the
  farthest neighbour, and **no robustness iterations** (`it = 0`, lead decision 2). Robust
  reweighting would treat an island as an outlier and keep it out of the trend, so the detrend
  would depend on the signal being tested and the calibration would have to cover that
  data-adaptive step. Plain LOESS is linear in y: one matrix per set of positions, so ADR-0020 can
  apply it to B permuted replicates. The implementation reproduces
  `statsmodels.nonparametric.smoothers_lowess.lowess(y, pos, frac=k/N, it=0, delta=0)` operation
  for operation (≤ 1e-10; neighbourhood ties keep the window on the left, fewer than two weights
  above 1e-12 give fit = y, weighted variance of x floored at 1e-12).
- **Span in points.** k = min(n, max(3, ⌈frac·n − 1e-9⌉)) (the 1e-9 keeps 0.1·30 at 3). `pos` are
  gene indices or TU ordinals, never base pairs; gaps are allowed and respected by the neighbour
  search. ADR-0020 passes ordinals 0 … n − 1 among the genes with values.
- **Circular padding.** On a circular replicon the last k points are copied in front at
  `pos − period` and the first k behind at `pos + period`; the fit at the original points is kept
  and the padded copies feed back into their original columns of H. This removes the end effect
  (a cosine over 100 points, span 10: max residual 0.0071 circular vs 0.0098 linear, the latter at
  the ends) and makes the fit equivariant under rotation of the origin.
- **Small n.** n = 0 gives an empty fit (k = 0). For 1 ≤ n < 3 no local line exists and the fit is
  mean(y) (residuals centred, k = n).
- **The 3-point floor is degenerate** (lead decision 4). With k = 3 the farthest neighbour sits at
  the radius and gets weight 0, and the line through the remaining two points (one is the point
  itself) returns y: fit = y and every residual is 0, for any strictly increasing positions. At the
  10 % span this covers n ≤ 30 values. The detrended series does not exist there, so M5-05 and
  M5-06 report NA with flag `few_values` for any replicon whose k < 4 instead of computing
  statistics over zero residuals. Power by replicon size is recorded in M5-07.

## Statistics and nulls

### Index and lags

Each replicon is a series of genes in gene-index order. Gene i and gene j form a lag-k pair when:

- **circular replicon:** j = (i + k) mod n, for i = 0 … n − 1. The lag must satisfy 1 ≤ k ≤ ⌊(n − 1)/2⌋.
  A larger k would count each pair twice (k and n − k give the same pairs) or pair a gene with itself.
- **linear contig:** j = i + k, with i + k < n. The lag must satisfy 1 ≤ k ≤ n − 1.

A statistic at an invalid lag does not exist and is NaN.

### Moran's I at lag k

Only genes with a value count. With n such genes, x̄ is the mean of their values and z_i = x_i − x̄.
W_k is the number of lag-k pairs in which both genes have a value.

    I_k = ((1/W_k) Σ_{pairs} z_i z_j) / ((1/n) Σ_i z_i²)

I_k does not exist (NaN) when:

- the lag is invalid;
- n < 2;
- W_k = 0;
- Σ z² ≤ 1e-24 · max(1, max x²), meaning there is no variance up to rounding.

Example: the alternating series +1, −1, … gives I_1 = −1 and I_2 = +1 (Appendix B.4).

- I_k is invariant to affine maps x → a·x + b with a ≠ 0.
- I_k is not bounded by [−1, 1] in small samples: the numerator averages over W_k pairs while the
  denominator averages over n genes.

### Join count: raw value and Pearson residual test statistic

**Raw count.** d_i marks whether eligible gene i is detected. The raw count is BB_k, the number of
lag-k pairs in which both genes are detected. The table reports it as `value`.

**Test statistic.** The test does not use BB_k, because detection clusters for reasons other than
expression: §8.2 requires detection tests to adjust for detectability. The test uses the Pearson
residual (ADR-0020 Amendment 1, lead decision 17):

    r_i = (d_i − p_i) / √(p_i (1 − p_i)),       S_k = Σ_{lag-k pairs} r_i r_j

- p is the cross-fitted `p_detect`, unclipped.
- A gene with p ∈ {0, 1}, with p NaN, or with d missing counts as missing: it has no residual
  (r = 0), so it joins no pair.
- A non-NaN p outside [0, 1] can only come from an upstream bug, so `join_count` and
  `join_count_matrix` raise `ValueError` rather than drop the gene silently (decision 21).
- A valid lag with no pair left gives S_k = 0.
- `effect` = S_k − its null mean.

**Assumption.** The test is valid when the d_i are independent Bernoulli(p_i) given a *calibrated*
p. It is therefore conditional on the cross-fitted `p_detect` being calibrated (M3 calibration). If
`p_detect` is biased in a region (e.g. too low in an island), r has a non-zero mean there, and the
test reports that miscalibration as detection clustering. A rejection means "clustering beyond the
detectability model" (decision 22).

**Why the earlier statistic was replaced.** ADR-0020 decision 3 first used the raw residual
Σ (d_i − p_i)(d_j − p_j). Its variance is Σ v_i v_j with v = p(1 − p). The permutation null averages v
over random pairs, so when low `p_detect` clusters along the replicon, the null is too narrow and the
test rejects too often. The statistician's attack set this up on a circular replicon of 1,000 genes:
30-gene islands with p ~ U(0.01, 0.05) among genes with p ~ U(0.3, 0.9), α = 0.05 (Slurm 1861073,
1861078). Rejection rates under both nulls:

| Island share | 10 % | 20 % | 30 % | 50 % |
|---|---|---|---|---|
| raw residual | 0.059 | 0.074 | 0.086 | 0.105 |

**Pearson residual results.** Dividing by √(p(1 − p)) gives every gene unit variance under the null.
- Level is 0.045–0.056 in every scenario above and in the original ADR-0020 attacks.
- Power costs about 10 %.
- At very small p (0.001–0.005), α = 0.01 gives a rejection rate of 0.019. This is recorded, not gated.
- A floor on p in the denominator was rejected because it brings the clustering failure back.

The engine attack file `tests/calibration/attacks/test_attack_autocorr_engine.py` gates the low-p
island arms. It also has an arm that must fail (raw residual, 50 % islands), which shows the gate can
detect the failure.

### The two nulls (B permutations each)

**Gene-label null** (`label_permutations`) asks whether there is any structure.
- Values are shuffled uniformly among the replicon's genes that have values (`movable`). Genes without
  a value stay where they are, so the pattern of missing values is fixed and W_k does not change.
- For the join count, p moves together with d.

**TU-block null** (`tu_block_permutations`) asks whether there is structure beyond operons.
- A block is a maximal run of consecutive genes with the same TU label. On a circular replicon, a run at
  the end and a run at the start with the same label form one block (a TU that crosses the origin).
- A TU interrupted by another TU's gene forms several blocks.
- Whole blocks are placed in a uniformly random order. Each block keeps its internal order and
  orientation; a wrapped block is emitted as its end part followed by its start part.
- Values, gaps and p move with their TU.
- Under this null, within-TU correlation is part of the null, so only structure between TUs counts as
  a signal.

### p-value, the constant-null rule and the null summaries

`lag_test(obs, null)` drops NaN null entries; `n_valid` counts the rest, including ±inf (decision 13).

**p-value.** p = `perm_p(obs, null, "two-sided")`. It is centred on the pooled mean of the observation
and the null and uses the +1 rule, so it is exactly valid under exchangeability
(docs/methods/permutation-and-shifts.md).

**Constant-null rule (decisions 6 and 14).** If the null is constant, p is NA and the row is flagged
`not_estimable`.
- "Constant" means the spread of a null of finite values is ≤ 1e-12 · max(1, max |v|), which is
  `perm_p`'s tie scale. A null whose entries are all the same inf is also constant.
- This applies even when the observation differs from the constant null. The observed arrangement is
  always in the permutation support (the identity, or a rotation for the circular statistics), so a
  constant sampled null with a different observation arises only by chance from a sparse null. The
  +1-rule p, 1/(B + 1), would then be valid, but GIPPA reports NA (decisions 20 and 22).
- The rule is conservative: it never adds a rejection. A p of 1 from a degenerate null would make BH
  more conservative, not inflate the FDR. Reporting NA avoids enlarging the BH family (m) without
  adding information.
- Examples of degenerate nulls:
  - a replicon that is one TU (the TU-block null is the identity);
  - a circular replicon with two TU blocks (both block orders are rotations of each other, and both
    statistics are invariant to rotation);
  - all genes detected or none detected, with every p equal.

**Null summaries.** The table reports:
- `null_mean`, the mean of the valid null;
- `null_lo` and `null_hi`, its linear 2.5 % and 97.5 % quantiles.

For a constant null all three equal the constant. If the observation is NaN or no null value is valid,
every field except `value` and `n_valid` is NaN.

**V10 envelope.** The envelope in V10 plots these quantiles. Two points about it:
- It is **pointwise**: each lag has its own 95 % band, not a simultaneous band over the 20 lags.
- It is on the **statistic's scale**: I_k for Moran's I and the Pearson statistic S_k for the join
  count, not the raw count BB_k.

Lag statistics get no confidence interval (ADR-0020 decision 4, flag `ci_unestimated`).

### Block length

`choose_block_len` implements ℓ̂ from §7.3.5(d) step 5:
- The Moran's I of the TU scores is tested at lags k = 1 … min(20, max_lag) under B gene-label
  permutations. One permutation matrix is shared by all lags, so the draws depend only on B and n.
- ℓ̂ is the first lag whose raw two-sided p is above 0.05 or NaN. This is a stopping rule, not a
  reported test, so the p is not BH-adjusted.
- If every lag is significant, ℓ̂ = min(20, max_lag + 1, n). With n < 3, ℓ̂ = 1.

The block length that callers actually use is ℓ' = `block_len_prime(ℓ̂, T)` = min(20, 2ℓ̂ + 1, T)
(ADR-0021); see docs/methods/bootstrap.md "Block length".

## Statistician sign-off

Signed off by the statistician (PR #116 re-review). The review covered "Statistics and nulls" and
the Pearson residual implementation, after ADR-0020 Amendment 1 and lead decisions 17–22.
