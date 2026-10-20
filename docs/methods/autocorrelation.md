# Autocorrelation along the replicon

This note covers `giap.stats.autocorr`: `morans_i`, `join_count`, `label_permutations`, `tu_blocks`,
`tu_block_permutations`, `lag_test`, `choose_block_len` and `block_len_prime`. Sources: GIAP_SPEC.md
§7.3.5(c), §7.3.5(d) step 5, §7.3.6, §8.3 and Appendix B.4; issue M5-02, lead decisions 4–20;
ADR-0020 with Amendment 1; ADR-0021. The M5-05 step handles the detrending (ADR-0020 decisions 1–2) and
the table layout. The statistician has not signed off yet; see the end of this note.

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
- A gene with p ∈ {0, 1}, with p NaN, or with d missing has no residual: r = 0, so it joins no pair.
- A valid lag with no pair left gives S_k = 0.
- `effect` = S_k − its null mean.

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
- This applies even when the observation differs from the null: `perm_p` would then return 1/(B + 1),
  a p that carries no information.
- A p of 1 from a degenerate null would inflate the BH family, so it is never reported.
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

**Pending.** The statistician reviewed PR #116; that review led to ADR-0020 Amendment 1 and lead
decisions 17–20. Sign-off on this note and on the Pearson residual implementation is still to come.
