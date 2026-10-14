# De novo precision, estimated from the genome

Implements `giap.mapping.denovo_precision` (GIAP_SPEC.md §7.2.5 "De novo precision, estimated from
the genome" and step 4, §8.1, §8.4, §9.3, §9.4; issue M3-11; ADR-0010).

**Status:** awaiting statistician sign-off (§14.1).

## Estimator

De novo peptides that enter mapping are split into score bins b and length classes l (M3-10 labels
`7`, `8`, `9`, `10-11`, `12-14`, `15-19`, `20-40`). Per bin and class:

- m_bl: share matching the six-frame target translation **exactly** (I = L);
- r_bl: share matching the shuffled decoy exactly;
- c_l: share of reference peptides (`--reference-peptides`, targets with q ≤ 0.01) matching the
  target exactly. Without reference peptides c_l = 0.99 and every row carries `c_assumed` (also a
  run flag). With reference peptides, a class that has none gets c_l = NA.

A peptide is correct with probability p_bl; a correct one matches the target with probability c_l,
a wrong one with the decoy-match probability r_bl. So m_bl = p c_l + (1 − p) r_bl, and

    prec_b = Σ_l w_bl (m_bl − r_bl) / (c_l − r_bl),   clipped to [0, 1],

with w_bl the share of bin-b peptides in class l.

## Bins

Deciles of `score`: edges `np.quantile(scores, [0.1, …, 0.9])`, right-closed; coinciding edges merge
bins (tied scores share a bin). Peptides without a score get bin −1 and no row.

## Reportability

The estimator is unstable when c_l − r_bl is small. A class with c_l − r_bl < 0.05 (or c_l = NA) is
dropped from that bin and the remaining weights renormalised (c − r = 0.05 is kept). A bin with no
class left has precision NA and flag `unreported`. `c_value` = Σ_l w_bl c_l over the classes used.

## Bootstrap

1,000 resamples of the bin's peptides with replacement (drawn as one multinomial over the
class × target-match × decoy-match cells, which is the same distribution and is invariant to row
order), c_l held fixed; the same reportability rule is applied to each resample. CI = 2.5 % and
97.5 % percentiles. RNG: `giap.rng.for_step("placement", seed)`. A bin with n = 1 yields a degenerate
CI and is reported (n is visible).

## π̂ for M3-10

π̂_wrong (m = 0, de novo) = 1 − ci_low of the peptide's bin: the conservative bound (ADR-0010; a
point-estimate π̂ that is too small breaks placement FDR). Unreported bins, bins without a row and
unscored peptides get π̂ = 1.

## Calibration

`tests/calibration/test_denovo_precision_calibration.py`: known p_b per bin, two classes; per-bin
|bias| ≤ `denovo_precision_max_abs_bias` (0.05), 95 % CI coverage, and P(π̂ ≥ π) ≥ 0.975 within a
binomial band. Genome-level validation with giapsim de novo errors is M3-21.

## Assumptions and limits

- c_l from reference peptides transfers to de novo peptides of the same length class.
- Wrong de novo peptides match the target no more often than the shuffled decoy (r_bl unbiased for
  the wrong-peptide target-match rate).
- c_l's own sampling error is ignored in the CI (small with hundreds of reference peptides).
