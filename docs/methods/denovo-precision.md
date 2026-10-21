# De novo precision, estimated from the genome

Implements `giap.mapping.denovo_precision` (GIAP_SPEC.md §7.2.5 "De novo precision, estimated from
the genome" and step 4, §8.1, §8.4, §9.3, §9.4; issue M3-11; ADR-0010).

**Status:** signed off by the statistician (M3-11 review, 2026-10-14; see Sign-off).

## Estimand (ADR-0019 decision 5)

**Status: statistician sign-off pending** (added with the ADR-0019 implementation, M3-25).

`prec_b` estimates the **correct share among the distinct de novo peptides of GIAP score bin b**.
Two consequences of "distinct" and "GIAP's bins":

- the unit is the distinct peptide string, after §5.3 collapses the rows of a peptide (`score` =
  max over its rows), not the spectrum or the read; a peptide is **correct** when its string
  (I = L) equals the `true_peptide` of at least one of the reads it was collapsed from;
- the bins are GIAP's own score deciles over those aggregated scores (see Bins), not the bins or
  the design parameters of whatever produced the input.

Calibration therefore gates against that realised share, computed per bin from the truth table.
giapsim's `truth/denovo_precision.tsv` is the **per-read design precision** of the error model —
a different quantity, over another unit and other bins — and is not the target of these gates.
The two drift apart whenever correct and wrong reads collapse into peptides at different rates
(a correct peptide is usually read several times; a wrong one once), which is the ordinary case.

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

Equivalence. The estimator depends on a bin's rows only through the 4k cell counts (k classes), so
one Multinomial(n, cell shares) draw has exactly the distribution of the cell counts of n rows
drawn with replacement; the two bootstraps are identical, not approximately so. The class weights
w_bl are resampled too (the bin's class mix is treated as random), which is right for a bin whose
peptides were not stratified by class.

What it ignores (and the direction):
- c_l's sampling error (held fixed). Matters only when a class has few reference peptides; an
  underestimate of c_l is anti-conservative (see Calibration and attacks, `hard_reference`).
- Dependence between rows: nested or overlapping de novo peptides of one locus, or several
  candidates from one spectrum, share correctness and match status. Rows are treated as
  independent, so the CI is too narrow (per-bin lower-bound coverage ≈ 0.70–0.84 with runs of 5
  dependent rows).
- The bin edges are data-dependent quantiles and are held fixed.
- Resamples that become unreportable (every class c_l − r_bl < 0.05) are dropped, so ci_low is
  taken over the reportable resamples only. Near the 0.05 boundary this conditions on a small
  denominator; counting such resamples as precision 0 would be the conservative choice (follow-up,
  not observed in calibration).
- The percentile interval at the clip: when every row in a small bin matches the target, the
  point estimate is clipped to 1 and the bootstrap distribution is degenerate at 1, so ci_low = 1
  and π̂ = 0. Per-bin coverage of the lower bound is therefore poor in small high-precision bins
  (P(1 − ci_low ≥ π) ≈ 0.75 / 0.64 at p = 0.98 / 0.99 with 50 peptides per bin, ≈ 0.5 with 20,
  versus 0.95–1.0 with 400).

## π̂ for M3-10

π̂_wrong (m = 0, de novo) = 1 − ci_low of the peptide's bin: the conservative bound (ADR-0010; a
point-estimate π̂ that is too small breaks placement FDR). Unreported bins, bins without a row and
unscored peptides get π̂ = 1.

Why the per-bin lower bound is enough. ADR-0010 needs W_m(e) = Σ π̂_j over decoy hits not to
underestimate the expected number of wrong target hits. A per-bin undercoverage in a small,
high-precision bin costs at most π_b per decoy hit from that bin (≤ 0.02 here), and the + c term
(c = max π̂ over searched peptides, ≈ 1 whenever any bin has low precision or is unreported) adds a
whole false hit; small lists accept nothing at q = 0.01. The attacks below confirm the placement
FDR holds with 5 peptides per bin.

## Calibration and attacks

`tests/calibration/test_denovo_precision_calibration.py`: known p_b per bin, two classes; per-bin
|bias| ≤ `denovo_precision_max_abs_bias` (0.05), 95 % CI coverage, and P(π̂ ≥ π) ≥ 0.975 within a
binomial band. Genome-level validation with giapsim de novo errors is M3-21.

`tests/calibration/attacks/test_attack_denovo_precision.py` (§8.5) feeds π̂ = 1 − ci_low into the
ADR-0010 rule written out directly, FDR(e) = (Σ π̂ over decoy hits with E ≤ e + max π̂) / T(e),
with the precision estimated from the same exact-match calls the placement uses. Ten bins with
p = 0.30 … 0.99, chance rate r = 0.1. Mean FDP at q = .01 / .05 / .10, R = 300 (SE ≤ .0015):

| Scenario | FDP | Gate |
|---|---|---|
| oracle π (reference) | .0098 / .032 / .032 | – |
| base, N = 4000 | .0083 / .032 / .032 | at level |
| small_bins, 20 per bin | .0033 / .025 / .032 | at level |
| tiny_bins, 5 per bin | .0000 / .013 / .027 | at level |
| chance_by_bin (r 0.25 → 0.05 with p) | .0084 / .043 / .064 | at level |
| clustered (runs of 5 dependent rows) | .0085 / .033 / .033 | at level |
| easy_reference (c_ref .99 > c_true .90) | .0064 / .033 / .035 | at level (conservative) |
| c_true_one (c_assumed .99, true 1.0) | .0086 / .032 / .032 | at level |
| hard_reference (c_ref .85 < c_true .99) | **.0124** / .032 / .032 | direction: anti-conservative |
| near_correct_wrong (wrong hit target 2r) | **.0186 / .062** / .062 | direction: fails with oracle π too (.0195 / .062) |

Directions to remember:
- Reference peptides **easier** than the de novo peptides (c_ref > c_true) bias precision DOWN and
  π̂ UP: conservative. The dangerous direction is reference peptides **harder** than correct de
  novo peptides (c_ref < c_true: the reference search covers variant or contaminant peptides, or a
  class has only a handful of reference peptides): precision is biased up by (c_true − r)/(c_ref − r)
  and π̂ is too low. A shortfall of 0.04 in c (.95 vs .99) stayed at level (.0094); 0.14 did not.
- Wrong de novo peptides that match the target more often than the shuffled decoy (near-correct
  errors that hit a paralog or the true locus exactly) inflate m − r and break the placement FDR
  even with the true π: a decoy-model violation owned by M3-09/M3-10 and measured in M3-21.

## Assumptions and limits

- c_l from reference peptides transfers to de novo peptides of the same length class.
- Wrong de novo peptides match the target no more often than the shuffled decoy (r_bl unbiased for
  the wrong-peptide target-match rate).
- c_l's own sampling error is ignored in the CI (small with hundreds of reference peptides).
- Rows are independent within a bin (see Bootstrap: dependence narrows the CI).

## Sign-off

Statistician (M3-11 review, 2026-10-14): **approved** for `giap.mapping.denovo_precision` as
implemented: the estimator is §7.2.5's (a de novo peptide counts as matching only on an exact
target hit, I = L; r_bl from exact shuffled-decoy hits; c_l per class from q ≤ 0.01 target
reference peptides, entering only the denominator), the multinomial bootstrap is exactly the row
bootstrap, and π̂ = 1 − ci_low keeps the ADR-0010 placement FDR at level in every at-level attack,
including 5-peptide bins. Forwarded to M3-21 (genome-level, giapsim): measure c_ref − c_true and
the wrong-peptide target-vs-decoy match ratio, the two directions that make π̂ anti-conservative;
consider a minimum reference count per class (or resampling c_l) and counting unreportable
resamples as precision 0.
