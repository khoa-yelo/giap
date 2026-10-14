# MAYU-style gene-level FDR (closed form)

Implements `giap.stats.mayu.mayu_gene_fdr` (GIAP_SPEC.md §7.3.2 "MAYU-style closed form (G7)",
§8.3, Appendix B.3; issue M4-02). Used by the `mayu` detection tier (M4-12): when the input has no
pairable decoy rows, or when more than 20 % of decoys are unpaired (ADR-0009 keeps this as an upper
guard; below it the picked tier rescales decoy wins by 1/f).

## What is estimated

Peptide-level q-values do not control gene-level error (§8.7): false peptides pile up on genes that
have no true peptide, and each such hit creates a whole false gene. Given F expected false peptides
among the accepted ones, MAYU (Reiter et al. 2009) estimates how many of the genes with ≥ 1 accepted
peptide owe their detection to false peptides alone.

## Model

1. **Where false peptides land.** Each false peptide independently picks one eligible gene, gene g
   with probability a_g / Σ_h a_h, where a_g = `n_observable_peptides` (floored at 1). This is
   "false peptides spread uniformly over the database's observable peptides": a random match is
   equally likely to be any observable peptide, so a gene with twice the observable peptides takes
   twice the wrong hits.
2. **Size bins.** Genes are split into `mayu_size_bins` equal-count bins by rank of a_g (ties broken
   by input order; fewer bins when there are fewer genes). Within bin b every gene gets the bin's
   mean share w_b = mean(a_g ∈ b) / Σ_h a_h, so Σ_b N_b w_b = 1. Rank bins are deterministic and
   always non-empty; with many ties (e.g. all sizes equal) they give equal w_b and the one-bin
   answer.
3. **Per-gene hit probability.** With F independent false peptides, a given gene in bin b is missed
   by all of them with probability (1 − w_b)^F, so p_b = 1 − (1 − w_b)^F (binomial, sampling with
   replacement). F is an expectation and need not be an integer; the power is evaluated directly so
   w_b = 1, F = 0 gives p_b = 0.
4. **Fixed point.** False genes can only arise among genes with no true peptide. Let M_b be that
   number in bin b. Then E[false genes] = M_b p_b and E[undetected genes] = N_b − n_b = M_b (1 − p_b).
   Eliminating M_b: FP_b = (N_b − n_b + FP_b) p_b, i.e.

       FP_b(t) = (N_b − n_b(t)) · p_b / (1 − p_b),     FDR(t) = Σ_b FP_b(t) / Σ_b n_b(t)

   (FP_b = N_b − n_b when p_b = 1; FDR clipped to 1; NaN when no gene passes). Under the model this
   is a ratio of unbiased estimates: E[FP_b] equals the expected number of false genes given M_b.
5. **F(t).** P(t) = accepted peptides with q ≤ t. With PEPs, F(t) = Σ_{P(t)} PEP (the expected
   number of false peptides if PEPs are calibrated). Without PEPs, F(t) = t·|P(t)|, the q-value's
   bound on the false count. When q is the running mean of sorted PEPs the two agree at every
   threshold (the `no_pep` attack below is numerically identical to the PEP run). PEPs are used for
   all peptides or none: all-NaN means "no PEPs"; partly-NaN raises (the tier decides once per run).
6. **Gene q.** Thresholds are the distinct finite best-peptide q values; q_g = min over t ≥ best_q_g
   of FDR(t) (`monotone_q` from `giap.stats.tda`). Detected at q_g ≤ 0.01. This is one test family
   per run (all eligible genes), matching §8.1's "Gene detected" row.
7. **≥ 2-peptide variant.** The chance that a gene gets ≥ 2 false peptides is approximated by
   Poisson(λ = F·w_b): p_b⁽²⁾ = 1 − e^(−λ)(1 + λ), plugged into the same FP formula with the same
   n_b, N_b. Reported as `fdr_2pep` in the curve and in `detection_summary.tsv` only; it is not a
   valid FDR for the set of ≥ 2-peptide genes (n_b still counts ≥ 1-peptide genes), only an
   indication of how much of the error a two-peptide rule would remove.

### Approximations

- **Binomial vs hypergeometric.** Reiter's MAYU draws false protein IDs from the database without
  replacement (hypergeometric within a size bin). Here false peptides choose genes with
  replacement; several false peptides on one gene are allowed and count once. For F ≪ N_b the two
  coincide; for large F the binomial form is the exact one for independent random matches.
- **Bin averaging.** p_g = 1 − (1 − w_g)^F is concave in w_g, so within a bin the mean of p_g is
  ≤ p(mean w). The bins also mix genes whose absence probabilities differ (short genes are more
  often absent). Both effects are small with 10 rank bins (simulation below: estimated / realised
  loose-threshold FDR within ±3.5 % for log-SD 0.6–1.2 sizes); a single bin overstates FDR by ≈ 22 %
  when absent genes are shorter than average (one bin: est. 0.0323 vs realised 0.0264).
- **Expectation, not a bound.** FDR(t) estimates E[false]/n; it carries no +1 and no finite-sample
  guarantee. Its calibration rests on the model and on F being right.

## Assumptions and when they fail (direction)

| Assumption | Fails when | Direction |
|---|---|---|
| False hits ∝ a_g | wrong hits concentrate on long genes beyond a_g (e.g. ∝ a_g²; long, well-covered proteins attract more near-miss matches) | **conservative** (realised / est. ≈ 0.78): more false peptides land on genes that are present anyway |
| | wrong hits spread more evenly than a_g (∝ a_g^0.5 or uniform; e.g. semi-specific or unspecific search while a_g counts tryptic peptides) | **anti-conservative**: uniform hits give realised / est. ≈ 1.24 and FDP 0.0122 at q = 0.01; a_g^0.5 ≈ 1.14 (R = 100) |
| a_g is the right weight | a_g is a noisy proxy of the true hit weight (hits ∝ length, a_g = Poisson(length/25), zeros floored) | mildly anti-conservative (≈ 1.05; FDP 0.0102 at q = 0.01, R = 300): genes with a_g overstated by noise are more often detected, so the undetected pool is under-weighted |
| a_g ≥ 1 floor | genes with no observable peptide | conservative (they cannot really be hit) |
| F is right | PEPs miscalibrated (too optimistic), or q-values anti-conservative upstream (Wen 2025, §8.7) | FDR scales roughly linearly with F; an underestimated F is the main real-world risk and is not detectable here (validation: entrapment, §15.5) |
| Independent hits | false peptides cluster (several wrong peptides from one homologous or contaminant protein onto one gene) | conservative: fewer distinct false genes than the independent model predicts |
| "No true peptide" = undetected | present genes whose true peptides were all rejected get hit by a false peptide | counted as false by the estimator, as true by an entrapment truth: conservative |
| Many single-peptide genes (one-hit wonders) | not a violation: n_b counts genes with ≥ 1 peptide and the model is per gene | at level (attack `one_hit_wonders`); the share of single-peptide genes is what `fdr_2pep` speaks to |

The anti-conservative cases need a hit distribution flatter than the weights. `n_observable_peptides`
comes from the configured enzyme (§7.3.3); when the input search was semi-specific or unspecific the
tier should say so in the manifest (recommendation for M4-12: warn and report `fdr` with weights
a_g ∝ length as a sensitivity row).

## Relation to the picked tier and ADR-0009

Picked target–decoy (Savitski 2015) estimates gene FDR from decoy wins, with a +1 and a finite-sample
argument; MAYU estimates it from a hit model and F. Picked is preferred whenever decoys are
available. ADR-0009 rescales picked decoy wins by 1/f for unpaired decoys and keeps the > 20 %
unpaired fallback to MAYU as an upper guard: beyond that share the paired decoys may no longer be a
random sample of all decoys, and MAYU's model-based estimate is the safer one. MAYU's own weak spot
(hit distribution flatter than a_g) is independent of the decoy-loss problem, so the fallback does
not compound errors. M4-12's `gene_fdr` suite calibrates both tiers (± PEP) and the decoy-loss
cases.

## Simulation evidence

`tests/calibration/test_mayu_calibration.py` (test-author): under the exact model (two size groups
aligned with bins), FP is unbiased, the pooled hit rate matches p_b, and mean FDP ≤ mean FDR + 2 SE.

### Attack nulls (`tests/calibration/attacks/test_attack_mayu.py`)

G = 3,000 genes, lengths log-normal (median 300 aa), a_g = round(length/25) ≥ 1, presence logistic
in log a_g (absent genes are short), Poisson(0.3·a_g) true peptides per present gene, calibrated
PEPs (each peptide false with probability PEP), accepted at peptide q ≤ 0.01, 10 bins. Criterion
§8.4: mean FDP ≤ q + 2 SE. Seed 20261014. Mean FDP at gene q ≤ 0.01 / 0.05 / 0.10 (the loosest
estimated FDR is ≈ 0.026–0.033, so 0.05 and 0.10 accept every gene):

| Attack | R = 50 | R = 300 | realised / est. (R = 300) |
|---|---|---|---|
| `detectability_bias` (base) | 0.0095 / 0.0267 / 0.0267 | 0.0098 / 0.0266 / 0.0266 | 1.005 |
| `heavy_tail` (log-SD 1.2) | 0.0092 / 0.0318 / 0.0318 | 0.0095 / 0.0321 / 0.0321 | 0.982 |
| `one_hit_wonders` | 0.0053 / 0.0053 / 0.0053 | 0.0053 / 0.0053 / 0.0053 | 0.997 |
| `long_gene_concentration` (∝ a_g²) | 0.0076 / 0.0202 / 0.0202 | 0.0076 / 0.0206 / 0.0206 | 0.782 |
| `no_pep` (F = t·\|P\|) | as base | as base | 1.005 |

All stay at level (SE ≤ 0.0004 at q = 0.01, 0.0007 at the looser levels). Note the base case sits
at 0.0098 for q = 0.01: the estimator is essentially unbiased, with no conservative margin.

Violations (direction tests, not level gates): `uniform_hits` realised / est. 1.24 (FDP 0.0122 at
q = 0.01), `noisy_size` 1.05 (FDP 0.0102 at q = 0.01; at the 2-SE edge).

## Numerical notes

- 1 − (1 − w)^F loses relative precision only for w ≲ 1e-9 (relative error ~1e-7 at w = 1e-9,
  F = 0.01); real genomes have w ≥ 1/Σa ≈ 1e-6. `-expm1(F·log1p(−w))` with a w = 1 guard would be
  exact; not needed.
- p⁽²⁾ is clipped to ≤ p (always true mathematically: Poisson P(≥ 2) ≤ Poisson P(≥ 1) ≤ binomial
  P(≥ 1)).
- Not validated by the function: PEPs outside [0, 1] (a negative PEP lowers F silently), non-finite
  sizes (NaN sizes give an all-NaN q), non-integer `n_bins`. Callers (M4-12) build these from
  validated tables; input validation would be cheap and is recommended.
- No effect size or CI: gene q is an error-rate estimate, not a test of an effect. The curve's `fp`
  and `n_genes` (the estimated false and total detected genes) are what `detection_summary.tsv`
  should report next to the q threshold; no CI is available for a plug-in model estimate, which is
  why the calibration suites carry the evidence.

## Sign-off

Statistician (M4-02 review, 2026-10-14): **approved** for `mayu_gene_fdr` as implemented in
`src/giap/stats/mayu.py`: formula matches §7.3.2 and B.3, bins and monotone q are correct, at level
under the model and the attacks above. Forwarded to M4-12: the flat-hit-distribution direction
(semi-/unspecific searches), the base case's zero margin at q = 0.01, and input validation of PEPs
and sizes.
