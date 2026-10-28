# Detectability model (cross-fitted p_detect)

Implements `gippa.genetable.detectability` (GIPPA_SPEC.md §7.3.3, §8.2 "Adjust for detectability",
§8.4 row "Detectability", §8.5 attack 2; G5; issue M4-13, docs/specs/M4-13-detectability.md).

**Status:** signed off by the statistician (M4-13 review, 2026-10-15; see Sign-off).

## Model

Response: `detected` of the eligible genes (non-partial protein-coding rows; others get NaN).
Covariates are physicochemical only: log1p(`n_observable_peptides`), log(`length_aa`), GRAVY, pI,
`tm_segments`. Expression determinants and tested categories are never covariates (an included
category would have O − E = 0 by construction). NaN or infinite values take the median of the
transformed column; a zero-variance column is dropped and reported (`dropped_zero_variance`).

Per covariate a cubic `SplineTransformer(n_knots=5, knots="quantile", extrapolation="linear")`
(`tm_segments` enters raw when it has fewer distinct values than knots), then an L2 logistic
regression (`l1_ratio=0`, lbfgs, `tol=1e-6`). Rows are processed in locus_tag order, so the result
does not depend on the input row order.

## Cross-fitting and the choice of C

- Outer folds: the sorted unique TU IDs are permuted with `gippa.rng.for_step("detectability",
  seed)`; fold = position mod 5. Every gene's `p_detect` comes from the model fitted without its TU
  (genes without a TU are their own TU).
- Inner search, per outer fold, on its training TUs only: 3 TU folds from the same generator, C ∈
  {0.01, 0.1, 1, 10} by the mean per-fold held-out log loss; ties go to the smaller C.
- Fewer than 50 detected or 50 undetected eligible genes, or fewer TUs than folds: `p_detect` = the
  overall rate, flag `detectability_constant` (`too_few_tus` is reported and logged too).

Why TU folds. Genes of one TU share abundance and therefore detection; a gene-level fold would let
the model learn a TU's detection state from its neighbours. With physicochemical covariates only
this leakage is small and **cannot be detected statistically** (gene folds gave the same slope as TU
folds in the review simulations); it is guarded structurally by unit tests (example 4; oracle checks
of the chosen C and of p_detect to 2e-3).

## Calibration diagnostics (`detection_summary.tsv`, section `detectability`)

- 10 equal-count bins of p̂ (stable rank order): mean p̂, observed rate, n.
- Brier score.
- Recalibration fit d ~ a + b·logit(clip(p̂, 1e-6)) by unpenalised Newton–Raphson; `ci_low` /
  `ci_high` are Wald 95 % intervals from the inverse observed information X'WX at the MLE. Other
  rows carry no CI. If the data are (quasi-)completely separated in either direction (every
  detected gene's logit ≥ every undetected one's, or every one ≤) the MLE does not exist:
  intercept, slope and CIs are NaN and the `recalibration_separated` row is 1 (0 otherwise; spec
  M4-13 decision 17). These are descriptive: there is no p-value and no test family (§8.1 has no
  detectability claim).

Do not test b = 1 with this Wald interval on cross-fitted p̂: cross-fitted predictions are shrunk
by estimation error, so b ≈ 0.97 ± 0.02 at n = 2,000 even when the model is right, while the
interval half-width is ≈ 0.11. The interval is a statement about the realised p̂, not about the
true model.

## Calibration evidence

Slurm jobs 1848615 (pre-code, spec simulation) and 1848717 (this implementation). Generator: the
M4-07 layout with M4-08 `covariates_only` detection (logistic truth in log length, log1p n_obs,
GRAVY), unless stated.

| Check (suite `detectability`) | Result on the implementation |
|---|---|
| D1 (§8.4 gate, mean cross-fitted slope ∈ [0.8, 1.2]) | slope 0.972 ± 0.019 per replicate at n = 2,000 (200 reps), 0.981 ± 0.011 at 4,500; gate replicates 0–19: mean 0.974. False alarm ≈ 0. Median fit time 4.7 s / 6.9 s single-threaded |
| D2 (bins χ²₁₀, slope, intercept at known p; Bonferroni × 3) | rejection 0.008 / 0.043 / 0.086 at α 0.01 / 0.05 / 0.10 (4,000 reps); components 0.049 / 0.045 / 0.052 at 0.05 |
| Wald CI coverage (known a, b; n 200–5,000; b ∈ {0.7, 1, 1.3}) | 0.944–0.960 for slope and intercept (2,000 reps per cell) |
| Attack 2 subgroup gate, correct model, prophage O − E z | mean −0.35, SD 0.89; rejection 0.005 / 0.039 / 0.089 (1,000 reps) |
| same, shortest-length decile | mean −0.66, SD 0.69; rejection 0.002 / 0.026 / 0.087 |
| Teeth: model without length and n_obs, prophage | z mean −3.60, SD 0.75; rejection at 0.05 is 0.99; gate blocks 20/20 and 20/20 |
| Implementation vs the statistician's independent sklearn simulation of the spec (8 replicates, same data and generator) | same C in 39 of 40 folds; mean |Δp| 0.0007–0.0017, max ≤ 0.028; one replicate with a near-tie inner C (1 vs 10) differs by up to 0.11 in that fold |

False-alarm budget of the suite (probability that some detectability check is a confirmed fail
under a correct implementation, two-block rule, levels 0.01 / 0.05 / 0.10 evaluated together):
fast (R = 20) 0.0066 (D2 0.0035, prophage 0.0023, short decile 0.0008, teeth and D1 ≈ 0); full
(R = 100) 0.0029. The pre-review design (separate D2–D4) spent ≈ 0.023.

## Power and what the gates do not see

- D1 catches gross mis-specification: C fixed at 0.01 (slope 2.5) and logit scaling by 0.7 or 1.4.
  It misses in-sample prediction (slope 1.03–1.05), gene folds, and dropping the size covariates
  (0.93): its marginal slope is blind to subgroup miscalibration.
- The subgroup gate (`tests/calibration/attacks/test_attack_detectability.py`) covers the §8.5
  attack 2 confounding for p_detect itself: dropping the size covariates is a confirmed fail.
- Short genes are over-predicted a little even by the correct model (decile z −0.66; prophage
  −0.35): the penalised spline flattens at the low-length edge. The gate holds because the subgroup
  z is under-dispersed, not because the bias is zero. Region tests and ORA (M4 region tests, M7)
  must run their own attack-2 gate on this p_detect; a region made only of very short genes will
  show a small negative O − E.

## Sign-off

Statistician (M4-13 review, 2026-10-15): **approved** for `gippa.genetable.detectability` as
implemented. Cross-fitting, the inner C search and the recalibration match an independent
simulation of the spec; D1 and the merged D2 hold their level; the Wald intervals cover at the
nominal rate; the attack-2 subgroup gate is at level for the correct model and fails the model
without size covariates; the suite's false-alarm budget is ≈ 0.007 per run. Forwarded: the
short-gene over-prediction to the region-test and ORA attack-2 suites.
