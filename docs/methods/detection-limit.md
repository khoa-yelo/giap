# Detection-limit model (P(detect | abundance) and P(miss))

Implements `gippa.operons.detection_limit` (GIPPA_SPEC.md §7.4.3 "One detected, one not", §7.3.3
conventions, §8.7; issue M6-02, docs/specs/M6-02-detection-limit-model.md). Consumer: M6-05H.

**Status:** method reviewed by the statistician before implementation (2026-10-21). The estimator
below replaces the naive logistic fit that the spec's wording suggests; the model, its features,
its training rows and its outputs are the spec's.

## Why the naive fit is wrong

Both training sets in §7.4.3 are **selected on the response**:

- the peptide model is trained on the observable peptides of *detected* genes, and a gene is
  detected exactly when at least `k_min` of those peptides were seen;
- the replicate model is trained on genes seen in *some but not all* runs, i.e. with
  1 ≤ S ≤ R − 1 per-run detections out of R.

An ordinary logistic regression on such rows estimates P(seen | features, gene selected), not
P(seen | features). At low abundance a detected gene has, almost by definition, one peptide seen
out of several, so the per-peptide rate looks flat in abundance: the abundance coefficient is
attenuated and the intercept is pushed up. The fitted P_det is then too high at the low end, which
is exactly where the undetected partner lives, so **P(miss) is too small (anti-conservative)**.

At R = 2 the replicate training set is degenerate: every retained gene has S = 1 of 2, the response
is 1/2 for every gene whatever its abundance, and the logistic MLE is the constant 1/2 (slope
exactly 0). No estimator can recover the slope from those rows: P(S = 1 | 1 ≤ S ≤ 1) = 1.

## Estimator: conditional (truncated) likelihood

Same model, same rows, the likelihood conditioned on the selection event (the conditional
likelihood of capture–recapture with covariates; Huggins 1989, Biometrika 76:133–140).

**Peptide model.** For a training gene g with observable peptides j, linear predictor
η_gj = β₀ + β_a·a_g + β_L·L_j + β_G·GRAVY_j, p_gj = σ(η_gj), and P_det,g = P(X_g ≥ k_min) the
Poisson-binomial survival over g's peptides (the quantity `p_detect` already returns):

    ℓ_c(β) = Σ_g [ Σ_j ( o_gj·η_gj − log(1 + e^η_gj) ) − log P_det,g(β) ]

Score: Σ_g Σ_j ( o_gj − P(O_gj = 1 | X_g ≥ k_min) )·x_gj, with
P(O_gj = 1 | X_g ≥ k) = p_gj · P(X_g,−j ≥ k − 1) / P(X_g ≥ k); for k_min = 1 this is p_gj / P_det,g.

- The conditioning event must be the event P_det models: **training genes = genes with finite
  `abundance_log2` and at least `k_min` distinct observable peptides seen in C11** (any route, any
  uniqueness, Q2). A detected gene whose evidence lies only outside the observable set (missed
  cleavages, length outside 7–30) contradicts the conditioning and is not a training gene; count it
  in the detection summary. For the detected genes of the spec that satisfy the event this is the
  spec's training set.
- Genes with exactly `k_min` observable peptides contribute 0 (all seen given detection); they are
  still counted in `n_fit`, which stays "rows of the training genes".
- Fit unpenalised (Newton or BFGS with the analytic score; centre and scale the features for the
  optimiser, report `coef` on the raw scale). **Do not use a cross-validated ridge on the naive
  log loss**: shrinking β_a towards 0 is the anti-conservative direction. If the optimiser fails
  to converge or the data are separated, `usable = False` (spec degenerate row "separated
  logistic"). A fixed, negligible ridge purely for numerical safety is acceptable only if the
  calibration criteria below still pass.
- Cross-fitting, folds and `model_for_tu` are unchanged (GroupKFold by TU, folds from
  `gippa.rng.for_step("operons.detection_limit", seed)`).

**Replicate model.** Per-run detection probability q_g = σ(γ·z_g) on the gene's mean
`abundance_log2` plus `detectability.covariates`; S_g of R runs. Training genes as in the spec
(1 ≤ S_g ≤ R − 1); the likelihood is the doubly truncated binomial

    ℓ_c(γ) = Σ_g [ S_g log q_g + (R − S_g) log(1 − q_g) − log(1 − (1 − q_g)^R − q_g^R) ]

(the binomial coefficient is constant). It is identified for R ≥ 3. **At R = 2 the replicate model
is not fitted and `usable = False`** with reason `replicate_two_runs`: under the spec's training
set its MLE is the constant model, which spec item 6 already declares unusable; the code must
decide this from R, not from a numerically small slope. When runs differ per gene (a run in which
the gene could not be quantified), R is the gene's own number of runs.

`p_detect` for the replicate model is the per-run probability q (spec item 4). The partner is
undetected in the pooled table, i.e. in all R runs, so 1 − q overstates P(miss) relative to
(1 − q)^R under run independence: the spec-literal choice is **conservative**. Changing it is an
ADR (see Open decisions).

## Evidence (statistician simulation, 2026-10-21, Slurm job 1864733)

Calibration design of `tests/calibration/test_detection_limit_calibration.py` (truth
logit = −10 + 0.5a − 0.03L + 0.8·GRAVY, a ~ N(15, 2.5), 5–20 observable peptides per gene),
R = 100 replicates per row, scored on an independent 4,000-gene genome. CITL = calibration in the
large, the intercept of d ~ offset(logit p̂) + c.

| genes | k_min | fit | β_a (0.5) | β₀ (−10) | slope | intercept | CITL | mean p̂ − mean d |
|---|---|---|---|---|---|---|---|---|
| 600 | 1 | naive | 0.349 | −7.15 | 1.067 | −0.951 | −0.889 | +0.138 |
| 600 | 1 | conditional | 0.501 | −10.01 | 1.000 | +0.024 | +0.026 | −0.004 |
| 600 | 2 | naive | 0.316 | −6.40 | 1.141 | −1.463 | −1.358 | +0.194 |
| 600 | 2 | conditional | 0.499 | −9.98 | 1.005 | +0.002 | +0.005 | +0.000 |
| 1500 | 1 | naive | 0.351 | −7.20 | 1.065 | −0.960 | −0.899 | +0.139 |
| 1500 | 1 | conditional | 0.502 | −10.04 | 1.001 | +0.009 | +0.011 | −0.002 |
| 1500 | 2 | naive | 0.320 | −6.47 | 1.134 | −1.447 | −1.348 | +0.191 |
| 1500 | 2 | conditional | 0.504 | −10.08 | 0.993 | +0.017 | +0.024 | −0.002 |

The true coefficients give slope 1.00 and |CITL| ≤ 0.01 on average. The naive fit **passes** the
slope gate [0.8, 1.2] in every row while over-predicting detection by 14–19 points.

Replicate model, truth per-run logit = −8 + 0.5a, a ~ N(16, 2), 600 genes, R = 100:

| runs | naive slope (0.5) | doubly truncated, mixed genes | zero-truncated, all seen genes |
|---|---|---|---|
| 2 | 0.000 (exactly) | not identified | 0.498 |
| 3 | 0.149 | 0.515 | 0.502 |
| 4 | 0.245 | 0.500 | 0.498 |
| 6 | 0.361 | 0.501 | 0.499 |
| 8 | 0.417 | 0.501 | 0.502 |

Propagation to P(miss) under the operon null (900 OP pairs, a_u = a_v + Δ, Δ ~ N(0, σ²), the r
sample = oriented Δ of both-detected pairs; opportunities = (v detected, u) ordered pairs):

| σ | k_min | fit | mean P(miss) | observed undetected share | P(miss) ≤ 0.05 among one-detected | P(miss) ≤ 0.10 among one-detected |
|---|---|---|---|---|---|---|
| 0.6 | 1 | true β | 0.301 | 0.300 | 0.013 | 0.035 |
| 0.6 | 1 | naive | 0.199 | 0.300 | 0.036 | 0.097 |
| 0.6 | 1 | conditional | 0.302 | 0.300 | 0.014 | 0.035 |
| 1.5 | 1 | true β | 0.304 | 0.310 | 0.018 | 0.047 |
| 1.5 | 1 | naive | 0.206 | 0.310 | 0.046 | 0.119 |
| 1.5 | 1 | conditional | 0.304 | 0.310 | 0.019 | 0.048 |

The naive fit shrinks P(miss) by a third and multiplies the rate of small P(miss) under the null
by 2.5–3. The selection of the r sample on both-detected pairs was negligible in this design (sd of
the r sample 0.60 / 1.48 against 0.6 / 1.5). Scripts and raw output:
`$GIPPA_ROOT/stat/m6-02/` (sim_pep.py, sim_rep.py, sim_pmiss.py, out/).

## Calibration criteria

The slope alone has no teeth here: unlike M4-13 (trained on the population it scores, so a
logistic MLE with an intercept is calibrated in the large by its score equation), this model is
trained on a selected subset and scored on everyone. Proposed for DL1 (Khoa's decision; the
tolerance file is protected): mean CITL over R within ±0.2 in addition to the slope in [0.8, 1.2].
The conditional fit gives |mean CITL| ≤ 0.03; the naive fit −0.89 to −1.36. Until a key exists,
DL1 records the CITL.

## P(miss) is a probability, not a p-value

M6-index decision 1: no q-value, no BH family. For a pair whose detected partner has abundance
a_v, P(u undetected | a_v, H0) = P(miss); so P(P(miss) ≤ α and u undetected) ≤ α per
opportunity **if P(miss) is correctly estimated** — the level is valid over opportunities, not
conditionally on the pair being in the one-detected arm. An underestimated P(miss) breaks this:
a pair with true P(miss) 0.2 estimated at 0.04 contributes 0.2 to the level at α = 0.05.

## Owed to M6-09 (attack C)

- Report the level with both denominators: per opportunity (v detected) and per one-detected
  pair. Only the first is a validity property; a 2 × α trigger on the second can fire under a
  correct model (or stay silent under a biased one, as above at σ = 0.6).
- Record the calibration in the large of P(miss): mean P(miss) over opportunities against the
  observed undetected share (the conditional fit matches to 0.01; the naive fit is 0.10 low).
- Run the arm with k_min = 2 (`two_peptide`) and with a replicate-only input at R = 3 and R = 2
  (the latter must give `no_detection_limit_model`).
- Measurement error with selection in `abundance_log2` of low-abundance detected genes (estimated
  from one or two observed peptides, biased upwards) is not in these simulations; it pushes the
  same direction and needs an M6-03 generator in which abundance is estimated, not given.

## Open decisions (ADRs for Khoa)

1. Replicate model at R = 2: accept `usable = False` (inside the spec), or widen the training set
   to every gene seen in at least one run with the zero-truncated binomial likelihood, which is
   identified at R = 2 (slope 0.498 above). The second changes the spec's training set.
2. Replicate `p_detect`: per-run q (spec, conservative) or pooled 1 − (1 − q)^R (correct under
   run independence, anti-conservative if detections are over-dispersed across runs).
3. The DL1 CITL tolerance (above).

## Statistical review of the implementation (PR #157, head cf164a5, 2026-10-21)

Simulations ran the code under review: Slurm jobs 1865238 and 1865261. The scripts and output are in
`$GIPPA_ROOT/stat/m6-02/review/` (review_sim.py, review_sel.py, out/). R = 100 replicates per arm,
60 for the selection arm. The values below are means, with the Monte Carlo SE in brackets.

**Verdict: blocked on one item (B1). Everything else matches this note.**

### Matches the note

- **Likelihoods.** The peptide conditional likelihood and its score are correct, including the
  leave-one-out factor P(X_−j ≥ k − 1) for k_min = 2. The same holds for the doubly truncated binomial and
  its score. Since cf164a5, `_newton` returns only at a stationary point with −H positive definite.
- **Peptide model on the operon null.** 900 OP pairs, opportunities = (v detected, u). The level is
  per opportunity.

| σ | k_min | β_a (0.5) | gene slope | CITL | mean P(miss) | undetected share | level α=0.05 | level α=0.10 |
|---|---|---|---|---|---|---|---|---|
| 0.6 | 1 | 0.502 (0.002) | 0.998 | +0.006 (0.008) | 0.302 | 0.304 | 0.004 | 0.011 |
| 0.6 | 2 | 0.504 (0.002) | 0.992 | +0.015 (0.014) | 0.434 | 0.443 | 0.003 | 0.009 |
| 1.5 | 1 | 0.503 (0.002) | 1.000 | +0.012 (0.009) | 0.305 | 0.313 | 0.006 | 0.016 |
| 1.5 | 2 | 0.504 (0.002) | 0.999 | +0.017 (0.014) | 0.420 | 0.459 | 0.005 | 0.014 |

  At σ = 1.5 with k_min = 2, P(miss) falls 0.04 below the observed share. This gap comes from the
  r sample, which is selected on both-detected pairs; the model itself has CITL +0.017. It belongs to
  M6-09 attack C, not to M6-02.
- **Cross-fitting is honest.** The folds are a seeded partition of the training TUs. `fold_coef[j]` never
  saw fold j, and a TU that trained nothing gets the full model, which never saw it either. An
  undetected partner is never a training gene in either model. M6-05H must call `model_for_tu`
  with the TU of the **detected** partner v. That is the only gene of the pair that can have trained.
- **Replicate model at R ≥ 3.** The slope is 0.504 (0.004) at R = 3 and 0.499 (0.003) at R = 4,
  against a truth of 0.5. At R = 2 the model is decided from the run count. With the per-run q, P(miss)
  is very conservative: mean 0.45 against an observed undetected share of 0.15, and level 0 at every α.
  This is as the note predicted.
- **Calibration tests.** The DL1 slope gate has no teeth against the naive fit. The teeth are elsewhere:
  - DL2 (halved slopes);
  - the every-commit β_a recovery test (±0.08; the naive fit gives 0.35);
  - the slow unit CITL test (|CITL| ≤ 0.1 on an independent genome; the naive fit gives −0.89);
  - the k_min = 2 recovery test (±0.04; the k = 1 score gives 0.44).

  These tests are adequate for the failure this note predicted. They cannot detect B1, because in
  `simulate_peptides` `detected` is exactly X_obs ≥ k_min.

### B1 (blocking): the training genes are selected on `detected`, not on the conditioning event

`fit_peptide_model` keeps `called & finite(abundance)` and then `sum(o) >= k`. The likelihood
conditions on X_obs ≥ k_min. In the `mayu` and `decoy_picked` tiers, `detected` is a gene-level q ≤ α.
That is a strict subset of the conditioning event: single-hit genes fail it at random. Genes that
fail are missing from the training set, yet the likelihood treats them as possible. Because they are
low-abundance genes with few peptides, the bias that this note's estimator was chosen to remove
comes back, in the same direction (anti-conservative).

Attack arm (review_sel.py): the pair design above with σ = 0.6. Every gene with ≥ 2 peptides seen is
called, and a gene with one peptide seen is called with probability π₁.

| π₁ | X ≥ 1 genes not called | fit | β_a | CITL vs X ≥ 1 | level α=0.05 (event: u not called) | level α=0.10 |
|---|---|---|---|---|---|---|
| 0.5 | 22 % | code (`detected`) | 0.447 | −0.40 (0.01) | 0.022 | 0.049 |
| 0.5 | 22 % | training = X ≥ 1 | 0.503 | +0.01 (0.01) | 0.014 | 0.033 |
| 0.2 | 36 % | code (`detected`) | 0.400 | −0.75 (0.01) | **0.046** | **0.092** |
| 0.2 | 36 % | training = X ≥ 1 | 0.503 | +0.01 (0.01) | 0.023 | 0.050 |

The fix is in the implementation and the tests, which are not mine to edit:

- Training genes = finite `abundance_log2` and ≥ k_min distinct observable peptides seen in C11, whatever
  `detected` says. This is decision 5's wording and this note's; spec item 2's "detected" is the special
  case. M4-11 gives finite abundances to undetected genes that have accepted peptides.
- The test-author adds a case in which a gene with ≥ k_min observable peptides seen but
  `detected = False` contributes its rows to `n_fit`, and a recovery test with π₁ < 1.

If the lead reads spec item 2 as binding, this needs a Khoa ADR. Until one exists, it stays blocking.

**B1 re-check (bac9046, Slurm 1865546, R = 60).** B1 is fixed and the block is lifted. The
selection arm against bac9046 reproduces the "training = X ≥ 1" rows to every printed digit:

| π₁ | β_a | CITL vs X ≥ 1 | level at α = 0.05 (event: u not called) | level at α = 0.10 |
|---|---|---|---|---|
| 0.5 | 0.503 | +0.011 | 0.014 | 0.033 |
| 0.2 | 0.503 | +0.011 | 0.023 | 0.050 |

The two fits, with `detected` as called and with `detected` = X ≥ 1, are identical in both arms.
N1 (now spec decision 9) stands for M6-05H.

### Non-blocking

- **N1: the event P(miss) models.** Even with B1 fixed, P_det models X_obs ≥ k_min, not `detected`.
  Against the event "u not called", P(miss) is too small in the large: CITL −0.84 (π₁ = 0.5) and
  −1.39 (π₁ = 0.2). The per-opportunity level still held in these arms: at most 0.023 at α = 0.05 and
  0.050 at α = 0.10. M6-05H must define the one-detected arm by u having < k_min observable peptides
  seen, or ADR the gap. This is owed to Khoa (M6-05H) and to M6-09 attack C.
- **N2: stray run label.** With 5 % of the run labels dropped to `""`:
  - 2 real runs give R = 3 and the replicate model is fitted instead of `replicate_two_runs`;
  - with 3 real runs the slope is 0.473 (0.003) against 0.5.

  P(miss) stays conservative under the per-run q: level 0. A C16 validation item, not blocking.
- **N3: unstored covariate median.** At scoring, a NaN covariate gives NaN, so the pair is not
  tested. This costs power, not validity. Store the median. If M6-05H passes `cov=None`, every
  replicate P(miss) is NaN, so M6-05H must pass the M4-10 covariates.
- **N4: replicate model under `two_peptide`.** A run counts as seen with ≥ 1 peptide, while detection
  needs ≥ 2 distinct peptides pooled. P(miss) is then 0.41 against an observed share of 0.46 (R = 3,
  1–4 peptides per gene), still level ≤ 0.003 at α = 0.10 because 1 − q is so conservative.
  Record it in the summary; it only matters if the pooled 1 − (1 − q)^R ADR is accepted.
- **N5: imputed abundances.** A detected gene with an imputed `abundance_log2` (no intensity) trains
  with a fill value. Excluding `abundance_imputed` rows from training would be cleaner; the current
  behaviour is measurement error in the direction this note already lists.
- **N6: the DL1 docstring** still describes the naive oracle as "recorded, not gated" with intercept
  −0.97. It should describe the conditional fit.
- The other backlog nits are numerical or test gaps, not statistical issues: NaN `tu_id` pooled into
  one fold group, rng consumption order, I/L, a constant feature giving a singular H.

### ADRs owed to Khoa

1. The DL1 CITL tolerance of ±0.2 (`tolerances.yaml` is protected). The measured |CITL| is ≤ 0.02 at
   k_min = 1 and 2. A conditional-fit CITL record should also be added at k_min = 2 (DL1 runs only
   k_min = 1).
2. If B1's fix is read as a change to spec item 2: the training set is the conditioning event, not
   `detected`.
3. The one-detected event for M6-05H (N1).
4. The two existing open decisions above: R = 2 zero-truncation, and the pooled replicate `p_detect`.

## References

Huggins RM (1989). On the statistical analysis of capture experiments. Biometrika 76:133–140.
