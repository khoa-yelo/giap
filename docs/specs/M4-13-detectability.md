# M4-13: Detectability model (cross-fitted p_detect) and its calibration diagnostics

Milestone M4 · track D (detectability) · depends on: M4-06, M4-10 (calibration needs M4-07, M4-08) ·
estimate ≈ 300 code lines (+ ≈ 300 test lines) · **can start now**

## Spec sections implemented

- §7.3.3 quoted: "Response. `detected` for eligible genes." Covariates (M4-10). "Expression
  determinants (CAI, promoter, RBS) are excluded so Module 8 is not circular. No tested category is
  ever a covariate." "Model (G5). Per covariate: `SplineTransformer(n_knots=5, degree=3,
  knots="quantile")`. Then: `LogisticRegression(penalty="l2")`, with C chosen from {0.01, 0.1, 1, 10}
  by inner `GroupKFold(3)`." "Cross-fitting. `GroupKFold(5)` with groups = TU (folds from seeded TU
  order). Each gene's `p_detect` comes from the model that did not see its TU." "Calibration check. 10
  equal-count bins, Brier score, and the calibration intercept and slope (logistic regression of d on
  logit p̂), written to `detection_summary.tsv` (section `detectability`) and plot V18." "Edge case.
  With fewer than 50 detected or 50 undetected eligible genes, `p_detect` = the overall rate and the
  flag is `detectability_constant`."
- §5.2 `detectability: {covariates: [n_observable_peptides, length_aa, gravy, pi, tm_segments],
  spline_knots: 5, folds: 5}`.
- §8.2 "Adjust for detectability. Detection tests compare against Σ p_detect from cross-fitted
  predictions, never against raw counts."
- §8.4 "Detectability | logistic truth on physicochemical covariates | calibration slope 1 ± 0.2 |
  20 | 100"; tolerances `detectability_slope: [0.8, 1.2]`.
- §2.4 G5; §6.6 seeds.

## Public interface

`src/giap/genetable/detectability.py`
```python
@dataclass(frozen=True, slots=True)
class DetectabilityResult:
    p_detect: pd.Series             # index locus_tag; NaN for non-eligible genes
    constant: bool                  # detectability_constant
    fold_of_tu: Mapping[str, int]
    chosen_C: tuple[float, ...]     # one per outer fold
    diagnostics: pd.DataFrame       # detection_summary rows, section "detectability"

def transform_covariates(cov: pd.DataFrame, names: Sequence[str]) -> pd.DataFrame: ...
    # log1p(n_observable_peptides), log(length_aa), gravy, pi, tm_segments; NaN → column median
def tu_folds(tu_ids: Sequence[str], k: int, rng: np.random.Generator) -> dict[str, int]: ...
def fit_detectability(detected: pd.Series, cov: pd.DataFrame, tu_of_gene: pd.Series,
                      cfg: DetectabilityConfig, rng: np.random.Generator) -> DetectabilityResult: ...
def calibration_diagnostics(d: ArrayLike, p: ArrayLike, n_bins: int = 10) -> pd.DataFrame: ...
```
- Pipeline per fit: `ColumnTransformer` of one `SplineTransformer(n_knots=cfg.spline_knots, degree=3,
  knots="quantile", extrapolation="linear")` per covariate → `LogisticRegression(penalty="l2",
  C=…, solver="lbfgs", max_iter=2000)`. Integer covariate `tm_segments` with fewer distinct values
  than knots: its spline is replaced by the raw value (prevents duplicate quantile knots).
- `tu_folds`: sort unique TU IDs, permute with `rng` (`giap.rng.for_step("detectability", seed)`),
  fold = position mod k. Outer folds use these; the inner 3-fold C search applies the same rule
  to the training TUs with the same stream. C selection metric: mean held-out log loss; ties → the
  smaller C.
- Constant rule: eligible detected < 50 or undetected < 50 → `p_detect` = detected share for every
  eligible gene; `constant = True`; flag `detectability_constant` on eligible genes (M4-15 sets it).
- `calibration_diagnostics`: 10 equal-count bins of p̂ (rank order, ties by input order) → rows
  `group_type = "bin"`, `group = 1..10`, metrics `mean_p`, `observed_rate`, `n`; plus
  `group_type = "overall"` metrics `brier`, `calibration_intercept`, `calibration_slope` (unpenalised
  logistic regression of d on logit p̂, p̂ clipped to [1e-6, 1 − 1e-6]), `n_detected`, `n_eligible`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 49 detected, 951 undetected eligible genes | `constant = True`; every eligible `p_detect = 0.049` |
| 2 | 50 detected, 50 undetected | model fitted (`constant = False`) |
| 3 | `tu_folds(["TU_1", …, "TU_10"], 5, default_rng(0))` twice | identical maps; each fold has 2 TUs |
| 4 | genes of one TU | always in the same outer fold (no TU straddles train and test) |
| 5 | `calibration_diagnostics(d=[0,0,1,1], p=[0.1,0.2,0.8,0.9], n_bins=2)` | bin 1: mean_p 0.15, observed 0.0; bin 2: mean_p 0.85, observed 1.0; brier = (0.01 + 0.04 + 0.04 + 0.01)/4 = 0.025 |
| 6 | perfectly calibrated synthetic p (d ~ Bernoulli(p), n = 20,000, seeded) | slope within 1 ± 0.05, intercept within 0 ± 0.05 |
| 7 | M4-08 layout with `detection = "covariates_only"`, 4,500 genes | `p_detect` correlates with `p_detect_true` (Spearman ≥ 0.8) |
| 8 | partial gene and non-eligible rows | `p_detect` NaN |
| 9 | same inputs and seed twice | identical `p_detect` to the last bit (TE7) |

## Edge cases

- A covariate with zero variance (e.g. all `tm_segments = 0`) is dropped from the fit (logged).
- Fewer TUs than folds (tiny genome) → constant rule applies with warning `too_few_tus`.
- A training fold with only one class (possible in tiny data) → that fold uses the fold's detected
  share as constant prediction; flag recorded in diagnostics (`metric = degenerate_fold`).
- sklearn convergence warnings are captured into `diagnostics` (`metric = convergence_warning`), not
  printed.

## Acceptance criteria

- [ ] `tests/unit/genetable/test_detectability.py`: examples 1–9.
- [ ] Property: permuting gene row order leaves `p_detect` unchanged (folds depend on TU IDs only).
- [ ] **Calibration** `tests/calibration/test_detectability_calibration.py` (suite `detectability`, R = 20 fast /
  100 full): M4-07 layout (2,000 genes fast, 4,500 full) + M4-08 `covariates_only`; fit; the mean
  calibration slope over R lies in `tolerances.detectability_slope` = [0.8, 1.2]; per-replicate
  slopes, intercepts and Brier scores written to the calibration log.
- [ ] Attack null 2 input check: with `short_genes_in="prophage"`, `p_detect` of prophage genes is
  lower on average (the model sees length) — consumed by M4-17's attack test.
- [ ] Fit on 4,500 genes in < 20 s on 4 CPUs (5 outer × 4 C × 3 inner fits).

## Methods note (statistician sign-off)

`docs/methods/detectability.md`: covariates (with M4-10 definitions), why expression determinants and
tested categories are excluded, spline + L2 logistic (G5), TU-grouped cross-fitting, C selection,
constant fallback, the calibration diagnostics and the simulation evidence (slope distribution).

## Questions for the lead

1. Calibration criterion: the **mean** slope over R replicates in [0.8, 1.2], or every replicate?
   **Proposed default:** the mean (per-replicate slopes on 2,000 genes are noisy, SD ≈ 0.1);
   per-replicate values are logged.
2. Inner C selection by log loss (not AUC)? **Proposed default:** log loss, since calibration is the goal.
3. Spline `extrapolation="linear"` for held-out genes outside the training range? **Proposed
   default:** yes (sklearn's default `constant` flattens tails, hurting short sORFs).

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (non-eligible genes `<NA>` in a boolean `detected`, NaN
   p_detect; fold and inner-split RNG order; C chosen by mean per-fold log loss; tolerances to the oracle
   2e-3 / 1e-3; stable-sort bins; NaN covariates → median of the transformed column).
2. `too_few_tus` and dropped covariates are reported in the returned diagnostics (machine-readable), and
   also logged at WARNING — not only as a Python warning.
3. D1 (§8.4 slope gate, mean cross-fitted slope in the tolerance interval) stays outside `calibrate()`
   (mean-in-interval is not a `calibrate()` criterion); D2–D4 go through `calibrate()`. The statistician
   signs off the claims before merge.

## Lead decisions after the statistician's review (changes requested)
4. D2–D4 are merged into one check (smallest of the three p-values × 3, Bonferroni) through `calibrate()`,
   to free false-alarm budget (they used ≈ 0.023 of the 0.04/night).
5. Blocking (§8.5): a subgroup calibration gate in tests/calibration/attacks/ — prophage and
   shortest-length-decile O−E z through `calibrate` (level check), plus a confirmed fail for a model
   without size covariates (statistician's draft:
   /home/classes/bios/270/khoa/giap/local/scratch/review-m4-13/draft_test_attack_detectability.py).
6. Structural leakage tests stay (example 4; oracle checks on chosen C and p_detect).
7. `detection_summary` gets Wald 95 % CIs for the slope and intercept in `ci_low`/`ci_high`.
8. D1 runtime is recorded in benchmarks/baseline.json.
9. Implementation deviations accepted: `LogisticRegression(l1_ratio=0.0, tol=1e-6)` (sklearn 1.9 deprecates
   `penalty="l2"`); fits under `threadpoolctl.threadpool_limits(1)` (31 s → 3 s; threadpoolctl is an
   sklearn dependency — declare it in pyproject at the next lock update, backlog); rows sorted by
   locus_tag internally.

## Lead decisions after code review (PR #83, changes requested)
10. Inputs are the real contracts: covariates as returned by M4-10 `covariate_table` and TUs as in C6
    (`locus_tag` and `tu_id` as COLUMNS) are accepted as-is (or the function takes a `locus_tag` column);
    duplicate locus_tags or genes missing from the covariates/TUs raise ValueError naming them (no
    silent median fill, no silent singleton TUs). Tests on the real M4-10/M4-11/C6 tables. The docstring
    says M4-11 `detected_2pep` needs NA masking of partial genes first.
11. `_recalibration` returns NaN slope/intercept/CIs when the predicted logits are constant.
12. (Corrected — the spec rule wins.) Fewer TUs than folds → constant fallback + `too_few_tus` (spec);
    ValueError only for `folds` < 2 and for too few TUs for the inner C search (not a bare `min()` crash); the constant
    fallback (< 50 genes) is logged at WARNING and recorded in diagnostics.
13. `diagnostics` columns have numeric dtypes for value/ci_low/ci_high.
14. Tests kill the surviving mutants: inner C search must not see the held-out fold; inner folds from
    the step stream (not `default_rng(0)`); degenerate fold predicts its detected share; `n_tus == k`
    boundary; convergence-warning filter and its diagnostics row; `dropped_zero_variance` row;
    missing-TU fallback.
15. Benchmark: budget set from the measurement (≈ 3× of 1.47 s/replicate); backlog: a CI check that
    compares benchmarks with baseline.json; add `scikit-learn>=1.8` at the next lock update.
16. docs/methods/detectability.md (statistician) is required before merge.
