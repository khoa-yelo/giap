# M8-08: Determinants models — univariate ρ / AUROC, repeated grouped CV, drop-one-block ΔR²

Milestone M8 · track R · **pass 1** · depends on: M8-07 · estimate ≈ 340 code lines
(+ ≈ 420 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.8.2 step 1 "Rows. Abundance models use detected genes with a non-imputed log2 abundance.
  Detection AUROC uses all eligible genes. Ribosomal-protein genes … are excluded from model
  fitting and from univariate CAI statistics, and are reported separately."
- §7.8.2 step 2 "Univariate. Spearman ρ of each feature with log2 abundance, and AUROC for
  detection, each with a TU-block bootstrap CI."
- §7.8.2 step 3 "Models: LightGBM regression (§5.2 parameters) and an elastic-net baseline
  (standardised features, one-hot categoricals, `l1_ratio` ∈ {0.1, 0.5, 0.9} with inner
  `GroupKFold(3)`). Cross-validation: `GroupKFold(5)` by TU, repeated 5 times with seeded TU
  shuffles. Report cross-validated R² as mean ± SD, and keep out-of-fold predictions averaged over
  repeats."
- §7.8.2 step 5 first bullet "Drop-one-block ΔR² per repeat, with a t-interval over the repeats."
- §5.2 `determinants:`; §6.6 "LightGBM uses `deterministic=true`, `force_row_wise=true`, a fixed
  `num_threads` and a seed"; §9.3 `determinants_model.tsv`.

## Public interface

`src/giap/determinants/models.py` (pure; file I/O none)
```python
@dataclass(frozen=True, slots=True)
class CvResult:
    model: Literal["lightgbm", "elastic_net"]
    r2: F64                          # per repeat
    oof: F64                         # out-of-fold prediction averaged over repeats, per row
    folds: I64                       # (repeats, n_rows) fold index
    delta_r2: Mapping[str, F64]      # block -> per-repeat ΔR²
    boosters: tuple[object, ...]     # repeat x fold fitted LightGBM models (M8-09 SHAP)

def model_rows(features: pd.DataFrame) -> B8: ...
def univariate(features: pd.DataFrame, B: int, rng: np.random.Generator) -> pd.DataFrame: ...
def tu_folds(tu: Sequence[str], n_folds: int, repeats: int, rng: np.random.Generator) -> I64: ...
def cross_validate(X: pd.DataFrame, y: F64, folds: I64, blocks: Mapping[str, str],
                   cfg: DeterminantsConfig, seed: int, threads: int) -> dict[str, CvResult]: ...
def model_rows_table(cv: Mapping[str, CvResult], uni: pd.DataFrame) -> pd.DataFrame: ...
```

1. **Rows.** `model_rows` = eligible, `detected`, `abundance_log2` finite, `abundance_imputed` false,
   `ribosomal_protein` false. AUROC rows = all eligible genes (ribosomal included, except for the
   `cai` row: ribosomal excluded from both ρ and AUROC of `cai` only, as §7.8.2 says).
2. **Univariate.** Per encoded feature column: Spearman ρ with `abundance_log2` on model rows;
   AUROC for `detected` on AUROC rows (Mann–Whitney form, ties ½). CIs: 95 % percentile TU-block
   bootstrap (`tu_block_bootstrap`, B = `stats.bootstrap`, stream `determinants.bootstrap`).
   Rows: `model = univariate`, `feature_or_block` = column, `metric` ∈ {spearman_rho, auroc},
   `n` = rows used; `n_perm` 0. No p-values (§9.3 metric enum has none for univariate; the
   statistical claim of Module 8 is the permutation null of M8-09).
3. **Folds.** `tu_folds`: per repeat, the TU IDs are shuffled with the stream `determinants.cv`
   and assigned round-robin into 5 groups balanced by gene count (largest-first greedy on the
   shuffled order), so `GroupKFold(5)` semantics hold — no TU in two folds — and repeats differ.
4. **LightGBM.** `LGBMRegressor` with §5.2 parameters, `force_row_wise=True`,
   `num_threads = threads`, `seed` from `SeedSequence([seed, repeat, fold])`, `verbose=-1`;
   NaN passed through (native missing handling); categoricals as their one-hot columns (same
   matrix as the elastic net, decision: one encoding for both models).
5. **Elastic net.** `Pipeline(SimpleImputer(median), StandardScaler, ElasticNetCV(l1_ratio=[0.1,
   0.5, 0.9], cv=<inner GroupKFold(3) splits by TU>, n_alphas=100, max_iter=10000,
   selection="cyclic"))`, fitted per outer training fold (imputer and scaler fitted on that fold
   only — decision 14).
6. **CV R².** Per repeat: `1 − Σ(y − ŷ_oof)² / Σ(y − ȳ)²` over all model rows (pooled OOF;
   decision 15). `cv_r2_mean`, `cv_r2_sd` (ddof 1) over repeats. `oof` = mean over repeats.
7. **Drop-one-block.** Per repeat and block present: refit with that block's columns removed on
   the same folds; ΔR² = R²_full − R²_drop. Row `metric = delta_r2`, `value` = mean over repeats,
   CI = `mean ± t_{0.975, R−1} · sd/√R` (R = repeats = 5 → t = 2.776).
8. **Minimum data** (decision 18): < 100 model rows or < 25 TUs → no CV; a single row per model
   with `metric = cv_r2_mean`, `value` NA, `flags = too_few_rows` (index decision 12 adds the
   column). Univariate rows still written when ≥ 10 rows.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| a constant feature on model rows | ρ NA (undefined), AUROC 0.5 only if the classes exist; LightGBM ignores it; elastic net coefficient 0 |
| no undetected genes (or none detected) | AUROC NA |
| a block with every column dropped (layer missing) | no ΔR² row for it |
| fewer TUs than folds in a repeat | cannot happen above the 25-TU minimum |
| y constant | R² undefined → NA rows, `flags = too_few_rows` is **not** used; `flags = unestimated` |
| `threads` changes | identical results (fixed `num_threads` per fit is the caller's `threads`; LightGBM `deterministic=true` makes the result independent of it, tested) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 400 synthetic rows, y = 2·x₁ + noise (SD 1), x₁ ~ N(0, 1) in block "codon", 5 noise columns | LightGBM and elastic-net `cv_r2_mean` within 0.05 of the true R² 0.8; ΔR²(codon) ≈ 0.8, other blocks ≈ 0 with CIs covering 0 |
| 2 | same seed twice; `threads` 1 vs 4 | identical tables |
| 3 | TU column with 80 TUs of mixed size | every TU in exactly one fold per repeat; fold gene counts within ±15 % of each other |
| 4 | feature x and y, Spearman ρ by hand on 6 rows | `spearman_rho` equals `scipy.stats.spearmanr` |
| 5 | AUROC with 3 detected (scores 3, 2, 1) and 2 undetected (scores 2, 0) | pairs: (3>2, 3>0, 2=2 ½, 2>0, 1<2, 1>0) → 4.5/6 = **0.75** |
| 6 | ribosomal genes with high CAI and high abundance | excluded from fitting and from the `cai` rows; included in the `gravy` AUROC |
| 7 | 60 model rows | `cv_r2_mean` NA with `flags too_few_rows`; univariate rows present |
| 8 | elastic-net imputer leakage check: a feature NA in the test fold only | the test-fold value is imputed with the **training** median |

## Acceptance criteria

- [ ] `tests/unit/determinants/test_models.py` covers examples 1–8 and every degenerate row.
- [ ] Property: `tu_folds` never splits a TU; repeats with different seeds differ; the same seed
  reproduces.
- [ ] Property: permuting row order leaves every metric unchanged (row order is canonicalised).
- [ ] Performance (benchmark): 2,500 rows × 35 columns, 5 × 5 CV for both models plus five
  drop-one-block refits each, in < 120 s on 8 CPUs; recorded in `benchmarks/baseline.json`.
- [ ] `docs/methods/determinants.md` §"Models": rows, folds, both models, pooled-OOF R², ΔR² and its
  t-interval (5 repeats share data, so the interval is optimistic — stated).

## Calibration or attack work owed

- Recovery of planted block effects with known R² (§7.8.3 "Calibration") is recorded in M8-15 on
  M8-14's `determinants_planted`; this issue's example 1 is its unit-sized version.

## Hardening owed (pass 2)

- §8.6 rerun of CV R² under `min`/`exclude` (M8-16).
- Check the `ElasticNetCV(l1_ratio=…)` call against the locked scikit-learn (the backlog notes
  ≥ 1.8 changes the `l1_ratio` API); no lock change in M8.

## Questions for the lead

1. Elastic-net missing values: fold-wise median imputation (proposed), or drop rows with any NA
   (loses most rows on runs without one layer)?
2. CV R² per repeat: pooled OOF (proposed) or the mean of the 5 fold R²s?
3. Minimum data: 100 model rows and 25 TUs (proposed)?
4. Univariate rows carry no p-value (§9.3's metric enum has none). Confirm — or should ρ and AUROC
   get bootstrap-based p and a BH family? **Proposed default:** no p; CIs only.

## Lead decisions

1. **Questions 1–4:** the proposed defaults are accepted.
   - Elastic-net missing values are imputed with the fold-wise median, fitted on the training fold only.
   - CV R² per repeat is the pooled out-of-fold R².
   - The minimum is 100 model rows and 25 TUs.
   - Univariate rows carry CIs and no p-value or BH family.
2. **From the M8-07 review (#162):** the model matrix has constant columns on real data (14 on sim100k: `start_codon_other`, `replicon_type_unknown`, 12 `nterm_residue_*`).
   - Drop zero-variance and all-NaN columns per training fold before standardising.
   - Record the dropped columns for the step's events.
   - LightGBM takes the same matrix.
3. **The test author's choices** (commit 6c2ad3f) are accepted as the contract.
   - `enough_data`, and `CvResult.dropped`.
   - `univariate` encodes with `model_matrix` and bootstraps by TU.
   - Rows are sorted internally.
   - `flags` is a `;`-joined string, and constant y gives `unestimated`.
   - The LightGBM tolerances in example 1 are ±0.15 with ΔR²(codon) > 0.7. The spec's ±0.05 holds only for the elastic net: LightGBM's drop model overfits the noise.
   - Examples 4 and 5 are doubled, so each univariate row has ≥ 10 rows.
4. **Inner CV grouping:** `cross_validate` takes the TU IDs as a keyword argument (`tu=`), so the inner `GroupKFold(3)` is by TU, as the spec says. A test for it is owed to the code-review round, and the lead adds it.
5. **The ΔR² t-interval under-covers.** Under the null it misses 0 in 0.30–0.57 of replicates for LightGBM and 0.04–0.15 for the elastic net, against a bound of about 0.09.
   - In pass 1, every ΔR² row carries `ci_unestimated`: the interval is reported, but it makes no coverage claim (R4).
   - **Owed to Khoa (ADR):** a calibrated interval, for example a TU-block bootstrap over repeats or nested CV, or relabelling it.
6. **Protected file:** the calibration key `determinants_delta_r2_max_abs_bias` (proposed 0.05) is missing from `tests/calibration/tolerances.yaml`.
   - The file is protected, so the key needs Khoa.
   - Until then the M8-08 calibration tests fail on the missing key. They are `calibration_nightly`/`calibration_full`, not every-commit.
