# M8-09: SHAP block attribution, the dominant block, residual flags, and the TU-size-preserving null

Milestone M8 · track R · **pass 1** · depends on: M8-08 · estimate ≈ 260 code lines
(+ ≈ 380 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.8.2 step 4 "Null. Permute the response at TU level, B = 1,000 times: response blocks move among
  TUs with the same number of included genes. TU sizes that occur only once stay in place, which is
  conservative. Refit with one CV repeat each time. p = (1 + b)/(1 + B), where b counts
  permutations with CV R² ≥ the observed mean CV R². Computed separately for each model."
- §7.8.2 step 5 second bullet "Per-gene SHAP values from LightGBM `pred_contrib=True` on the
  out-of-fold models, summed by block. No `shap` package is needed."
- §7.8.2 step 6 "The dominant block for each top-decile protein = the block with the largest
  positive SHAP sum. Out-of-fold residuals > 2 SD are flagged `abundance_not_explained_by_sequence`:
  candidates for regulation or stability."
- §7.8.3 "Unit: TU-size-preserving permutation blocks; the p formula; SHAP block sums equal the
  prediction minus the expected value."
- §8.1 row "Determinant model | model | abundance permuted across TUs of equal size | permutation p
  | 0.05"; §8.2 permutations (1,000; 10,000 with `--publication`; p = (1 + b)/(1 + B), seeds
  recorded); §9.3 `determinants_attribution.tsv`; §9.4 gene flag
  `abundance_not_explained_by_sequence`.

## Public interface

`src/giap/determinants/attribution.py` and `src/giap/determinants/null.py` (pure)
```python
def shap_blocks(cv: CvResult, X: pd.DataFrame, blocks: Mapping[str, str]) -> pd.DataFrame: ...
def attribution_table(features: pd.DataFrame, rows: B8, cv: CvResult,
                      shap: pd.DataFrame) -> pd.DataFrame: ...
def tu_size_blocks(tu: Sequence[str]) -> list[I64]: ...          # equal-size groups of TUs
def permute_response(y: F64, tu: Sequence[str], rng: np.random.Generator) -> F64: ...
def null_p(X, y, tu, folds_rep1: I64, observed_mean_r2: Mapping[str, float],
           cfg: DeterminantsConfig, B: int, seed: int, threads: int) -> dict[str, NullResult]: ...
```

1. **SHAP.** For every repeat r and fold f, `booster.predict(X[test_rf], pred_contrib=True)`
   gives one column per encoded feature plus the expected value (last column). Per gene, average
   the contribution vectors over the 5 repeats (each repeat holds the gene out exactly once), then
   sum columns by block → `shap_transcription`, `shap_translation`, `shap_codon`, `shap_protein`,
   `shap_context` (NA for a block with no column). Invariant (§7.8.3): Σ blocks + mean expected
   value = mean OOF LightGBM prediction, per gene, to 10⁻⁹ (tested on every row).
2. **Attribution rows** = model rows only (detected, non-imputed, non-ribosomal), plus the
   ribosomal genes "reported separately": they get a row with `predicted_log2` from a model that
   never saw them (the fold models applied to them, averaged), SHAP likewise, and the flag
   `cai_reference` — they are never in the residual SD or the top decile (decision in Q3).
   `predicted_log2` = mean OOF LightGBM prediction; `residual` = observed − predicted.
3. **Top decile** = model rows with `abundance_log2` ≥ the 90th percentile (linear interpolation);
   `top_decile` true/false. `dominant_block` (top decile only) = argmax of the positive block SHAP
   sums; NA when none is positive; ties by `BLOCKS` order.
4. **Residual flag** (decision 17): |residual| > 2 · SD(residuals over model rows, ddof 1) →
   `flags` gains `abundance_not_explained_by_sequence`.
5. **Permutation blocks.** TUs are grouped by their number of **model rows** (included genes);
   within each size group with ≥ 2 TUs the TU-level response vectors (genes in TU order) are
   permuted among the TUs; groups of one TU stay put (the spec's conservative rule).
6. **Null loop.** For b = 1 … B: permute y (stream `determinants.null`, child seed b from
   `SeedSequence.spawn`), refit **one** CV repeat on the repeat-1 fold assignment (decision 16) for
   each model, record its pooled-OOF R². `b_count` = #{R²_b ≥ observed `cv_r2_mean`};
   `p = (1 + b_count)/(1 + B)` via `giap.stats.perm_p` with `alternative="greater"`. B =
   `determinants.permutations` (10,000 with `--publication`, §8.2). Permutations run in parallel
   over `threads` workers, each fit single-threaded, results collected in b order: identical for
   any thread count.
7. **Model rows** (to M8-10): per model, `metric = null_p`, `value = p`, `n_perm = B`; the
   observed and null mean R² in the methods-note diagnostics (V18 consumes the null vector, M9).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| every TU size unique | no movable block: every permutation is the identity → b = B, p = 1 (conservative, stated); flag `null_degenerate` on the row |
| < 100 model rows (M8-08 minimum) | no SHAP, no null; attribution table has the rows with `predicted_log2` NA and `flags too_few_rows` |
| no block positive for a top-decile gene | `dominant_block` NA |
| residual SD 0 | no flag set |
| elastic net | null p computed; no SHAP (§9.3 has one SHAP column set; LightGBM only) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | TUs A (2 genes), B (2), C (3), D (1), E (2) | size groups {A, B, E} (size 2) movable, {C} and {D} fixed; a permutation never maps a 2-gene response block onto C or D |
| 2 | y = [a1 a2 \| b1 b2 \| c1 c2 c3], A→B, B→A | permuted y = [b1 b2 \| a1 a2 \| c1 c2 c3] (gene order inside a block kept) |
| 3 | B = 9, null R² = [0.1, 0.2, 0.3, 0.5, 0.05, 0.6, 0.0, 0.1, 0.2], observed mean 0.45 | b = 2 → p = 3/10 = **0.3** |
| 4 | B = 999, no null ≥ observed | p = 1/1000 |
| 5 | a 3-gene toy LightGBM on 200 rows | per row Σ block SHAP + expected = prediction to 10⁻⁹ |
| 6 | residuals with SD 1, one gene at +2.5, one at −2.3 | both flagged (two-sided) |
| 7 | top decile of 100 rows = 10 genes; a gene with block sums (+0.4, −0.1, +0.6, 0, −0.2) | `dominant_block` = **codon** (the largest positive) |
| 8 | same seed, threads 1 vs 8 | identical p, SHAP and flags |

## Acceptance criteria

- [ ] `tests/unit/determinants/test_null.py`: examples 1–4, 8 (the §7.8.3 "TU-size-preserving
  permutation blocks" and "the p formula" units).
- [ ] `tests/unit/determinants/test_attribution.py`: examples 5–7 (the §7.8.3 SHAP identity).
- [ ] Property (Hypothesis): `permute_response` preserves the multiset of TU sizes, each TU's
  internal order, and fixes every TU whose size is unique.
- [ ] Performance: measure and record in `benchmarks/baseline.json` the null at B = 1,000 on 2,500
  rows on 16 CPUs (the index estimates 30–40 CPU-minutes per model; the measured figure replaces
  the estimate and sets the `process_high` request in M8-10).
- [ ] `docs/methods/determinants.md` §"Null and attribution": why one repeat per permutation is
  conservative against a 5-repeat observed mean, the unique-size rule, SHAP block sums.

## Calibration or attack work owed

- `determinants` suite (null at level) and the planted-R² recovery: M8-15 on M8-14's generators.
  The statistician decides there whether decision 16 (fixed repeat-1 folds) is acceptable.

## Hardening owed (pass 2)

- Besag–Clifford early stopping (h = 10, `stats.besag_clifford_h`) is **not** in the spec for this
  null; if runtime is prohibitive at `--publication`, it needs an ADR, not a silent change.

## Questions for the lead

1. Null refits: reuse the repeat-1 folds for every permutation (proposed) or draw fresh folds per
   permutation (adds fold noise to the null, making it more conservative still)?
2. Residual flag two-sided |r| > 2 SD (proposed) or only r > 2 SD (abundant beyond prediction)?
   The spec's "candidates for regulation or stability" fits both directions.
3. Ribosomal genes "reported separately": rows in `determinants_attribution.tsv` with
   `cai_reference` and predictions from models that never saw them (proposed), or no rows at all?
