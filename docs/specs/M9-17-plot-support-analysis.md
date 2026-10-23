# M9-17: Plot support from the analysis steps: running scores, SHAP, partial dependence, null draws

Milestone M9 · **pass 2** · depends on: M9-03 (ADR-M9-A); PR #178 (determinants step) and #181
(function step) merged · estimate ≈ 340 code lines (+ ≈ 350 test lines) · **stats (light)**:
partial dependence is a model-explanation computation; null draws are written as computed

## The gap this closes

V14(a), V15(b,c) and V18(c) need things that only a fitted model or a running null has: the ranked
TU scores and set positions behind each GSEA running sum, per-feature SHAP values (the table keeps
block sums only), partial-dependence curves, and the null statistics of regions, islands and the
determinants permutation test. §7.9.4 forbids `giap report` from recomputing statistics, so the
producing steps write them once, as support files (ADR-M9-A).

## Spec sections implemented

- §10.2 V14(a) "GSEA running-score plots for the top 6 sets", V15(b) "partial dependence for CAI,
  tAI, OSTIR score and promoter score", V15(c) "SHAP beeswarm (top 20 features)", V18 "null
  distributions against the observed statistic (regions, islands, determinants)".
- §7.7.4 (TU-collapsed fgsea; ADR-0028 median-centred scores), §7.8.2 (25 OOF boosters:
  5 folds × 5 repeats; `pred_contrib`; the permutation null), §7.3.5(b) (circular shifts),
  §7.3.5(d) (island max-statistic null), §6.6 (determinism).

## Support files (all under `plots/support/`, ADR-M9-A rules)

| File | Written by | Columns | Size bound |
|---|---|---|---|
| `gsea_running.tsv.gz` | `function` step | `collection`, `term`, `rank_position` (1…n TUs), `tu_id`, `score` (the centred TU score fgsea ranked), `in_set` (bool), `running_es` | top 6 sets × n TUs (≈ 6 × 2,500 rows on K-12) |
| `determinants_shap.tsv.gz` | `determinants` step | `locus_tag`, `feature`, `feature_value`, `shap` | model rows × top-20 features (by mean \|SHAP\|) |
| `determinants_pd.tsv` | `determinants` step | `feature` ∈ {cai, tai, tir_ostir, promotech_max}, `grid_value`, `pd_mean`, `pd_low`, `pd_high` (min/max over boosters) | 4 × 20 rows |
| `null_regions.tsv.gz` | `context_stats` step | `region_id`, `stat`, `shift_k`, `null_value`, `observed` | regions × K ≤ 2 × N per region; capped (Q2) |
| `null_islands.tsv` | `context_stats` step | `contig`, `direction`, `replicate`, `max_scan_s`, `observed_max` | per contig × direction × B |
| `null_determinants.tsv` | `determinants` step | `model`, `permutation`, `cv_r2`, `observed_cv_r2_mean` | 2 × B |

Definitions:

- **Top 6 GSEA sets**: smallest `q_value`, then largest |`nes`|, then collection, term; rows with
  `p_value` NA excluded. `running_es` is fgsea's running sum for that set, computed by the same
  port function the `es` came from; the max-deviation of `running_es` equals the row's `es` to
  1e-9 (an assertion in the step).
- **SHAP per feature**: the mean over the 25 out-of-fold boosters of `pred_contrib` (the same
  matrix M8-09 sums by block); block sums of this file equal `determinants_attribution` block
  columns to 1e-9.
- **Partial dependence**: grid = the 5 %, 10 %, …, 95 % quantiles (19 points) plus min and max of
  the feature among model rows; for each booster, the mean prediction over model rows with the
  feature set to the grid value; `pd_mean` the mean over boosters, `pd_low`/`pd_high` their
  min/max (a spread, labelled so, not a CI). PD is computed only for LightGBM (the nonlinear model;
  ADR-0033 keeps it secondary).
- **Null draws** store the statistic each test compares, not the reported effect: for region
  `o_minus_e` the studentised statistic of ADR-0018 (and its observed value), for `hl_shift` the
  statistic of ADR-0018's amendments; the column header names it. They are the values the step already computes; no extra replicates are drawn, so no
  RNG stream changes (§6.6).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| function off / no abundance | no `gsea_running` file |
| determinants below 100 rows / 25 TUs (`too_few_rows`) | no SHAP, PD or null files |
| a PD feature missing (`layer_missing:ostir`) | that feature's rows absent |
| a region with K = 0 shifts (linear contig, region too wide) | no rows for it |
| `--publication` (B = 10,000) | null files grow tenfold; still < 5 MB (gz) on K-12 (to be measured) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 6 TUs scores (centred) 2, 1, 0.5, −0.5, −1, −2; set {TU_1, TU_3} | `running_es` rows match the fgsea port's running sum; max deviation = table `es` |
| 2 | determinants on the M8-09 fixture | `determinants_shap` summed by block = attribution block columns (1e-9) |
| 3 | feature cai, 25 boosters | 21 grid rows; `pd_low ≤ pd_mean ≤ pd_high` |
| 4 | region with N = 100 genes on a circular contig | 99 `shift_k` rows per stat |
| 5 | the same run twice | byte-identical support files |

## Acceptance criteria

- [ ] `tests/unit/steps/test_support_function.py`, `test_support_determinants.py`,
  `test_support_context.py`: examples 1–5; files absent in the degenerate cases.
- [ ] No change to any `tables/` file bytes (a test diffs the tables with and without the support
  writer).
- [ ] The statistician signs the PD definition and the "spread, not CI" wording.

## Questions for the lead

1. This issue edits the M5 `context_stats` and M8 `function`/`determinants` steps (new outputs).
   **Proposed default:** acceptable as M9 work after #178/#181 merge, with the M8 owners reviewing.
2. Region null size: exhaustive shifts give N − 1 values per region (≈ 4,500 on K-12, × 2 stats ×
   ≈ 30 regions ≈ 270,000 rows, ≈ 3 MB gz). **Proposed default:** write them all if < 5 MB (R11),
   else the null's 1st–99th percentile grid (99 values) per region, recorded in the file header.
3. PD for LightGBM only. **Proposed default:** yes; the elastic net is linear (its PD is a line
   whose slope is the coefficient).

## Lead decisions (2026-10-23)
1. **Q1 yes:** M9-17 edits the M5 `context_stats` and the M8 `function` and `determinants` steps
   (new outputs). Stacked on #178 (M8-10). #181 is merged. The code review covers the M8 owners'
   concerns.
2. **Q2 yes:** write every exhaustive-shift null value if the file is under 5 MB (R11); otherwise
   write the 1st–99th percentile grid (99 values) per region, and record which in the file header.
3. **Q3 yes:** PD for LightGBM only; for the elastic net, PD is the coefficient line.
4. **Stacking:** after #178 squash-merges, merge the old base tip, then main, and retarget.
5. **Test-author (a3224acf): gaps 1–10 accepted as proposed.**
   - One optional directory output per step (`function_support`, `determinants_support`,
     `context_stats_support`). Local `publish` places the files in `OUT/plots/support/`, and the
     gzip rule follows ADR-0043.
   - The Q2 grid file's third column is `percentile`; the size bound is
     `context_stats.NULL_REGIONS_MAX_BYTES`.
   - The block-sum identity holds for the full per-column SHAP matrix (`shap_features`); the
     top-20 file is cut from it.
   - The PD grid is on the feature's own scale (`tir_ostir` is logged before predicting).
   - `null_islands` is per scanned contig, using its own per-replicon maximum, with `direction` =
     `abs`. `null_regions` gains a `contig` column. No `gsea_running` file with zero testable sets.
   - The lead updates the existing step tests the new outputs contradict
     (test_determinants_step.py:122-123; test_function_step.py:221-222 and 256-266) and regenerates
     `workflow/tests/data/annotate/step_registry.json` after the implementer.
   - Next: the implementer.
