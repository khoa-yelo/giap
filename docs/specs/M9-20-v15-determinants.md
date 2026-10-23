# M9-20: V15 `determinants`: what explains abundance?

Milestone M9 · **pass 2** · depends on: M9-17 (SHAP, PD support); PR #178 merged · estimate
≈ 340 code lines (+ ≈ 320 test lines) · **stats (light)**: the statistician checks that the panels
present ADR-0032/0033/0034's claims correctly

## Spec sections implemented

- §10.2 V15: "(a) block ΔR² with CIs; (b) partial dependence for CAI, tAI, OSTIR score and
  promoter score; (c) SHAP beeswarm (top 20 features); (d) observed vs. out-of-fold predicted";
  Plotly twin.
- §7.8.2 (models, CV R² mean ± SD, drop-one-block ΔR² with a t-interval over repeats, SHAP
  `pred_contrib`, ribosomal-protein genes "reported separately", residual flag
  `abundance_not_explained_by_sequence` > 2 SD), §9.3 `determinants_model.tsv`,
  `determinants_attribution.tsv` (branch `origin/m8/determinants-step`), ADR-0031 (ΔR² tolerance),
  ADR-0032 (univariate CIs `ci_unestimated`), ADR-0033 (elastic-net ΔR² primary; LightGBM and SHAP
  secondary), ADR-0034 (`null_uncalibrated` on `null_p`).

## Panels (double width, 2 × 2)

- **(a) Block ΔR²:** `determinants_model` rows `metric = delta_r2`, one row per block
  (transcription, translation, codon, protein, context), elastic net (filled circle, primary) and
  LightGBM (open diamond, secondary) side by side, `ci_low`–`ci_high` bars; vertical line at 0;
  legend "elastic net (primary, ADR-0033)", "LightGBM (secondary check)".
- **(b) Partial dependence** (`determinants_pd.tsv`): four small panels cai, tai, tir_ostir,
  promotech_max ("promoter score" = `promotech_max`, which M8-07 builds from `promoter_score`);
  `pd_mean` line with the `pd_low`–`pd_high` spread shaded; a rug of feature values; label "spread
  over 25 CV models (not a CI)". Feature absent → "not available (`layer_missing:<layer>`)".
- **(c) SHAP beeswarm** (`determinants_shap.tsv.gz`): top 20 features by mean |SHAP| top to bottom;
  points jittered deterministically by rank within the feature; colour viridis by the feature value
  rank within the feature (low → high); x "SHAP (log2 abundance)"; header "LightGBM (secondary)".
- **(d) Observed vs predicted:** `determinants_attribution` rows: x = `predicted_log2`, y =
  `abundance_log2`; rows flagged `abundance_not_explained_by_sequence` in `#D55E00`;
  `cai_reference` (ribosomal) genes drawn as `#E69F00` crosses and excluded from the R² text (they
  are not in the model rows, §7.8.2); identity line; text "CV R² `<m>` ± `<sd>` (elastic net)".

Footnotes: "null p `<p>` (uncalibrated: valid only for exchangeable TUs, ADR-0034)"; "univariate
CIs uncalibrated (ADR-0032)" when univariate rows appear in the table section; the §7.8.2 caveat
"expect a lower R² on native genes" (from `giap.determinants.CAVEATS`).

## Data TSV `V15_determinants.tsv`

`panel` enum{a, b, c, d}, `model`, `block`, `feature`, `x`, `y`, `ci_low`, `ci_high`, `locus_tag`,
`category`, `flags`. Sorted by panel, model, block/feature (drawing order), x.

## Title

"Elastic net explains CV R² `<m>` ± `<sd>`; largest block ΔR²: `<block>` `<d>` (95 % CI `<lo>`–
`<hi>`)"; if no block CI excludes 0: "… no block's ΔR² CI excludes 0". With `too_few_rows`: "Too
few genes for the determinants model (`<n>` rows, `<t>` TUs)".

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `too_few_rows` (sim100k) | (a)–(c) empty; (d) drawn if attribution rows exist; title as above |
| support files absent | (b), (c) "not available"; (a), (d) drawn |
| a block with no feature left (all its layers missing) | no row for it (M8 rule) |
| ΔR² CI NA | point without bar, "CI unestimated" |
| `predicted_log2` NA on some rows (non-model rows) | not drawn in (d) |
| > 20,000 genes | (d) binned 2-D counts (hexbin) per M9-01 rule; (c) capped by the same rule per feature |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | elastic_net cv_r2_mean 0.21, sd 0.03; delta_r2 transcription 0.08 (0.05, 0.11), others smaller | title "Elastic net explains CV R² 0.21 ± 0.03; largest block ΔR²: transcription 0.08 (95 % CI 0.05–0.11)" |
| 2 | all elastic-net ΔR² CIs include 0 | title ends "no block's ΔR² CI excludes 0" |
| 3 | null_p row flags `null_uncalibrated` | footnote contains "(uncalibrated" |
| 4 | gene with flag `abundance_not_explained_by_sequence` | (d) row category `unexplained`, colour `#D55E00` |
| 5 | `cai_reference` gene | (d) category `ribosomal_reference`; not in the R² text |
| 6 | model table with `too_few_rows` | title "Too few genes for the determinants model (…)" |

## Acceptance criteria

- [ ] `tests/unit/viz/test_v15_determinants.py`: examples 1–6; the primary/secondary legend text;
  PD "not a CI" label.
- [ ] Plotly twin registered.
- [ ] Coded against the branch schemas; re-checked when #178 merges.

## Questions for the lead

1. Which model's ΔR² leads the title: elastic net per ADR-0033. **Proposed default:** elastic net.
2. "promoter score" = `promotech_max` (M8-07 maps it from `promoter_score`). **Proposed default:**
   `promotech_max`.

## Lead decisions (2026-10-23)
1. **Q1 yes:** the elastic net's ΔR² leads the title (ADR-0033).
2. **Q2 yes:** "promoter score" is `promotech_max`.
3. **Stacking:** on m9/plot-support-analysis (M9-17 SHAP and PD support files, decision 5).
4. **Test-author (a85b325a): gaps 1–9 accepted as proposed.**
   - The title names the largest block whose ΔR² CI is wholly above 0; a CI wholly below 0 does
     not count.
   - The title quotes a finite CI even though M8 flags every ΔR² row `ci_unestimated`. The
     statistician reviews the wording in the PR (backlog: the meaning of that flag on finite CIs).
   - Too-few-rows counts are model rows (detected, non-imputed, not ribosomal) and their distinct
     `tu_id`.
   - The missing-layer name is read from the run flags, falling back to the map tir_ostir→ostir,
     promotech_max→promoters, cai/tai→codon.
   - V15 declares `layers=("abundance",)`.
   - The pinned axes labels and category names are the interface.
5. **Review round (code review of 4efbd957). Decision 4 amended.**
   - **B1, ΔR² interval:** keep the numbers but name the interval honestly.
     - Title: "… largest block ΔR²: <block> <v> (repeat t-interval <lo>–<hi>, uncalibrated)".
     - Panel (a) header: "(a) block ΔR² (t-interval over CV repeats, uncalibrated)".
     - Footnote: "ΔR² intervals uncalibrated: CV repeats share data (ADR-0031)".
     - "CI unestimated" stays for NaN intervals. The naming rule (interval wholly above 0) is unchanged.
   - **B3:** the sentence for no named block reads "no block's ΔR² interval lies above 0". Add an example
     with an interval wholly below 0.
   - **B2, panel (d):**
     - Retitle it "(d) observed vs. LightGBM out-of-fold".
     - Its R² text quotes both models: "CV R² <v> ± <sd> (LightGBM, plotted); elastic net <v> ± <sd> (primary)".
     - The flag label is "unexplained (> 2 SD of LightGBM residuals)".
     - Elastic-net out-of-fold predictions in M8 go to the backlog.
   - **B4, hexbin:**
     - Bin only `category == "model"`, with `bins="log"`.
     - Add a colourbar "genes per bin (log count)" and a legend proxy "model genes (binned)".
     - The big-genome test runs the M9-24 checklist.
   - **S1:** label panel (b) "partial dependence, LightGBM (secondary); spread over <n> CV models (not a CI)".
   - **S2:** neutral title "Elastic net CV R² …", with a Unicode minus and "unestimated" for NA.
   - **S3, twin:**
     - Dodged numeric y, with circle vs. diamond-open markers.
     - Reversed y axis.
     - The footnotes as an annotation, and "CI unestimated" marks.
   - **Nits:**
     - N1: legend at loc="best".
     - N3: tests for the four surviving mutants (low end exactly 0, imputed rows, a feature's points at its own tick, a half-NA CI).
     - N4: break ties in BLOCKS order.
     - Delete the unreachable `few and d_rows` branch.
     - N2 and N5 are backlogged.
