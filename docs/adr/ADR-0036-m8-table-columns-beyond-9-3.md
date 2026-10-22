# ADR-0036: M8 function tables carry four columns and five run flags beyond GIAP_SPEC.md §9.3 and §9.4
- Status: proposed
- Date: 2026-10-22
- Spec sections: §7.7.2, §7.7.3, §9.2, §9.3 (`category_summary.tsv`, `enrichment_detection.tsv`), §9.4 (Run scope); amends ADR-0029; docs/specs/M8-06-function-step-and-run.md (questions 4 and 9, lead decisions 1 and 2); docs/specs/M8-04-category-summary-ora.md (lead decisions 3 and 13); M8-index "ADR-A" (index decision 12)
- Needs Khoa: yes. It adds columns to two §9.3 dictionaries and five tokens to the §9.4 run vocabulary. The lead's default is option (a). M8-06a goes ahead on that default (M8-06 lead decision 2). If Khoa rejects this ADR, only the two table schemas and the `COLUMNS` tuples change.

## Context
- **§9.3 lists these columns:**
  - `category_summary.tsv`: `collection`, `category`, `name`, `n_genes`, `n_detected`, `expected`, `o_minus_e`, `ci_low`, `ci_high`, `median_abundance_log2`, `coverage`. There is no `flags` column.
  - `enrichment_detection.tsv`: `collection`, `term`, `name`, `direction`, `n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set`, `odds_weight`, `p_classic`, `parent_group`, `coverage`, plus the §9.2 test columns.
- **ADR-0029 changed what these columns mean.** The ORA unit is the TU and the null is the conditional Poisson-binomial (CPB). The four count columns are now counted in TUs. `p_classic` is now the gene-level Wallenius tail, kept only for comparison and goseq parity. In §7.7.3, `p_classic` meant "the non-elim p of the reported test". Under ADR-0029 that quantity has no column.
- **M8-index ADR-A (index decision 12)** added a `flags` column to `enrichment_detection`, `enrichment_gsea` and `determinants_model`, because R4 needs an `unestimated` token somewhere. It does not cover `category_summary`. ADR-A was recorded as a lead-level decision in the M8 index, not as an ADR file.
- **The code in PR #179 (origin/m8/ora) and the schemas on origin/m8/function-step** already contain the four extra columns:
  - `src/giap/function/summary.py` `COLUMNS` ends in `flags`;
  - `src/giap/function/ora.py` `COLUMNS` has `n_in_set_genes` and `k_detected_in_set_genes` after `k_detected_in_set`, and `p_noelim` after `p_classic`;
  - `src/giap/schemas/tables/category_summary.schema.json` and `enrichment_detection.schema.json` match those tuples. A schema-equals-`COLUMNS` test pins them (M8-06 lead decision 3).
- **Five run-level tokens** come out of the `function` step, through `df.attrs["run_flags"]` and the M8-06 step guard. §9.4's Run scope lists only `anti_sd_default`, `c_assumed`, `layer_missing:<layer>`, `bakta_db_light` and `operon_source_distance_model`.

## Decision
Option (a): accept all four columns and all five run tokens as written on origin/m8/function-step. Every addition is non-breaking under §9.2's semantic versioning: a new column goes last or next to its partner, and no existing column is renamed or retyped.

### 1. Columns

| Table | Column | Position | Type and constraints (schema) | Meaning | Why it is needed |
|---|---|---|---|---|---|
| `category_summary` | `flags` | last | `;`-joined tokens; pattern allows only `ci_few_units` today | Row flags. `ci_few_units`: the category has fewer than `CI_MIN_UNITS` = 10 TUs (`summary.py`), so the TU block-bootstrap CI is not computed and `ci_low`/`ci_high` are `NA`. | §9.2 requires flags in a `flags` column. R4 requires every `NA` to say why. Without the column, an `NA` CI is unexplained. ADR-A gave the other M8 tables this column but missed this one. The column also gives a home to any later CI-coverage token (statistician, #179: percentile CIs cover 0.886 at 10 TUs; backlog). |
| `enrichment_detection` | `n_in_set_genes` | after `k_detected_in_set` | int ≥ 0, required | Eligible genes in the term: the gene-level n that `p_classic` (Wallenius on genes) uses. | ADR-0029 moved `n_in_set` to TUs. Without this column the inputs of `p_classic` are not in the table, so it cannot be checked or compared with goseq. The M8-13 hardening golden test compares these counts with goseq's per-category counts (M8-04 "Hardening owed"). Readers also see partial-membership dilution as the gap between `n_in_set` (TUs) and `n_in_set_genes`. |
| `enrichment_detection` | `k_detected_in_set_genes` | after `n_in_set_genes` | int ≥ 0, required | Detected eligible genes in the term: the gene-level k behind `p_classic`. | Same as `n_in_set_genes`. The gene-level k cannot be recovered from the TU counts, because a detected TU can hold one detected member gene or several. |
| `enrichment_detection` | `p_noelim` | after `p_classic` | number in [0, 1]; `NA` when the law is undefined | The CPB tail, in TUs, for the row's direction, on the term's full TU set before elim (M8-04 §6). For rows without elim (`under`, `cog`, `kegg`, `go_slim`) it equals `p_value`. For GO `over` rows it is the pre-elim p, and `p_value` is the elim p. | It restores §7.7.3's "classic" (non-elim) p of the reported test, which ADR-0029 took out of `p_classic`. It shows when a GO parent is significant only through a child elim removed. **It is descriptive only: it is never adjusted and never enters BH** (M8-04 §5; methods note `docs/methods/function.md` on origin/m8/ora, "Tests per term"). This matters because it is not FDR-safe. With one planted child per parent it rejects 0.0093 / 0.0556 / 0.1172 at α 0.01 / 0.05 / 0.10, against bounds 0.0118 / 0.0540 / 0.1055, while elim stays at 0.0008 / 0.0073 / 0.0189 (`$GIAP_ROOT/stat/m804b/probe_elim.py`, `elim.log`, 12,000 parents, R = 600; the job ID is not recorded in the log: owed). |

The schema descriptions say which columns are in TUs (`n_universe`, `n_in_set`, `n_detected`, `k_detected_in_set`) and which are in genes (the two `_genes` columns). §9.2 puts units in column names only for physical quantities.

### 2. Run-scope flag tokens
These tokens go in the `function.run_flags` event now (M8-06 question 9, lead decision 1). They join the §9.4 Run vocabulary through the manifest schema during hardening, as `operon_source_distance_model` did.

| Token | Emitted when | Source | Also a row token? |
|---|---|---|---|
| `tu_recal_skipped` | ADR-0035's cross-fitted TU recalibration did not run: `folds` is None or empty (M4-13's constant model; `tu_folds.parquet` has only a header), a fold's training TUs have fewer than `MIN_CLASS` = 10 detected or undetected TUs, or the Newton fit failed. q falls back to q_raw. | `ora.tu_universe` | yes: every `enrichment_detection` row of that run (M8-04 decision 13) |
| `too_few_ranked_tus` | Fewer than 2 TUs have an abundance score, so `enrichment_gsea.tsv` is written with only a header. | `gsea.gsea_table` | no |
| `no_eligible_genes` | The gene table has no eligible (non-partial) gene, or is empty. All three tables are written with only a header (M8-06 step 2 guard). | `steps/function` | no |
| `p_detect_clipped` | At least one M4-13 `p_detect` was exactly 0 or 1 and was clipped to [1e-12, 1 − 1e-12] for the null (ADR-0035 decision 2). | `ora.tu_universe` | yes: rows whose full TU set holds a clipped TU |
| `null_degenerate` | The observed detection has probability 0 under the TU null (an undetected TU with q = 1). Every ORA p is `NA` and a `UserWarning` is raised. | `ora.detection_enrichment` | yes: every row, with `unestimated` |

The run-only tokens (`layer_missing:*`, `bakta_db_light`, `tu_recal_skipped`, `too_few_ranked_tus`, `no_eligible_genes`) are not added to the `category_summary` or `enrichment_gsea` row patterns.

### 3. ADR-0029 amendment line
ADR-0029 now reads with this addition: "`enrichment_detection.tsv` counts are in TUs. `p_classic` is gene-level Wallenius, with its gene counts in `n_in_set_genes` and `k_detected_in_set_genes`. `p_noelim` is the full-set TU CPB p. Only `p_value` enters BH."

## Options
- **(a) Accept all four columns and the five run tokens (lead default; chosen above).** Non-breaking. It matches the code and schemas already under review. Nothing is lost.
- **(b) Drop some columns and fold them into `flags` or the methods note.** Rejected.
  - `p_noelim` and the gene counts are numbers, and a flag cannot carry a number.
  - The methods note is per method, not per row. Without the gene counts, `p_classic` cannot be checked against goseq for a single term. Without `p_noelim`, a reader cannot tell an elim-reduced GO parent from one with no signal.
  - The gene counts cannot be re-derived from the published tables: that needs the gene table and the collection resources.
  - Dropping `category_summary.flags` leaves `NA` CIs unexplained (R4).
- **(c) Put the extra columns in a sidecar table** (for example `tables/enrichment_detection_detail.tsv`). Rejected.
  - It adds a fourth function table, a schema, a join key (`collection`, `term`, `direction`) and a publishing rule. M6-08 decision 7 publishes only the three `tables/*.tsv` files.
  - It splits one row's statistics across two files, so a reader of `p_classic` would need a join to see its inputs.
  - It solves no compatibility problem, because the added columns are already non-breaking.

## Evidence (before → after: metrics, tests, simulations)
- Before: §9.3 has no home for the gene-level inputs of `p_classic` after ADR-0029, for the non-elim CPB p, or for the reason a `category_summary` CI is `NA`.
- After:
  - the columns are in `summary.COLUMNS` and `ora.COLUMNS` (origin/m8/ora) and in both schemas (origin/m8/function-step), pinned by the schema-equals-`COLUMNS` test (M8-06 lead decision 3);
  - the run tokens are covered by M8-06 examples 1 and 4–7 (`tu_recal_skipped`, `no_eligible_genes`, `layer_missing:*`) and by #179's skip-path and degenerate-run tests;
  - the elim-versus-`p_noelim` rates above are the evidence that `p_noelim` must stay out of BH.

## Consequences
- §9.3 reads as amended for these two tables. The §9.4 Run row gains five tokens, registered through the manifest schema during hardening.
- Report and V7 legends label `p_noelim` and `p_classic` as "not adjusted; not a discovery claim".
- Rejecting this ADR, wholly or in part, changes only the two schemas, the `COLUMNS` tuples and their pinning test. No statistic changes.
- **Open for the lead:**
  - The function row tokens (`ci_unestimated`, `unestimated`, `odds_out_of_range`, `elim_empty`, `ora_uncalibrated`, `ci_few_units` and GSEA's `size_filtered`, `genomic_cluster`, `detection_confounded`) are not in any §9.4 scope either. They are registered only in the table schemas, which §9.4's title allows ("extend only through the schema"). Should §9.4 gain a "Function" scope row in hardening?
  - ADR-A (index decision 12) exists only as a decision in the M8 index. Should it become an ADR file, or should this ADR be read as covering it for the function tables?
