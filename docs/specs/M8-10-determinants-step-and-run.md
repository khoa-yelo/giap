# M8-10: The `determinants` step, its three tables, and wiring them into `gippa run`

Milestone M8 · track T (tables and integration) · **pass 1 — the Module 8 integration milestone** ·
depends on: M8-07, M8-08, M8-09 (and merged MVP-01 wiring) · estimate ≈ 300 code lines
(+ ≈ 330 test lines, + nf-test) · not statistical

## Spec sections implemented

- §6.2 step row: "`determinants` | py | gene table, features → `determinants_*` | core".
- §7.8.3 outputs: "`tables/determinants_features.tsv`; `determinants_model.tsv`: per model and
  metric: CV R² mean and SD, null p, n_perm, block ΔR² with CI, univariate rows;
  `determinants_attribution.tsv`: per-gene SHAP by block, dominant block, residual, flags."
- §7.8.2 "Caveats (printed in the report)" — this issue stores them as constant text in
  `gippa.determinants.CAVEATS` for M9's report; it does not render them.
- §9.1, §9.2, §9.3 (`determinants_model.tsv`, `determinants_attribution.tsv`), §9.4
  (`abundance_not_explained_by_sequence`, `cai_reference`), §9.5; §6.3 process conventions
  (labels), §6.4 executor equality, §6.6 (LightGBM determinism, fixed `num_threads`).

## Public interface

`src/gippa/steps/determinants.py`
```python
@step("determinants",
      inputs=["gene_table", "genes", "tus", "layers", "replicons"],
      optional_inputs=["regions", "region_genes", "gene_context"],
      outputs=["determinants_features", "determinants_model", "determinants_attribution"])
def determinants_step(ctx: StepContext) -> None: ...
```
`src/gippa/io/determinants_tables.py` — three writers; three Frictionless resources in
`src/gippa/schemas/tables/`.

1. **Inputs.** `gene_table` (detectability output), C2 `genes` (`aa`, `ribosomal_protein`), C5
   `tus`, C6 `layers` (`n_promoter_windows_tu`), C1 `replicons`; optional C7/C8 and M5-09's
   `context/gene_context.parquet` (published by CONTEXT, so DETERMINANTS runs after CONTEXT in
   Nextflow and in topological order locally).
2. **Order.** features (M8-07) → univariate + CV + ΔR² (M8-08) → SHAP + null (M8-09) → writers.
   Streams `determinants.cv`, `determinants.null`, `determinants.bootstrap`. `threads` =
   `cfg.threads` capped by the process CPUs.
3. **`determinants_model.tsv`** rows, sorted by model (lightgbm, elastic_net, univariate), then
   metric (cv_r2_mean, cv_r2_sd, null_p, delta_r2, spearman_rho, auroc), then `feature_or_block`:
   - per model: `cv_r2_mean` (CI NA), `cv_r2_sd`, `null_p` (`n_perm` = B), one `delta_r2` per block
     (t-interval CI); `feature_or_block = all` for the first three;
   - univariate: one `spearman_rho` and one `auroc` per encoded feature column, with CIs;
   - the `flags` column (index decision 12, ADR-A) carries `too_few_rows`, `unestimated`,
     `null_degenerate`.
4. **`determinants_attribution.tsv`**: §9.3 columns in order; `flags` from §9.4
   (`abundance_not_explained_by_sequence`, `cai_reference`).
5. **No abundance** → features written; model table = univariate `auroc` rows only;
   attribution header-only; run flag `layer_missing:abundance`. Module 8 is core (§6.2), so it is
   never switched off by `modules`.
6. **Nextflow.** `workflow/subworkflows/local/determinants.nf` (`GIPPA_STEP` alias `DETERMINANTS`,
   `stepIn`, `kv`, `pick`), called from `main.nf` after `CONTEXT` (it takes
   `CONTEXT.out`'s gene context). `GIPPA_STEP` is `process_single`; `workflow/conf/base.config`
   gains `withName: 'DETERMINANTS' { label … }`-equivalent resources (`process_high`: 16 CPU,
   48 GB) so the null's permutations run in parallel, and `params.threads` reaches the step.
7. **publishDir / local publish**: `tables/determinants_*.tsv` only (M6-08 decision 7).
8. **Caveats** (`CAVEATS`, three strings, verbatim from §7.8.2) and the feature/layer drop list go
   to the step's events (`determinants.caveats`, `determinants.features_dropped`) for the M9 report.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| sim100k (≈ 100 genes, ≈ 50 model rows) | below the 100-row minimum: features table full, univariate rows, one `too_few_rows` row per model, attribution rows with NA predictions — a **passing** pass-1 run |
| no abundance | step 5 |
| no Promotech / OSTIR / Promoter Calculator layer | those features NA and dropped from the models (M8-07); events list them; no block row for an empty block |
| `--publication` | B = 10,000 (§8.2); the run takes ~10× longer, stated in the CLI help |
| CONTEXT skipped (no C7/C8) | `gene_context` absent → `replichore_position` NA; step runs |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `gippa run --executor local` on sim100k | three files under `tables/`; `determinants_model.tsv` has univariate rows and `too_few_rows` rows; exit 0 |
| 2 | the same through Nextflow (`-profile test`) | byte-identical tables |
| 3 | a 1,000-gene simulated genome (M4-07 layout, M4-08 expression) through `--executor local` | CV, null (B from config), ΔR² and SHAP rows present; same seed → identical tables |
| 4 | the same with `threads` 1 vs 8 | identical tables |
| 5 | run without `--abundance` | features + AUROC rows; attribution header-only |
| 6 | `frictionless validate` | valid; `model` enum {lightgbm, elastic_net, univariate}; `metric` enum as §9.3 |
| 7 | nf-test stub of `DETERMINANTS` | three keys emitted, §9.1 names |
| 8 | `sbatch scripts/slurm/gippa_run.sbatch` on sim100k | completes; three tables exist |

## Acceptance criteria

- [ ] `tests/unit/steps/test_determinants_step.py`: examples 1, 3 (with B = 49 to stay in the unit
  budget), 4, 5, 6; `tests/integration` adds example 2.
- [ ] `workflow/tests/determinants.nf.test` (stub) + the registry key test.
- [ ] **Evidence that counts:** one real `sbatch scripts/slurm/gippa_run.sbatch` on sim100k reaching
  the end with all six M8 tables (after M8-06) — job ID in PROGRESS.md; plus the measured runtime
  of the step on the 4,500-gene simulated genome (`benchmarks/baseline.json`).

## Hardening owed (pass 2)

- A golden `determinants_model.tsv` / `determinants_attribution.tsv` on a 1,000-gene fixture
  (sim100k is below the minimum).
- `-resume` across `DETERMINANTS`; contract assertions in the pipeline; the real-tool arm.
- `tables/genes.tsv` gains `abundance_not_explained_by_sequence` when `finalize` exists (Q1).

## Questions for the lead

1. `abundance_not_explained_by_sequence` is a §9.4 **gene** flag, but `tables/genes.tsv` is written
   by `context_stats`, which runs before `determinants`. **Proposed default:** the flag lives in
   `determinants_attribution.tsv` now; a backlog line adds it to `genes.tsv` when the §6.2
   `finalize` step is built (M9/M11). Alternative: `determinants` rewrites `genes.tsv` (two writers
   of one file — rejected).
2. Resources: `withName: 'DETERMINANTS'` at `process_high` (proposed), or a second `GIPPA_STEP_HIGH`
   module with `label 'process_high'`?
3. `--publication` multiplies the null's cost by 10 (B = 10,000, §8.2) — roughly 5–7 CPU-hours per
   model on a 4,500-gene genome if M8-09's estimate holds. **Proposed default:** honour §8.2 as
   written and print the expected runtime; a Besag–Clifford early stop would need an ADR.

## Lead decisions (2026-10-22)

1. **Q1:** `abundance_not_explained_by_sequence` goes in `determinants_attribution.tsv` for now. A backlog line adds it to `genes.tsv` when `finalize` (§6.2) is built. One writer per file.
2. **Q2:** `withName: 'DETERMINANTS'` at `process_high`. No second module.
3. **Q3:** `--publication` follows §8.2 as written (B = 10,000) and logs the expected runtime. A Besag–Clifford early stop needs an ADR and goes to hardening.
4. **Carry-overs from M8-09** (statistician #173, ADR-0034):
   - the `null_p` row carries `null_uncalibrated` on every model, and `unestimated` when the observed R² is NaN;
   - `cv_r2_mean`, the null p's effect, gets `ci_unestimated` (ADR-0031 pattern);
   - `n_perm` counts the non-NaN null values, as `perm_p` does.

   The univariate rows keep `ci_unestimated` (ADR-0032).
5. **Branch:** m8/determinants-step is stacked on m8/attribution (#173). MVP: wire the step end to end through `gippa run` before any polish.
