# M8-16: Sensitivity families "enriched categories" and "determinants"

Milestone M8 · track R · **pass 2** · depends on: M8-06, M8-10 (and merged M4-19 framework) ·
estimate ≈ 160 code lines (+ ≈ 250 test lines) · statistical (light)

## Spec sections implemented

- §8.6 "`basic` (default) reruns the context statistics, operons, GSEA and the determinants CV R²
  (without its permutation null) under `abundance.missing ∈ {min, exclude}`. Under the
  `two_peptide` tier it also counts single-peptide genes as detected." "`full` adds UniOP
  thresholds 0.3 and 0.7, and Top3 instead of iBAQ." "Permutations. Sensitivity reruns use
  B = 199." "Columns: `family`, `variant`, `n_primary`, `n_variant`, `n_shared`, `jaccard`; plus
  `agreement` for the dominant determinant block." "Families: … enriched categories (q ≤ 0.05),
  determinants."
- §6.2 `sensitivity` row ("reruns the `context_stats`, `operons`, `function` and `determinants`
  logic under `min` and `exclude`"); §2.4 G9.

## Public interface

`src/gippa/function/sensitivity.py`, `src/gippa/determinants/sensitivity.py` — two `Family`
callables registered in M4-19's `FAMILIES` (`gippa.context.sensitivity`), keyword-only tables bound
by name (M4-19 decisions 5, 13).
```python
def enriched_categories_family(genes: pd.DataFrame, cfg: GippaConfig, rng: np.random.Generator,
                               *, collections: Sequence[GeneSets], tus: pd.DataFrame) -> set[str]: ...
def determinants_family(genes: pd.DataFrame, cfg: GippaConfig, rng: np.random.Generator,
                        *, features: pd.DataFrame) -> DeterminantsSensitivity: ...
```

1. **Enriched categories.** Flagged IDs = `"{collection}:{term}:{kind}"` with q ≤
   `stats.alpha.enrichment` (0.05), where kind ∈ {`gsea`, `ora_over`, `ora_under`}. Variants
   (decision 25): `min` and `exclude` rerun **GSEA** (abundance-driven) with `n_perm_simple = 199`
   and the multilevel stage at `sampleSize = 101`; the `two_peptide` single-peptide variant reruns
   **ORA** (detection-driven; analytical, no B) and GSEA; `top3` (full) reruns GSEA. ORA under
   `min`/`exclude` is unchanged by construction (detection does not change) and is not rerun.
   **Amendment (M8-17b, ADR-0038; statistician, PR #214 N3):** every GSEA rerun goes through
   `gsea_table`'s engine, which is R fgsea (`RFgseaEngine`, a separate process in the `gippa-fgsea`
   env). Each variant therefore makes its own R calls (a main call and, with GO re-tests, a second
   call) at `n_perm_simple = 199`, and the `sensitivity` step needs the fgsea env whenever it
   reruns GSEA. Its tests inject a fake or `SIMPLE_ENGINE`; a `requires_r` test on Slurm should cover the
   real engine.
2. **Determinants.** Per variant, LightGBM CV R² (5 repeats, **no** permutation null) and the
   dominant block of each top-decile gene. Rows: `family = determinants`, `n_primary` /
   `n_variant` = counts of genes flagged `abundance_not_explained_by_sequence`, `n_shared`,
   `jaccard` of those sets; `agreement` = share of genes in the top decile of **both** runs whose
   dominant block is the same. The variant's CV R² goes to the step's events (the §8.6 column set
   has no R² column; the manifest's sensitivity summary holds it, §6.7).
3. **Wiring.** The `sensitivity` step gains optional inputs (the three function resources it
   needs through `collections`, `tus`, and `determinants_features`); a family whose inputs are
   absent (function module off, no abundance) emits **no** rows and an event — never a row of
   zeros. CONTEXT's `SENSITIVITY` call in `context.nf` gains the channels; `main.nf` order is
   unchanged because `sensitivity` already runs last in CONTEXT (if needed, the call moves after
   FUNCTION and DETERMINANTS — named in the PR).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| `output.sensitivity = off` | header-only table (M4-19 behaviour) |
| function module off | no enriched-categories rows |
| no category significant in primary or variant | row with n = 0, `jaccard` NA (M4-19's convention for 0/0) |
| determinants below the 100-row minimum | no determinants rows; event |
| top decile empty in either run | `agreement` NA |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | primary enriched {cog:K:gsea, go_bp:X:gsea}, `min` {cog:K:gsea} | `n_primary` 2, `n_variant` 1, `n_shared` 1, `jaccard` **0.5** |
| 2 | `two_peptide` tier, a variant counting single-peptide genes | ORA rerun; `ora_over` IDs can change; GSEA rerun |
| 3 | dominant blocks primary (codon, codon, transcription) vs `exclude` (codon, protein, transcription) on the same 3 top-decile genes | `agreement` = **2/3** |
| 4 | same seed twice | identical rows |

## Acceptance criteria

- [ ] `tests/unit/function/test_sensitivity_family.py`, `tests/unit/determinants/test_sensitivity_family.py`:
  examples 1–4 and the degenerate rows.
- [ ] `tests/integration/test_sensitivity_step.py` extended: a run with function and determinants
  on writes the two families.
- [ ] `docs/methods/sensitivity.md` §"Enriched categories", §"Determinants".

## Questions for the lead

1. Variants per kind (proposed above: GSEA under `min`/`exclude`/`top3`/`two_peptide`; ORA under
   `two_peptide` only)?
2. The determinants family's `n_*`/`jaccard` on the `abundance_not_explained_by_sequence` gene
   sets (proposed; §8.6 does not say what the determinants "set" is), with `agreement` for the
   dominant block?

## Lead decisions (2026-10-23)

1. **Q1 yes:** GSEA variants under `min`, `exclude`, `top3` and `two_peptide`; ORA under
   `two_peptide` only.
2. **Q2 yes:** the determinants family's `n_*` and `jaccard` are computed on the
   `abundance_not_explained_by_sequence` gene sets, with `agreement` for the dominant block.
   Reversible; listed for Khoa in PROGRESS questions.
3. **Stacking:** this branch starts from #178 (M8-10, queued). After #178 squash-merges, merge the
   old base tip, then main, and retarget.
4. **Test-author (26687185): gaps 1–12 accepted as proposed.**
   - Engine: keyword-only `engine=` with a default; the step binds the lazy R engine, the tests
     force `SIMPLE_ENGINE`.
   - ORA folds (gap 2): the family takes an optional `folds` keyword and the step an optional
     `tu_folds` input, so n_primary matches `enrichment_detection.tsv`; the implementer adds both.
   - `DeterminantsSensitivity(flagged, dominant, cv_r2)`; agreement is the share of genes in both
     top deciles with the same block, NA when none is shared; a block of None in both runs does not
     count as agreeing.
   - SENSITIVITY moves after DETERMINANTS in the workflow (gap 12).
   - Gap 11 (a finding, for Khoa): under `basic` on the 1 % FDR tier the determinants rerun is a
     no-op (M8-08 rule 1 keeps imputed genes out), so only `single_peptide_detected` and `top3`
     can change its rows. The family stays as specified; listed in PROGRESS questions.
   - Next: the implementer.

## Lead decision 5 (2026-10-23, review round: statistician 6371958269, code review 6372795283)
- **Statistician B1:** `agreement` leaves out genes whose block is `None` in both runs, and is NA
  when none remain. A gene with `None` in only one run still counts as a disagreement. A run
  compared with itself gives 1. This reverses part of decision 4.
- **Statistician B2:** the methods sections "Enriched categories" and "Determinants" come from the
  statistician's draft (lead commit).
- **Code-review B1:** `workflow/conf/base.config` gets a `withName: 'SENSITIVITY_RERUNS:SENSITIVITY'`
  entry sized like `FUNCTION_STEP` (6 GB, 4 h, process_low CPUs). The step runs R fgsea and the
  LightGBM CV.
- **Code-review B2:** new tests for
  - (a) ORA folds;
  - (b) operon_flags and determinants together (the `features` collision);
  - (c) detection calls taken from the variant gene table.

  Also new tests for N2 (`None` in both runs) and N3: a variant below the minimum gets an event
  `sensitivity.determinants_variant_skipped`, and no row and no crash.
- **N1:** threads from `os.sched_getaffinity`, as `steps/determinants._threads` does.
- **N4:** the skip message names its cause.
- **Backlog:** N5–N9. The size is accepted, because the PR is stacked and its split was agreed.
- The basic-mode no-op rows are documented in the methods note (pending Khoa).

## Lead decision 6 (Khoa, 2026-10-24)
- `agreement` covers the genes the features explain, not the top decile. A gene is explained in a
  run when it is a model row not flagged `abundance_not_explained_by_sequence` in that run.
- `DeterminantsSensitivity.dominant` gives a block for every explained model row:
  `agreement_block` on that row's LightGBM SHAP block sums, i.e. the block with the largest
  |sum| (NaN ignored, ties to the first block in `BLOCKS` order, `None` when no |sum| > 0). This
  is Khoa's option (a) (2026-10-24, PR #299 review): M8-08 rule 3 (`_dominant`, the largest
  positive sum) picks a noise block for genes predicted below the mean, and `agreement` then
  floors at about 0.7 with or without signal. `_dominant` itself is unchanged.
- `agreement` is the share of the genes explained in both runs whose dominant block is the same.
  A gene flagged in either run is excluded; a gene with `None` in both runs is excluded; `None` in
  one run only is a disagreement; NA when no gene is left. A run against itself gives 1.
- Genes that are model rows in only one run are outside the denominator. Each variant's
  denominator is `n_agreement` in the event `sensitivity.determinants_cv_r2`.
- `determinants_attribution.tsv` (`dominant_block` still rule 3, top-decile only) and the
  `sensitivity.tsv` columns are unchanged. This supersedes the top-decile wording of item 2,
  decision 4 and example 3 (whose three genes are now read as explained in both runs).
