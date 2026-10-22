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

`src/giap/function/sensitivity.py`, `src/giap/determinants/sensitivity.py` — two `Family`
callables registered in M4-19's `FAMILIES` (`giap.context.sensitivity`), keyword-only tables bound
by name (M4-19 decisions 5, 13).
```python
def enriched_categories_family(genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator,
                               *, collections: Sequence[GeneSets], tus: pd.DataFrame) -> set[str]: ...
def determinants_family(genes: pd.DataFrame, cfg: GiapConfig, rng: np.random.Generator,
                        *, features: pd.DataFrame) -> DeterminantsSensitivity: ...
```

1. **Enriched categories.** Flagged IDs = `"{collection}:{term}:{kind}"` with q ≤
   `stats.alpha.enrichment` (0.05), where kind ∈ {`gsea`, `ora_over`, `ora_under`}. Variants
   (decision 25): `min` and `exclude` rerun **GSEA** (abundance-driven) with `n_perm_simple = 199`
   and the multilevel stage at `sampleSize = 101`; the `two_peptide` single-peptide variant reruns
   **ORA** (detection-driven; analytical, no B) and GSEA; `top3` (full) reruns GSEA. ORA under
   `min`/`exclude` is unchanged by construction (detection does not change) and is not rerun.
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
