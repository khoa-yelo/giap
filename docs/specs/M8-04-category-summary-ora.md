# M8-04: `category_summary` and detection enrichment (Wallenius ORA, BH families, elim, `parent_group`)

Milestone M8 · track R · **pass 1** · depends on: M8-01 (Wallenius, elim), M8-03 (collections) ·
estimate ≈ 300 code lines (+ ≈ 400 test lines) · **statistical** (label `stats`)

## Spec sections implemented

- §7.7.2 "One row per COG category, KEGG pathway and GO slim term, with: n_genes, n_detected;
  expected (Σ p_detect), O − E with a TU-bootstrap CI; median log2 abundance of detected genes;
  coverage. Mass fraction is deferred (improvement 11)."
- §7.7.3 in full (setup, p-values, "Filters. Category sizes are limited to `function.ora_size`
  ([5, 500]). BH runs within each ontology and direction.", "GO redundancy (elim …) Report
  `p_classic` and `p_value` (elim); BH applies to the elim p-values. Group surviving terms under
  their parents (`parent_group`).").
- §8.1 row "Category detection | category | Wallenius with p_detect weights | BH per ontology and
  direction | 0.05"; §8.2 "Effect sizes with 95 % CIs accompany every p or q"; §9.2 test columns;
  §9.3 `category_summary.tsv` and `enrichment_detection.tsv`.

## Public interface

`src/giap/function/ora.py` and `src/giap/function/summary.py` (pure)
```python
def category_summary(genes: pd.DataFrame, colls: Sequence[GeneSets], B: int,
                     rng: np.random.Generator) -> pd.DataFrame: ...
def detection_enrichment(genes: pd.DataFrame, colls: Sequence[GeneSets],
                         size: tuple[int, int], alpha: float) -> pd.DataFrame: ...
```

1. **`category_summary`** — collections `cog`, `kegg`, `go_slim` only (§9.3 enum). Per category:
   `n_genes` (eligible members), `n_detected`, `expected` = Σ p_detect, `o_minus_e` with a 95 %
   TU-block bootstrap CI (`giap.stats.tu_block_bootstrap` over the category's TUs, B =
   `stats.bootstrap`, the M4-16 `basic_rows` construction: per-TU sums of d − p resampled),
   `median_abundance_log2` over detected members with a non-NA abundance (imputed values included,
   as everywhere under the primary policy), `coverage` of the collection. Every category is
   reported, whatever its size (descriptive table, no test). Rows sorted by collection
   (cog, kegg, go_slim), then category.
2. **`detection_enrichment`** — all six collections. Per term and direction (`over`, `under`):
   M, n, N, k over the eligible universe (decision 2: annotated or not); `odds_weight` = ω (M8-01);
   `p_classic` from `wallenius_tails`; `effect` = O − E = k − Σ_{C} p_detect, with CI NA in pass 1
   (decision below, Q6); `stat` = k; `method = wallenius_ora`, `null_model = wallenius_pdetect`,
   `n_perm = 0`.
3. **Size filter.** Terms with n outside `function.ora_size` ([5, 500]) are not tested: no rows.
   (Under elim, the filter is applied to the term's **full** n; a term whose elim-reduced set
   falls below 5 is still tested on what remains, decision Q3.)
4. **elim** (`go_bp`, `go_mf`, `go_cc`, direction `over` only — decision 4): `p_value` = the elim p
   from M8-01's walker with threshold 0.01 and the Wallenius test re-evaluated on the reduced set
   (ω recomputed on the reduced set; M, N fixed). For `under`, `cog`, `kegg` and `go_slim`:
   `p_value = p_classic`.
5. **BH families** = (collection, direction), on `p_value`, with `giap.stats.bh`; NaN p stay out
   and carry `flags = unestimated`. `q_value` per row.
6. **`parent_group`** (GO namespaces only; decision 11): for a row with `q_value ≤
   stats.alpha.enrichment`, the most general **significant** ancestor in the same direction and
   namespace (by depth, ties by ID); the term's own ID if it has none. NA otherwise and for
   non-GO collections.
7. **`coverage`** on every row (the collection's); `n_universe` = M; `n_in_set` = n;
   `n_detected` = N; `k_detected_in_set` = k (§9.3 names).

## Degenerate cases

| Situation | Behaviour |
|---|---|
| ω undefined (C = universe, or Σ p outside = 0) | row with p NaN, `flags = unestimated`, out of BH |
| no detected gene (N = 0) | every p = 1 (M8-01); rows written; nothing significant |
| `detectability_constant` (p_detect = overall rate) | ω = 1 exactly → central hypergeometric; the run flag from M4-13 is carried in the manifest, not repeated per row |
| an empty collection | no rows; coverage still in `category_summary` for cog/kegg/go_slim |
| elim reduces a term to 0 genes | `p_value` NaN, `flags = unestimated`; `p_classic` kept |
| category of 1 TU (operon-sized) | O − E CI from 1 TU is degenerate → CI NA (the M4-16 `ci_few_units` rule: < 10 TUs gives NA and flag) |
| sim100k | most terms below n = 5; `enrichment_detection.tsv` may be nearly empty — a **passing** pass-1 result |

## Examples

Universe of 10 eligible genes, p_detect 0.5 each except g1–g2 at 0.9; detected = g1, g2, g3, g5.

| # | Input | Expected |
|---|---|---|
| 1 | C = {g1, g2, g3, g4, g6} | n = 5, k = 3, N = 4, M = 10; ω = mean(0.9, 0.9, 0.5, 0.5, 0.5)/0.5 = 0.66/0.5 = **1.32**; expected = 3.3, O − E = **−0.3** |
| 2 | same with `ora_size` [6, 500] | no row |
| 3 | `category_summary` of COG letter K = C of example 1 | n_genes 5, n_detected 3, expected 3.3, o_minus_e −0.3, CI from the TU bootstrap (NA here: < 10 TUs, flag) |
| 4 | toy GO DAG of M8-01 example 8 with A1 elim-significant | A's `p_value` from {g4, g5}, `p_classic` from its full set; A1 and A rows in one BH family (go_bp, over) |
| 5 | A1 q ≤ 0.05 and A q ≤ 0.05 | `parent_group(A1) = A`'s most general significant ancestor = **A** (R not significant); `parent_group(A) = A` |
| 6 | the same term, direction `under` | `p_value = p_classic`, separate BH family |
| 7 | ω undefined for one term | its row has `p_value` NA, `flags unestimated`; the other q-values are those of BH without it |
| 8 | the same input twice with the same seed | identical tables (the bootstrap draws from `function.bootstrap`; ORA is deterministic) |

## Acceptance criteria

- [ ] `tests/unit/function/test_ora.py` and `test_summary.py` cover examples 1–8 and every
  degenerate row.
- [ ] Property: relabelling genes outside C leaves C's row unchanged; adding a non-member gene
  with p_detect equal to the outside mean leaves ω unchanged.
- [ ] Property: BH families are exactly (collection, direction); no NaN p enters one.
- [ ] Property: a GO term with no elim-significant descendant has `p_value == p_classic` (elim can
  move p either way otherwise, so no inequality is asserted).
- [ ] `docs/methods/function.md` §"Detection enrichment": the families, elim, the effect size.

## Calibration or attack work owed

- `ora` suite and attack 2 (M8-15). This issue makes both runnable on a frame without the step:
  the functions take a gene frame and `GeneSets`, nothing else.

## Hardening owed (pass 2)

- The O − E CI on enrichment rows (Q6) once the statistician confirms the TU bootstrap is the
  right interval there (it is the §7.7.2 one).
- Golden Wallenius p vs goseq through this frame (M8-13 compares the engine; a second check that
  the frame's counts equal goseq's per-category counts on the golden inputs).

## Questions for the lead

1. Universe: §7.7.3 says "the eligible genes" (annotated or not); goseq and topGO default to
   genes with ≥ 1 annotation in the ontology. **Proposed default:** the spec; the golden goseq call
   uses `use_genes_without_cat = TRUE`; M8-15 records the annotated-universe arm as a comparison.
2. elim for `under`? **Proposed default:** no (topGO's elim is an over-representation method);
   `under` rows carry `p_value = p_classic`.
3. Size filter under elim: on the full n (proposed) or on the reduced n?
4. Collection values: `cog`, `kegg`, `go_bp`, `go_mf`, `go_cc`, `go_slim` in
   `enrichment_detection.tsv`. **Proposed default:** yes (§9.3 gives no enum there).
5. `parent_group` = most general significant ancestor (proposed), or the nearest significant
   parent?
6. `effect`/`ci_low`/`ci_high` on enrichment rows: O − E with CI NA in pass 1 (proposed), or the
   TU-bootstrap CI now (it doubles the bootstrap cost over ~3,000 GO terms)?
