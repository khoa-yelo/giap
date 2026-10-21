# M6-10: Sensitivity family `operon_flags`

Milestone M6 · track R · **pass 2** · depends on: M6-07, M6-08 (+ merged M4-19 `FAMILIES`) ·
estimate ≈ 110 code lines (+ ≈ 150 test lines) · statistical (light)

## Spec sections implemented

- §8.6: "`basic` (default) reruns the context statistics, **operons**, GSEA and the determinants CV
  R² … under `abundance.missing ∈ {min, exclude}`. Under the `two_peptide` tier it also counts
  single-peptide genes as detected"; "`full` adds **UniOP thresholds 0.3 and 0.7**, and Top3
  instead of iBAQ"; "Permutations. Sensitivity reruns use B = 199"; "Families: … **operon flags
  (pair IDs)** …"; the `sensitivity.tsv` columns.
- §2.4 G9; M4-19's merged framework (`Family`, `Variant`, `flagged_ids`, `jaccard`,
  `sensitivity_rows`, `FAMILIES`).

## Behaviour

1. Register `FAMILIES["operon_flags"]` with the M4-19 signature `(genes, cfg, rng)` plus the
   keyword-only tables the family needs (`pairs_raw`, `tus`, `replicons`, `genes` and the three
   evidence layers), bound by the step exactly as `regions_family` is.
2. The family's flagged set = the `pair_id` values present in `operon_flags.tsv` for that variant.
   Because M6-04 decision 3 makes `pair_id` independent of the UniOP threshold, the sets are
   comparable across every variant, including the `full` UniOP variants.
3. Variants:
   - `basic`: `abundance.missing ∈ {min, exclude}`, and `two_peptide` when that tier is in force —
     all three already exist in M4-19's `variants()`;
   - `full`: `annotation.uniop.pair_threshold ∈ {0.3, 0.7}` (new variant kind: it re-derives the
     classes, not the abundances) and `abundance.method = top3`.
4. Reruns use `stats.permutations = 199` (M4-19's `SENSITIVITY_B`). Module 4 has no permutation
   test, so the only effect is on anything downstream that does; the value is still set, for
   consistency.
5. The rerun **must not** rebuild the whole pair geometry per variant: build the pair frame once
   per distinct `pair_threshold` and reuse it across the abundance variants.
6. Output: one `sensitivity.tsv` row per (family, variant) with `n_primary`, `n_variant`,
   `n_shared`, `jaccard`; `agreement` NA (it belongs to the determinants family).
7. **Reaching `giap run`:** the merged `sensitivity` step (registered by M4-19 decision 16 and
   wired by MVP-01's CONTEXT subworkflow) gains the new inputs; no new subworkflow is needed, but
   `context.nf` must pass the pair-related channels through. That wiring is part of this issue, not
   a follow-up.

## Degenerate cases

| Situation | Behaviour |
|---|---|
| no flagged pair in the primary run (sim100k) | `n_primary` 0, `n_variant` 0, `jaccard` NA (M4-19's convention for an empty union) |
| variant that cannot run (relative scale) | row written with NA `jaccard` and the reason in the log |
| UniOP threshold 0.7 leaves < 20 OP pairs | the variant's discordance is `unestimated`; `n_variant` counts only the missed-operon flags |
| `--sensitivity off` | the step does not run at all (existing behaviour) |

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | primary run with 40 flagged pairs; `missing = min` gives 38, 36 shared | `jaccard` = 36/42 = **0.857** |
| 2 | `missing = exclude` drops every imputed-abundance pair | `n_variant` ≤ `n_primary`; `jaccard` ≤ 1 |
| 3 | `pair_threshold = 0.3` (more OP pairs, fewer NO pairs) | `novel_operon_candidate` flags fall, discordance flags rise; `pair_id`s are the same objects, so the Jaccard is meaningful |
| 4 | `pair_threshold = 0.7` | the mirror of example 3 |
| 5 | `two_peptide` tier variant | more detected genes → more both-detected tests |
| 6 | sim100k | all four rows with `jaccard` NA and no error |
| 7 | same seed twice | identical rows |

## Acceptance criteria

- [ ] `tests/unit/context/test_sensitivity_operons.py`: examples 1–7 on hand-built frames.
- [ ] Property: the family is pure and order-independent; `jaccard` ∈ [0, 1] or NA.
- [ ] Property: the primary set equals `operon_flags.tsv` of the main run (no recomputation drift).
- [ ] Integration: `tables/sensitivity.tsv` from a real `giap run` contains the `operon_flags`
  rows — the "it must reach `giap run`" criterion.
- [ ] Runtime: the four reruns add < 2× the single-run cost of Module 4 on 4,500 genes.

## Hardening owed

- The `full` variant set is only exercised when `--sensitivity full` is requested; add it to the
  nightly end-to-end so it does not rot.
- If the Jaccard against the UniOP variants turns out to be very low, the report must say that the
  flag set is threshold-sensitive — a text item for M9.

## Questions for the lead

1. The UniOP-threshold variants change the **classes**, not the abundances, so they are a new kind
   of variant in M4-19's framework. **Proposed default:** extend `Variant` with an optional config
   patch applied before the pair frame is built, rather than special-casing operons.
2. `jaccard` when both sets are empty: NA (proposed, M4-19's convention) or 1.0? **Proposed
   default:** NA — "no flags either way" is not stability evidence.
