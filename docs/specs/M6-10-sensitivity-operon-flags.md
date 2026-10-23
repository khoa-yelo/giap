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

## Lead decisions (2026-10-22)

1. **UniOP-threshold variants:** the proposed default. `Variant` gains an optional config patch
   that is applied before the pair frame is built; operons are not special-cased. Existing variants
   (no patch) are byte-identical.
2. **`jaccard` of two empty sets:** NA, per the M4-19 convention, flagged as in M4-19.
3. **Test-author resolutions (588b219c) accepted:**
   - the C2 gene model binds as `gene_models`;
   - C3 `features` is required by the family;
   - the variant names are `pair_threshold=<t:g>` in the shared `variants()`;
   - the `regions` and `islands` rows under the threshold variants are left unpinned;
   - example 6 has three rows;
   - a relative scale logs "relative".
4. **Step inputs:** C2 `genes`, `pairs_raw`, `features` and the three evidence layers are
   **optional** step inputs. When `genes`, `pairs_raw` or `features` is absent, the `operon_flags`
   family is skipped and the step logs `sensitivity.operon_flags_skipped`, so the M4-19/M5 step
   tests that pass only gene_table, C1, C5, C7 and C8 keep passing.
5. **Lead edits to existing tests** in tests/unit/context/test_sensitivity.py: `Variant` gains
   `config_patch`, `FAMILIES` may hold `operon_flags`, and `full` lists the two
   `pair_threshold` variants.
6. **Runtime** ("four reruns < 2× a single run") is not pinned; this is a backlog item.
7. **Implementer round 1 (a1c218f7):**
   - Without the pair tables, the `pair_threshold` variants are dropped along with the family,
     because their only patch is the UniOP threshold. Accepted: it keeps
     `test_full_adds_top3_with_peptide_intensities`.
   - CONTEXT gains a 12th input, `layers_for_pairs`; the lead edits context.nf.test to add
     `input[11] = channel.empty()`.
   - Behaviour 5 (one pair frame per threshold, reused) is deferred to the decision-6 backlog item.
8. **Review #243 (comment, 2026-10-23):**
   - **B1 (lead, done):** `test_m6_10_sensitivity_has_operon_flags_rows` in
     tests/integration/test_m6_08_operons_run.py checks that a real `giap run` writes
     `operon_flags` rows (n 0, J NA on sim100k). The backlog line is dropped.
   - **B2 (accepted, option a):** a variant that only patches config runs only for the families
     that read the patched key.
     - `Variant` gains `families: tuple[str, ...] | None`, where `None` means every family. The
       `pair_threshold=*` variants carry `("operon_flags",)`.
     - `regions` and `islands` therefore get no `pair_threshold` rows: J = 1 there would be true
       by construction, and those reruns cost four permutation/bootstrap runs each.
     - A `Variant.patch_only` property replaces the two hand-written conditions (N6).
   - **B3:** if any one of `genes`, `pairs_raw` or `features` is absent, the family is skipped
     with one `sensitivity.operon_flags_skipped` event. It never raises a KeyError. The test is
     parametrised over each single missing key and over all of them.
   - Nits N1–N5 go to the backlog.
9. **Test-author B2/B3 tests (9f77b30d):** the B3 skip path already held, and the new tests now pin
   it (mutant M6 is killed). For B2 the lead edited two existing tests:
   - the `Variant` field list gains `families`;
   - `test_patch_variants_reach_the_family` registers its recorder as `operon_flags`.
