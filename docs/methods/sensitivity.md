# Sensitivity analysis: stability of flagged sets under analysis choices

Implements `giap.context.sensitivity` and `giap.steps.sensitivity` (GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.2, §5.2; issue M4-19, docs/specs/M4-19-sensitivity-regions.md, lead decisions 1–11).
Statistician sign-off: pending.

## Variants

`output.sensitivity = off` runs nothing and writes no `sensitivity.tsv`. `basic` reruns under
`missing=<p>` for each policy in `abundance.sensitivity` other than the primary `abundance.missing`,
and, under the `two_peptide` tier, `single_peptide_detected`. `full` adds `top3` (unless
`abundance.method` is already `top3`; skipped with a log line when the gene table has no Top3
values, i.e. no peptide intensities). UniOP thresholds are M6's family and are not run here.

A variant's gene table is built from the primary one:

- `missing=<p>`: the values the primary imputed are reset to missing and M4-14 `impute` is rerun
  with policy p on the observed log2 abundances of the detected non-partial genes. Detection and
  `p_detect` are the primary's, so O − E statistics are unchanged; only abundance-based ones move.
- `top3`: the same imputation applied to `top3_log2` instead of the iBAQ values.
- `single_peptide_detected`: genes with status `single_peptide_unestimated` count as detected, the
  primary policy fills their missing abundance, and `p_detect` is refitted (M4-13) on the new
  calls with the primary's generator `for_step("detectability", seed)`, so the TU folds are the
  primary's and only the detection calls differ. The primary is never refitted.

## Permutations and common random numbers

Every rerun uses B = 199 (`stats.permutations`), and the primary flagged set is recomputed at the
same B, so differences reflect the choice and not Monte Carlo error at different B. All calls of a
family (primary and every variant) start from one generator state, `for_step("sensitivity/<family>",
seed)` (common random numbers; adding a family or a variant changes no other family's draws).

For the `regions` family the flagged set is the IDs `<region_id>:<contig>:<stat>` of the M4-17
tests with `null_model = circular_shift`, a p-value, and q ≤ `alpha.context` (a `split_contigs`
region gives one ID per contig). With `context.exhaustive_shift_max_genes = 20000`, every replicon
of up to 20,000 genes uses all shifts, so the recomputed primary equals the set reported in
`regions_stats.tsv`. Above 20,000 genes the shifts are sampled, the p-value floor is 1/200 = 0.005,
and `n_primary` may differ from the count in `regions_stats.tsv`.

## Output

One row per (family, variant): `n_primary` = |P|, `n_variant` = |V|, `n_shared` = |P ∩ V| and
`jaccard` = |P ∩ V| / |P ∪ V| (`NA` when both sets are empty). `agreement` is reserved for the
dominant determinant block (M8) and is `NA` for other families.
