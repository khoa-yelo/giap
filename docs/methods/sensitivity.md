# Sensitivity analysis: stability of flagged sets under analysis choices

Implements `giap.context.sensitivity` and `giap.steps.sensitivity` (GIAP_SPEC.md §8.6, §2.4 G9,
§7.3.2, §5.2; issue M4-19, docs/specs/M4-19-sensitivity-regions.md, lead decisions 1–18).

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
- `single_peptide_detected`: genes with status `single_peptide_unestimated` count as detected.
  M4-14 keeps observed abundances on undetected genes, so these genes enter HL with their observed
  one-peptide abundance (noisier than a multi-peptide value); the primary policy fills only those
  with no intensity at all. `p_detect` is refitted (M4-13) on the new calls with the primary's
  generator `for_step("detectability", seed)`, so the TU folds are the primary's and only the
  detection calls differ. The primary is never refitted.

## Permutations and common random numbers

Every rerun uses B = 199 (`stats.permutations`), and the primary flagged set is recomputed at the
same B, so differences do not come from Monte Carlo error at different B. All calls of a family
(primary and every variant) start from one generator state, `for_step("sensitivity/<family>",
seed)`; adding a family or a variant changes no other family's draws.

These common random numbers matter only when shifts are sampled:

- **Up to `context.exhaustive_shift_max_genes` (20,000) genes per replicon**, which covers every
  bacterial replicon, the region tests use all shifts. The flagged sets then use no randomness at
  all (the generator feeds only the bootstrap CIs, which the flagged sets ignore), so every
  difference between primary and variant comes from the variant.
- **Above 20,000 genes** the shift sample (`choose_shifts`) and both TU-block bootstrap CIs of every
  region are drawn from the one family stream, in region order, and the HL bootstrap runs only when
  HL is finite. Once one region's data differ between primary and variant, every later region's
  draws differ, so the common random numbers hold only up to that region. Differences then include
  Monte Carlo error. In the statistician's probe (giapsim fixture with sampled shifts; the variant
  changed one region's abundances), the O − E p of unchanged regions moved by up to 0.075. A
  per-region shift stream in M4-17 `test_region` would remove this (backlog).

For the `regions` family the flagged set is the IDs `<region_id>:<contig>:<stat>` of the M4-17
tests with `null_model = circular_shift`, a p-value, and q ≤ `alpha.context` (a `split_contigs`
region gives one ID per contig). With `context.exhaustive_shift_max_genes = 20000`, every replicon
of up to 20,000 genes uses all shifts, so the recomputed primary equals the set reported in
`regions_stats.tsv`. Above 20,000 genes the shifts are sampled, the p-value floor is 1/200 = 0.005,
and `n_primary` may differ from the count in `regions_stats.tsv`.

## Interpretation

- Flagged sets are small (at most 2 × regions × contigs IDs, usually fewer than 20), so a single
  flip moves Jaccard by a large step. Read `jaccard` together with `n_primary`, `n_variant` and
  `n_shared`.
- Changes concentrate in flags whose q is near `alpha.context` (0.05). A flag far below it that
  disappears under a variant is the signal worth inspecting.
- BH is recomputed within each rerun's (region_type, stat) family, so each set is at its own FDR.
  A variant can move a region's q, and so flip its flag, even when that region's p is unchanged.
- What each variant holds fixed: every variant keeps the region definitions (C7 / C8), the replicon
  topology (C1), the detectability covariates and the primary's B. `missing=<p>` and `top3` also
  keep detection and `p_detect`, so their O − E flags can change only through BH or sampled-shift
  noise. `single_peptide_detected` keeps the observed abundance values and the TU folds of the
  detectability fit.

## Output

One row per (family, variant): `n_primary` = |P|, `n_variant` = |V|, `n_shared` = |P ∩ V| and
`jaccard` = |P ∩ V| / |P ∪ V| (`NA` when both sets are empty). `agreement` is reserved for the
dominant determinant block (M8) and is `NA` for other families.

Signed off by the statistician (PR #117 review).
