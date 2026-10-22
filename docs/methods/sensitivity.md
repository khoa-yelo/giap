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

## Islands

For the `islands` family (M5-10) the flagged set is `<direction>:<locus_tag>` for every gene of
every island that M5-06 reports (an island is reported only when its FWER p ≤ `alpha.context`).
The genes are the island's best-window gene list (M5-06 Q3), so a change in how clusters merge
shows up as a changed set. The direction is part of the ID: a gene that moves from a high to a low
island counts as a change, and a gene in both a high and a low island gives two IDs. The family
reruns the M5-06 pipeline (`tu_series` per C1 replicon on the C5 TUs, then `scan_replicons` and
`islands_frame`) with both `stats.bootstrap` and `stats.permutations` set to 199 for the primary
and every variant (the `regions` family keeps the user's `stats.bootstrap`). All calls start from
the one family stream `for_step("sensitivity/islands", seed)`. The FWER null is a block
bootstrap, so the islands sets always use randomness. With B = 199 the smallest attainable p is
1/200 = 0.005.

The recomputed primary is an independent Monte Carlo realisation, not `islands.tsv` at a coarser
p floor. It uses its own stream and B = 199 for both the block bootstrap and the l-hat test, so an
island with p near `alpha.context` can differ from `islands.tsv` in either direction, not only
through the 0.005 floor.

Common random numbers hold only in part. `scan_replicons` draws the l-hat permutations of every
scanned replicon before any bootstrap draw, and how many it draws depends on each replicon's
number of scored TUs T. A variant that changes T, or whether a replicon is scanned, anywhere in
the genome therefore changes the bootstrap draws of every replicon, including earlier ones. A
change of the block length l' at replicon r changes the draws of r and of later replicons only.
The triggers are `missing=exclude` (a TU losing every scored gene, since detected genes without an
abundance are not scored, M5-06 step 1) and the `single_peptide_detected` refit of `p_detect`.
With `context.islands.block_len` fixed, only T matters. When the draws differ, the set difference
includes Monte Carlo error (see Interpretation).

## Interpretation

- Flagged sets are small (at most 2 × regions × contigs IDs, usually fewer than 20), so a single
  flip moves Jaccard by a large step. Read `jaccard` together with `n_primary`, `n_variant` and
  `n_shared`.
- Changes concentrate in flags whose q is near `alpha.context` (0.05). A flag far below it that
  disappears under a variant is the signal worth inspecting.
- BH is recomputed within each rerun's (region_type, stat) family, so each set is at its own FDR.
  A variant can move a region's q, and so flip its flag, even when that region's p is unchanged.
- What each variant holds fixed: every variant keeps the region definitions (C7 / C8), the replicon
  topology (C1), the detectability covariates and B = 199 (the primary is recomputed at the same
  B). `missing=<p>` and `top3` also keep detection and `p_detect`, so every O − E p, q and flag is
  identical to the primary's when shifts are exhaustive; above 20,000 genes they can differ only
  through sampled-shift noise. `single_peptide_detected` keeps the observed abundance values and the TU folds of the
  detectability fit.

- Islands (M5-10) are flagged by one two-sided genome-wide FWER (both directions), not
  by BH. Their IDs are
  genes, so one island flipping moves the Jaccard by that island's share of the genes.
- Islands whose p lies within about ±0.03 of `alpha.context` can flip on Monte Carlo noise alone.
  In the statistician's probe (identical data, 30 independent streams; PR #184 review), a planted
  island at 2.2 SD gave Jaccard 0 against stream 0 in 59 % of streams; at 3.0 SD it was always 1.
  At B = 199 the Monte Carlo SD of p at 0.05 is about 0.015. An island lost under a variant whose
  primary p ≤ 0.01 is the real signal to inspect. The Jaccard has no confidence interval.

## Output

One row per (family, variant): `n_primary` = |P|, `n_variant` = |V|, `n_shared` = |P ∩ V| and
`jaccard` = |P ∩ V| / |P ∪ V| (`NA` when both sets are empty). `agreement` is reserved for the
dominant determinant block (M8) and is `NA` for other families.

Signed off by the statistician, 2026-10-20 (PR #117 reviews
https://github.com/khoa-yelo/giap/pull/117#issuecomment-6304191087 and
https://github.com/khoa-yelo/giap/pull/117#issuecomment-6304219014), conditional on the
shared-stream limitation above sampled-shift sizes (backlog: per-region shift streams).
