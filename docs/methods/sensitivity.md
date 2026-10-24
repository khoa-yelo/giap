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

## Enriched categories

For the `enriched_categories` family (M8-16) the flagged set is the IDs `<collection>:<term>:<kind>`
of every row with q ≤ `stats.alpha.enrichment` (0.05). The kind is one of:
- `gsea`: M8-05 abundance GSEA, BH per collection;
- `ora_over` or `ora_under`: M8-04 detection ORA, BH per ontology and direction (§8.1).

Rows flagged `redundant`, `genomic_cluster` or `detection_confounded` count like any other row. A
`gsea` ID has no direction, so a set whose ES changes sign under a variant still counts as shared.

The collections are built once from the primary gene table, and every run reuses them. The ORA
uses the primary's TU folds (`tu_folds`) when they are given, as the function step does. ORA
p-values are analytical (Wallenius), so the ORA IDs change only when the detection calls change,
i.e. under `single_peptide_detected`. Under `missing=<p>` and `top3` they equal the primary's by
construction.

GSEA is rerun under every variant (`missing=<p>`, `top3`, `single_peptide_detected`). It goes
through R fgsea (`fgseaMultilevel`, ADR-0038) at `n_perm_simple = 199`, `sample_size = 101` and
`eps = 0`, whatever `function.gsea` says. Each variant makes its own R calls, so the step needs the
`giap-fgsea` env.

The family is skipped, with the event `sensitivity.enriched_categories_skipped`, when `function` is
not in `modules` or no function resource is given. It then writes no rows.

**Permutation depth.** With `eps = 0` the multilevel stage estimates p with no 1/(1 + B) floor. A
row can reach q ≤ 0.05 only with p ≤ 0.05·k/m. fgsea refines such small p-values in the multilevel
stage, whose precision depends on `sample_size` and not on `n_perm_simple`. A p above 0.05 cannot
lower any BH q to 0.05 or below. So 199 simple permutations hardly change which rows are flagged.

In the statistician's probe (PR #265 review, Slurm 1876268), the design was:
- 1,000 TUs;
- 150 sets of 15 TUs, 15 of them shifted;
- 20 sets of 40 TUs, 4 of them shifted;
- 6 worlds × 21 streams on identical data.

The mean pairwise Jaccard between streams was the same at `n_perm_simple` 199 and 1,000 in 4 of 6
worlds, and within 0.03 in the other two.

**Monte Carlo noise.** The flagged set still varies between streams. A set near the BH cut can
flip, and one flip moves the cut for the whole collection. In the same probe, two streams on
identical data gave different sets in 26–72 % of pairs, and the mean pairwise Jaccard per world was
0.57–0.90.

The recomputed primary uses `for_step("sensitivity/enriched_categories", seed)` at B = 199. It is
therefore an independent realisation, not the set in `enrichment_gsea.tsv` (stream
`function.gsea`, `n_perm_simple` 1,000), and `n_primary` can differ from that table's count.

Common random numbers hold only while a variant ranks the same number of TUs: the tie order, the R
seed and fgsea's sampled index sets are then the same. This holds for `missing=min`, because imputed
genes stay ranked. It fails for:
- `missing=exclude`, when a TU loses all its scoring genes;
- `top3` and `single_peptide_detected`, when a TU enters or leaves the ranking.

The GSEA differences of those variants include Monte Carlo noise of the size above.

## Determinants

For the `determinants` family (M8-16; §8.6 "the determinants CV R² (without its permutation null)")
each run repeats the M8-10 chain without its TU-permutation null:
- model rows: detected genes with a finite, non-imputed log2 abundance, ribosomal proteins excluded
  (M8-08 rule 1);
- outer `GroupKFold(5)` TU folds, `determinants.repeats` repeats;
- LightGBM and elastic-net CV;
- LightGBM SHAP summed by block;
- the attribution table.

Features come from `determinants_features.tsv`. `detected`, `abundance_log2` and
`abundance_imputed` come from the run's gene table. The folds are drawn from
`for_step("sensitivity/determinants", seed)`, and LightGBM is seeded with `seed`. The recomputed
primary is therefore an independent realisation of the fold assignment, not
`determinants_attribution.tsv` (stream `determinants.cv`).

- **Set.** The flagged set is the locus tags flagged `abundance_not_explained_by_sequence`: model
  rows whose out-of-fold residual is more than 2 SD from 0. The test is two-sided, uses the run's
  own residual SD, and flags about 4 % of model rows. A variant that changes the CV R² moves the
  threshold for every gene. This set is a reversible default (M8-16 Q2) because §8.6 does not
  define the determinants set.
- **Agreement.** A gene is *explained* in a run when it is a model row and is not flagged
  `abundance_not_explained_by_sequence` in that run. "Explained" means only "not flagged at 2
  residual SD": it is not a claim that the features predict the gene well. Its block is the block
  with the largest **absolute** LightGBM SHAP sum (ties go to the first block in the order
  transcription, translation, codon, protein, context); it is `None` when every block sum is 0.
  The absolute value matters: for a gene predicted below the mean, the block that drives the
  prediction carries a large negative sum, and M8-08 rule 3 (the largest positive sum) would pick
  whichever noise block happens to sit just above 0. Under rule 3, `agreement` floored at about 0.7
  whether or not the features carried signal (PR #299 review). `agreement` is computed over the
  genes explained in **both** runs: a gene flagged in either run is left out, and so is a gene
  whose block is `None` in both runs. A gene that is a model row in only one run is left out too.
  This matters for `single_peptide_detected` and `top3`, which move genes in or out of the model
  rows; those genes never enter the denominator. `agreement` is the share of the remaining genes
  whose block is the same. A gene with a block in one run and `None` in the other counts as a
  disagreement. `agreement` is NA when no such gene remains. A run compared with itself gives 1
  (NA if none of its explained genes has a block). This follows Khoa's decisions of 2026-10-24
  (M8-16 lead decision 6, option (a)); before them, `agreement` covered the top decile only and
  used rule 3. `determinants_attribution.tsv` is unchanged: its `dominant_block` keeps rule 3 and
  is still reported for the top-decile genes only.
- **Denominator.** The number of genes behind each variant's `agreement` (`n_agreement`) goes to
  the event `sensitivity.determinants_cv_r2`, next to the CV R². §8.6 has no column for it.
- **CV R².** The mean LightGBM CV R² of the primary and of every variant goes to the event
  `sensitivity.determinants_cv_r2` (`data.cv_r2`; `data.n_agreement` holds the denominators).
  §8.6 has no R² column, and the manifest's sensitivity block is built from `sensitivity.tsv`
  only, so the event is the one record of both.
- **Minimum.** A run with fewer than 100 model rows or 25 TUs is not fitted (M8-10 rule 8). If the
  primary is below the minimum, the family writes no rows and emits
  `sensitivity.determinants_skipped`. If a variant is below it, that variant has no row.

**What the variants can change.** Model rows exclude imputed values. So `missing=min` and
`missing=exclude` keep the model rows, the TU list, the folds and the LightGBM seed of the primary.
Their rows are exact copies of the primary (`jaccard` 1, or NA when nothing is flagged;
`agreement` 1) and **say nothing about robustness**. Only two variants can change the determinants
rows:
- `single_peptide_detected`, which adds the single-peptide genes, with their one-peptide
  abundance, as model rows;
- `top3` under `full`, which changes every abundance and can move genes in or out of the model
  rows.

**Noise.** A variant that changes the model rows also changes the TU list. `tu_folds` then
reassigns the folds, so common random numbers are lost, and the difference includes
fold-assignment noise.

**Fold-noise reference.** The determinants calibration test
(`tests/calibration/test_determinants_agreement_calibration.py`; PR #299 review) reruns planted
worlds of 400 genes (log2 abundance = signal × CAI + noise) under fold streams 1–5 against
stream 0 on identical data, so only the fold assignment differs:
- with negligible signal (mean LightGBM CV R² 0.01–0.08), `agreement` was 0.73–0.86;
- with strong signal (CV R² about 0.55–0.59), it was 0.89–0.99.

On every world, each strong-signal pair beat each no-signal pair; the test gates that separation,
not the levels. An `agreement` of about 0.8 is therefore what fold noise alone gives when the
features explain nothing; values near 1 need signal. The Jaccard figures of the earlier probe
(PR #265 review, Slurm 1876267: 11 fold streams, CV R² 0.05–0.77) still hold for the flagged set:
0.67–1.0 against stream 0, median about 0.77.

A determinants Jaccard of about 0.7 or more is therefore within fold-assignment noise. Read it
together with `n_primary`, `n_variant` and `n_shared`. The denominator of `agreement` (the genes
explained in both runs, most of the model rows) is not in the table; it is `n_agreement` in the
event `sensitivity.determinants_cv_r2`. Neither `jaccard` nor `agreement` has a confidence
interval.

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
  through sampled-shift noise. `single_peptide_detected` keeps the observed abundance values and
  the TU folds of the detectability fit.

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
