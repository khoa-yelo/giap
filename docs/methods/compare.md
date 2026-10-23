# `giap compare`: detection, abundance agreement and result stability

Implements `giap.compare.{detection,abundance,stability}` and `giap.io.compare_io`
(GIAP_SPEC.md §7.9.3, §8.2, §8.3; issue M9-22; layout ADR-0046). Statistician sign-off: block B1 fixed (decision 11), 2026-10-22.

## Inputs and checks
Two or more `OUT/` trees on the same genome (genome SHA-256 from `run_manifest.json`, else
`preflight/preflight.json`) with identical locus-tag sets. Labels are `--labels` or the directory
basenames (repeats suffixed `_2`, `_3`).

## Detection
A gene x run matrix of `detected` / `detection_status` (genes in §6.5 order of the first run) and
the exact intersections of the detected sets: a gene detected in exactly A1 and A3 counts in
`A1&A3` only, as an UpSet plot draws it; `none` holds the genes detected in no run.

## Abundance agreement
Per unordered pair, shared genes are those detected with a measured (`abundance_imputed = false`)
log2 abundance in both runs.

- **Spearman ρ** with average ranks for ties (equal to `scipy.stats.spearmanr`); NA for a
  constant side.
- **CI**: the 95 % percentile CI of `tu_block_bootstrap` with `B = stats.bootstrap` resamples of
  the TUs of the pair's first run (a gene without a TU is its own unit), drawn from
  `giap.rng.for_step("compare", seed)`. Genes of one TU share promoter and transcript, so their
  measurement errors in both runs are positively correlated; resampling genes would treat them
  as independent and give a CI too narrow by about the square root of the design effect
  (docs/methods/bootstrap.md). The calibration test's teeth arm (gene-level units) shows the
  under-coverage.
- **Bland–Altman** on log2: bias = mean(a − b), limits = bias ± 1.96 SD(a − b), ddof = 1.

Flags: `too_few_shared` (< 3 shared genes: every statistic NA), `ci_few_units` (fewer than 10 resampling units, reported in `n_units`: the CI is kept from 2 units but is not a 95 % interval; below 2 it is NA),
`ci_unestimated` (fewer than B/2 finite replicates), `scale_differs` (`abundance_type` or
`abundance_scale` differ: the bias is not interpretable but is written), `missing_in:<label>` (a
run without abundance: the row is NA).

## Stability
Jaccard |A ∩ B| / |A ∪ B| of five sets per run: genes of islands with `p_value ≤
stats.alpha.context`; (`gene_a`, `gene_b`) of `operon_flags` with `q_value ≤
stats.alpha.operon_flags`; (table, collection, term, direction) of `enrichment_detection` and
`enrichment_gsea` rows with `q_value ≤ stats.alpha.enrichment` (GSEA direction = sign of NES);
called novel ORFs; (contig, codon_start, strand, alt_aa) of `variants` rows with `call = ALT`.
Two empty sets give NA with `both_empty`; an absent table (or the `function` module off) gives NA
with `missing_in:<label>`. Jaccard is descriptive: no test, no CI.

## Calibration
`tests/calibration/test_compare_rho_ci.py` (suite `compare`, R = 100 / 1,000 per ADR-0042):
coverage of the ρ CI on simulated TU-correlated pairs (gated, anti-conservative side), recorded
coverage, width and bias, and the gene-level teeth arm.
