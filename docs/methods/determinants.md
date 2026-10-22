# Determinants: what explains abundance

Module 8 track R (GIAP_SPEC.md §7.8). Issues M8-07 onward; this note grows with each.

## Features

Issue M8-07, module `giap.determinants.features` (pure). `feature_frame` builds
`tables/determinants_features.tsv`; `model_matrix` encodes it for M8-08.

**Rows.** One per eligible gene (§7.3: not partial, finite `p_detect`), in §6.5 order (C1 contig
order, start, end, strand, locus tag). Ribosomal-protein genes, the CAI reference set, stay in the
table with `ribosomal_protein = true`; M8-08 excludes them from fitting and from univariate CAI
statistics. `ribosomal_protein` is the C2 column when the gene table carries it, otherwise the
§7.1.2 product / gene-name rule or the `cai_reference` gene flag.

**Columns.** Identity (`locus_tag`, `tu_id`, `detected`, `abundance_log2`, `abundance_imputed`,
`ribosomal_protein`), then the 31 features of §7.8.1 in table order. Name map (feature ← source):

| Feature | Source |
|---|---|
| `promotech_max`, `promotech_tu` | gene table `promoter_score`, `promoter_score_tu` |
| `n_promoters_tu_to_gene` | C6 `n_promoter_windows_tu` (windows ≥ `promotech_yes`, lead decision 1) |
| `mw` | gene table `mw_da` |
| `replichore_position` | M5-09 `context/gene_context.parquet` (NA without it or without an oriC) |
| `in_prophage`, `in_bgc` | C8 `region_genes` with C7 `region_type` `prophage` / `bgc`; NA without C7/C8 |
| `dist_to_tu_start_bp` | computed: bp from the TU's first start codon to this gene's, in the direction of transcription, modulo the length on circular replicons; NA for a gene in no TU |
| `nterm_residue` | computed from C2 `aa` (the Met rule below) |
| every other feature | the gene-table column of the same name |

**The Met rule (decision 19).** `nterm_residue` is residue 2 when residue 1 is M and residue 2 ∈
{A, C, G, P, S, T, V} (methionine aminopeptidase removes the initiator Met before small residues),
else residue 1; NA for proteins shorter than 2 aa or with an empty `aa` (pseudogenes). It is a
covariate only: nothing in mapping or annotation uses it (§16 excludes a recoding catalogue).

**Transforms.** None in the TSV (raw values). In `model_matrix`: `tir_ostir` and `length_aa` →
natural log (lead decision 1; values ≤ 0 → NaN); `start_codon` one-hot over {ATG, GTG, TTG, other},
`replicon_type` over {chromosome, plasmid, unknown}, `nterm_residue` over the 20 residues (every
level kept), named `<feature>_<level>`, an NA value giving NaN in every level; booleans 0/1. A
feature that is NA on every row (a missing layer) stays in the TSV and is dropped from the matrix
with all its encoded columns; the returned column → block map omits it, and a block left without
columns is absent. Run flags (`layer_missing:<layer>`) are the step's (M8-10).

## Models

Owed by M8-08 (`giap.determinants.models`); not written yet.

## Null and attribution

Issue M8-09, modules `giap.determinants.null` and `giap.determinants.attribution` (pure). Outputs:
one `NullResult` per model and `tables/determinants_attribution.tsv` (§9.3).

### The null
- **Statistic:** the observed mean CV R² over the 5 repeats (`cv_r2_mean`); one p per model, no BH
  (§8.1 "Determinant model"). `p = (1 + b)/(1 + B)`, b = #{null R² ≥ observed} (`perm_p`,
  `alternative="greater"`); B = `determinants.permutations` (10,000 with `--publication`).
- **Permutation unit:** TUs are grouped by their number of model rows (detected, non-imputed,
  non-ribosomal genes). Within a group of ≥ 2 TUs, whole TU response vectors move among the TUs,
  gene order inside a TU kept. A TU whose size is unique stays put (conservative).
- **Degenerate:** if every TU size is unique, every permutation is the identity: no refit, b = B,
  p = 1, flag `null_degenerate`.
- **Refit:** each permutation refits one CV repeat on the fixed repeat-1 folds (decision Q1). The
  folds depend only on the TUs, not on y, so conditioning on them is sound and the permutation test
  of the repeat-1 R² is exact.
- **Seeds:** permutation b uses child b of `SeedSequence([seed, crc32("determinants.null")]).spawn(B)`,
  so permutation b is the same for any B, and results are identical for any thread count.

### Conservativeness of the 5-repeat statistic
The 5-repeat observed mean against a 1-repeat null is not exact, but it is conservative in practice.
Statistician runs (Slurm 1867304, 1867334, 1867348, 1867349; B = 99, 50 trees, 400 genes,
R = 200 replicates), rejection at α 0.01 / 0.05 / 0.10 (binomial bound α + 2 SE: 0.024 / 0.081 /
0.142):

| World | LightGBM | Elastic net |
|---|---|---|
| clustered | 0.010 / 0.040 / 0.105 | 0.000 / 0.005 / 0.020 (mean p 0.55) |
| size | 0.000 / 0.020 / 0.075 | 0.005 / 0.025 / 0.045 |
| clustered, exact (observed = repeat-1 R²) | 0.010 / 0.060 / 0.115 | 0.015 / 0.060 / 0.140 |

- LightGBM is at nominal level; the elastic net is clearly conservative (0.005 against 0.060 for
  the exact version at α 0.05).
- Swapping the observed statistic for the repeat-1 R² would restore power; this is an option, not
  a requirement.
- The §8.4 gate tier (R = 20) is not a level claim: 0 of 20 rejections happens 36 % of the time
  for a level-0.05 test.

### Pass-1 limit (ADR-0034, proposed)
The p is valid only when TU response blocks are exchangeable within size groups. It fails the
statistician's attacks (R = 100, same jobs; bound at α 0.05: 0.094):

| Attack | LightGBM at α 0.05 | Elastic net at α 0.05 |
|---|---|---|
| inter-TU AR(1) abundance (ρ 0.9) with an independent AR(1) `gc_local_dev` | 0.12 | 0.12 |
| + `replichore_position` | 0.95 | 0.41 |
| detection selection (y ⊥ X; detection depends on y and length) | 0.53 | 0.83 |

- **Why:** permuting TUs destroys the response's autocorrelation, so the null is too narrow;
  scattered TU folds let a positional feature interpolate a gene's response from its neighbours;
  detection is a collider of abundance and the detectability covariates in the protein block.
- **Consequence:** every `NullResult` carries `null_uncalibrated` (M8-10 copies it onto the
  `null_p` row). The p is not a discovery claim.
- **Committed arms:** tests/calibration/attacks/test_attack_determinants_null.py; they gate the
  level unless every result of an arm carries `null_uncalibrated` (the ADR-0032 pattern).

### Hardening owed (pass 2)
- Spatially contiguous CV folds (chunks of consecutive TUs) for both the observed CV and the null.
  In the statistician's probe this brings the AR(1) attack to level (EN 0.06, LightGBM 0.02 at α
  0.05) with unchanged planted power; circular shifts within size groups did not.
- A decision on positional features (`replichore_position`): exclude them from the null claim, or
  use a spatial null. Contiguous folds alone leave LightGBM at 0.70.
- A detection-selection model. Stratifying the permutation by TU-mean detectability still fails
  (EN 0.46 at α 0.05).
- Fresh folds per permutation stay a sensitivity option (decision Q1).

### SHAP attribution
- **Scope:** LightGBM only (`pred_contrib=True`, no `shap` package); the elastic net has no SHAP.
  The primary attribution is the elastic-net ΔR² (ADR-0033); SHAP describes the LightGBM fit only.
- **Per gene:** in each repeat, the fold model that held the gene out gives its contribution
  vector; the vectors are averaged over the 5 repeats and summed by block into `shap_transcription`,
  `shap_translation`, `shap_codon`, `shap_protein`, `shap_context` (NA for a block with no column).
- **Additivity (§7.8.3):** Σ block SHAP + mean expected value = the OOF LightGBM prediction, per
  gene. A guard raises if any row misses by more than 10⁻⁶ (the rows of X are not in fitted order).
- **Dominant block:** for top-decile genes (model rows with `abundance_log2` ≥ the linear 90th
  percentile), the block with the largest positive SHAP sum; NA when none is positive; ties by
  block order.

### Ribosomal rows
- Ribosomal-protein genes are the CAI reference, so their CAI is self-referential and lies outside
  the training range. They are excluded from fitting.
- Detected, non-imputed ribosomal genes get rows predicted by the mean of all repeat × fold models
  (none saw them), with SHAP likewise averaged; additivity holds for the average.
- These predictions are extrapolations (trees cannot extrapolate), labelled `cai_reference`. The
  rows are kept out of the residual SD and the top decile and are never residual-flagged.
- Undetected ribosomal genes get no rows in pass 1.

### Residual flag
- `residual` = observed − `predicted_log2` (the mean OOF LightGBM prediction).
- |residual| > 2 SD (ddof 1, over model rows; two-sided, decision Q2) adds
  `abundance_not_explained_by_sequence`; SD exactly 0 sets no flag.
- A descriptive screen, not a claim: about 4.6 % of genes are flagged under any model (statistician
  review of #173).
- Fewer than 100 model rows: no SHAP and no null; rows have NA predictions and `too_few_rows`.
