# Gene-level abundance: iBAQ, Top3, protein-level route, missing values

Implements `giap.genetable.abundance` (`peptide_abundance`, `protein_abundance`,
`resolve_c10_scale`, `impute`, `gene_abundance`); GIAP_SPEC.md §2.1 decisions 2–3, §7.3.4, §7.2.6,
§9.3–9.4, §5.2 `abundance`; issue M4-14 (docs/specs/M4-14-abundance.md, including lead decisions
1–6 on PR #82). The output is one point estimate per gene in one sample. This module runs no test,
so it has no p-values, no test family (§8.1) and no CI. Uncertainty enters downstream through
rank-based statistics, TU-block resampling and the `basic`/`full` sensitivity reruns (§8.6).

## What is estimated

The aim is the relative molar amount of each gene's protein in a single sample, known only up to
one common constant per run. Only ratios between genes and ranks carry meaning. Values are not
absolute copy numbers. Values from different runs or inputs cannot be compared without
normalisation, which GIAP does not do.

## Which route is used (and when abundance is NA)

| Condition (checked in this order) | `source` | `abundance_type` | `abundance_scale` |
|---|---|---|---|
| some `abundance_eligible` accepted peptide has a non-NA intensity | `peptide_intensity` | `abundance.method` (`ibaq` default, or `top3`) | `molar` |
| otherwise, C10 (`abundance.tsv`) given | `protein_table` | `abundance.type` (declared; default `unknown`) | `molar` iff type ∈ {ibaq, top3}, else `relative` |
| neither | `none` | `none` | NA; every abundance column NA, `imputation_policy` NA |

The two routes are never mixed. When peptide intensities exist, C10 is ignored. One edge case
follows from this: if every eligible intensity is exactly 0, the run takes the peptide route, every
gene comes out missing, and C10 is still not used. A `no_observed_abundance` warning is logged.

Under the `relative` scale (`maxlfq`, `intensity`, `unknown`), Module 4 marks its abundance tests
`unestimated (relative_scale)` and skips stoichiometry (§7.3.4). The default C10 type is `unknown`.
A protein table therefore gives relative abundance unless the user declares `--abundance-type ibaq`
or `top3`.

## Peptide route

Inputs:
- The M4-11 accepted peptides with `abundance_eligible`. These are unique exact hits plus tolerant
  hits under `related`. `multi_identical`, `shared`, N- and C-terminal extension evidence are
  excluded.
- A finite intensity > 0. An intensity of 0 or NA means the peptide is skipped, not counted as 0.
- Each peptide's intensity is already summed over charges and modified forms within a run, then the
  median across runs is taken (§5.3, C9).

Estimates:
- **iBAQ** (Schwanhäusser et al. 2011): iBAQ_g = Σ_quantified intensity / max(1, n_obs,g), where
  n_obs,g = `n_observable_peptides`. This count comes from M4-10: configured enzyme, 0 missed
  cleavages, length 7–30, counted on the strain gene sequence (on the partial sequence for partial
  genes).
  - The sum of observed peptide intensities scales with the gene's molar amount times the number of
    its peptides the instrument can see. Dividing by the number of theoretically observable peptides
    removes the length (peptide-count) factor.
- **Top3** (Silva et al. 2006): the arithmetic mean, on the linear scale, of the (up to) 3 largest
  intensities, then log2 (lead decision on spec question 1).
  - With 1–2 quantified peptides it is the mean of those, flagged `top3_lt3`.
  - `top3_log2` is always reported on the peptide route. `abundance_log2` follows
    `abundance.method`.
- `abundance_raw` is the linear iBAQ (or Top3) value, and `abundance_log2 = log2(abundance_raw)`.
  - Sums are taken over sorted values with `math.fsum`, so results do not depend on row order.
  - Multiplying every intensity by c > 0 shifts every observed `abundance_log2` by exactly log2 c
    (a unit property test).

### Assumptions and when they fail (direction)

| Assumption | Fails when | Effect |
|---|---|---|
| Peptide intensity ∝ molar amount × a peptide response factor independent of the gene | response factors vary a lot (ionisation, sequence) | noise only. Ranks degrade smoothly (Spearman with truth 0.98 → 0.91 as the log-SD of peptide noise goes 0.5 → 1.5; table below) |
| Missing peptides are missing at random (MAR) | low-abundance genes lose more peptides (MNAR, the usual case) | iBAQ of low-abundance genes is pulled further down, because unseen observable peptides count as 0 in the numerator. The scale is stretched at the low end. Ranks are mostly kept because the distortion is monotone in abundance |
| Numerator and denominator count the same peptides | accepted peptides include missed-cleavage or out-of-range-length peptides that n_obs does not count. Overlapping missed-cleavage forms count the same sequence twice | iBAQ is somewhat inflated for genes rich in missed cleavages. This matches MaxQuant's iBAQ convention and is not corrected |
| n_obs reflects observability | genes with n_obs = 0 that still have quantified peptides (e.g. only missed-cleavage peptides) | the denominator is floored at 1. For these genes iBAQ is the raw sum, which is too high relative to other genes |
| Top3: the 3 best-responding peptides are similar across proteins | genes with many quantified peptides | the mean of the top 3 of n order statistics grows with n. Top3 is biased upward for long, well-covered genes. iBAQ has no such bias, which is why it is the primary method |
| Tolerant (`related`) hits measure the strain protein | strain differences change ionisation | small per-peptide noise. These hits are labelled `strain_difference` upstream |

## Protein route (`abundance.tsv`, C10)

C10 values are joined to genes through the M4-12 protein → gene map (`protein_id, gene, method`).

**Measurement scale** (lead decision 1). The scale is the explicit `abundance.scale` if it is not
`auto`, otherwise the preflight scale (`M4Inputs.abundance_scale`), otherwise a `ValueError`. Two
details:
- C10 normally has no `scale` column. If a `scale` column carries more than one value under `auto`,
  the run errors.
- Preflight detects the scale as `log2` if any value is negative or every value is ≤ 50, else
  `linear` (§5.3).

**Transform.** Linear values become log2. Values ≤ 0 and ±inf become missing, because log2 is
undefined for them. On a log2 scale, finite values are kept, including negatives. `abundance_raw`
keeps the input value on the input's own scale, so for a log2 input it is the log2 number.

**Several proteins on one gene** (§7.2.6). The gene takes the largest log2 value (NA values lose;
ties go to the smallest `protein_id`), with flag `two_proteins_one_gene`.
- Taking the maximum rather than summing is the spec's choice. Two input proteins on one locus are
  usually redundant entries (isoform or accession duplicates) quantified from overlapping peptides,
  so summing would double count.
- The maximum is biased upward when the entries are independent noisy measurements of one protein.
  The flag lets users filter these genes.

**Cross-strain proteins** (lead decision 2). A protein mapped by `protein_hit` whose C12 hit row
carries `cross_strain_protein` passes its value to the strain gene, and the gene gets the flag.

## Missing values and imputation (§2.1 decision 2; §7.0 "two kinds of missing")

Two kinds of missingness are kept apart:
1. **Never identified** (not detected): the gene is never imputed. It is a detection outcome,
   handled by the detection and detectability models (left-censored ranks in §7.3.5 islands, p_detect
   weighting in ORA). Imputing it would invent a quantity for a protein with no evidence and would
   double count the detection-limit information those models already use.
2. **Identified but not quantified** (detected, no value): for example, every peptide of the gene
   is `multi_identical`, an extension, or has NA intensity; or the gene has no C10 row.
   `abundance.missing` applies to these genes.

**Fill statistic.** It is computed on the observed `abundance_log2` of detected, eligible
(non-partial) genes, pooled over all replicons:

| Policy | Fill value | Role |
|---|---|---|
| `median` (default) | median of observed log2 | primary: the genes are detected, so a typical value is used |
| `mean` | mean of observed log2 | accepted |
| `min` | minimum observed log2 | sensitivity (`basic`): treats unquantified as near the limit of quantification |
| `zero` | minimum − 1 (half the minimum on the linear scale) | accepted; "zero" cannot be represented on a log scale |
| `exclude` | none: stays NA, stays detected, leaves abundance statistics | sensitivity (`basic`) |

Rules:
- Filled genes get `abundance_imputed = true` and flag `imputed_abundance`. `abundance_raw` stays NA.
- `imputation_policy` records the policy on every row whenever the source is not `none`.
- If no detected eligible gene has an observed value, nothing is filled, `fill_value` is None, and
  `no_observed_abundance` is logged at WARNING.
- Detected **partial** genes with a missing value are filled (lead decision 5). They do not enter
  the fill statistic, because their n_obs and their quantified peptides cover only part of the
  protein.
- Undetected genes that have an observed value keep it (spec question 3; e.g. single-peptide genes
  under `two_peptide`). This is reporting only. Abundance statistics use detected genes.

**Direction of the default.** Detected-but-unquantified genes are most likely low in abundance (MNAR).
Filling them at the median pulls them toward the centre. This attenuates region and island shifts
driven by such genes and adds a block of tied values at the median. The `basic` sensitivity reruns
under `min` (the opposite direction) and `exclude` (no fill) bracket this effect, and
`tables/sensitivity.tsv` reports the Jaccard overlap of flagged sets.

Two places react to imputed values:
- Operon discordance excludes imputed pairs from its t fit (§7.4.2).
- Preflight warns when more than 20 % of values are missing (`abundance_missing_warn_frac`).

The fill statistic is genome-wide, not per replicon. On a high-copy plasmid the imputed value is the
chromosome-dominated median. Region tests treat plasmids as replicon-level effects (§8.5 attack 4),
so this affects only within-plasmid comparisons.

## Simulation evidence

Two claims are registered in `tests/calibration/test_abundance_calibration.py` (suite `regions`,
criterion `fwer` at 0.05 on per-replicate 0/1 indicators, ADR-0015 two-block confirmation):
- `abundance.spearman.{ibaq,top3}`: the indicator 1{Spearman(estimate, truth) < 0.7} over detected
  genes with ≥ 3 quantified peptides. This is a **sanity** claim, not an error-rate claim.
- `abundance.ibaq_denominator`: the indicator 1{the 5 % paired gene-bootstrap lower bound of
  ρ_iBAQ − ρ_rawsum ≤ 0}.

The generator is the test oracle (`tests/helpers/abundance_oracle.py`):
- giapsim layout with 400 chromosome genes and M4-08 expression;
- observable peptides quantified independently with probability p_obs = 0.6;
- peptide noise: ln I = a·ln 2 + N(0, 0.5²).

Statistician check (Slurm job 1848760; exact nightly streams, replicates 0–1999 = both blocks of
the full run; head 0fcc695):

| Quantity | min | 1 % | median | max |
|---|---|---|---|---|
| ρ iBAQ | 0.969 | 0.974 | 0.984 | 0.993 |
| ρ Top3 | 0.968 | 0.976 | 0.984 | 0.992 |
| ρ raw sum | 0.911 | 0.924 | 0.951 | 0.976 |
| bootstrap lower bound of Δρ | 0.0087 | 0.0136 | 0.024 | 0.049 |
| genes per replicate (≥ 3 peptides) | 186 | 194 | 216 | 248 |

Results:
- 0 failing replicates of 2,000 for each of the three claims. The fast (R = 100) and full
  (R = 1,000) checks pass in both blocks.
- The per-replicate failure rate is < 0.0015 (95 % upper bound). The checks' added confirmed false
  alarm on the nightly run is therefore ≈ 0, well inside the ≈ 0.04 family budget.
- Teeth: the permuted-n_obs denominator check failed in 40/40 replicates (largest lower bound
  −0.042), and the shuffled-estimate Spearman check fails. Each teeth check needs ≥ 3/20 failures
  in both blocks.
- Cost: about 0.6 s per replicate per claim, so about 30 CPU-min per nightly run for the three
  claims and about 3 CPU-min in the fast suite.

### Stress (outside the gate; R = 300 per cell, seed [777, i])

| p_obs, peptide log-SD | ρ iBAQ median (min) | ρ Top3 median (min) | ρ raw sum | denominator indicator rate |
|---|---|---|---|---|
| 0.6, 0.5 (oracle default) | 0.983 (0.971) | 0.984 (0.975) | 0.949 | 0.000 |
| 0.6, 1.0 | 0.959 (0.931) | 0.941 (0.904) | 0.923 | 0.000 |
| 0.6, 1.5 | 0.909 (0.851) | 0.877 (0.801) | 0.873 | 0.003 |
| 0.3, 0.5 | 0.969 (0.946) | 0.981 (0.965) | 0.954 | **0.103** |
| 0.3, 1.0 | 0.933 (0.883) | 0.929 (0.874) | 0.916 | **0.270** |
| 0.3, 1.5 | 0.871 (0.788) | 0.853 (0.761) | 0.850 | **0.403** |

How to read the stress table:
- The Spearman sanity claim never failed: the 0.7 threshold is far below every cell.
- The claim also cannot tell iBAQ from the raw sum, the max, or Top3 (reviewer finding 4). It
  guards only against estimators that ignore the data.
- iBAQ ranks better than Top3 once peptide noise is large, as expected: Top3 discards information
  and has the order-statistic bias.
- The denominator claim measures a property of the generator, not of the code. At low peptide
  coverage (p_obs = 0.3), binomial sampling of the quantified peptides swamps the n_obs
  correction. The per-replicate indicator rate is then 0.10–0.40, which would be a confirmed
  nightly fail at R = 1,000 (bound 0.064).

**Before swapping the generator to M4-09** (backlog; logistic peptide detectability, missed
cleavages), re-measure this rate. If the M4-09 coverage is low, restate the claim as a pooled
effect (for example, the mean Δρ over replicates > 0 with a replicate-level CI) rather than
per-replicate power.

The bootstrap resamples genes rather than TUs. Here the estimation errors are independent per gene,
but the true abundances are TU-correlated (M4-08 AR(1)). The gene bootstrap therefore somewhat
understates the variance of Δρ, which makes the claim easier to pass, not harder. This is
acceptable for a teeth-style check. It is not an inference reported to users, so §8.2 "resample
TUs" does not bind here.

## Sign-off

Statistician (M4-14, PR #82, 2026-10-15): **approved**. The following match §7.3.4/§7.2.6 and the
lead decisions:
- the iBAQ and Top3 definitions;
- the eligibility and zero/NA rules;
- the protein-route scale resolution and the maximum rule;
- the fill rule (statistic over non-partial detected genes, filling only detected genes, never
  undetected ones).

The two calibration claims are sound as sanity and teeth checks, with zero false-alarm cost on the
current generator. Forwarded:
1. The denominator claim's fragility at low peptide coverage, to re-check with the M4-09 generator
   before swapping. Suggested attack: `p_obs = 0.3`.
2. The MNAR direction of median imputation, for the Module 4 methods notes.
3. Recording the abundance claims' fast-suite time in `benchmarks/baseline.json`.

## References

- Schwanhäusser B. et al. (2011) Global quantification of mammalian gene expression control.
  Nature 473:337–342 (iBAQ).
- Silva J.C. et al. (2006) Absolute quantification of proteins by LCMS^E. Mol Cell Proteomics
  5:144–156 (Top3).
