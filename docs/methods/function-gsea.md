# Abundance GSEA (M8-05)

`giap.function.gsea.gsea_table` tests whether a gene set's genes rank high or low by protein abundance (§7.7.4). It runs the fgsea port of M8-02 (`docs/methods/fgsea-port.md`) on TUs, not genes.

## Ranking unit: TUs
- **Scoring genes:** detected, not partial, with a finite `p_detect` and an abundance (`abundance_log2` not NA). Imputed abundances count (ADR-0028); the §8.6 sensitivity reruns vary the imputation.
- **TU score:** the mean `abundance_log2` of the TU's scoring genes. Only TUs with at least one scoring gene are ranked.
- **Centring:** scores are centred on their median over ranked TUs (ADR-0028). Raw log2 abundances are all positive, so `scoreType = std` would make the test effectively one-sided. Centring changes the weights only, so validity under exchangeability is unchanged.
- **Ties:** ties enter fgsea in one seeded random order per run, drawn from the `function.gsea` stream. Re-tests share that order.
- **Membership:** a TU is in a set when any of its scoring genes is. A detected gene without an abundance gives its TU no membership; §7.7.4's "detected genes" is read as scoring genes.
- **Size:** set sizes are counted in TUs, and the `function.gsea_size` filter applies to them.

## Family and multiple testing
- **BH family:** BH runs per collection over the size-filtered sets, with both ES signs in one family (§8.1). Unestimated sets (NaN p) are left out of m.
- **Engine:** pass 1 uses `fgsea_simple`, whose p-value floor is 1/(1 + modeFraction), about 1/500 at the default `n_perm = 1000` (M8-05 decision 4; statistician review). This departs temporarily from §7.7.4 and §8.1, which ask for multilevel fgsea with no floor. M8-11 removes the floor, and `enrichment_gsea.tsv` is not release-ready until it does.
- **`log2err`:** a Monte Carlo precision label. It slightly understates the error of the conditional p.

## GO redundancy
- **The re-test:** a GO set that is significant together with one of its significant descendants is re-tested without the TUs of all its significant descendants. The re-test uses min_size 1.
- **Deciding `redundant`:** the re-test p replaces the set's main p in a copy of the collection's p vector, and `redundant` is true when that substituted BH q is above α.
- **Reported q:** `q_value` stays the main-run BH q.
- **Degenerate re-tests:** a 0-TU re-test gives `redundant = true` plus `size_filtered`. A NaN re-test p gives `redundant = NA`.
- **`parent_group`:** the most general significant ancestor, meaning the one with the shortest path (in `parents` steps) to a term without parents, ties broken by ID. A set with no significant ancestor is its own group.

## What a row claims (restricted claim)
- p-values are valid under TU exchangeability.
- There is no level claim for genomically clustered sets under positional autocorrelation (ADR-0025).
- Results are conditional on detection, with no `p_detect` adjustment (ADR-0026).
- NES has no CI; every row carries `ci_unestimated` (ADR-0027, a deviation from §8.2).
- "Significant and non-redundant" is not an FDR-controlled list.
- "No significant set" never means "no enrichment", because the floor is about 1/500 and BH needs several sets at the floor before it rejects any.

Pass-1 rows do not yet carry the `genomic_cluster` and `detection_confounded` flags (M8-05b), so the restricted claim applies to every row. M8-06 (the writer) does not merge before M8-05b.

## Evidence
- **Single-set level at R = 4,000** (statistician, Slurm 1867122): rejection 0.0508 at α = 0.05 under the null. Attack 3, the ties, median imputation and scattered AR(1) are all at level.
- **Family probe** (GO-like nested sets, R = 2,000): FWER 0.0015 at q 0.05 under the complete null, and FDR 0.017–0.053 with shifted children.
- **Committed arms:** tests/calibration/test_gsea_tu_calibration.py and tests/calibration/attacks/test_attack_gsea_structure.py (Slurm 1867070: 14 passed, 2 strict xfails).
