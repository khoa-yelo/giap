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
- **The re-test:** a significant GO set with at least one significant descendant is re-tested without the genes of all its significant descendants. TU membership is rebuilt, so a TU stays when another of its scoring genes is still in the reduced set. The re-test uses min_size 1.
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

The two position flags below mark the rows where the last two limits are most likely to matter. An unflagged row is not a guarantee (see Limits).

## Position flags (M8-05b; spec decision 9, ADR-0025, ADR-0026)
Both flags are computed for every main-run row, are pure and use no RNG. Redundancy re-tests compute no flags; a re-tested parent keeps the flags of its main run. They follow the existing flags in the order `ci_unestimated`, `unestimated`, `size_filtered`, `genomic_cluster`, `detection_confounded`.

- **`genomic_cluster` (δ ≥ τ = 0.15).** δ estimates how much the set's member TUs inflate the variance of its mean score, relative to the TU-label permutation null, if TU scores were AR(1) with φ₀ = 0.5 along the genome.
  - *Order:* ranked TUs are split by replicon (`contig`) and kept in natural TU-ID order (genome order, from `build_tus`); unranked TUs are skipped, so members with only unranked TUs between them are adjacent. On a circular replicon (C1 `topology == "circular"`) d = min(|i − j|, n_r − |i − j|); every other topology is linear. Members on different replicons are never paired.
  - *Statistic:* with w(d) = 0.5^d for 1 ≤ d ≤ 10 (else 0), k member TUs and n ranked TUs, C = (2/k) Σ_{member pairs} w(d), E[C] = 2(k − 1) P_w / (n(n − 1)) with P_w = Σ w(d) over all pairs of ranked TUs, and δ = (C − E[C]) / (1 − k/n); δ = 0 when k < 2 or k ≥ n. δ is centred at 0 for randomly placed sets of any size.
- **`detection_confounded` (|z| ≥ z₀ = 1 with the sign that produces the observed ES).** Counted in genes over the eligible universe (non-partial, finite `p_detect`): M eligible genes, N scoring genes, n eligible genes of the gene-level set (undetected members included), x of them scoring. With f = N/M and v = n f (1 − f)(M − n)/(M − 1), z = (x − n f)/√v (z = 0 when v = 0). The row is flagged when ES is finite and ≠ 0 (`unestimated` rows included, since their ES is a real statistic) and z · ES < 0 with |z| ≥ 1; when ES = 0, when |z| ≥ 1. A NaN ES, which `fgsea_simple` does not return, is treated as 0.

### Limits
- `genomic_cluster` is calibrated only for φ ≤ 0.5: the kernel fixes φ₀ = 0.5, and at φ = 0.7 unflagged rows are anti-conservative at α = 0.01 (0.0165 against 0.0143). It uses ranked-TU order, not bp distance. The claim is per-row level, not FDR among unflagged rows.
- `detection_confounded` is a screen, not a guarantee. In mixed-confounding worlds (s ~ U(0, 1.5), β ∈ {1, 2}) unflagged rows reject at about 1.3α / 1.15α / 1.07α at α = 0.01 / 0.05 / 0.10 (pooled R ≈ 12,600: 0.0128 / 0.0567 / 0.1073), above the §8.4 bound; unscreened rows in the same worlds reject at about 0.25 at α = 0.05. Under an unconfounded null, unflagged rows reject at about 1.1α (0.0551 / 0.1088 at 0.05 / 0.10), because flagged null rows are conservative. Rows with |z| < 1 keep a small bias that no detected-fraction threshold can screen out.
- q values of unflagged rows come from a BH family that includes flagged rows; with 10 flagged sets among 40 the FWER of unflagged null rows rose from 0.004–0.02 to 0.022–0.057 (still ≤ α under the 1/500 floor; M8-11 re-checks). It reacts to observed detection, so a real effect with an opposite-sign detection difference is flagged too; about 13–16 % of rows are flagged under an unconfounded null.
- The hardening ADRs (a rotation null; censored or `p_detect`-weighted ranks) remove both limits.

## Evidence
- **Single-set level at R = 4,000** (statistician, Slurm 1867122): rejection 0.0508 at α = 0.05 under the null. Attack 3, the ties, median imputation and scattered AR(1) are all at level.
- **Family probe** (GO-like nested sets, R = 2,000): FWER 0.0015 at q 0.05 under the complete null. With shifted children, FDR is 0.017 at q 0.05 and 0.053 at q 0.10 (shift +1.5), and 0.0035 / 0.024 (shift +0.8).
- **Committed arms:** tests/calibration/test_gsea_tu_calibration.py and tests/calibration/attacks/test_attack_gsea_structure.py (Slurm 1867070: 14 passed, 2 strict xfails).
- **Position flags** (statistician, #174): full tier R = 300 (Slurm 1867323: 33 passed, 2 strict xfails); R = 4,000 per arm (Slurm 1867332; cluster φ 0.5 mix at level, detection as in Limits); BH-family spillover (Slurm 1867330). Committed arms: tests/calibration/attacks/test_attack_gsea_flags.py (detection arms gated at ≤ 1.5α).
