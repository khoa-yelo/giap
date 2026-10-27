# Abundance GSEA (M8-05, engine M8-17b)

`gippa.function.gsea.gsea_table` tests whether a gene set's genes rank high or low by protein abundance (§7.7.4). It ranks TUs, not genes, and hands the ranking to **R fgsea 1.36.2** (`fgsea::fgseaMultilevel`), run as a separate process (ADR-0038, M8-17). The Python/numba port of M8-02 (`docs/methods/fgsea-port.md`) is superseded; it survives only as `SIMPLE_ENGINE`, the in-process stand-in that unit tests inject.

## Engine (ADR-0038)
- **Environment.** `workflow/envs/fgsea.yml` (`gippa-fgsea`): `r-base` 4.5.3, `bioconductor-fgsea` 1.36.2 and `r-data.table` 1.18.6.1 from conda-forge and bioconda, installed from the explicit lock `workflow/envs/locks/fgsea.explicit.txt`. No Bioconductor data package is needed, so nothing is fetched from bioconductor.org. fgsea is conda-only under every profile for now (M8-17 lead decision 9).
- **Call.** `workflow/bin/gippa_fgsea.R` runs under `Rscript --vanilla` and calls `fgseaMultilevel(sampleSize = 101, eps = 0, nPermSimple = 1000, minSize, maxSize, scoreType = "std", nproc = 1)` with the `function.gsea` defaults. `eps = 0` means the multilevel stage estimates p with no lower bound, so there is no 1/(1 + nPermSimple) floor. The engine's `method` is `fgsea_multilevel`.
- **Seed.** `gsea_table` keeps one generator, `for_step("function.gsea", seed)`. It draws the tie order (`perm = rng.permutation(n)`) and then one integer `r_seed = rng.integers(1, 2**31 - 1)`. R runs `set.seed(r_seed)` before every job, and `nproc = 1` (SerialParam) makes a rerun byte-identical (acceptance A2, Slurm 1868835).
- **At most two R calls per step.** The main call holds one job per collection with an in-range set (`main:<collection>`). The second call holds every GO redundancy re-test of every GO collection (`retest:<collection>:<term>`) with the same ranking and the same `r_seed`; it is made only when at least one re-test exists. With no job at all, R is not called and the table is header-only.
- **Shared seed.** The re-test's Monte Carlo draws are common random numbers with the main call's. Each p is still a valid marginal estimate, and p_re only enters its own collection's BH vector, whose main p values already share one permutation set inside fgsea; so the shared seed adds only dependence of the kind BH already faces (statistician, PR #214). BH never crosses collections.
- **Failure.** Any R error fails the `function` step with exit 4 and writes no partial GSEA table (M8-17 lead decision 5). R's warnings (for example fgsea's ties warning) go into the `function.gsea_engine` event, never into errors.
- **Licence.** fgsea is MIT (NOTICE.md); R is GPL-2. Both run only as a separate program started by the step; nothing under `src/` imports or links them (R5; `scripts/ci/check_imports.py`, `check_no_rscript.py`).

## Ranking unit: TUs
- **Scoring genes:** detected, not partial, with a finite `p_detect` and an abundance (`abundance_log2` not NA). Imputed abundances count (ADR-0028); the §8.6 sensitivity reruns vary the imputation.
- **TU score:** the mean `abundance_log2` of the TU's scoring genes. Only TUs with at least one scoring gene are ranked.
- **Centring:** scores are centred on their median over ranked TUs (ADR-0028). Raw log2 abundances are all positive, so `scoreType = std` would make the test effectively one-sided. Centring changes the weights only, so validity under exchangeability is unchanged.
- **Ties:** ties enter fgsea in one seeded random order per run, drawn from the `function.gsea` stream: `stats.tsv` lists the TUs in that order, and fgsea's stable sort keeps it. Re-tests share that order.
- **Membership:** a TU is in a set when any of its scoring genes is. A detected gene without an abundance gives its TU no membership; §7.7.4's "detected genes" is read as scoring genes.
- **Size:** set sizes are counted in TUs, and the `function.gsea_size` filter applies to them. R's reported `size` must equal Python's count, or the call fails.
- **A set covering every ranked TU** has no complement: its running sum has no miss step, and every permutation gives the same set, so ES and its null are degenerate. It is never sent to R. When its size is within `function.gsea_size` and its collection has a job, it gets an NA row flagged `unestimated`; a collection with no job adds no row, and no job anywhere gives a header-only table (M8-17 lead decision 6a). This mirrors the ORA whole-universe rule.

## Family and multiple testing
- **BH family:** BH runs per collection over the size-filtered sets, with both ES signs in one family (§8.1). Unestimated sets (NaN p, including the whole-ranking row) are left out of m. R must return a row for every set it was sent (`read_results`), so it cannot silently shrink m.
- **p-values:** fgsea multilevel p with `eps = 0`, so small p go below 1/1001 (minimum 1.2e-4 in the R = 300 null run) and BH can reject single sets in large collections. This is the engine §7.7.4 and §8.1 ask for.
- **`unestimated`:** fgsea's unbalanced-signs case returns NA p and NES; the row keeps its ES, and its q is NA.
- **`log2err`:** fgsea's Monte Carlo precision label for p. When the multilevel stage cannot estimate it (`isCpGeHalf` FALSE), fgsea sets `log2err` to NA and warns that the p is likely overestimated. GIPPA keeps that p in BH, which is level-safe because the p is conservative, and flags the row `log2err_unestimated`. The flag excludes `unestimated` (an NA p implies an NA `log2err`) and sits where `unestimated` would in the flag order.

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
- "No significant set" never means "no enrichment": power is limited by set sizes in TUs, by the number of ranked TUs and by detection, and the result is conditional on detection (above).

The two position flags below mark the rows where the last two limits are most likely to matter. An unflagged row is not a guarantee (see Limits).

## Position flags (M8-05b; spec decision 9, ADR-0025, ADR-0026)
Both flags are computed for every main-run row, are pure and use no RNG. Redundancy re-tests compute no flags; a re-tested parent keeps the flags of its main run. They follow the existing flags in the order `ci_unestimated`, `unestimated`, `size_filtered`, `genomic_cluster`, `detection_confounded`.

- **`genomic_cluster` (δ ≥ τ = 0.15).** δ estimates how much the set's member TUs inflate the variance of its mean score, relative to the TU-label permutation null, if TU scores were AR(1) with φ₀ = 0.5 along the genome.
  - *Order:* ranked TUs are split by replicon (`contig`) and kept in natural TU-ID order (genome order, from `build_tus`); unranked TUs are skipped, so members with only unranked TUs between them are adjacent. On a circular replicon (C1 `topology == "circular"`) d = min(|i − j|, n_r − |i − j|); every other topology is linear. Members on different replicons are never paired.
  - *Statistic:* with w(d) = 0.5^d for 1 ≤ d ≤ 10 (else 0), k member TUs and n ranked TUs, C = (2/k) Σ_{member pairs} w(d), E[C] = 2(k − 1) P_w / (n(n − 1)) with P_w = Σ w(d) over all pairs of ranked TUs, and δ = (C − E[C]) / (1 − k/n); δ = 0 when k < 2 or k ≥ n. δ is centred at 0 for randomly placed sets of any size.
- **`detection_confounded` (|z| ≥ z₀ = 1 with the sign that produces the observed ES).** Counted in genes over the eligible universe (non-partial, finite `p_detect`): M eligible genes, N scoring genes, n eligible genes of the gene-level set (undetected members included), x of them scoring. With f = N/M and v = n f (1 − f)(M − n)/(M − 1), z = (x − n f)/√v (z = 0 when v = 0). The row is flagged when ES is finite and ≠ 0 (`unestimated` rows included, since their ES is a real statistic) and z · ES < 0 with |z| ≥ 1; when ES = 0, when |z| ≥ 1. A NaN ES (the whole-ranking row) is treated as 0.

### Limits
- `genomic_cluster` is calibrated only for φ ≤ 0.5: the kernel fixes φ₀ = 0.5, and at φ = 0.7 unflagged rows are anti-conservative at α = 0.01 (0.0165 against 0.0143). It uses ranked-TU order, not bp distance. The claim is per-row level, not FDR among unflagged rows.
- `detection_confounded` is a screen, not a guarantee. In mixed-confounding worlds (s ~ U(0, 1.5), β ∈ {1, 2}) unflagged rows reject at about 1.3α / 1.15α / 1.07α at α = 0.01 / 0.05 / 0.10 (pooled R ≈ 12,600: 0.0128 / 0.0567 / 0.1073), above the §8.4 bound; unscreened rows in the same worlds reject at about 0.25 at α = 0.05. Under an unconfounded null, unflagged rows reject at about 1.1α (0.0551 / 0.1088 at 0.05 / 0.10), because flagged null rows are conservative. Rows with |z| < 1 keep a small bias that no detected-fraction threshold can screen out. It reacts to observed detection, so a real effect with an opposite-sign detection difference is flagged too; about 13–16 % of rows are flagged under an unconfounded null.
- q values of unflagged rows come from a BH family that includes flagged rows. A flagged row whose p is anti-conservative acts like a signal in BH and lowers the cut for every row of the collection, so unflagged null rows can be rejected more often. The only measurement (10 flagged sets among 40: FWER of unflagged null rows 0.004–0.02 without them, 0.022–0.057 with them; Slurm 1867330) was made on the port, whose 1/500 floor kept BH from rejecting, so it does not carry over to the R engine. **Open (statistician):** the spillover limit is not re-justified for the R engine and has to be re-measured on R; until then no FWER or FDR claim is made for unflagged rows in a collection with flagged rows.
- The detection-screen figures above and the `genomic_cluster` φ = 0.7 figure were measured on the port engine. The flag screens themselves use no RNG and do not depend on the engine, but their at-level numbers must be re-measured on R before the hardening pass removes the port (statistician, PR #214 N4).
- The hardening ADRs (a rotation null; censored or `p_detect`-weighted ranks) remove both limits.

## Evidence
### R engine (shipped)
Rates are rejection rates at α = 0.01 / 0.05 / 0.10 with the §8.4 bounds from `tolerances.yaml`. The R-engine arms run on Slurm only (`slow`, `requires_r`), in the weekly `gsea` task of the full-calibration array (`calsuite gsea`); that task needs the `gippa-fgsea` env built. In every arm the fast tier (R = 50) reuses replicates 0..49 of the full arm's stream, so the fast arm is a prefix of the full arm and not independent evidence.
- **Single-set level** (Slurm 1868835, `tests/calibration/test_gsea_r_engine_calibration.py`, R = 300, all pass):
  - null: 0.017 / 0.057 / 0.127 (bounds 0.021 / 0.075 / 0.135); KS against U(0, 1) p = 0.78, so the slight excess over nominal is Monte Carlo noise, not a shape problem;
  - attack 3, operon-packed sets after the TU collapse: 0.003 / 0.020 / 0.047 (KS p = 0.08, slightly conservative);
  - ties (40 % of TUs tied at the minimum score, as with `abundance.missing = min`; contiguous 15-TU sets) in the seeded order: 0.007 / 0.030 / 0.093 (KS p = 0.19).
- **Reproducibility:** two step runs with one seed are byte-identical (A2, Slurm 1868835).
- **BH per collection, without the floor** (statistician's attacks, Slurm 1869053, `tests/calibration/attacks/test_attack_gsea_r_engine.py`, code at e3c4614; 12 passed; R = 300, records in `$GIPPA_ROOT/logs/calibration/pr214/gsea.jsonl`). GO-like nested collection: 1,000 TUs (E. coli operon lengths, ρ 0.8), 10 parents each the union of 3 children of 10–25 TUs, 20 standalone sets of 10–60 TUs.
  - complete null, FWER (= FDR here) at q = 0.01 / 0.05 / 0.10: 0.007 / 0.050 / 0.123 (bounds 0.021 / 0.075 / 0.135);
  - first child of 3 parents shifted by +1.5, mean FDP: 0.006 / 0.037 / 0.073 (bounds 0.015 / 0.062 / 0.115). Here parents and children are significant together, so the second R call (the GO re-test with the shared seed) runs live; the test asserts that it ran.
- **Structure arms on R** (same job, R = 300, all at level):
  - median imputation: 0.007 / 0.043 / 0.087;
  - size-biased membership: 0.000 / 0.017 / 0.067;
  - scattered between-TU AR(1) (φ 0.5): 0.010 / 0.050 / 0.120.
- **Teeth on R** (same job): attack 3 ranked at gene level (R = 50) and ties sent to R in TU-ID order (R = 300) are both confirmed fails through `calibrate`, so the R arms above can fail. Teeth log to a private directory, so their rates are not in the record file. The port's teeth stay in `tests/calibration/test_gsea_tu_calibration.py`.

### Port (superseded, history only)
These were measured on `fgsea_simple`, whose 1/500 floor kept BH from rejecting; they do not carry over to the R engine.
- **Single-set level at R = 4,000** (statistician, Slurm 1867122): rejection 0.0508 at α = 0.05 under the null. Attack 3, the ties, median imputation and scattered AR(1) are all at level.
- **Family probe** (GO-like nested sets, R = 2,000): FWER 0.0015 at q 0.05 under the complete null. With shifted children, FDR is 0.017 at q 0.05 and 0.053 at q 0.10 (shift +1.5), and 0.0035 / 0.024 (shift +0.8).
- **Committed arms:** tests/calibration/test_gsea_tu_calibration.py (which also holds the port's attack 3 and ties teeth) and tests/calibration/attacks/test_attack_gsea_structure.py (Slurm 1867070: 14 passed, 2 strict xfails).
- **Position flags** (statistician, #174): full tier R = 300 (Slurm 1867323: 33 passed, 2 strict xfails); R = 4,000 per arm (Slurm 1867332; cluster φ 0.5 mix at level, detection as in Limits); BH-family spillover (Slurm 1867330). Committed arms: tests/calibration/attacks/test_attack_gsea_flags.py (detection arms gated at ≤ 1.5α).
