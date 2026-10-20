# ADR-0019: Stratum-0 π̂_wrong for database input from exact-match shares per length class

- Status: accepted-by-khoa (2026-10-20)
- Date: 2026-10-20
- Issue: M3-21 (branch m3/calibration-gate); statistician validation, Slurm 1860123, 1860138, 1860180
  (scripts and outputs under local/scratch/stat-m321-m418-m419/b2/)
- Spec sections: §7.2.5 step 4 (database, m = 0), §5.3 (run aggregation, unchanged), §8.1 "Peptide placed
  at a locus", §8.4 "Placement FDR" and "De novo precision", §8.5; extends ADR-0010

## Context
§7.2.5 step 4 sets the stratum-0 π̂ for database input to the mean PEP of the input peptides. M3-21 found
two defects:

1. **B2: cross-strain input.** With `search_db = other`, some correct identifications are database-strain
   sequences that span a strain difference. They have no exact true locus, so an exact chance hit elsewhere
   is a false stratum-0 placement. Their PEP measures identification error only, so it is small. Mean PEP
   counts only the identification errors among peptides without an exact origin, which under-estimates
   false placements by about (π_wrong + δ)/π_wrong; δ is the spanning share of correct identifications.
2. **Run collapse.** §5.3 collapses a peptide's runs with `pep = min`. The mean of min-collapsed PEP is not
   a peptide-level wrong share: with 3 runs M3-21 measured π̂ ≈ 0.0026 against a true wrong share of
   ≈ 0.024 (Slurm 1860161). No PSM-level aggregation (min, mean, max, product) is calibrated in general.

A single pooled π̂ is also invalid when the no-origin share correlates with the chance-hit rate, which
changes by about ×14 per residue of length (short wrong identifications; an input padded with long exact
peptides).

## Decision
1. **Estimator.** Stratum 0, database input. The population P is the stratum-0 database population: all
   database rows except `same_db` exact peptides (unchanged from M3-10). For each length class ℓ of
   `mapping.length_classes`:
   - n_ℓ = number of peptides of P in class ℓ;
   - k_ℓ = number of them with no exact target hit (m_T ≠ 0, or no hit at all);
   - r_ℓ = share of them with an exact decoy hit;
   - π̂_ℓ = 1 if r_ℓ = 1; otherwise π̂_ℓ = min(1, U_ℓ / (1 − r_ℓ)), where U_ℓ is the one-sided 95 %
     Clopper–Pearson upper bound of k_ℓ/n_ℓ (U_ℓ = 1 when k_ℓ = n_ℓ).

   This is the database analogue of the de novo precision estimator with c_l = 1: a peptide with an exact
   origin always matches; one without matches at the decoy rate. The upper bound mirrors step 4's "1 −
   lower 95 % bound". Edge cases: r_ℓ = 0 gives π̂ = U_ℓ; m < r is clipped to 1; an empty class has no
   weight and does not enter c.
2. **Entry into the FDR.** Every peptide of P carries the π̂ of its class. The stratum-0 FDR is ADR-0010's
   FDR_0(e) = (W_0(e) + c)/T_0(e), with W summing π̂ over decoy hits and c = max π̂_ℓ over classes with
   n_ℓ > 0. A peptide outside every class gets π̂ = 1. PEP and the q-value fallback are not used for the
   stratum-0 database π̂. m ≥ 1 is unchanged (π̂ = 1).
3. **Run aggregation.** §5.3's `pep`/`q_value` = min stays for ranking, detection and evidence; placement
   no longer reads it.
4. **Outputs.** `mapping_fdr.tsv` columns are unchanged; `pi_wrong` = (W + c)/(D + 1) per row now reflects
   the class π̂. `detection_summary.tsv` gains QC rows: the stratum-0 unmapped share k/n (the §7.5.1
   genome-inconsistent rate), overall and per class. Run warning `input_prefiltered` when k falls below the
   0.001 binomial quantile of n × (row-level mean PEP) × (1 − r): an input already filtered to genome
   matches cannot be calibrated by any genome-based π̂.
5. **De novo precision estimand (§8.4 "known precision per bin").** The estimand is the correct share among
   distinct de novo peptides per GIAP score bin (bins on the aggregated, max score); a distinct peptide is
   correct when its string (I = L) equals the `true_peptide` of at least one of its reads. Calibration gates
   against that realised share; giapsim's per-read `truth/denovo_precision.tsv` is the per-spectrum design
   precision, a different quantity.

## Evidence
FDP = mean over replicates; failing = some level with mean FDP > q + 2 SE.

- **List level through `placement`** (Slurm 1860138, 1860180; N = 4,000, R = 300 and N = 300, R = 1,000;
  δ ∈ {0, .02, .05, .1, .2, .3} × π ∈ {.005, .01, .05} × q ∈ {.01, .05, .10}). Failing cells of 18:

  | Scenario | mean PEP (spec) | pooled exact share (± bound) | per class, no bound | **this ADR** |
  |---|---|---|---|---|
  | flat chance rate p = 0.25 | 13 (9.97×) | 0 | 0 | **0** |
  | flat p = 0.05 | 6 (2.44×) | 0 | 0 | **0** |
  | long peptides span more | 6 (2.54×) | 0 | 0 | **0** |
  | short wrong / spanning peptides | 13 (8.56×) | 13 (2.4×) | 2 (1.03×) | **0** |
  | padding with long exact peptides | 10 (4.81×) | 10 (2.4×) | 10 (1.09×) | **0** |
  | N = 300 | 13 (9.96×) | 0–1 | 0 | **0** |
  | PEP deflated 10× (min-collapse) | 14 (11.2×) | 1–13 | 1–5 | **0** |

  Power at δ = 0 relative to mean PEP: 0.984–1.0 (N = 4,000), 0.85–0.97 (N = 300); π̂ = 1 instead gives
  0.43–0.84 and 0.25. An oracle using the true unweighted no-origin share fails the short-peptide and
  padding scenarios (2.19×, 2.02×), so no pooled π̂ can work there.
- **Genome level, real chain** (Slurm 1860123; giapsim 100 kb / 1 Mb, R = 300; 5 Mb, R = 100; tryptic,
  7-mer and M3-18 strain arms with realised δ 0.03–0.24). Stratum-0 FDP ≤ 0.0022 for every rule, but the
  estimated/true false count is 0.03–0.46 under mean PEP (q-values 6–30× too small) and 0.96–2.4 under this
  estimator. FDP exceeds q once the list binds: short-peptide-heavy input, very large genomes, or q < ~0.005.
- **Run collapse.** With input PEP × 0.1 this ADR's results are unchanged.
- **Limitation.** If the upstream pipeline keeps only genome-matching peptides, every genome-based rule
  fails, mean PEP included (FDP 0.061 at q = 0.01, δ = 0.2); hence the warning in decision 4.

## Consequences
- `placement_fdr.placement`: the database stratum-0 π̂ comes from the class estimator, not
  `pi_wrong_database`. Stratum-0 q-values and `fdr_stratum` become honest under cross-strain input.
  Genome-scale acceptance at q = 0.01 is unchanged in simulation; small or binding lists lose up to ~15 %
  power at δ = 0.
- ADR-0010's forwarded item (mean PEP correlated with chance-hit rate) is resolved along length.
- Tests to update with the implementation: database examples in tests/unit/test_placement_fdr.py; the
  `known_pi` / `calibrated_pep` / `max_q` expectations in tests/calibration/test_placement_fdr_lists.py; the
  reference rule in tests/calibration/attacks/test_attack_placement_same_db.py; re-verify
  test_placement_fdr_genome.py. New attack: tests/calibration/attacks/test_attack_placement_strain_span.py. Strict xfails to flip: the 1-Mb related_dbstrain estimate check and `test_database_stressed_stratum0_estimate_covers_false_placements` (1-Mb stressed arm; W0/V0 ≈ 0.2 under mean PEP, ≈ 1.04 under this rule, statistician review of #113 round 2).
- GIAP_SPEC.md §7.2.5 step 4 is superseded by this ADR; docs/methods/placement-fdr.md (new section
  "Cross-strain input and run collapse") and docs/methods/denovo-precision.md (estimand) are updated.

## Alternatives
- **Pooled π̂ = max(mean PEP, 1 − (m − r)/(1 − r))** (the first candidate): fails length-correlated input
  and padding (2.4×) and is marginal on small lists.
- **π̂ = 1 for `search_db = other`**: valid everywhere, including pre-filtered input, and free at genome
  scale in simulation, but loses 16–75 % power on binding lists. Kept as the fallback if Khoa prefers
  simplicity.
- **Change the run aggregation of PEP** (mean or max across runs): no aggregation is calibrated when run
  counts differ between correct and wrong peptides, and it does not address B2.
