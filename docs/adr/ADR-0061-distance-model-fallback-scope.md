# ADR-0061: Scope of the TU-blocked claims under the distance-model operon fallback
- Status: proposed
- Date: 2026-10-26
- Spec sections: §7.1.3, §7.3.5, §7.4.2, §8.1, §8.5; R3, R4; `docs/specs/M11-14-container-install.md`
  14s1 (outcome rule), 14d, 14d2; ADR-0022 and its 14s2 addendum
- Needs Khoa: yes. It adds a row flag and a caveat, it bears on the 14d2 default, and one
  finding concerns the UniOP path too (question 2).

## Context

When UniOP is absent or fails, the ADR-0022 distance model (d₀ 102.69, s 47.16, linked when
d ≤ 102 bp) builds the pair classes and the TUs. Most §8 nulls use those TUs as their blocks.
14s1 attacked every TU-blocked claim with the fallback as the TU source, in
`tests/calibration/attacks/test_attack_operon_distance_fallback.py`:
- arm A re-derives the classes of the M6-09 operon scenarios;
- arm B uses K-12-like gaps per class, so true TUs split at 8.4 % and merge at 32.5 % per pair.

Both arms ran at the existing α, bounds and R. The spec's outcome rule applies: a failing check
does not move its bound. Instead, the failing rows are flagged, caveated and scoped here.

## Evidence

Full tier, Slurm 1895399; per-statistic re-run and the K-12-rate control, Slurm 1895506. The
fast tier passed (1895398, 1895507), but its R and genome size do not resolve either failure.
Logs are under `$GIAP_ROOT/stat/m11-14s1/`. Rates are given at α 0.01 / 0.05 / 0.10, followed
by the bound.

| Claim (gated, distance model) | Result |
|---|---|
| Missed-operon decoy FDR, 3 scenarios × shipped / non-vacuous, R 300 | pass (non-vacuous mean FDP at q 0.1: 0.037 null, 0.048 AR(1), 0.031 detection; shipped vacuous, as with UniOP) |
| Discordance pooled level, 3 scenarios, R 1,000 | pass, **conservative**: 0.0002 / 0.0078 / 0.033 on `operon_null` (bounds 0.0104 / 0.051 / 0.101) |
| Discordance BH-family FDR, R 1,000 | pass (0 false calls) |
| Discordance KS uniformity of null p, R 1,000 | **fail** (KS p 0 in both blocks), conservative direction |
| GSEA (a) median imputation, (b) size-biased, (c) AR(1) scattered, R 300 | pass: 0.0067 / 0.043 / 0.100; 0.0067 / 0.010 / 0.040; 0 / 0.060 / 0.100 (bounds 0.0215 / 0.075 / 0.135) |
| ORA 4b (c 1.5, 2.5), R 1,000 | rates fail on large-TU categories, under tail (c 2.5: 0.022 / 0.084 / 0.143 against 0.012 / 0.054 / 0.106). This passes by the existing rule, since every ORA row already carries `ora_uncalibrated` (methods "Pass-1 limit"). UniOP-like errors: 0.016 / 0.065 / 0.119 |
| Determinants TU-block CI, R 500 | rows flagged `ci_unestimated` (existing rule). Miss rates ρ 0.054, AUROC 0.038 (bound 0.069). With the feature constant within *true* TUs (recorded): 0.074 / 0.072 |
| Autocorrelation TU-block, join count, lags 1 and 20, R 1,000, 4,500 genes | pass |
| Autocorrelation TU-block, **Moran's I lag 1** | **fail**: 0.152 / 0.342 / 0.445 (bounds 0.016 / 0.064 / 0.119), confirmed in block 2 |
| Autocorrelation TU-block, Moran's I lag 20 | pass: 0.012 / 0.050 / 0.101 |

Two controls on the same worlds were recorded and not gated:
- UniOP-like Beta errors (2 % splits): Moran's I lag 1 gave 0.017 / **0.068** / 0.110, a
  confirmed fail at α 0.05.
- UniOP's own K-12 error rates at 0.5 (14.7 % splits, 17.0 % merges, as independent per-pair
  errors): **0.431 / 0.668 / 0.763**.

**Mechanism.**
- *Moran's I.* A split true TU becomes two analysis blocks whose abundances share the TU
  effect. The TU-block permutation moves them apart, so the observed lag-1 similarity counts as
  signal. Merges only make the null conservative.
- *Discordance.* Merged boundary pairs enter the OP family. Their Δ comes from two independent
  TUs, which widens the fitted t, so the p of true within-TU pairs pile up near 1. This costs
  power and is not anti-conservative. Calls on merged pairs are true discoveries, because a
  promoter lies between the genes.

## Decision (proposed)

1. **Fallback rows are flagged.** When `operon_source = distance_model`, the following rows
   carry `tu_source_distance_model`, and the report states the caveat:
   - every `morans_i` row with `null_model = tu_block` (all lags; only lags 1 and 20 are
     gated, so the short lags in between are not cleared). Caveat: "Moran's I under the TU-block
     null is anti-conservative at short lags when TUs come from the distance model; rejection
     at α = 0.05 was 0.34 in simulation";
   - every discordance `both_detected_t` row. Caveat: "conservative under the distance-model
     fallback (merged TU boundaries widen the within-operon fit); fewer discordant pairs are
     found".

   Implementation belongs to 14d (implementer). The flag joins `RUN_FLAGS`/row-flag schemas
   there. The two failing checks are strict xfails naming this ADR (ADR-0049 / ADR-0057
   pattern).
2. Every other TU-blocked claim attacked in 14s1 is cleared under the fallback at its existing
   bound. ORA and determinants stay under their existing flags.
3. No bound, α, R or tolerance moves (R3). No threshold is re-chosen (14s2: precision 0.828
   [0.805, 0.847], recall 0.916 [0.901, 0.931] at 0.5).

## Questions for Khoa

1. Accept decision 1? If accepted, 14d2 (`enabled: auto` as the default) can ship with the flag
   and the caveats. The fallback then delivers at-level results everywhere 14s1 looked, except
   Moran's I at short lags (anti-conservative, flagged) and discordance (conservative).
2. **The Moran's I finding is not specific to the fallback.** With UniOP's K-12 error rates, the
   TU-block Moran's I at lag 1 rejects 0.67 at α 0.05 under the null, about twice the fallback's
   0.34, because UniOP splits more TUs than the fallback does (14.7 % vs 8.4 % on K-12). The
   M5-05 calibration assumed exact TUs. The statistician recommends a separate issue: either
   flag short-lag TU-block Moran's I rows on every run, or make the TU-block null robust to
   split TUs, for example by joining adjacent analysis TUs whose pair probability is in an
   uncertain band. Until then, the claim stays as M5-05 shipped it. The independent-error model
   may overstate the effect if UniOP's real errors cluster.

## Consequences

- 14d emits the flag and the caveat. 14d2 needs this ADR accepted and the statistician's
  sign-off.
- The attack file re-runs nightly (fast) and in `calibration_full`. If a future TU model
  passes, the strict xfails XPASS and force this ADR to be revisited.
