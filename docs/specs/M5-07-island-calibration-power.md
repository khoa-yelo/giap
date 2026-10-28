# M5-07: Island FWER calibration (the M5 gate), attack null 1 and the recorded power curve

Milestone M5 · track C (calibration) · depends on: M5-06 (islands), M5-04 (scenarios); M4-05 harness
(merged: `calibrate`, `check_fwer`, two-block confirmation ADR-0015) · estimate ≈ 40 code lines (power
recorder script under `scripts/dev/`) + ≈ 350 test lines (test-author owns the calibration and
statistician owns the attack) · **statistical** (label `stats`) · **needs Slurm**: full mode (R = 500 × 3
backgrounds) and the power grid run as `gippa-` batch jobs (`scripts/slurm/calibration_full.sbatch`,
existing array task `islands`)

## Spec sections implemented

- §13 row M5 gate quoted: "island FWER ≤ 5 % under simulated nulls; power curve recorded".
- §7.3.6 "Calibration: island FWER ≤ 5 % under AR(1) backgrounds without islands"; "Power (recorded, not
  gated): islands of 5, 10 and 20 TUs with effects of 1, 1.5 and 2 SD".
- §8.4 "FWER. The share of replicates with ≥ 1 false discovery is ≤ 0.05 + 2 SE"; table row "Islands |
  AR(1) backgrounds without islands | FWER | 50 | 500".
- §8.5 attack 1 "Autocorrelated abundance with no islands → islands, autocorrelation".

## Public interface (tests and a recorder; no `src/` change)

- `tests/calibration/test_island_fwer.py` (suite key `islands`, so `-k calsuite_islands` selects it and
  the existing nightly array task runs it):
  - replicate i: `ar1_background(n_genes, ρ)` from M5-04 with `replicate_rng("islands", i)`; gene table
    from the simulator truth (detected = `detected_true`, `p_detect = p_detect_true`, abundance
    `abundance_log2_true` for detected genes, TUs from the layout); run M5-06 with
    `stats.bootstrap = 1000`, `stats.permutations = 1000`, default `context.islands`;
  - indicator = 1 if ≥ 1 island is reported (every island is false: no island is planted);
  - criterion `check_fwer(indicators, 0.05)` through `calibrate(..., criterion="fwer")`.
  - Fast (per PR): 1,000-gene chromosome, ρ = 0.5, R = 50 (`replicates_fast.islands`).
  - Full (nightly / gate): 4,500-gene chromosome + one 60-gene plasmid, ρ ∈ {0, 0.5, 0.9}, R = 500 each
    (`replicates_full.islands`); one stream per ρ (`islands.rho0`, `islands.rho05`, `islands.rho09`).
- `tests/calibration/attacks/test_attack_islands.py` (statistician): attack 1 variants, FWER check, full
  mode only: (a) ρ = 0.95 (strong background, block length capped at 20); (b) `oriter_gradient` with
  2 log2 cosine at ρ = 0.5 (LOESS must absorb it); (c) `tu_null` with mean TU size 3 (large TUs, no
  inter-TU structure); (d) `detection = "abundance"` at ρ = 0.5 (detection follows abundance, so
  censoring carries signal).
- `scripts/dev/island_power.py` (implementer; runs inside a Slurm job): grid n_tus ∈ {5, 10, 20} × effect
  ∈ {1, 1.5, 2} SD × direction ∈ {high, low} on the 4,500-gene layout at ρ = 0.5, R = 100 per cell
  (`island_power` scenario, seeds `calibration_base` + cell index); per cell writes one JSON line to
  `local/logs/calibration/islands_power.jsonl`: `n_tus`, `effect_sd`, `direction`, `R`, `power` (share of
  replicates with a same-direction island sharing ≥ 1 TU with the planted one), `mean_jaccard` (TU
  Jaccard of the best overlapping island with the planted island, over hits), `false_island_rate`
  (share with an island not touching the planted one), commit and run ID (lead decision 14 of M4-17).

## Examples

| # | Scenario | Expected |
|---|---|---|
| 1 | fast, R = 50, ρ = 0.5 | ≤ 5 of 50 replicates with an island (bound 0.05 + 2·√(0.05·0.95/50) = 0.1116) |
| 2 | full, R = 500 per ρ | rate ≤ 0.0695 for each ρ ∈ {0, 0.5, 0.9} (two-block confirmation applies) |
| 3 | attack (a)–(d) | each ≤ 0.0695 at R = 500 |
| 4 | power cell 20 TUs, +2 SD | recorded; expected high (the planning target is ≥ 0.8; not gated) |
| 5 | power cell 5 TUs, ±1 SD | recorded; expected low |
| 6 | the same seeds twice | identical indicator vectors |

## Edge cases

- A replicate whose chromosome happens to give ℓ = 20 at ρ = 0 is legal; the check is over replicates.
- Planted islands are kept away from regions (`island_power` guarantees no overlap) so region effects do
  not add power.
- The fast test must stay < 60 s in the 8-CPU session (else it is marked `calibration_nightly`, as the
  harness allows); measured time recorded in `benchmarks/baseline.json`.

## Acceptance criteria

- [ ] Fast calibration green in CI (`calsuite_islands`), under 60 s.
- [ ] Full calibration: one green run of the `islands` array task on main as a Slurm job (job ID and the
  three FWER rates in PROGRESS.md) — this is the **M5 gate** (same practice as the M4 gate, until the
  nightly scrontab entry exists).
- [ ] Attacks (a)–(d) at level (statistician sign-off; a failing attack blocks the merge, §8.5).
- [ ] Power grid run once as a Slurm job; the 18 cells summarised as a table in
  `docs/methods/islands.md` and PROGRESS.md (the "power curve recorded" half of the gate).
- [ ] `tests/calibration/tolerances.yaml` unchanged (protected).

## Methods note (statistician sign-off)

`docs/methods/islands.md` §"Calibration and power": the FWER table per ρ and attack, the power table, the
block-length distribution observed per ρ, and the limitation that blocks underrepresent dependence at
joins.

## Questions for the lead

1. Which AR(1) strengths make up the gate? **Proposed default:** ρ ∈ {0, 0.5, 0.9} full, ρ = 0.5 fast;
   ρ = 0.95 as an attack (pass required, but it is the statistician's).
2. If FWER fails at ρ = 0.9 (blocks of ≤ 20 TUs too short for the dependence): the fix changes a method
   detail (for example the 20-lag bound or a stationary bootstrap) and needs an ADR; tolerances never
   change (R3). **Proposed default:** escalate with numbers, do not tune.
3. Power definition: "same-direction island sharing ≥ 1 TU with the planted island". **Proposed default:**
   as stated, plus the Jaccard for location accuracy.

## Lead decisions (2026-10-20, after the statistician's attack check)

ADR-0021 accepted by Khoa 2026-10-20. When accepted:
28. Q1: full gate ρ ∈ {0, 0.5, 0.9} at 4,500 genes + plasmid; fast gate ρ 0.5 at 1,000 genes; attack ρ 0.95
    required; gated attacks: draft assembly, short-gene prophage, gradient, large TUs, abundance-driven
    detection. (ADR-0021 **Amendment 1** supersedes decision 2: the size rule and the `uncalibrated`
    flag are gone. Gate cells stay at 4,500 genes; add **recorded, not gated** cells at 1,000 and 1,800
    genes x rho 0.9, and label the low-direction row a **structural zero** with the ceiling and the
    critical value beside it, per M5-06 decisions 47/50.)
29. Q2: escalate with numbers, never tune (done: ADR-0021).
30. Q3: power = same-direction island sharing ≥ 1 TU with the planted one, plus Jaccard and false-island
    rate; report the null FWER next to each power cell. Power is recorded, not gated; run the grid also with
    `detection = "abundance"`. The ≥ 0.8 planning target is not met under gippasim defaults (≈ 0.18 at
    20 TUs × 2 SD).

## Lead note (2026-10-21, from the statistician's #131 review)
The short-gene prophage arm of M5-06 decision 27 is **missing** from
tests/calibration/attacks/test_attack_islands.py: `scenarios.tu_null(1800, short_genes_in="prophage")`
through `idat.scenario_contigs`, one gated arm. It is the same mechanism class as the ori–ter breach
(clustered censoring → local-scale heterogeneity), so it must be gated rather than assumed. **It must
land before M5-07 claims the island gate.**
