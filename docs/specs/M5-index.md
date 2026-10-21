# M5 "Autocorrelation and islands": issue index

Spec: GIAP_SPEC.md §13 row M5. Deliverables: "Moran's I and join counts with both nulls, islands,
ori–ter". **Gate:** "island FWER ≤ 5 % under simulated nulls; power curve recorded", i.e. the `islands`
calibration suite (R = 50 per PR / 500 nightly, `tests/calibration/tolerances.yaml`, unchanged) passes the
§8.4 FWER check under AR(1) backgrounds without islands, and the power grid of §7.3.6 (islands of 5, 10,
20 TUs × 1, 1.5, 2 SD) is recorded. The `autocorrelation` suite (R = 100 / 1,000, level) ships with it
(§7.3.6 "autocorrelation tests at level under the TU-block null"), but is not part of the M5 gate text.

Sources: §7.3.5(c)–(e), §7.3.6, §8.1–§8.6 (claims, rules, engine functions `morans_i`, `join_count`,
`scan_max`, `circular_block_bootstrap`, `blom_scores`; calibration rows "Autocorrelation" and "Islands";
attack 1; sensitivity family "islands"), §9.2–§9.4 (`autocorrelation.tsv`, `islands.tsv`, genes.tsv
context columns, flags), §6.2 (`context_stats`), §6.5 (`ISL_<n>`), §6.6, §11.3 (AR(1), planted islands,
oriC), §10.2 V2/V10/V18 (consumers in M9), Appendix B.4. ADRs read: ADR-0015 (two-block confirmation,
applies to both suites), ADR-0018 (region tests; its lessons — detectability confounding, label-dependent
variances, fail-closed on non-positive variance — are applied to the M5 questions below), ADR-0007
(size counts code).

One issue = one PR; size counts **code** only (ADR-0007; simulator lines under `tests/sim/` are tooling and
listed separately). "Stats" = statistical change: label `stats`, the statistician reviews the methods note
and simulations before the code (§14.1).

| Order | Issue | Track | Scope | Depends on | Code lines | Stats | Tools / data / Slurm |
|---|---|---|---|---|---|---|---|
| 1 | [M5-01](M5-01-stats-blom-loess.md) | S | `stats/`: `blom_scores` (left-censored), `loess_detrend` (circular padding) | none | 200 | yes | none |
| 2 | [M5-02](M5-02-stats-moran-joincount.md) | S | `stats/`: `morans_i`, `join_count`, matrices, label and TU-block permutations, `lag_test`, `choose_block_len` | none | 300 | yes | none |
| 3 | [M5-03](M5-03-stats-scan.md) | S | `stats/`: `scan_max`, `scan_null_max`, `fwer_p`, `merge_windows` (numba allowed) | none | 240 | yes | none |
| 4 | [M5-04](M5-04-sim-oric-gradient.md) | F | C15 `features` contract entry; giapsim oriC, ori–ter gradient, replichore truth, null/power scenarios | none | 50 (+200 sim) | no | none |
| 5 | [M5-05](M5-05-autocorrelation.md) | R | Moran + join count per replicon, both nulls, BH across lags; **autocorrelation calibration**, attack 1 | 01, 02 (cal: 04) | 330 | yes | full calibration on Slurm |
| 6 | [M5-06](M5-06-islands.md) | R | scores → TU means → detrend → ℓ → bootstrap max → FWER → merged islands, effect CI, gene map | 01, 02, 03 | 380 | yes | none |
| 7 | [M5-07](M5-07-island-calibration-power.md) | C | **island FWER gate** (fast + full), attack 1 variants, **power grid** recorder | 04, 06 | 40 (+tests) | yes | **Slurm** (R = 500 × 3 ρ; 18-cell power grid) |
| 8 | [M5-08](M5-08-oriter.md) | R | oriC from C15, replichore columns, OLS per replichore + pooled, block-bootstrap CI, `no_oric` | 02, 04 | 260 | yes | recorded Bakta JSONs (K-12, Nissle, sim100k; already in the repo) |
| 9 | [M5-09](M5-09-context-tables-step.md) | T | writers + schemas for `autocorrelation.tsv`, `islands.tsv`; flags; `gene_context.parquet`; `context_stats` step extension | 05, 06, 08 + **M4-15, M4-18** (+ M4-20 for its integration test) | 320 | no | none |
| 10 | [M5-10](M5-10-sensitivity-islands.md) | R | sensitivity family `islands` (Jaccard of island gene sets) | 06, 09 + **M4-19** | 100 | yes (light) | none |
| 11 | [M5-11](M5-11-integration-golden.md) | T | end to end on fixtures, golden `sim50k` + `sim1000_island`, TE7, benchmark, gate evidence | 01–10 + **M4-20** | 100 | no | Slurm for the full array tasks |

Total ≈ 2,320 code lines in `src/` (+ ≈ 200 simulator lines, ≈ 3,300 test lines). Effort fits the
spec's 1 week if the M4 tail (M4-15, 18, 19, 20) lands first or in parallel.

## Dependency DAG

```text
Track S (engine)       01        02        03                (mutually independent; 3 parallel PRs)
Track F (simulator)    04
Track R (statistics)   01, 02 ─► 05 (cal ◄─ 04)
                       01, 02, 03 ─► 06 ─► 07 ◄─ 04           (07 = the gate)
                       02, 04 ─► 08
Track T (tables)       05, 06, 08 + M4-15, M4-18 ─► 09 ─► 10 (◄─ M4-19) ─► 11 (◄─ M4-20)
```
Unit tests of 05, 06 and 08 use small hand-built frames, so the code can start as soon as the engine
issues merge; only calibration files wait for 04. Only 09–11 wait for unmerged M4 issues.

## Parallel tracks and waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01, 02, 03, 04 | four independent PRs; 04 is tooling (implementer) plus a 50-line contract entry |
| B | 05, 06, 08 | after A; 06 is the critical path to the gate |
| C | 07 (gate, Slurm), 05 full calibration (Slurm) | test-author + statistician; submit `scripts/slurm/calibration_suite.sbatch` (the M4-05 harness array, renamed from `calibration_full.sbatch` by CI-01) tasks `islands`, `autocorrelation` (≤ 2 jobs, ≤ 32 CPUs) |
| D | 09, 10, 11 | after M4-15, M4-18, M4-19, M4-20 merge |

## What M5 consumes (checked against merged code on main, 05ba500)

| ID | Artefact | Producer (status) | What M5 relies on |
|---|---|---|---|
| G | gene table `work/genetable/genes.parquet` | M4-15 (**not merged**) | `locus_tag`, `contig`, `start`, `end`, `partial`, `tu_id`, `detected`, `p_detect`, `abundance_log2` — exactly the columns `giap.context.regions.replicon_index` (merged, M4-17) reads; M5 reuses that index |
| C1 | `replicons.parquet` | M2-04 `parse_bakta` (merged) | `contig`, `length`, `topology`, `kind` (row order = FASTA order) |
| C5 | `tus.parquet` | M2-12 (merged) | `tu_id`, `contig`, `start`, `end`, `genes`; every gene in exactly one TU (M4-17 decision 13: missing `tu_id` raises) |
| **C15** | `features.parquet` | M2-04 `parse_bakta` (merged; **new in the contract**, M5-04) | rows `type = oric`: `contig`, `start`, `end`. K-12 has 2 oriC calls, Nissle 1, sim100k 0 (recorded Bakta 1.12.2 JSONs) |
| E | `giap.stats` | M4-01/03/04 (merged) | `perm_p` (two-sided centred, "greater"), `bh`, `circular_block_bootstrap`, `tu_block_bootstrap`, `percentile_ci` |
| S | `context_stats` step, `run_context_stats` | M4-18 (**not merged**) | M5-09 extends the step; regions/detection-summary outputs must stay byte-identical |
| Z | sensitivity `FAMILIES` | M4-19 (**not merged**) | M5-10 registers `islands` |
| X | fixture builder, golden `sim50k` | M4-20 (**not merged**) | M5-09/M5-11 integration and golden |
| K | config | merged `ContextConfig`: `moran_lags [1, 20]`, `islands {min_tus 2, max_tus 20, detrend_span_frac 0.10, block_len auto}`; `stats.permutations/bootstrap/alpha.context` | no `defaults.yaml` change needed |
| T | calibration harness | M4-05 (merged): `calibrate`, `check_fwer`, `check_level_pooled`, `replicates_fast/full.{autocorrelation, islands}`, nightly array tasks | `tolerances.yaml` unchanged (protected) |
| Sim | giapsim | M4-07/08 (merged): AR(1) per replicon (restarts per replicon, does not wrap at the origin), `IslandSpec`, island truth, `tu_sd`, `detection` modes | M5-04 adds oriC, gradient, scenarios |

Columns M5 fills that M4-15 leaves NA: genes.tsv `island_id`, `replichore`, `replichore_position`.

## What happens when an estimator does not exist (summary; details per issue)

| Estimator | Does not exist when | Output |
|---|---|---|
| Blom score | detected gene without abundance (`exclude`) | gene not scored; TU dropped if empty |
| LOESS | n < 3 points | fit = mean; residuals centred |
| Moran's I_k | invalid lag (circular k > ⌊(N−1)/2⌋), < 2 values, no lag-k pair, zero variance | value, p, q NA; flag `not_estimable` |
| join count BB_k | invalid lag; null constant (no or all genes detected; one TU for the TU null) | value kept when defined; p, q NA; `not_estimable` |
| block length ℓ | < 3 TUs; constant scores | ℓ = 1 |
| scan / island | < 2·min_tus scored TUs on a replicon; no replicon scanned | replicon not scanned; empty `islands.tsv` |
| island effect CI | < B/2 valid replicates | CI NA; < 10 TUs → `ci_few_units` |
| ori–ter slope | no oriC (`no_oric`); < 3 TUs with values or zero variance of position on a replichore | NA row with `descriptive;no_oric` or `descriptive;not_estimable` |
| ori–ter CI | < B/2 valid bootstrap fits | CI NA, `not_estimable` |

No silent fallbacks (the M4-17 rule): every NA carries a flag, and no NA p enters a BH family.

## Calibration suites owned by M5 (replicates from `tolerances.yaml`)

| Suite (key) | Issue | Generator (M5-04 scenario) | Criterion | R fast / full |
|---|---|---|---|---|
| `autocorrelation` | M5-05 (+ engine check in M5-02) | `label_null` (both nulls), `tu_null` (TU-block null) | level, pooled over 20 lags per (statistic, null) | 100 / 1,000 |
| `islands` | M5-07 (**gate**) | `ar1_background(ρ)`, ρ = 0.5 fast; ρ ∈ {0, 0.5, 0.9} full | FWER ≤ 0.05 + 2 SE | 50 / 500 |
| attack 1 | M5-05, M5-07 (statistician) | TU noise + ori–ter gradient; ρ = 0.95; large TUs; abundance-driven detection | level / FWER | as the suites |
| power (recorded) | M5-07 | `island_power`: 5/10/20 TUs × 1/1.5/2 SD × high/low, ρ = 0.5, 4,500 genes | none (recorded) | 100 per cell, Slurm |

## Methods notes (statistician sign-off, §14.1, §14.9)

`docs/methods/autocorrelation.md` (01 detrending, 02, 05), `docs/methods/islands.md` (01 scores, 03, 06,
07 evidence and power table), `docs/methods/oriter.md` (08); updates to `bootstrap.md` §"Block length"
(02) and `sensitivity.md` §"Islands" (10); `multiple-testing.md` already lists the autocorrelation family.

## Risks

- **FWER at strong autocorrelation.** Blocks are capped at 20 TUs; at ρ ≥ 0.9 the bootstrap under-represents
  dependence across block joins, so FWER may exceed 5 % (docs/methods/bootstrap.md). Mitigation: the full
  gate includes ρ = 0.9 and attack ρ = 0.95; a failure is escalated with numbers and fixed by an ADR on the
  method (bound or bootstrap type), never by tolerances (R3).
- **Detrend vs island length.** LOESS span = 10 % of TUs; on a 1,000-gene chromosome (≈ 555 TUs, span ≈
  55 TUs) a 20-TU island is partly absorbed; on small plasmids nearly everything is. Valid, but the power
  curve may look poor on small genomes — recorded on 4,500 genes as the spec intends.
- **Detectability confounding (ADR-0018 lesson).** Islands rank undetected genes lowest; a run of short or
  hydrophobic genes is a low "island" of detection, not of abundance. The spec's estimand is the combined
  detection/abundance score, so this is by design, but the methods note and report text must say so; the
  join count has the same issue (M5-05 Q4). Attack (d) in M5-07 measures it.
- **Unmerged M4 tail.** 09–11 depend on M4-15/18/19/20; 01–08 do not.
- **Runtime.** All CPU-light: per replicon ≤ 2·10⁸ vectorised operations at B = 1,000. Fast suites < 60 s
  each; full suites ≈ 10–30 min on 16 CPUs (Slurm).

## Feasibility on FarmShare

Nothing in the M5 spec is infeasible: no external tool, database or network access; memory < 4 GB;
the full calibration (≈ 1,500 + 2,000 replicates) and the 18-cell power grid fit one or two
`giap-calibration` array tasks within R14. The only operational dependency is the nightly scrontab
entry (Khoa, §14.8 item 7); until it exists the gate is one green full run submitted by hand, as in M4.

## Decisions taken in these specs (for the lead to confirm)

- Engine helpers added beyond §8.3: `loess_detrend` (M5-01), `choose_block_len`, `lag_test` (M5-02),
  `fwer_p`, `merge_windows` (M5-03).
- Moran/join-count index = M4-17's eligible-gene index; lags valid up to ⌊(N−1)/2⌋ on circular replicons
  (M5-02).
- Gene-label null for Moran keeps the missingness pattern fixed (M5-05 Q1).
- Moran/join-count rows: `effect` = value − null mean, CI NA (deviation from §8.2 for this table; ADR if
  accepted) (M5-05 Q2).
- BH family per (contig, statistic, null_model) (M5-05 Q3).
- Island FWER over the whole genome (max over replicons per bootstrap replicate) (M5-06 Q1).
- Detected genes without abundance are not scored (M5-06 Q2); islands reported by their best window
  (M5-06 Q3); threshold from `stats.alpha.context` (M5-06 Q4).
- Multiple oriC calls → circular mean (K-12 case) (M5-08 Q1); replichore `right`/`left`; slope rows
  labelled in `stat` (M5-08 Q2); linear chromosomes skipped (M5-08 Q3).
- New flags `not_estimable`, `few_values` through the schema (M5-09 Q1); C15 `features` added to the
  input contract as optional (M5-04 Q1).
- `context_stats` writes `gene_context.parquet`; the genes.tsv writer joins it (M5-09 Q2).
- Separate RNG streams `context_stats.autocorrelation`, `.islands`, `.oriter`, so M4 outputs stay
  byte-identical (M5-05, M5-06, M5-08).
