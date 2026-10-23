# M6-09: The M6 gate — missed-operon FDR, discordance level, the attacks, the power grid

Milestone M6 · track C (calibration) · **pass 2** · depends on: M6-03, M6-05, M6-06 (+ M6-02/M6-05H
for the one-detected arm) · estimate ≈ 60 code lines (+ ≈ 450 test lines) · **statistical**
(statistician owns the attacks) · **Slurm**

## What the gate says

§13 row M6: "decoy FDR calibrated on simulations; UniOP and TerminatorNet cut-offs recalibrated
(ADR); UniOP AUROC recorded". This issue owns the **first** clause. The other two are M6-11/M6-12
and are blocked on RegulonDB, so M6's gate is claimed in two parts and PROGRESS.md must say which
part is green.

§7.4.7 "Calibration" lists exactly three criteria, and they are this issue's checks:
1. "missed-operon FDR ≤ 0.1 + Monte Carlo error on simulations with planted co-abundant NO pairs";
2. "the same on the attack null (autocorrelated abundance with no operon effect)";
3. "discordance p-values uniform for OP pairs simulated without discordance".

## Suites and criteria

| Suite | Generator (M6-03) | Criterion (§8.4) | R fast / full |
|---|---|---|---|
| `missed_operons` (**gate**) | `missed_operon_power(1800, 50, 0.1)` | mean FDP ≤ 0.1 + 2 SE at q ≤ 0.1 | 50 / 300 |
| `discordance` | `operon_null(1800)` | level at α ∈ {0.01, 0.05, 0.10}, pooled; KS uniformity p > 0.001 | 100 / 1,000 |
| attack A (gated) | `operon_attack_ar1(1800, 0.5)` | FDR as above | as `missed_operons` |
| attack A′ (recorded) | `operon_attack_ar1(1800, 0.9)` | FDR, **recorded**; escalated if > 0.15 | 20 / 100 |
| attack B (gated) | `operon_attack_detection(1800, "prophage")` | FDR as above | as `missed_operons` |
| attack C (recorded) | `operon_null` with the one-detected arm on (M6-02) | level of `one_detected_pmiss` **alone** | 50 / 300 |
| power (recorded) | `missed_operon_power` over \|Δ\| ∈ {0.25, 0.5, 1.0} SD × 20/50/100 planted | recorded table | 50 per cell |

`tolerances.yaml` already carries `discordance: 100 / 1,000` and `missed_operons: 50 / 300`
(protected file, **unchanged**). Everything goes through M4-05's `calibrate` so ADR-0015's
two-block confirmation and sharding apply (the standing backlog item about migrating suites).

## Behaviour

1. `tests/calibration/test_operon_calibration.py` registers the two suites; the attacks live in
   `tests/calibration/attacks/test_attack_operons.py` as M5-07 does.
2. The fast subset must finish inside 10 minutes on 4 vCPUs (§8.4); the full suites run as
   `scripts/slurm/calibration_full.sbatch` tasks `discordance` and `missed_operons` (≤ 2 jobs,
   ≤ 32 CPUs; sum the running `giap-*` CPUs before submitting, R14).
3. **Evidence recorded in PROGRESS.md:** test names, the measured FDP and level with their SEs, the
   job IDs, the power table, and the explicit note that the §15.4 half of the gate is outstanding.
4. A failure is escalated with numbers and fixed by an ADR on the method — never by editing
   `tolerances.yaml` (R3). The candidate ADR contents are already named: π̂₀ estimation, a third
   decoy class, or separate BH families for the two discordance tests.

## Examples (what a passing run looks like)

| # | Check | Expected |
|---|---|---|
| 1 | `missed_operons`, R = 50 fast | mean FDP ≈ 0.02–0.08, ≤ 0.1 + 2 SE; ≥ 60 % of planted pairs recovered at \|Δ\| = 0.1 SD |
| 2 | `discordance`, R = 100 fast | rejection rate at α = 0.05 within 0.05 ± 2√(0.05·0.95/100) = 0.05 ± 0.044 |
| 3 | KS of the both-detected p | p > 0.001 |
| 4 | attack A (ρ = 0.5) | FDP ≤ 0.1 + 2 SE |
| 5 | attack A′ (ρ = 0.9) | recorded; a value > 0.15 opens an ADR, it does not fail the merge silently |
| 6 | attack C | recorded level of the one-detected arm; if it exceeds 2× α the M6-05 Q3 decision is revisited |
| 7 | the power grid | 9 cells recorded in `docs/methods/operons.md` with R and seeds |
| 8 | reruns with the same seed | identical numbers (fixed seeds in `tolerances.yaml`) |

## Acceptance criteria

- [ ] Both suites registered through `calibrate`, with fast and full variants, and listed in the
  nightly array.
- [ ] Criteria 1–3 of §7.4.7 pass at full R on Slurm, with the output committed as evidence.
- [ ] Attacks A and B gated; A′ and C recorded with their numbers in `docs/methods/operons.md`.
- [ ] The power grid recorded (9 cells × R = 50).
- [ ] `benchmarks/baseline.json` gains the fast-suite runtimes.
- [ ] PROGRESS.md records the gate as **partially** met, naming M6-11/M6-12 as the remainder.

## Hardening owed

- Re-run the full suites after any change to M6-01's FDR formula, and after M6-12 changes the
  default UniOP threshold (the simulations plant their own probabilities, so this is a check, not a
  dependency).
- Add attack D (plasmid copy-number differences inside a run) if the M10 pilot shows it matters.

## Questions for the lead

1. Is a **recorded** (not gated) attack A′ at ρ = 0.9 acceptable, as ADR-0021 did for islands?
   **Proposed default:** yes, with an escalation threshold of 0.15 written into the test.
2. The gate text says "decoy FDR calibrated on simulations" — does attack B (detection confounded
   with distance) belong to the gate or to the record? **Proposed default:** gate it; it is cheap
   and it is the ADR-0018 failure mode in a new place.

## Khoa's decision (2026-10-22): missed-operon power floor

Option 1 of docs/methods/operons.md "Fix space", **keep the method and label it** (no method ADR):
- The missed-operon calls are reported as they are, and lfdr is the ranking.
- The report and M6-09's power table state that this family is "near-powerless by construction":
  a call needs near-identical abundance of the two genes, and the absence of calls is not evidence
  against missed operons.
- The `missed_operons` FDR gate is recorded as a **vacuous pass** when it makes no calls, and the
  record names that reason. It is never read as power evidence.
- Option 2 (a stronger combined score) may be revisited after the M10 pilot.
