# M8-15: Gate, half 2 — `ora`, `gsea` and `determinants` at level; attacks 2 and 3; recovery

Milestone M8 · track C (calibration) · **pass 2** · depends on: M8-04, M8-05, M8-09, M8-11, M8-14 ·
estimate ≈ 60 code lines (+ ≈ 500 test lines) · **statistical** (the statistician owns the attacks)
· **Slurm**

## What the gate says

§13 row M8 "null p-values uniform"; §7.7.5 "Calibration: ORA p-values uniform when detection
depends only on p_detect; the attack null 'detectability confounded with category membership' …
stays at level with weights and fails without them; the GSEA attack null 'gene sets packed into one
operon' stays at level after TU collapse"; §7.8.3 "Calibration: null p-values uniform when the
simulated response is independent of the features but keeps TU structure; recovery of planted
block effects with known R²"; §8.4 rows ORA, GSEA, Determinants null — all with criterion
**level**. §8.4 also says "Uniformity (continuous p only) … Discrete permutation p-values use the
level check only": Wallenius p are hypergeometric-type (discrete), the determinants p is a
permutation p, and fgsea p are permutation-based with a multilevel tail. So "uniform" is gated as
the §8.4 level criterion for all three; KS is computed and **recorded** (decision 24).

## Suites and criteria

| Suite | Generator (M8-14) | Criterion (§8.4) | R fast / full |
|---|---|---|---|
| `ora` (**gate**) | `ora_null(4000)`, all six collections of the synthetic ontology | level at α ∈ {0.01, 0.05, 0.10}, **pooled** over tested categories (`check_level_pooled`, m = categories per replicate), per direction | 100 / 1,000 |
| attack 2, weighted (gated) | `ora_attack_confounded(4000, −1.5)` with the **fitted** M4-13 p_detect (backlog: never the truth) | level of the confounded categories' over **and** under p | as `ora` |
| attack 2, teeth (gated) | the same, ω forced to 1 (unweighted hypergeometric) | the under-representation rejection rate at α = 0.05 is **> 2 ×** α (the test must be able to fail) | as `ora` |
| attack 2, truth arm (recorded) | the same with p_true | level, recorded beside the fitted arm | as `ora` |
| `gsea` (**gate**) | `gsea_null(4000)`, ≥ 100 sets per replicate | level, pooled over sets | 50 / 300 |
| attack 3 (gated) | `gsea_attack_operon_packed(4000, 20)` | level after TU collapse | as `gsea` |
| attack 3, teeth (recorded) | the same with gene-level ranking | anti-conservative (rate > 2 × α at 0.05) — recorded to show the attack bites | 20 / 100 |
| `determinants` (**gate**) | `determinants_null(400)` fast (decision 23) / `(1500)` full | level of `null_p` per model, R replicates (1 p per model per replicate) | 20 / 100 |
| recovery (recorded) | `determinants_planted(1500, …)` at true R² ∈ {0.2, 0.5} | LightGBM `cv_r2_mean` within ±0.1 of the truth; dominant block = the planted largest in ≥ 70 % of top-decile genes; recorded | 20 |
| centring (recorded, ADR-B) | `gsea_null` and a planted `gsea` power arm, centred vs raw scores | level of both; power difference recorded | 50 |

`tolerances.yaml` already carries `ora: 100 / 1,000`, `gsea: 50 / 300`, `determinants: 20 / 100`
(protected, **unchanged**). Everything goes through M4-05's `calibrate` so ADR-0015's two-block
confirmation and sharding apply.

## Behaviour

1. `tests/calibration/test_function_calibration.py` and `test_determinants_calibration.py`
   register the three suites; attacks in `tests/calibration/attacks/test_attack_function.py`.
2. **Fast tier** (per PR, < 10 min on 4 vCPU, §8.4): `ora` and `gsea` at the stated sizes;
   `determinants` at 400 genes, B = 99, LightGBM `n_estimators = 50`, elastic net as configured,
   1 CV repeat for the observed value (decision 23 — model-size settings are generator settings,
   not tolerances, so R3 is untouched). **Full tier** (nightly Slurm array): §5.2 parameters,
   B = 999, 5 repeats, 1,500 genes; tasks `ora`, `gsea`, `determinants`; sum running `giap-*` CPUs
   before submitting (R14); `determinants` is sharded (`--shard K/N`) so ≤ 48 CPUs run at once.
3. **Evidence** in PROGRESS.md: test names, the level per α with its bound, R, job IDs, the
   recorded arms, the recovery table.
4. A failure is escalated with numbers and fixed by an ADR on the method — never by editing
   `tolerances.yaml` (R3). Named candidates: attack 2 fails → ADR "ORA null = Poisson-binomial of
   Σ_C d given p_detect" (`giap.stats.poisson_binomial_sf`); `determinants` anti-conservative →
   ADR on the null (fresh folds per permutation, or repeats in the null).

## Examples (what a passing run looks like)

| # | Check | Expected |
|---|---|---|
| 1 | `ora`, R = 100 fast | pooled rejection at α = 0.05 ≤ bound (Wallenius is discrete, so typically below α) |
| 2 | attack 2 weighted | at level; unweighted under-representation rate ≈ 0.3–0.6 at α = 0.05 (teeth) |
| 3 | `gsea`, R = 50 | pooled rate at α = 0.05 ≤ bound |
| 4 | attack 3 | TU-collapsed at level; gene-level rate > 0.1 (recorded) |
| 5 | `determinants`, R = 20 fast | ≤ 20 · (0.05 + 2√(0.05·0.95/20)) ≈ 2.95 → at most 2 rejections at α = 0.05 |
| 6 | recovery at true R² 0.5 | LightGBM CV R² in [0.4, 0.6]; planted top block dominant in ≥ 70 % of top-decile genes |
| 7 | reruns, same seeds | identical numbers |

## Acceptance criteria

- [ ] Three suites registered through `calibrate`, fast and full variants, listed in the nightly
  array (`calibration_full.sbatch` task list).
- [ ] Full-tier runs pass on Slurm, outputs committed as evidence; job IDs in PROGRESS.md.
- [ ] Attacks 2 and 3 gated (with teeth); truth, centring and gene-level arms recorded in
  `docs/methods/function.md` / `determinants.md` / `calibration.md`.
- [ ] `benchmarks/baseline.json` gains the fast-suite runtimes; the fast tier is < 10 min on 4 CPUs.
- [ ] PROGRESS.md records M8's gate as met only when **both** M8-13 and this issue are green.

## Hardening owed

- The annotated-universe ORA arm (M8-04 Q1) recorded beside the spec's universe.
- A real-data sanity check (K-12 recording, M8-06 hardening): ORA p-value histogram is near-uniform
  for collections unrelated to detectability (descriptive).

## Questions for the lead

1. Determinants fast tier at reduced model size (400 genes, B = 99, 50 trees, 1 repeat) with the
   full tier at the §5.2 parameters (proposed)? The alternative — the real parameters per PR — does
   not fit 10 minutes.
2. "Uniform" = the §8.4 level criterion for all three suites, KS recorded (proposed)?
3. Teeth thresholds (unweighted attack 2 must reject > 2 × α): acceptable as a gated check?
