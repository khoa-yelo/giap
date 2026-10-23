# M5-11: M5 integration: fixtures end to end, golden tables, determinism, benchmarks, gate evidence

Milestone M5 · track T · depends on: M5-01 … M5-10; **M4-20** (fixture builder, golden `sim50k`) ·
estimate ≈ 100 code lines (fixture extension) (+ ≈ 250 test lines + committed golden TSVs) · **not
statistical** · **Slurm:** the full calibration array tasks `autocorrelation` and `islands` and the power
grid (M5-07) run as `giap-` jobs; nothing heavy in the session allocation

## Spec sections implemented

- §13 row M5: "Deliverables: Moran's I and join counts with both nulls, islands, ori–ter. Gate: island
  FWER ≤ 5 % under simulated nulls; power curve recorded."
- §11.4 golden fixture `tests/data/golden/sim50k/`; §11.6 TE7 ("two runs with the same seed give
  byte-identical tables"); §9.5 Frictionless validation in every integration test.
- §14.9 definition of done (tests at every layer, coverage ≥ 90 %, methods notes signed off, schema
  updated and validated, benchmark recorded, reviews, PROGRESS.md).

## Interface

- `tests/sim/giapsim/fixtures.py` (tooling): `build_m4_fixture(..., oric_frac: float | None = 0.8,
  islands: Sequence[IslandSpec] = ())` — the M4-20 builder also writes `features.parquet` (C15) and can
  plant islands; `sim50k` keeps its M4 golden files byte-identical (the oriC is a new file; the gradient
  stays 0).
- Call order: `run_gene_table` → `run_detectability` → `run_context_stats` (now with M5) →
  `run_sensitivity` → genes.tsv writer.
- Golden: `tests/data/golden/sim50k/autocorrelation.tsv`, `islands.tsv` (header only is acceptable: a
  50-kb genome has ≈ 25 TUs) and a second golden `sim1000_island` (1,000 genes, one planted +3 SD island of
  20 TUs, oriC; strong enough to be found at any seed) with its `autocorrelation.tsv`, `islands.tsv` and the three genes.tsv columns
  (compressed, < 1 MB, R11).

## Examples

| # | Scenario | Expected |
|---|---|---|
| 1 | `sim50k` end to end | `autocorrelation.tsv`, `islands.tsv` written and Frictionless-valid; M4 tables byte-identical to their M4 golden files |
| 2 | `sim1000_island` | one `high` island overlapping the planted TUs; oriter rows with `descriptive`; tables equal the golden files byte for byte |
| 3 | either fixture twice with the same seed | byte-identical `tables/*` (TE7) |
| 4 | `sim1000_island` with `--set seed=…` changed | the planted island is still reported (shares ≥ 1 TU), other numbers may change |
| 5 | the recorded K-12 Bakta JSON through `parse_bakta` + M5-08 | oriC resolved (M5-08 example 3); no proteomics needed |

## Acceptance criteria

- [ ] `tests/integration/test_m5_pipeline.py`: examples 1–4; `tests/golden/test_m5_golden.py` compares
  the golden files.
- [ ] Benchmarks (`benchmarks/test_context_m5.py`): `context_stats` M5 part on the 4,500-gene fixture
  (B = 1,000): wall time and peak memory recorded in `benchmarks/baseline.json`; budget < 2 min on 4 CPUs.
- [ ] Gate evidence in PROGRESS.md: Slurm job IDs and results of the full `islands` and `autocorrelation`
  array tasks on main, the power table (M5-07), coverage of `giap.context.{autocorrelation,islands,oriter}`
  and `giap.stats.{scores,smooth,autocorr,scan}` ≥ 90 %, methods notes `autocorrelation.md`,
  `islands.md`, `oriter.md` signed off.
- [ ] `docs/backlog.md` "M5 autocorrelation and islands" ticked with the PR links.

## Questions for the lead

1. Golden `sim50k` islands table will be header-only (too few TUs). **Proposed default:** keep it (it
   tests the empty path) and add `sim1000_island` for a non-empty golden.
2. Gate run: as for M4, one green full calibration run on main submitted as a Slurm job until Khoa
   installs the nightly scrontab entry. **Proposed default:** yes.

## Lead decisions (2026-10-20, after the statistician's attack check)

39. Q1, Q2: accept; check the `sim1000_island` golden (+3 SD, 20 TUs) against ℓ' before freezing it.

## Lead decisions on the tests (2026-10-22, after fc34e31)
40. **Island golden.** Decision 39's `sim1000_island` reports no island under giapsim's default
    `covariates_only` detection (even at +10 SD; ~45 % of genes censored regardless of abundance, the
    ADR-0018 lesson) and is marginal under `detection="abundance"` (FWER p ≥ 0.035). The golden is
    `sim4500_island` (oriC 0.8, one +3 SD island of 20 TUs, `detection="abundance"`), where the planted
    island is found at every seed tried (FWER p 0.003–0.008, ℓ' 9–13 < 20 TUs). `sim1000_island` stays in
    `tests/helpers/m5_chain.py` for reference. Backlog (statistician, hardening): island power under
    `covariates_only` detection.
41. **oriC default.** As the spec's interface says, `build_m4_fixture(oric_frac=0.8)` by default, so
    sim50k gets an oriC; the implementer regenerates `tests/data/golden/sim50k/m4/genes.tsv` (only the
    replichore columns change) and any existing test that pinned NA replichores on a default fixture
    is reported to the lead. The `oric_frac` keyword wins over a `layout={"oric_frac": …}` pass-through.
42. `expression=` (ExpressionConfig pass-through, as `layout=`) is accepted as an interface addition.
43. The benchmark times the whole `context_stats` step (the M5 part has no public entry point); the
    120 s budget holds a fortiori.

## Lead decisions on the review round (2026-10-22, review of 1e72caf)

44. **Join-count null columns compared numerically.** On `join_count_bb` rows of the island golden,
    `null_mean`, `null_lo`, `null_hi` and `effect` are compared with `math.isclose` (rel 1e-3, abs
    1e-3); every other cell stays byte for byte. Cause: Pearson residuals of `p_detect` from the
    sklearn lbfgs detectability fit, which differs across BLAS builds in the 4th–6th significant
    digit (CI ubuntu 3.11 / macOS 3.12 vs Slurm); observed counts, p and q agree everywhere. Same
    pattern as M4-20 decision 9. The golden is never regenerated on a CI runner.
45. NITs N1, N2, N6 and the N3 fixture pins are fixed by the lead in this round; N4 (the `giap-bench`
    script) and the gate evidence (the full `autocorrelation` record, `calibration_full.sbatch
    --array=4`) go to the backlog and PROGRESS.md.

46. **`stat` joins the numeric join-count columns.** The standardised `stat` of `join_count_bb` rows
    is derived from the same BLAS-dependent null, so it drifts like `null_mean` (621.875 vs 621.876
    on CI ubuntu 3.11, #219 run 39890993815); it is compared with the decision 44 tolerance too.
47. **`value` joins them too.** CI (#229, run 39904322170, ubuntu 3.11 and macOS) shows the observed
    `join_count_bb` value drifting in the fourth decimal (−2.45409 vs −2.45423), because it is
    computed from the same Pearson residuals as the null; it is compared with the same tolerances.
48. **Permutation p and BH q of `join_count_bb` rows** are compared within 5 permutation counts
    (5 / (n_perm + 1)). CI #229 (run 39911158467, ubuntu 3.11) shows line 60 p 0.682318 vs
    0.68331, one count apart, because the drifted observed value moves the permutation count.
