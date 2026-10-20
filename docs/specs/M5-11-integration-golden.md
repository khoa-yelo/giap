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
