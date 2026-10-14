# M4-20: M4 integration on simulated fixtures, golden tables, determinism and gate evidence

Milestone M4 · track T · depends on: M4-01 … M4-19 · estimate ≈ 120 code lines (a fixture-building
script) (+ ≈ 250 test lines + committed golden TSVs) · **can start before M2/M3 merge** (fixture path)

## Spec sections implemented

- §13 row M4: "Deliverables: gene table, detection tiers, detectability, abundance, basic statistics,
  region tests. Gate: circular-shift tests calibrated; gene-level FDR calibrated (picked, MAYU)."
- §11.4: "The 50-kb golden fixture (seed `GOLDEN_SEED` in the tolerances file) has its expected tables
  committed under `tests/data/golden/sim50k/`."
- §11.6 TE7: "two runs with the same seed give byte-identical tables".
- §9.5: "Every integration and end-to-end test validates outputs with `frictionless validate`".
- §14.9 definition of done (tests at every layer, coverage ≥ 90 %, methods notes signed off, schema
  updated and validated, benchmark recorded, reviews, `PROGRESS.md`).

## Interface

`tests/sim/giapsim/fixtures.py` (tooling): `build_m4_fixture(work_dir, seed, size: Literal["sim50k",
"sim1000", "sim4500"], evidence: EvidenceConfig) -> None` runs M4-07 → 08 → 09 and writes the contract
tables (C1–C14) into `work_dir`. `sim50k` uses `layout_from_simulation(simulate(SimConfig(), GOLDEN_SEED))`
plus synthetic regions (one prophage block, one AMR set) so every M4 output is non-empty.

M4 end-to-end call order (as the future step runner will do): `run_gene_table` → `run_detectability`
→ `run_context_stats` → `run_sensitivity`.

## Examples

| # | Scenario | Expected |
|---|---|---|
| 1 | `sim50k`, default evidence (decoys emitted) | tier `decoy_picked`; `genes.tsv`, `detection_summary.tsv`, `regions_stats.tsv`, `sensitivity.tsv` written; each validates with frictionless against its resource |
| 2 | example 1 vs committed `tests/data/golden/sim50k/*.tsv` | byte-identical |
| 3 | example 1 run twice in two directories | byte-identical outputs (TE7) |
| 4 | `sim1000` with `emit={"q_value", "pep"}` | tier `mayu`; `detection_summary` has `fdr_2pep` curve rows |
| 5 | `sim1000` with `emit=set()` | tier `two_peptide`; `single_peptide_unestimated` genes present; sensitivity has `single_peptide_detected` |
| 6 | `sim1000` with `proteins_only=True`, `abundance_table="maxlfq"` | tier `protein_list`; `abundance_scale = relative` |
| 7 | `sim1000` with `linear_contigs=(30,)` and a plasmid | partial genes `partial_excluded`; plasmid rows `dosage_plus_regulation` |
| 8 | contract violation (a gene missing from TUs) | `run_gene_table` stops with the contract error; no table written |

## Edge cases

- Golden files are regenerated only by an explicit script (`scripts/golden/make_m4_golden.py`, run by
  the implementer) and reviewed as data; a test diff on them needs the statistician's note in the PR.
- The golden fixture does not depend on M2/M3 code; M4-21 re-records it through the real path.

## Acceptance criteria

- [ ] `tests/integration/test_m4_end_to_end.py`: examples 1, 3–8; `tests/golden/test_m4_golden.py`:
  example 2.
- [ ] Gate evidence collected in `PROGRESS.md`: calibration logs for `gene_fdr` (picked; MAYU with and
  without PEP) and `regions` (both stats, each type, α ∈ {0.01, 0.05, 0.10}) from the nightly full
  suite (R = 300 / 1,000), plus the fast-suite runs on the PR.
- [ ] Coverage of `src/giap/stats`, `genetable`, `context` ≥ 90 %.
- [ ] Methods notes signed off: `multiple-testing`, `target-decoy`, `mayu-gene-fdr`,
  `permutation-and-shifts`, `bootstrap`, `calibration`, `detection-tiers`, `detectability`,
  `abundance`, `detection-summary`, `region-tests`, `sensitivity`.
- [ ] Benchmarks recorded (`benchmarks/baseline.json`): gene table, detectability, context stats on
  `sim4500`.

## Methods note

None new (collects the sign-offs above).

## Questions for the lead

1. Gate on the nightly full suite (R = 300 / 1,000) or accept the per-PR fast suite for the M4 gate?
   **Proposed default:** the gate needs one green nightly full run on `main` after M4-17 and M4-12 merge.
