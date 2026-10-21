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

## Lead decisions after the test-author's report (2026-10-21)
1. **Golden alpha.** The sim50k golden runs at `stats.alpha.identification = 0.05` (33 of 47 genes
   detected, real p-values), documented in the module docstrings. The vacuous default stays pinned by
   `test_default_alpha_detects_nothing_on_50_kb`, which also asserts `len(genes) < 100` so the pin shouts
   if the fixture grows. A golden of zeros would test nothing; moving to sim1000 is rejected because
   §11.4 names the 50-kb fixture. Record in the methods/README text that the 50-kb fixture cannot reach
   the default alpha (picked q floor 1/29 = 0.0345) — this is a property of the fixture, not a defect.
2. **`tables/genes.tsv` (gap 1): write it at the end of the M4 chain.** `write_genes_tsv` has no caller
   anywhere, so the pipeline's headline gene table is never produced — that is an end-to-end hole and the
   MVP push is the right time to close it. §6.2 gives genes.tsv to `finalize`; when `finalize` arrives it
   takes the call over. Note that in the backlog beside the `genes_pre.parquet` publishDir item.
3. **detection_summary tier rows (gap 2): yes.** `gene_table` persists `det.summary` as
   `genetable/detection_summary_parts/tier.parquet`, an input to `context_stats`, so the `fdr_2pep`
   q-threshold curve reaches `tables/detection_summary.tsv`. Today `steps/gene_table.py` discards it.
4. **`build_m4_fixture` takes an extra `layout={...}` kwarg** (LayoutConfig fields) for example 7's
   linear contigs and plasmid; `size` alone cannot express them. Added to the spec's Interface line.
5. **Fixture robustness confirmed:** the builder plants the prophage and AMR regions relative to the
   actual gene count (it varies ~40–50 by seed) and gives AMR rows NA coordinates
   (`contract.regions.amr_coordinates`). The property tests build at arbitrary seeds and must stay that way.
6. **Gate scope confirmed** as the spec's proposed default: one nightly full run (R ≥ `replicates_full`;
   gene_fdr 300, regions 1000) for `gene_fdr.{picked,mayu_pep,mayu_nopep}.pi{0.01,0.05}` and the regions
   level/FDR streams. Generating the PROGRESS.md evidence block from it is **deferred** to the hardening
   pass — useful, not integration.
7. **`sensitivity = off` stays pinned in M4-19's own tests**, not here; M4-20 keeps its scope.
8. The test-author was right not to write `tests/sim/giapsim/fixtures.py`: giapsim is implementer-owned
   tooling (AGENTS.md), and writing it there would have made the tests pass before the implementation.
   The implementer writes it to the interface `tests/helpers/m4_chain.build` and the property tests pin.
