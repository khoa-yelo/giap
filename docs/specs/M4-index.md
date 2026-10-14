# M4 "Gene table, detectability, region tests": issue index

Spec: GIAP_SPEC.md §13 row M4. Deliverables: "gene table, detection tiers, detectability, abundance,
basic statistics, region tests". **Gate:** "circular-shift tests calibrated; gene-level FDR calibrated
(picked, MAYU)" — i.e. the `regions` suite (R = 100 per PR / 1,000 nightly) passes the §8.4 level
check for every region type and statistic, and the `gene_fdr` suite (R = 50 / 300) passes the FDR
check for the `decoy_picked` and `mayu` tiers (`tests/calibration/tolerances.yaml`).

Sources: §7.3 (Module 3), §8 (claims, rules, shared engine `src/giap/stats/`, calibration, attack
nulls, sensitivity), §9.3/§9.4 (genes.tsv, detection_summary, regions_stats, flags), §11.3 (simulator
expression and proteomics), Appendix B.1, B.3, B.4; decisions 2, 3, 4, 13, 14; G5, G6, G7, G9.

One issue = one PR; size counts **code** only (ADR-0007: ≤ ~400 lines under `src/`, `workflow/`,
`scripts/`; simulator lines under `tests/sim/` are tooling and listed separately).

| Order | Issue | Track | Scope | Depends on | Code lines | Before M2/M3? |
|---|---|---|---|---|---|---|
| 1 | [M4-01](M4-01-stats-bh-tda.md) | S | `stats/`: `bh`, `tda_qvalues`, `picked_fdr`, `monotone_q` | none | 180 | yes |
| 2 | [M4-02](M4-02-stats-mayu.md) | S | `stats/`: `mayu_gene_fdr` (+ ≥ 2-pep variant); `poisson_binomial_*` is M3-09's | 01 | 140 | yes |
| 3 | [M4-03](M4-03-stats-permutation-shifts.md) | S | `stats/`: `perm_p`, `circular_shifts`, `choose_shifts`, `besag_clifford` | none | 220 | yes |
| 4 | [M4-04](M4-04-stats-bootstrap-hl.md) | S | `stats/`: `tu_block_bootstrap`, `circular_block_bootstrap`, `hodges_lehmann` | none | 260 | yes |
| 5 | [M4-05](M4-05-calibration-harness.md) | S | `tests/calibration/_criteria.py` (test-author only) + engine calibration tests | 01, 03 (tests) | 0 (+150 test infra) | yes |
| 6 | [M4-06](M4-06-input-contract.md) | C | **input contract** C1–C14 (M1/M2/M3 → M4), validators, `load_m4_inputs` | none | 330 | yes |
| 7 | [M4-07](M4-07-sim-layout-tables.md) | F | giapsim sequence-free layout, regions, operon pairs; contract-table writer | 06 | (380 sim) | yes |
| 8 | [M4-08](M4-08-sim-expression.md) | F | giapsim expression (AR(1), region/plasmid/island effects) and detection truth | 07 | (300 sim) | yes |
| 9 | [M4-09](M4-09-sim-evidence.md) | F | giapsim peptides, false IDs ∝ size, decoys, q/PEP, C8–C11 tables | 08 | (380 sim) | yes |
| 10 | [M4-10](M4-10-protein-covariates.md) | G | covariates: observable peptides, GRAVY, pI, MW, TM segments (G6) | none | 220 | yes |
| 11 | [M4-11](M4-11-accepted-evidence.md) | G | accepted evidence (§7.3.1), per-gene counts, coverage, evidence flags | 06 | 300 | yes |
| 12 | [M4-12](M4-12-detection-tiers.md) | G | tiers: picked, MAYU, two-peptide, protein list; protein→gene map; **gene-FDR calibration** | 01, 02, 11 (cal: 09) | 380 | yes |
| 13 | [M4-13](M4-13-detectability.md) | D | cross-fitted spline-logistic `p_detect` (G5), diagnostics; slope calibration | 06, 10 (cal: 08) | 300 | yes |
| 14 | [M4-14](M4-14-abundance.md) | A | iBAQ/Top3, protein-level route, log2, missing policies | 06, 10, 11, 12 | 320 | yes |
| 15 | [M4-15](M4-15-gene-table-assembly.md) | T | `genes.tsv` assembly (context from M2 `layers.parquet`), function pivot, writer, schema, steps | 06, 10–14 | 340 | yes |
| 16 | [M4-16](M4-16-detection-summary.md) | R | basic O − E summary, QC rows, `detection_summary.tsv` writer + schema | 04 | 330 | yes |
| 17 | [M4-17](M4-17-region-tests.md) | R | region O − E and HL under circular shifts, CIs, BH; **region calibration**, attacks 2/4 | 01, 03, 04 (cal: 05, 07, 08) | 380 | yes |
| 18 | [M4-18](M4-18-regions-stats-table.md) | R | plasmid replicon-level rows, prophage subsets, depth offset, `regions_stats.tsv`, `context_stats` step | 15, 16, 17 | 330 | yes |
| 19 | [M4-19](M4-19-sensitivity-regions.md) | R | sensitivity framework + `regions` family, `sensitivity.tsv` (off the gate path) | 12, 13, 14, 18 | 250 | yes |
| 20 | [M4-20](M4-20-integration-golden.md) | T | end-to-end on simulated fixtures, golden `sim50k`, TE7, gate evidence | 01–19 | 120 | yes |
| 21 | [M4-21](M4-21-m2-m3-conformance.md) | T | contract conformance on real M2/M3 outputs; golden through the real path | 20 + **M2 + M3 merged** | 150 | **no** |

Total ≈ 4,550 code lines in `src/` + ≈ 1,060 simulator lines + ≈ 6,000 test lines.

## Dependency DAG

```text
Track S (engine)     01 ─► 02            03            04            05 (needs 01, 03 for its tests)
Track C (contract)   06
Track F (simulator)  06 ─► 07 ─► 08 ─► 09
Track G (evidence)   10                  06 ─► 11 ─► 12 ◄─ 01, 02         (12 calibration ◄─ 09)
Track D (detect.)    06, 10 ─► 13                                          (13 calibration ◄─ 08)
Track A (abundance)  06, 10, 11, 12 ─► 14
Track R (context)    04 ─► 16            01, 03, 04 ─► 17                  (17 calibration ◄─ 05, 07, 08)
Track T (assembly)   06, 10–14 ─► 15
Joins                15, 16, 17 ─► 18 ─► 19 (◄─ 12, 13, 14) ─► 20 ─► 21 ◄─ M2 + M3 merged
```
Calibration tests reference the simulator: 12 ← 09, 13 ← 08, 17 ← 07 + 08 (+ 05 harness). Unit tests
of 11–19 use small hand-built frames, so code can start before the simulator issues land; only the
calibration test files wait for them.

## Parallel tracks

| Track | Issues | Owner pattern | Starts | Notes |
|---|---|---|---|---|
| **S — stats engine** | 01 → 02; 03; 04; 05 | one worktree; 01, 03, 04 are mutually independent (3 parallel PRs), 02 after 01, 05 by the test-author in parallel | immediately | **Highest priority:** 01 should land before M3-10 so placement FDR reuses `tda_qvalues` (M4-01 Q3); M5–M8 reuse 03/04. `poisson_binomial_*` is built by M3-09. |
| **C — input contract** | 06 | one PR, reviewed by the M2 and M3 leads | immediately | Freezes what M2/M3 write; send its table (below) to the M2/M3 spec-writers now. |
| **F — simulator** | 07 → 08 → 09 | implementer (tooling), test-author tests it | after 06 (07 compares column lists with `TABLES`) | Feeds every calibration suite; sequential within the track. |
| **G — evidence & detection** | 10; 11 → 12 | 10 is independent; 11 after 06; 12 after 01, 02, 11 | 10 immediately; 11 after 06 | 12 carries the gene-FDR gate. |
| **D — detectability** | 13 | after 06, 10 | early | slope calibration needs 08. |
| **A — abundance** | 14 | after 06, 10, 11, 12 interfaces | mid | can stub `Detection.table` to start earlier. |
| **R — context statistics** | 16; 17 → 18 → 19 | 16 and 17 need only the engine (04; 01, 03, 04) | as soon as S lands | 17 carries the circular-shift gate; can be developed in parallel with all of G/D/A. |
| **T — assembly & integration** | 15 → 20 → 21 | joins G, D, A, R | last | 21 is the only issue that waits for M2 and M3. |

Suggested waves (≤ 6 concurrent subagents, PROGRESS.md):
- **Wave A (now, parallel):** 01, 03, 04, 06, 10, 05 (test-author).
- **Wave B:** 02, 07, 11, 16, 17 (code + unit tests), 13 (code + unit tests).
- **Wave C:** 08, 09, then the calibration tests of 12, 13, 17; 12, 14.
- **Wave D:** 15, 18, 19, 20. **After M2 + M3 merge:** 21.

**Before vs after M2/M3 merge.** Issues 01–20 are developed and gated entirely on simulator fixtures
written in the M4-06 contract format. Only **M4-21** needs the merged M2 and M3 code. Risk: if M2/M3
deviate from the contract, M4-21 absorbs the fix; keep the contract table below in the M2/M3 specs.

## The input contract in one table (details, invariants and alignment notes in M4-06)

Names aligned with the M2/M3 drafts now in `docs/specs/` (M2-04, M2-11, M2-12, M2-19, M2-20, M3-16).
Parquet = internal 0-based coordinates; M3's final TSVs are converted on read.

| ID | File | Producer | Key columns M4 relies on |
|---|---|---|---|
| C1 | `replicons.parquet` | M2-04 `parse_bakta` | `contig`, `length`, `topology`, `kind` (row order = FASTA order) |
| C2 | `genes.parquet` | M2-04 | `locus_tag`, `feature_type`, `contig`, `start`, `end`, `strand`, `partial`, `crosses_origin`, `hypothetical`, `aa`, `nt`, `start_codon` |
| C3 | `xrefs.parquet` | M2-04 | `locus_tag`, `db`, `id` → `cog_id`, `go_terms`, `kegg_ko`, `ec` |
| C4 | `annotation_meta.json` | M2-04 | `annotation_source` |
| C5 | `tus.parquet` | M2-12 `operon_pairs` | `tu_id`, `genes` (transcription order); **every gene in exactly one TU** |
| C6 | `layers.parquet` | M2-20 `context_layers` | `tu_id`, `tu_position`, `tu_size`, `uniop_p_prev/next`, `region_ids`, regulatory/codon columns, `flags` |
| C7 | `regions.parquet` | M2-19 `region_layer` | `region_id`, `region_type`, coordinates (AMR NA) |
| C8 | `region_genes.parquet` | M2-19 | `region_id`, `locus_tag`, `overlap_frac`, `subset`, `core_biosynthetic` |
| C9 | `peptides.norm.parquet` | M1-12 | M1-06 `AGG_COLUMNS` (incl. decoys, intensity) |
| C10 | `abundance.norm.parquet` (opt) | M1-12 | `protein_id`, `value`, `scale` |
| C11 | `mapping/peptide_hits.tsv` | M3-16 `loci` | `placement_q`, `locus_class`, `uniqueness`, `route`, `gene`, `flags`, locus |
| C12 | `mapping/protein_hits.tsv` | M3-16 | `protein_id`, `rank`, `class`, `gene`, `flags` |
| C13 | `mapping/unmapped.tsv` | M3-16 | QC only |
| C14 | `annotation/id_map.tsv` (opt) | M2-11 | `source_id` |

Open alignment items: pseudogene assignment (M3-14 Q3 says no, M4-06 Q3 says yes);
`gc_local_dev` and `cog_category` have no producer yet (M4-06 Q6).

## Calibration suites owned by M4 (replicates from `tolerances.yaml`)

| Suite (key) | Issue | Generator | Criterion | R fast / full |
|---|---|---|---|---|
| `gene_fdr` | M4-12 (+ BH in M4-05) | M4-07/08/09, false peptides ∝ n_observable | FDR at q ∈ {0.01, 0.05, 0.10}, picked and MAYU (± PEP) | 50 / 300 |
| `detectability` | M4-13 | M4-08 logistic truth on covariates | mean calibration slope ∈ [0.8, 1.2] | 20 / 100 |
| `regions` | M4-17 (+ perm/shift/BC in M4-05) | M4-08 AR(1) abundance, detection from covariates only | level at α ∈ {0.01, 0.05, 0.10} per type × stat | 100 / 1,000 |
| attacks 2 and 4 | M4-17/18 (statistician) | short genes in prophages; 20-copy plasmid | level | as `regions` |

## Methods notes (statistician sign-off, §14.1, §14.9)

`docs/methods/`: `multiple-testing.md` (01), `target-decoy.md` (01), `mayu-gene-fdr.md` (02),
`permutation-and-shifts.md` (03), `bootstrap.md` (04), `calibration.md` (05), `detection-tiers.md`
(11, 12), `detectability.md` (10, 13), `abundance.md` (14), `detection-summary.md` (16),
`region-tests.md` (17, 18), `sensitivity.md` (19).

## Out of M4 (named so nobody builds them twice)

`morans_i`, `join_count`, `scan_max`, `blom_scores` (M5, autocorrelation and islands; `island_id`,
`replichore*` stay NA until then); `wallenius_ora`, `fgsea_multilevel` (M8); operon pair tests (M6).

## Decisions taken in these specs (for the lead to confirm)

- M4 reads M2's Parquet intermediates and M3's final `mapping/*.tsv` (C1–C14); every gene has a TU
  (M4-06 Q1); pseudogene placements assign to the pseudogene (M4-06 Q3, conflicts with M3-14 Q3).
- Picked FDR uses a rank-transformed lexicographic key (q, PEP, score) instead of −log₁₀ q (M4-12).
- Decoys pass the same engine q threshold as targets before picking (M4-12 Q1).
- A peptide-vote protein→gene map (not in the spec) handles cross-strain accessions without
  `proteins.faa` (M4-12 Q4).
- Region shifts index **eligible** genes; BH family = (region type, statistic) (M4-17 Q1–Q2).
- Engine calibration reuses existing `tolerances.yaml` keys, no protected-file change (M4-05 Q1).
- Proposed schema vocabulary additions needing approval: `few_shifts` (M4-17 Q3),
  `no_reference_replicon` (M4-18 Q2).
- The sensitivity framework is built in M4 (regions family only), off the gate path (M4-19 Q1).
- Steps register with M2-02's `@step` registry if merged, else in M4-20 (M4-15 Q1); schemas go to
  `src/giap/schemas/tables/` like M3-16 (M4-15 Q2).
- `poisson_binomial_*` is owned by M3-09 (`stats/poisson_binomial.py`); M3-10 should reuse `tda_qvalues` (M4-01 Q3);
  M3-19 and M4-09 share one score/PEP/q helper (M4-09 Q3).

## Lead decisions on the questions (2026-10-14)
1. Pseudogene peptides: no conflict — M3-14 Q3 was accepted as "placements on pseudogenes are assigned
   to the pseudogene", so M4-06 Q3 stands and pseudogenes stay eligible.
2. `cog_category`: from Bakta's COG cross-references mapped through the COG 2020 category table that
   `giap setup` downloads; NA when the table is not installed. `gc_local_dev`: produced by M2-20's
   layers table (added to the M2 index).
3. Shared code: M3-10 uses `giap.stats` `tda_qvalues` (M4-01) — M4-01 is therefore a foundation issue
   built early; M4-02 reuses M3-09's Poisson-binomial functions; giapsim gets one shared
   score/PEP/q helper used by M3-19 and M4-09.
4. Decoys pass the same engine q threshold as targets before picking. The decoy-pairing cut-off is NOT a
   code constant: it reuses `mapping.proteins.min_identity` / `min_coverage` (decision 6, 0.80).
5. Region tests: shifts index eligible genes only; BH within each (region type, statistic). New flags
   `few_shifts` and `no_reference_replicon` are added through the schema (§9.4). Accepted.
6. Calibration reuses the existing `gene_fdr` and `regions` replicate counts; tolerances.yaml unchanged.
   Detectability gate uses the mean calibration slope. Accepted.
7. M4 gate: one green full calibration run on main, submitted as a Slurm job until Khoa installs the
   nightly scrontab entry; sensitivity built but off the gate path. Accepted.

8. (Khoa, 2026-10-14; ADR-0009) Picked gene-level FDR rescales decoy wins by the paired share
   (1/f, f = paired decoys / all decoys); the 20 % MAYU fallback stays as an upper guard. M4-12 adds
   calibration under 10/20/30 % random decoy loss.
9. From the M4-01 statistician review: M3-10's π̂ for wrong placements needs its own calibration in the
   placement suite (an underestimated π inflates FDR).
