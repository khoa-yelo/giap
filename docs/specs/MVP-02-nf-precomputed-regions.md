# MVP-02: Nextflow parity for a precomputed regions layer (region_layer and CONTEXT run)

Track: **MVP / end-to-end** · depends on: **M4-21 (#145)**, whose `region_layer` builds C7/C8 from the
precomputed TSV · estimate ≈ 60 workflow lines (+ ≈ 80 nf-test lines) · not statistical · no tools

## The gap
After M4-21, `giap run --executor local` with `annotation.precomputed.regions` set plans `region_layer`
(it converts the M2-02 TSV into C7 `regions` / C8 `region_genes`) and then `context_stats` and
`sensitivity`. The Nextflow path still skips `REGION_LAYER` and therefore `CONTEXT`
(`workflow/subworkflows/local/annotate.nf` around line 178; `workflow/main.nf` line 6). The two
executors disagree on the same input, against §6.4 (one plan, two executors).

## Goal
With a precomputed regions layer, Nextflow:
- still skips `GENOMAD` and `ANTISMASH` (the layer replaces them);
- runs `REGION_LAYER` on the precomputed TSV (the step reads `annotation.precomputed.regions` from
  the config, as in M4-21) and emits the `c7`/`c8` channels;
- runs `CONTEXT`, so `tables/regions_stats.tsv`, `detection_summary.tsv`, `autocorrelation.tsv`,
  `islands.tsv` and `genes.tsv` are published.

## Examples
| # | Input | Expected |
|---|---|---|
| 1 | `annotate.nf.test` stub, precomputed regions | GENOMAD and ANTISMASH absent; REGION_LAYER succeeds; c7/c8 emitted |
| 2 | `main_genetable_context.nf.test` stub, precomputed regions layer | GENETABLE **and CONTEXT** run |
| 3 | no precomputed regions (unchanged path) | identical process list to today |
| 4 | local vs Nextflow on the M4-21 SIM100K real-path config | the same set of published tables |

## Acceptance criteria
- [ ] The two nf-tests that pin the old skip are updated by the test-author (superseded by this spec):
  `annotate.nf.test` "precomputed regions: no GENOMAD, no ANTISMASH, no region_layer" and
  `main_genetable_context.nf.test` "precomputed regions layer: GENETABLE still runs, CONTEXT does not".
- [ ] nf-test stubs green in CI (`nf-test-stub`).
- [ ] Example 4 as an integration test if the existing local-vs-Nextflow harness can run it in the
  stub profile; otherwise one Slurm `giap_run.sbatch` run, job ID in PROGRESS.md.

## Hardening owed (pass 2)
- The precomputed TSV carries no marker descriptions or core flags, so prophage subsets come from the
  Bakta product only and `core_biosynthetic` is always false (M4-21). Record as a run flag or extend
  the M2-02 TSV format; decide when M2-02's format is revisited.
