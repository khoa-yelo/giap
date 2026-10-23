# M2-14b: batch Promoter Calculator tasks

Milestone M2 (follow-up; Khoa 2026-10-23, pulled forward from the post-M9 performance pass) ·
depends on: M2-14 · estimate ≈ 40 code lines (+ ≈ 120 test lines)

## Why
`PROMOTER_CALCULATOR` runs one Nextflow task per TU upstream region. On FarmShare each task is its
own Slurm job, and `executor.submitRateLimit = '10/1min'` (farmshare.config) paces them. A K-12 run
spends about 5.5 h submitting about 2,600 jobs that each compute for seconds (Slurm 1873779, 1878443).

## What changes
1. **Batching.** annotate.nf sorts the `tu_upstream` regions by TU id and groups them into
   consecutive batches of at most `PC_BATCH = 100` regions (`collate(100, remainder = true)` over the
   sorted list). Each batch is one `PROMOTER_CALCULATOR` task. A K-12 run then has about 26 tasks.
2. **Per-region command unchanged.** Inside a task, each region runs exactly as before, one after
   another in batch order: `promoter-calculator -i <region> -o <TU>.csv -j 1`, plus the same
   "no promoter → header-only CSV" fallback. The process's output becomes `path "*.csv"` (one CSV
   per TU, the same file names), so `TX_LAYER`'s input is unchanged.
3. **Outputs byte-identical.** Every `<TU>.csv`, and therefore `tx_rates.parquet`, is byte-identical to the
   unbatched pipeline.
4. **Determinism and resume.** Batch membership depends only on the sorted TU ids, so `-resume`
   caches a batch whenever its regions are unchanged. The task tag is `pc_<first TU>_<last TU>`
   (or `pc_<n>` if that is simpler; tests pin whichever is chosen).
5. **Failure.** A failing region fails its batch. The existing 137/143 retry applies per batch.
   Every other exit terminates, as before.
6. **Resources.** The label stays `process_single`, and the time limit is
   `{ 2.h * task.attempt }`.
7. **Not a config key.** `PC_BATCH` is a constant in annotate.nf (defaults.yaml is protected; a key
   needs an ADR). Recorded in the backlog: make it configurable if a genome needs a different size.
8. **Stub.** The stub writes the header-only CSV for each region in the batch.

## Examples
| # | Input | Expected |
|---|---|---|
| 1 | 250 regions | 3 tasks of 100, 100 and 50, and 250 CSVs named `<TU>.csv` |
| 2 | 100 regions | exactly 1 task |
| 3 | 0 regions | no task. TX_LAYER behaves as with zero regions today |
| 4 | the same input run batched and unbatched | byte-identical `tx_rates.parquet` and per-TU CSVs |
| 5 | rerun with `-resume` after one region changes | only that region's batch reruns |
| 6 | one region whose tool exits 1 | its batch fails; the run terminates with the tool's error |

## Tests
- **nf-test:** stub runs for examples 1–3 and 6, and the `-resume` behaviour of example 5 with stub tasks.
- **Real:** the recorded-input real test (workflow/tests/real/promoter_calculator.nf) runs batched and compares the CSVs.
- **Static:** a test that annotate.nf collates sorted regions in batches of 100.
- **Byte identity:** example 4 on sim100k (288 TUs → 3 tasks), integration or slow.

## Acceptance
- Examples 1–6 hold.
- The fast suite and nf-test pass.
- A K-12 real run (MVP-04 or M6-12C) records the PROMOTER_CALCULATOR wall time in PROGRESS.md.
