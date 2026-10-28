# MVP-01: wire the registered steps into `gippa run` (gene table, detectability, context)

Track: **MVP / end-to-end** (Khoa, 2026-10-20: "achieve end to end integration of the tool first")
· depends on: nothing (every step is merged and registered)

## The gap
31 steps are registered, but `workflow/main.nf` calls only `ANNOTATE` and `MAPPING`. Everything M4 and
M5 built — `gene_table`, `detectability`, `context_stats`, `run_sensitivity`, and soon the
autocorrelation and island tables — is reachable from Python only. `gippa run` stops after mapping.
This issue closes that, and supersedes the three backlog wiring items (M4-15, M4-18, M4-19 decision 16).

## Goal
`gippa run` on the simulated 100 kb genome produces, in one invocation: `genetable/genes.parquet`,
the detectability summary, `tables/regions_stats.tsv` and `tables/sensitivity.tsv`.

## How this issue is sequenced (Khoa, 2026-10-20)
"After we integrate end to end, we can pick up the rigour… write the spec as usual, but prioritize end
to end integration first, then work on the rigour after." So this spec carries the whole job, and the
work is split into two passes:
- **Pass 1 (this issue, now):** wiring, a real run that reaches the end, nf-test stub coverage.
- **Pass 2 (MVP-01H, after the MVP runs):** the hardening the spec names below — contract assertions on
  the wired outputs, a golden end-to-end comparison, failure-path behaviour (a step erroring mid-graph,
  `-resume` correctness), and the per-step runtime budgets. Nothing here is dropped; it is scheduled.

## Hardening owed (pass 2, written down now so it is not lost)
- Contract validation of every wired output inside the pipeline, not only in the Python tests.
- A golden end-to-end fixture for sim100k so a wiring regression is caught by diff, not by eye.
- `-resume` correctness across the new subworkflows, and the behaviour when one step fails mid-graph.
- Per-step runtime recorded via scripts/measure_runtimes.py and a budget in benchmarks/baseline.json.
- The real-tool arm: the same run under `-profile farmshare,apptainer` rather than the local executor.

## Lead decisions
1. Two new subworkflows beside annotate.nf and mapping.nf: `genetable.nf` (`gene_table`,
   `detectability`) and `context.nf` (`context_stats`, plus sensitivity). Called from `main.nf` after
   MAPPING, each skipped cleanly when its inputs are absent, as MAPPING already is.
2. `run_sensitivity` is a plain function today (M4-19 decision 16): register it as a step so the
   subworkflow can call it like the others.
3. `cog2020_def` resolves from the setup manifest (M4-15 decision 5). `stratum0_classes` from #124 is an
   optional input to `context_stats`; wire it when #124 has merged, otherwise leave the channel empty.
4. **Scope discipline:** wiring only. Do not change step behaviour, do not add calibration, and if a step
   needs a small fix to run under Nextflow, make the smallest one and say so in the PR.
5. Evidence that counts: a real `sbatch scripts/slurm/gippa_run.sbatch` on sim100k that reaches the end
   and writes the four outputs. nf-test stub coverage for the new subworkflows, in the existing style.

## Definition of done
The real run above completes; the four outputs exist and are non-empty; `nf-test test workflow/tests`
passes; the three backlog wiring entries are deleted.
