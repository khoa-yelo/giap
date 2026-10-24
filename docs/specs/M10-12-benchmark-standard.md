# M10-12: the standard-run benchmark (`scripts/slurm/benchmark_standard.sbatch`)

Milestone M10 · **pass 2** · depends on: M10-04 (a real `giap run` on a 5-Mb genome has finished
end to end on FarmShare; soft: the MVP-04 K-12 inputs under `$GIAP_ROOT/mvp04/20261023/`) ·
estimate ≈ 230 code lines + ≈ 260 test lines · not `stats` (a median and a range; no inference) ·
roles: implementer (`scripts/slurm/benchmark_standard.sbatch`, `scripts/dev/benchmark_summary.py`,
the `SIM5M` giapsim scenario), test-author (`tests/unit/test_m10_12_*.py`), performance (runs it and
writes ADR-M10-D), docs (the user-guide runtime table) · Khoa: none (the numbers are reported, not
gated; §14.8 has no item for them).

## Spec sections implemented

- §13.2: "A standard run on 16 threads and 64 GB takes about 35 min (range 16–55; measured by
  `benchmark_standard.sbatch` on FarmShare, node class recorded). Core modules alone take about
  25 min." This issue **measures** both numbers; it does not have to reproduce them.
- §15.2 step 9: "Run `scripts/slurm/benchmark_standard.sbatch` (16 CPU, 64 GB, one pinned node
  class, median of 3); the measured runtimes replace the estimates in the docs."
- §15.6 bullet 4: "Measured runtimes replace the estimates." Index evidence row: "`benchmarks/
  baseline.json` entries with node class, ADR-M10-D, user-guide table updated."
- §3.7 row "16-thread / 64-GB runtime benchmark … `benchmark_standard.sbatch`, pinned to one node
  class"; §3.8 (job header, `sacct` Elapsed and MaxRSS into PROGRESS › Measured numbers, node
  class); §3.1 (partitions, `MaxMemPerCPU`).
- §14.4 speed loop: "wall time and peak memory per step (`sacct` MaxRSS, same node class)" on "the
  5-Mb simulated genome and K-12"; "each accepted improvement updates `benchmarks/baseline.json`
  through an ADR with before-and-after numbers". This issue writes the **first** entry M10-13
  improves on.
- R14 (≤ 128 CPUs, ≤ 20 queued or running jobs, names `giap-*`, look only at GIAP's own jobs),
  R3 (no tolerance is tuned on these numbers), R9 (3 attempts), R11 (nothing large committed).
- Index lead decision 1, L7 and L8: "the standard and core benchmark variants on wheat"; "runtimes
  go into baseline.json, the guide, PROGRESS and ADR-M10-D".

## What exists (main 1f567901)

| Piece | State | Used for |
|---|---|---|
| `scripts/slurm/weekly.sbatch` | 16 CPU, `--mem=60G`, `farmshare_local,apptainer` (the `farmshare` Slurm-executor profile is refused with exit 2, decisions 12, 22), a detached worktree per ref under `$GIAP_ROOT/ci/` with `flock`, a fresh node-local `NUMBA_CACHE_DIR` per task, `PYTHONPATH=$PWD/src:$PWD/tests/sim`, `scripts/measure_runtimes.py` on each trace | the model for checkout, profile guard and traps; no code is shared (copy the pattern) |
| `workflow/conf/farmshare_local.config` | local executor; `executor.cpus = SLURM_CPUS_PER_TASK`, `executor.memory = SLURM_MEM_PER_NODE − 2048 MB`; `resourceLimits` = min(allocation, `--threads`/`--memory`) | **the profile of the benchmark**: every task inside the one 16-CPU allocation |
| `workflow/conf/farmshare.config` | Slurm executor, `queueSize = 8`, per-task CPUs ≤ min(16, `--threads`) | refused here: with it, `--threads 16` caps each task, not the run (8 × 16 = 128 CPUs on other nodes) |
| `giap run` | `--threads`, `--memory` (default `64.GB`), `--modules`, `--annotation`, `--profile`, `--seed`, `--resume` (off by default), `-work-dir OUT/work`, `-with-trace OUT/logs/trace.tsv` | the measured command |
| `cli._bakta_json` | `--annotation` is a Bakta JSON only when its suffix is `.json`. **A `.json.gz` is taken for GFF3/GenBank and Bakta runs live** (M10-04 Q3; backlog chore "accept `.json.gz`") | the core variant decompresses the JSON first and asserts that no Bakta task ran |
| `run_manifest.json` | `git_commit`, `git_dirty`, `resume`, `profile`, `started_at`, `finished_at`, `exit_code`, `steps[]` with `status` (`completed`/`cached`/`failed`), `wall_time_s`, `peak_rss_bytes` | commit, resume and cache checks; per-step numbers |
| `scripts/measure_runtimes.py` | trace → markdown rows `step \| genome \| wall \| peak RSS \| node class (job)`; `CACHED`/`FAILED` left out | the per-step rows for PROGRESS, unchanged |
| `benchmarks/baseline.json` | 31 top-level entries of micro-benchmarks; `annotation_steps` is keyed by node class (`wheat`), each entry with `median_s`, `node`, `cpu`, `slurm_job`, `commit`, `date`, `rounds`. Not a protected file | gets a new top-level key `standard_run` in the same per-class shape |
| Inputs on disk | `tests/data/genomes/k12_mg1655.fna.gz` (U00096.3, `[topology=circular]`, committed); `tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz` (committed); `$GIAP_ROOT/mvp04/20261023/{peptides.tsv,config.yaml}` (giapsim real-genome peptides at `K12_EVIDENCE_SEED`, written by `mvp04_inputs.sbatch`); Bakta full DB 6.0, geNOMAD, Promotech, antiSMASH image, GO/KEGG/COG, contaminants under `$GIAP_DB_DIR` | inputs |
| giapsim | `simulate(SimConfig, seed)` writes `genome.fasta`, `bakta_like.json` (Bakta-shaped when a v2 option is on), `peptides.tsv` with `proteomics`; a 5-Mb genome in < 30 s (`benchmarks/test_giapsim.py`). **No named 5-Mb scenario exists** | `SIM5M` (new, below) |
| User guide | **does not exist**: no MkDocs site yet (M11), no runtime table anywhere under `docs/` | question L9 |
| Measured so far | Bakta full DB on K-12: 9.9 min, 20.6 GB (barley-01, 1849096); TerminatorNet on K-12: 60 min at `-n 17` (wheat-04, 1847778), ~3 h at 6 threads (M7-15); weekly annotate run with live Bakta: 2 h 08 min (wheat-05, 1852010) | the §13.2 range (16–55 min) is probably too low for the standard variant; this issue reports what it measures |

## Cluster facts that fix the job shape (checked 2026-10-23)

1. **Memory per CPU.** `normal` has `DefMemPerCPU=2400`, `MaxMemPerCPU=4000` (MB, Slurm's `M` is
   MiB) and `SelectTypeParameters=CR_CORE_MEMORY`. `--cpus-per-task=16 --mem=64G` asks for
   65,536 MB = 4,096 MB per CPU, so Slurm raises the request to **17 CPUs and allocates 18** (whole
   cores, `ThreadsPerCore=2`). Observed: jobs 1846846, 1848590, 1848630, 1852010 (`ReqCPUS 17`,
   `AllocCPUS 18`); jobs at `--mem=48G` or `60G` got exactly 16. The benchmark therefore asks
   **`--cpus-per-task=16 --mem=64000M`** (4,000 MB per CPU, the maximum): 62.5 GiB, which is
   67.1 GB in SI units, so it satisfies "64 GB". 16 CPUs at 4 GB each is the most memory 16 CPUs can
   have (question L1).
2. **What "16 threads" is.** On wheat (Intel Xeon Gold 6148, 2 sockets × 20 cores × 2 hardware
   threads = 80 CPUs, 191 GB) a Slurm CPU is a hardware thread: 16 CPUs = 8 physical cores. The
   report says "16 hardware threads (8 cores)".
3. **`--threads` does not cap the run's CPUs by itself.** Under `farmshare` it caps each Slurm task
   (up to 8 at once, on other nodes); under `farmshare_local` the cap on the sum is
   `executor.cpus = SLURM_CPUS_PER_TASK`, and the allocation's cgroup (`TaskPlugin=task/cgroup,
   task/affinity`, `ProctrackType=proctrack/cgroup`) confines every process, including tools that
   ignore their thread flag, to the 16 allocated CPUs. So the profile is load-bearing: only
   `farmshare_local,apptainer` measures "a run on 16 threads".
4. **Memory inside the allocation.** The cgroup holds the Nextflow JVM (`NXF_OPTS=-Xms512m -Xmx4g`,
   so up to ~5 GB resident), the `giap` Python head and the tasks. `farmshare_local` sets the task
   pool to 64,000 − 2,048 = 61,952 MB, which leaves too little for the JVM. The run therefore passes
   **`--memory 56.GB`** (57,344 MB of tasks + ~5 GB JVM + head < 64,000 MB), as `weekly.sbatch` does
   (question L2).
5. **Accounting.** `JobAcctGatherType=jobacct_gather/cgroup`, `JobAcctGatherFrequency=60`: `sacct`
   MaxRSS is sampled every 60 s and can miss a short peak; Nextflow's per-task `peak_rss` is the
   second, independent number. Each run is its own `srun` job step, so `sacct` gives Elapsed,
   MaxRSS and TotalCPU **per run**, not per job.
6. **Node class.** `--constraint=CPU_SKU:6148` selects wheat-01..08 (all `CPU_SKU:6148`; wheat-03
   was `drained` on 2026-10-23). Barley (`CPU_SKU:9384X`), oat (`6426Y`, GPU nodes) and rye
   (`9754`) are other classes. `OverSubscribe=NO`: allocated CPUs are never shared, but memory
   bandwidth, L3 cache, the node's NFS client and page cache are.

## Definition

**A benchmark** is one `BENCH_COMMIT` measured as a **matrix of 2 genomes × 2 variants × 3
repeats = 12 runs**, all on one node class.

| Genome key | Genome | Peptides | Bakta JSON (core variant) |
|---|---|---|---|
| `k12` | `tests/data/genomes/k12_mg1655.fna.gz` (4.64 Mb, one circular replicon) | `$GIAP_ROOT/mvp04/20261023/peptides.tsv` with its `config.yaml` (as MVP-04) | the committed `k12_mg1655.json.gz`, decompressed into the run's input directory |
| `sim5m` | giapsim `SIM5M` at `SIM5M_SEED`, written by the job (below) | the scenario's `peptides.tsv` | the scenario's `bakta_like.json` |

`SIM5M` (new in `tests/sim/giapsim/simulate.py`, next to `SIM100K`): `SIM100K`'s options with
`genome_length=5_000_000`, `n_plasmids=1`, `plasmid_length=(80_000, 80_000)`,
`n_rrna_operons=7`, `n_trna=86`, `prophage_bp=60_000`, `bgc_bp=40_000`, `n_pseudogenes=40`,
`n_sorfs=30`, `write_precomputed=False` (the real tools must run), `proteomics=ProteomicsConfig()`;
`SIM5M_SEED = 20261014`. The job writes it into its own input directory before any timing; the
summary checks that the three repeats' input sha256 values are equal.

| Variant | Modules | Bakta | Extra flags |
|---|---|---|---|
| `standard` | the default `core,mobile,bgc,function` (no `--modules` flag) | live, full DB 6.0 (`--set annotation.bakta.db=full`, no `--annotation`) | none |
| `core` | `--modules core` | not run: `--annotation <decompressed JSON>` | none |

(Index L7 as decided. The two variants differ in two things, modules and Bakta, so "core" is not
"standard minus modules"; the Bakta step's own median is reported so that both readings can be
derived, question L7.)

**Every run's command** (from a fresh launch directory, the pinned worktree on `PYTHONPATH`):

```text
srun --ntasks=1 --cpus-per-task=16 --mem=64000M --export=ALL \
  timeout --signal=TERM --kill-after=5m "$BENCH_RUN_TIMEOUT" \
  giap run --genome G --peptides P [--config C] [--annotation J --modules core] \
    [--set annotation.bakta.db=full] --db-dir "$GIAP_DB_DIR" --threads 16 --memory 56.GB \
    --profile farmshare_local,apptainer --seed 20261014 --outdir RUN_DIR/out
```

**Timed quantity.** `wall_s` = the run's `srun` step Elapsed from `sacct` (whole seconds), with the
script's own `date +%s.%N` around the step as a cross-check (they must agree within 5 s, else a
warning). It includes preflight, Nextflow start-up, every task and the report; it excludes queue
wait, the checkout, input generation, decompression, and the clean-up after the run.

**The benchmark number** of a (node class, genome, variant) cell is the **median of its 3 valid
repeats' `wall_s`**. Companions: the three values, `spread = (max − min) / median`, `noisy =
spread > 0.15` (the `benchmark_slowdown_max` of `ci/gates.yaml`, read, not changed: a cell noisier
than the gate's tolerance cannot back a later slowdown check), `peak_rss_sacct_max` (max over the
3 repeats of step MaxRSS), `task_peak_rss_max` (max over repeats and tasks of the trace's
`peak_rss`), `cpu_efficiency_median` (median of TotalCPU / (Elapsed × 16)), and per step (the
trace process name as `measure_runtimes.py` names it) the median over repeats of the step's summed
task `realtime` and the max of its `peak_rss`.

**Comparison with §13.2** (reported, never a gate): `standard` against 2,100 s (range 960–3,300);
`core` against 1,500 s (no range given). `ratio_to_estimate = median / estimate` and
`within_estimate_range` (standard only). The user-facing comparison is the `k12` row; `sim5m` is
the speed-loop baseline (question L8).

## Scheduling, isolation and noise

**One job per repeat, never two at once.** `#SBATCH --array=0-2%1`: array task *r* is repeat
*r* + 1 and runs the four cells one after another in a rotated order, so a warm page cache or a
busy hour does not always favour the same cell:

| Task | Order |
|---|---|
| 0 | k12/standard, k12/core, sim5m/standard, sim5m/core |
| 1 | k12/core, sim5m/standard, sim5m/core, k12/standard |
| 2 | sim5m/standard, sim5m/core, k12/standard, k12/core |

Job header: `--job-name=giap-bench`, `--partition=normal`, `--qos=normal`, `--cpus-per-task=16`,
`--mem=64000M`, `--constraint=CPU_SKU:6148`, `--time=22:00:00` (4 runs × the 5 h
`BENCH_RUN_TIMEOUT` + 2 h set-up), `--output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out`
(an absolute path, as in the other giap sbatch scripts: `#SBATCH` lines expand no variables). R14 cost: 16 CPUs and
at most 3 queued or running jobs (the pending tasks count toward the 20-job cap). The three repeats
land on whichever wheat nodes are free, so node-to-node variation inside the class is part of the
spread, which is what a user sees.

**Not exclusive** (question L4). `--exclusive` would allocate a whole wheat node (80 CPUs, 62 % of
the R14 cap), wait for an idle wheat node (none was idle on 2026-10-23) and measure conditions no
user has. Instead the isolation is guaranteed for GIAP's own work and **recorded** for the rest:

- **Before submission (lead's procedure, in the runbook):** no other `giap-*` job with ≥ 8 CPUs, and
  no `giap-*` Nextflow task jobs, are running or pending; PROGRESS › In progress gets a
  "benchmark window" line with the array ID; no heavy GIAP job is submitted until the array ends.
  The usual R14 sum (running + pending `giap-*` CPUs + 16 ≤ 128) applies.
- **During each run:** a background sampler writes, at the run's start, every 5 min and at its end,
  one row to `RUN_DIR/cotenancy.tsv`: time, the node's aggregate `CPUAlloc`, `CPULoad`, `AllocMem`
  (from `scontrol show node $SLURMD_NODENAME`; node totals only, no other user's job is looked
  at, R14), and this user's other running `giap-*` jobs (`squeue -u $USER -t R -h -o
  '%A|%j|%C|%N'`, the benchmark's own job excluded).
- **Validity:** another `giap-*` job on **the same node** at any sample makes the run invalid
  (`same_node_giap:<jobid>`); one on another node with ≥ 8 CPUs adds `contended_other_nodes`
  (valid; NFS is shared). Non-GIAP load is reported (median and max `CPUAlloc` of the node) and
  never invalidates.

**Nothing resumes, nothing is warm by accident.**

- Every run gets a new `RUN_DIR = $BENCH_DIR/<bench_id>/r<k>/<genome>_<variant>`; the job exits 2
  if it exists. `--resume` appears nowhere in the script (static test). The run's manifest must say
  `"resume": false` and list no step with `status = cached`; the trace must have no `CACHED` row.
- A fresh empty `NUMBA_CACHE_DIR` under the run's node-local `TMPDIR` for each run: numba
  compilation is part of a user's first run (question L6).
- Images and tool environments must already be present (a pull is set-up, not runtime): the file
  lists of `$NXF_APPTAINER_CACHEDIR` and `$NXF_CONDA_CACHEDIR` are taken before and after each run;
  any new entry makes the run invalid (`image_or_env_created`).
- The Bakta full DB and other databases sit on NFS and the node page cache cannot be dropped (no
  root): the rotation spreads this effect, and the report names it as a noise source.
- After the run is measured, `RUN_DIR/out/work` is deleted (outside the timed step); the tables,
  `logs/`, `run_manifest.json` and `report/` stay until M10-14's clean-up.

**Run validity** (first match decides `status`): `failed` (exit ≠ 0, other than the timeout's 124
or 137; and any step whose `sacct` State is `OUT_OF_MEMORY`, even with exit code 137: only a
non-OOM 124 or 137 is a timeout) → `timeout` → `invalid` (any reason below) → `ok`. Reasons for
`invalid`, all recorded: `cached_tasks:<n>` (trace `CACHED` rows), `cached_steps:<n>` (manifest
steps with `status = cached`), `manifest_resume`, `commit_mismatch` (manifest `git_commit` ≠ `BENCH_COMMIT` or
`git_dirty` true), `host_mismatch` (a trace `hostname` ≠ the job's node: a task left the
allocation), `node_class:<sku>` (the node's features lack `CPU_SKU:6148`, only possible under an
override), `bakta_ran_in_core` (a `BAKTA` process in a core run's trace), `same_node_giap:<id>`,
`image_or_env_created`, `alloc_cpus:<n>` (`sacct` AllocCPUS of the step ≠ 16). A run with Nextflow
task retries (the memory × 2 retry) is valid; `n_retried` is reported.

**Replacement.** A failed, timed-out or invalid repeat is replaced by resubmitting one array task
with `BENCH_ONLY=<genome>/<variant>` and `BENCH_ATTEMPT=b` (then `c`): it writes
`r<k>.<attempt>/<genome>_<variant>`. The summary uses, per cell and repeat, the latest valid
attempt and lists the superseded ones. After 3 attempts at one cell, R9 applies (a lead finding,
not a fourth attempt).

## Interface

**`scripts/slurm/benchmark_standard.sbatch`** (bash 3.2-compatible like the other scripts; starts
from the §3.8 header and `weekly.sbatch`'s traps):

| Variable | Default | Meaning |
|---|---|---|
| `BENCH_COMMIT` | **required** | full 40-hex SHA reachable from `origin`; checked out detached in `$GIAP_ROOT/ci/bench-<sha7>` under `flock`. Unset, short or unknown → exit 2. All three tasks use the same value, so `origin/main` moving between tasks changes nothing |
| `BENCH_ID` | `<YYYY-MM-DD of submission>_<sha7>` | the benchmark directory name; the lead exports it so the 3 tasks share it |
| `BENCH_DIR` | `$GIAP_ROOT/bench` | root of all runs |
| `BENCH_GENOMES`, `BENCH_VARIANTS` | `k12 sim5m`, `standard core` | subsets; an empty or unknown value → exit 2 |
| `BENCH_ONLY`, `BENCH_ATTEMPT` | unset, unset | one cell, for a replacement |
| `BENCH_RUN_TIMEOUT` | `5h` | per-run `timeout` |
| `BENCH_PROFILE` | `farmshare_local,apptainer` | any list containing `farmshare` → exit 2 (R14, weekly decisions 12, 22) |
| `BENCH_DRY_RUN` | unset | `1`: the preconditions and the planned `srun … giap run …` lines are printed, nothing runs (the CI tests use it) |

Preconditions before any run (exit 2, nothing written but the log): the profile guard;
`BENCH_COMMIT`; `giap setup --check --db-dir $GIAP_DB_DIR` for Bakta full, geNOMAD, Promotech,
antiSMASH, GO, KEGG, COG, contaminants; the K-12 inputs (the peptides and `config.yaml` exist; the
genome and JSON sha256 values equal the committed files'); a `RUN_DIR` that already exists.

Exit code of a task: 0 when its four runs are `ok`; 1 when any is `failed`, `timeout` or
`invalid` (the remaining runs still run); 2 a precondition; 143 on SIGTERM (trap writes the
unfinished run as `timeout`).

Per run, `RUN_DIR/bench_run.json` (sorted keys, 2-space indent, final newline):

```text
{schema: "giap-bench-run/1", bench_id, commit, genome, variant, repeat, attempt, order_index,
 job_id, array_task, step_id, node, cpu_sku, argv[], inputs: {name: {path, sha256}},
 sacct: {elapsed_s, max_rss_bytes, total_cpu_s, alloc_cpus, state},
 script_wall_s, wall_s, exit_code, status, invalid_reasons[], warnings[],
 trace: {n_tasks, n_cached, n_failed, n_retried, hosts[], task_peak_rss_max_bytes, failed_processes[]},
 steps: {<step>: {realtime_s, peak_rss_bytes, order}},
 cotenancy: {samples, node_cpualloc_median, node_cpualloc_max, giap_other: [{job, name, cpus, node}]}}
```

plus `RUN_DIR/runtimes.md` (`measure_runtimes.py` on the trace, unchanged) and `cotenancy.tsv`.
`trace.failed_processes` lists the step names (as `measure_runtimes.py` names them) of the trace's
`FAILED` rows that no later attempt of the same task replaced; `n_failed` is its length and the
summary names these steps for a repeat without a valid attempt (example 5). `steps.<step>.order`
is the step's first-start position in the trace (the JSON keys are sorted, so the order is kept
here). The job's inputs live outside the repeat directories, in `$BENCH_DIR/<bench_id>/inputs/`:
the decompressed `k12_mg1655.json` (identical bytes for every task) and `sim5m_r<k>[.<attempt>]/`
(SIM5M written by that task).

**`scripts/dev/benchmark_summary.py BENCH_DIR/<bench_id> [--baseline benchmarks/baseline.json
--write]`**: pure functions (`select_attempts(runs)`, `cell_summary(runs, estimate)`,
`baseline_entry(cell)`, `markdown_rows(cells)`) and a thin CLI. Writes `summary.json` and
`summary.md` (the PROGRESS › Measured numbers rows and the user-guide table) next to the runs.
With `--write`, it sets `standard_run.<node class>.<genome>_<variant>` in `baseline.json` for
**complete** cells only, keeps every other key byte-for-byte, and refuses (exit 1) when the
existing entry has a different commit and no `--replace` (an update is the speed loop's, by ADR).
Node class key: `wheat` for `CPU_SKU:6148` (the `annotation_steps` convention).

Baseline entry: `{median_s, runs_s[3], spread, noisy, nodes[3], cpu: "Intel Xeon Gold 6148
(CPU_SKU:6148)", threads: 16, mem_mb: 64000, task_memory: "56.GB", profile, modules, bakta,
peak_rss_sacct_max_bytes, task_peak_rss_max_bytes, cpu_efficiency_median, steps: {<step>:
{median_s, peak_rss_max_bytes}}, estimate_s, ratio_to_estimate, within_estimate_range,
slurm_jobs[], commit, date, rounds: 3, inputs_sha256, adr}`.

Cell status: `complete` (3 valid repeats), `incomplete` (1–2; `median_s` null, `mean_of_valid_s`
reported and labelled, no baseline entry), `no_data` (0). Exit 0 always when the files are
readable; 2 on an unreadable or schema-invalid `bench_run.json`.

## Outputs and where the numbers go

| Destination | What | Who |
|---|---|---|
| `$GIAP_ROOT/bench/<bench_id>/` | runs, `summary.json`, `summary.md` (never committed) | the job, the summary script |
| `benchmarks/baseline.json` › `standard_run` | the 4 cell entries | performance role, in the evidence PR |
| ADR-M10-D "Benchmark definition and first baseline" (next free number at merge) | this definition, the 4 medians against §13.2, the noise sources, the decision that §13.2's numbers are estimates superseded by the table | performance role |
| user-guide runtime table (question L9) | the `k12` rows: variant, median h:mm, range of the 3 runs, peak memory GiB, node class, date, commit | docs role |
| PROGRESS › Measured numbers | one row per cell (`step = giap run <variant>`), and the per-step rows of the median repeat of each `k12` cell | lead |
| M10-05's checklist | bullet 4 becomes `pass` when the 4 entries exist and ADR-M10-D is accepted | (M10-15 reads it) |

## Examples

Numbers are hand-made, not measurements.

| # | Input | Expected |
|---|---|---|
| 1 | Header `--cpus-per-task=16 --mem=64G` | rejected by the static test: 65,536 / 16 = 4,096 MB > `MaxMemPerCPU` 4,000, so Slurm would allocate 18 CPUs (observed 1846846, 1848590). `--mem=64000M` → 4,000 MB per CPU, `AllocCPUS` 16 |
| 2 | `SLURM_ARRAY_TASK_ID=1`, `BENCH_ID=2026-10-30_1f56790`, dry run | four planned runs in the order k12/core, sim5m/standard, sim5m/core, k12/standard, into `…/2026-10-30_1f56790/r2/k12_core`, …; the k12/core line has `--annotation …/inputs/k12_mg1655.json --modules core` and no `annotation.bakta.db`; the k12/standard line has `--set annotation.bakta.db=full`, no `--annotation`, no `--modules`; every line has `--threads 16 --memory 56.GB --profile farmshare_local,apptainer` and no `--resume` |
| 3 | k12/standard valid repeats `wall_s` 3,912, 3,700, 4,106 | `median_s` 3,912 (1:05:12); `spread` (4,106 − 3,700) / 3,912 = 0.1038, `noisy` false; `ratio_to_estimate` 3,912 / 2,100 = 1.863; `within_estimate_range` false (> 3,300). Status `complete`; a baseline entry is written |
| 4 | k12/core repeats 1,500, 1,520, 1,800 | `median_s` 1,520; `spread` 300 / 1,520 = 0.1974, `noisy` true; ratio 1,520 / 1,500 = 1.013; `within_estimate_range` null (no range for core); entry written with `noisy: true` and the summary prints "noisy (spread 19.7 % > 15 %)" |
| 5 | sim5m/core: r1 `ok` 1,490 s; r2 exit 4 (`failed`, trace FAILED row in `OPERONS`); r3 `ok` 1,510 s | cell `incomplete`, `median_s` null, `mean_of_valid_s` 1,500 labelled "2 of 3 repeats", no baseline entry; the summary names r2's failing step. After a replacement `r2.b` `ok` 1,505 s: `complete`, median 1,505, `superseded: ["r2"]` |
| 6 | a run whose trace has 3 `CACHED` rows (an `OUT/work` was reused) | `status: invalid`, `invalid_reasons: ["cached_tasks:3"]`; excluded from the median |
| 7 | cotenancy samples: (a) `giap-run-pilot-a1` (2 CPUs) on the same node at minute 15; (b) `giap-diann-a2` with 16 CPUs on wheat-07 while the run is on wheat-05; (c) node `CPUAlloc` 64–78 of 80 from other users | (a) `invalid`, `same_node_giap:<jobid>`; (b) `ok`, `contended_other_nodes`; (c) `ok`, `node_cpualloc_median`/`max` reported, no flag |
| 8 | k12/core trace contains `ANNOTATE:BAKTA (k12)` (the `.json.gz` was passed by mistake, so `cli._bakta_json` sent it down the GFF/GenBank path) | `invalid`, `bakta_ran_in_core` |
| 9 | `BENCH_PROFILE=farmshare,apptainer`; separately `BENCH_COMMIT` unset; separately `r1/k12_standard` already exists | each exits 2 before any `srun`, with a one-line reason; no `bench_run.json` is written |
| 10 | Per-step medians: `bakta` realtime 600, 590, 610 s and peak RSS 20.1, 20.6, 20.4 GiB over the three k12/standard repeats | `steps.bakta.median_s` 600, `peak_rss_max_bytes` = 20.6 GiB; so "standard with a recorded JSON" ≈ 3,912 − 600 s can be stated in ADR-M10-D as an approximation, not a measurement |

## Edge cases

- **Failed run.** Kept with its exit code, `sacct` state and the trace's FAILED process; the job
  goes on to the next cell; the cell is `incomplete` until a replacement is valid (example 5).
  `giap` failing before Nextflow starts gives `trace.n_tasks = 0` and `steps = {}`, still a
  `bench_run.json`. A cgroup out-of-memory kill of the step (`sacct` state `OUT_OF_MEMORY`) is
  `failed` with that state; it is a finding about `--memory 56.GB` (question L2), never retried with
  more memory inside the benchmark.
- **Timeout.** `timeout` returns 124 (or 137 after `--kill-after`): `status: timeout`, `wall_s` =
  the elapsed time, excluded from the median, reported as "> 5 h". A job-level SIGTERM (time limit)
  marks the current run `timeout` and the not-started ones absent.
- **Cached work directories.** Fresh `RUN_DIR`, no `--resume`, manifest `resume: false`, zero
  `CACHED` rows, no `cached` step (example 6). Nextflow's launch-directory history (`.nextflow/`)
  is per run (a fresh launch directory), so no session ID is shared.
- **Warm caches that are not work directories.** Numba: fresh per run. Images and conda envs:
  present before, unchanged after (else `invalid`). Page cache: uncontrolled, rotated, reported.
- **Node heterogeneity.** Inside `CPU_SKU:6148` the nodes are the same model; the three repeats may
  land on 1–3 different nodes (`nodes[]` records them; `single_node: true` when all are the same,
  not a flaw). A different class is possible only by overriding the constraint
  (`sbatch --constraint=…`); such runs are invalid (`node_class:<sku>`) and are never
  summarised or mixed with wheat (lead, 2026-10-24, PR #300 review nit 20). Every task must run on the job's node (`host_mismatch` otherwise).
- **No wheat node available** (all drained or the queue waits days): the array pends; after 48 h
  pending the lead logs it under PROGRESS › Blocked and decides on another class by a lead decision
  recorded in ADR-M10-D. Never a silent fallback.
- **Commit drift.** Tasks pinned to `BENCH_COMMIT`; a manifest commit that differs, or a dirty
  worktree, is `commit_mismatch`.
- **Input drift.** The summary compares `inputs.*.sha256` across the 3 repeats of a cell; a
  difference (for example `SIM5M` generation not deterministic) makes the cell `incomplete` with
  `inputs_differ`, and is a giapsim bug.
- **Empty input.** `BENCH_GENOMES` or `BENCH_VARIANTS` empty or unknown → exit 2. A benchmark
  directory with no `bench_run.json` → `summary.json` with every cell `no_data`, no baseline write.
- **Ties.** The median of 3 is the middle value, ties included (`[1500, 1500, 1800]` → 1,500,
  spread 0.2). Cells and steps are written in a fixed order (genome `k12`, `sim5m`; variant
  `standard`, `core`; steps in first-start order of the median repeat's trace). Between attempts
  of one repeat the latest valid attempt wins (letter order `a` < `b` < `c`).
- **Circular and linear contigs.** Both genomes are circular (`[topology=circular]` on K-12;
  giapsim's chromosome and plasmid circular); no linear contig is benchmarked. The topology does
  not change what is timed; the inputs' sha256 values pin it. A linear-contig runtime is not
  claimed anywhere.
- **Minus strand.** Not applicable: the benchmark measures time and memory, not coordinates; no
  coordinate is read or converted by either script.
- **Nextflow task retries.** Counted (`n_retried`), valid: the retry is part of what a user waits
  for. A run with retries is named in ADR-M10-D.
- **sacct gaps.** MaxRSS empty (job too short for one 60-s sample) → `null`, never 0; the trace's
  per-task peak is still reported. `sacct` unreachable after the run → retried 3 times over 2 min,
  then `wall_s` falls back to `script_wall_s` with a warning `wall_from_script`.

## Acceptance criteria

On CI (nothing submitted; stubs for `srun`, `sacct`, `squeue`, `scontrol`, `giap`, `timeout` on
`PATH` in `tmp_path`):

1. Static: `bash -n` passes; the header has `--job-name=giap-bench`, `--cpus-per-task=16`,
   `--mem=64000M` (and `mem / cpus ≤ 4000`), `--constraint=CPU_SKU:6148`, `--array=0-2%1`,
   `--partition=normal`, `--qos=normal`, a `--time` ≥ 4 × `BENCH_RUN_TIMEOUT`; the string
   `--resume` (and `-resume`) does not occur; the profile default is `farmshare_local,apptainer`.
2. Dry run (example 2): for each of the 3 task IDs the planned order equals the rotation table, and
   every planned line has the flags listed there for its variant; the K-12 core line points at a
   decompressed `.json`, never `.json.gz`.
3. Preconditions (example 9) exit 2 with nothing written: profile containing `farmshare`; unset or
   short `BENCH_COMMIT`; existing `RUN_DIR`; empty or unknown genome/variant; missing K-12
   peptides; K-12 genome or JSON sha256 ≠ the committed file's; a failing `giap setup --check` stub.
4. Stubbed full run: with `giap` stubbed to write a minimal manifest and trace, the task writes 4
   `bench_run.json` files that validate against the documented keys, with `wall_s` from the `sacct`
   stub, `status: ok`, exit 0; a second task in another `BENCH_DIR` gives identical files except
   for times and IDs.
5. Validity rules: each of `cached_tasks`, `manifest_resume`, `commit_mismatch`, `host_mismatch`,
   `bakta_ran_in_core`, `same_node_giap`, `image_or_env_created`, `alloc_cpus` is produced by a
   stub built for it (examples 6–8) and makes the task exit 1 while the other runs still run;
   `contended_other_nodes` and node load do not invalidate; precedence `failed` > `timeout` >
   `invalid` > `ok` holds on a run with two causes.
6. `benchmark_summary.py` pure functions reproduce examples 3, 4, 5 (with and without `r2.b`) and
   10 exactly, the tie case `[1500, 1500, 1800]`, `no_data`, and `inputs_differ`.
7. `--write` on a copy of the real `baseline.json`: adds `standard_run.wheat.<cell>` for complete
   cells only; every other top-level key serialises byte-identically before and after; an existing
   entry with another commit and no `--replace` → exit 1, file unchanged.
8. `SIM5M` exists, is exported next to `SIM100K`, has `write_precomputed=False` and a 5,000,000-bp
   chromosome, and two `simulate(SIM5M, SIM5M_SEED)` calls give byte-identical `genome.fasta`,
   `peptides.tsv` and `bakta_like.json` (marked `slow` if it exceeds the unit budget; the spec's
   bound is < 30 s).

On FarmShare (evidence; PR description, ADR-M10-D and PROGRESS, not CI):

9. The 3 array tasks finish; every cell is `complete` (replacements allowed, ≤ 3 attempts per
   cell); each run's `sacct` step shows `AllocCPUS 16` and a node with `CPU_SKU:6148`, recorded in
   `bench_run.json` (`alloc_cpus`, checked as invalid reason `alloc_cpus:<n>`). Lead, 2026-10-24:
   the step's command is `timeout` itself, so the in-step `SLURM_*` print is dropped; the `sacct`
   record is the evidence.
10. `benchmarks/baseline.json` › `standard_run.wheat` holds the 4 entries written by `--write`;
    ADR-M10-D states each median against §13.2 (including when outside 16–55 min or above 25 min,
    which is a finding handed to M10-13, not a failure of this issue); the user-guide table and
    PROGRESS › Measured numbers carry the same numbers as `summary.md` (a test compares the
    committed table rows with the baseline entries).
11. Disk: after the benchmark, no `RUN_DIR/out/work` remains; `$GIAP_ROOT/bench/<bench_id>` is
    reported in PROGRESS with its size.

## Test plan

Test-author, tests first (R2), new files only:

- `tests/unit/test_m10_12_sbatch.py`: criteria 1–5 (static checks as in
  `tests/unit/test_m2_21_nightly_weekly_sbatch.py` and `test_mvp04_sbatch.py`; the stub `PATH`
  pattern of the nightly tests).
- `tests/unit/test_m10_12_summary.py`: criteria 6–7.
- `tests/unit/test_m10_12_sim5m.py`: criterion 8.
- `tests/unit/test_m10_12_docs_sync.py`: the table-vs-baseline check of criterion 10 (skips while
  `standard_run` is absent from `baseline.json`).

Targeted local runs only (venv Python, `-n 2` at most); nothing submits a job from a test.

## Estimate

| Piece | Code lines |
|---|---|
| `benchmark_standard.sbatch`: header, guards, checkout, inputs, rotation, `srun` + timeout, co-tenancy sampler, validity checks, `bench_run.json`, traps | 130 |
| `benchmark_summary.py`: attempt selection, cell summary, baseline entry and merge, markdown | 90 |
| `SIM5M` scenario | 10 |
| **Total** | **≈ 230** (index: 160; +70 for the co-tenancy record, the validity checks and the replacement attempts) |

Tests ≈ 260 lines; user-guide page and ADR-M10-D are prose. One code PR; the evidence (baseline
entries, ADR, guide, PROGRESS) is a second, small PR after the runs. Wall time: 12 runs, probably
1–3 h each (the measured K-12 annotation run took 2 h 08 min), so 12–36 h of one 16-CPU job at a
time.

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | "64 GB" with 16 CPUs: `--mem=64000M` (exactly `MaxMemPerCPU` × 16; 62.5 GiB = 67.1 GB SI) or `--mem=64G` (Slurm then allocates 18 CPUs, so the benchmark would no longer be "16 threads")? | `--mem=64000M`; ADR-M10-D states the reading |
| L2 | The task pool inside the allocation: `--memory 56.GB` (room for the 4-GB-heap JVM and the head), or the profile's 61,952 MB (risk of a cgroup OOM of the JVM)? | `56.GB`, as `weekly.sbatch`; an `OUT_OF_MEMORY` run is a finding, not a retry with more |
| L3 | Concurrency: one repeat at a time (`%1`, 12–36 h, 16 CPUs) or 3 repeats in parallel (`%3`, 48 CPUs, the index's budget, but two repeats can share a node and its memory bandwidth)? The index row planned 6 jobs; this spec has 3 jobs of 4 runs each | `%1`, 3 jobs |
| L4 | `--exclusive` (a whole 80-CPU wheat node, 62 % of the R14 cap, long waits) for a cleaner number, or shared nodes with co-tenancy recorded and same-node GIAP jobs invalidating a run? | shared, recorded |
| L5 | A noisy cell (spread > 0.15): enter it in `baseline.json` with `noisy: true`, or re-run the cell's 3 repeats first? | enter it flagged; the lead may order one re-run of the cell (≤ R9) |
| L6 | Numba cold for every run (a user's first run) or warm (repeat runs)? | cold; the numba compile share is visible in the per-step times |
| L7 | Index L7 makes the variants differ in modules **and** Bakta. §13.2's "core modules alone" reads as "same input, fewer modules". Add a third variant `standard_json` (all modules, recorded JSON), or derive it from the Bakta step median (example 10)? | derive it, labelled an approximation; no third variant (+50 % wall time) |
| L8 | `sim5m/standard` runs live Bakta on simulated sequence, which few DB hits make unrepresentative. Keep it (the speed loop's §14.4 genome), or run only `sim5m/core`? | keep it, labelled "speed-loop baseline, not a user number"; the user guide shows only the `k12` rows |
| L9 | There is no user guide yet (MkDocs is M11). Create `docs/user-guide/runtimes.md` now with the table (docs role), for M11 to link into the site? | yes; the page holds only the table, the node class and the ADR link |
| L10 | Work directories on NFS (`$GIAP_ROOT/bench`, as a user runs) or on node-local `/tmp` (less NFS noise, but not what users do)? | NFS; ADR-M10-D names NFS as a noise source |
| L11 | §12's nightly benchmarks use a median of 5 with a 15 % slowdown fail; §15.2 says median of 3 for this one. Should the nightly ever re-run the standard benchmark against `standard_run`? | no: median of 3, run on demand by the speed loop (M10-13) and before release; the nightly is unchanged |
| L12 | `BENCH_COMMIT`: the main HEAD at submission, after M10-04's fixes have merged, and the report-freeze runs (M10-14) on the same or a later commit? | main HEAD at submission; ADR-M10-D names it; a re-benchmark is the speed loop's |

## Lead decisions (2026-10-24)

Every lead question L1–L12 takes its default. L9: docs/user-guide/runtimes.md is created by the docs role in the evidence PR.
