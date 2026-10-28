# ADR-0062: TerminatorNet gets the process_high CPUs; Promotech gets 2 CPUs and 10 GB
- Status: proposed
- Date: 2026-10-26
- Spec sections: §6.3 (process labels), §7.1.1, §14.4 (speed loop), R14
- Needs Khoa: no (no default, threshold, interface or dependency change; outputs byte-identical;
  the module labels stay as in the M2-07 table, only `withName` resource requests change)

## Context
Speed loop on the real runs (Nextflow traces and `sacct` of pilot pass 1–3 and MVP-04; 5-Mb
E. coli genomes, `-profile farmshare`, `--threads 14`):

- `ANNOTATE:TERMINATORNET` is the critical path of every run: 65–68 min on barley, 107–119 min
  on rye-01, 175 min on wheat-01, at 8 CPUs (`process_medium`) and 740–790 % CPU. In the
  pass-3 A1 run (232 min) it ended at minute 124, and nothing downstream of the context layers
  could start before it. TerminatorNet's scan is a `multiprocessing.Pool(-n)` `pool.map` over
  windows (order-preserving, each window independent), so its output does not depend on `-n`.
- `PROMOTECH`/`PROMOTECH_NOVEL` ask for 8 CPUs and 16 GB per 100,000-window chunk but RF-HOT
  predicts on one core (≈ 99 % CPU on an unloaded node, 7.4 GB peak RSS). Across the seven runs,
  8-CPU tasks waited a mean 9.5 min in the Slurm queue (max 118 min; pass-3 A2's Promotech
  chunks sat pending for 100 min and held the `queueSize = 8` slots, delaying OSTIR and the
  mapping steps by about 45 min), whereas tasks of at most 3 CPUs waited a mean ≤ 0.2 min.

## Decision
In `workflow/conf/base.config`:
- `withName: 'TERMINATORNET' { cpus = 16; memory = { 16.GB * task.attempt } }` — the
  `process_high` CPUs at the `process_medium` memory; `resourceLimits` clamp it to `--threads`
  (and to 16 on FarmShare), so 14 on the pilot runs.
- `withName: 'PROMOTECH|PROMOTECH_NOVEL' { cpus = 2; memory = { 10.GB * task.attempt } }`.
  `farmshare.config` repeats the `max(cpus, ceil(memory / 4 GB))` rule for this selector
  (`withName` outranks the label closure), so a FarmShare task asks for 3 CPUs and 10 GB (5 CPUs,
  20 GB on a 137 retry). The module still exports `OMP_NUM_THREADS=${task.cpus}`.

## Evidence (before → after: metrics, tests, simulations)
- TerminatorNet, pass-3 genome (GCF_021559835.1, 3 replicons, `-s`), rye-01 (AMD EPYC 9754):
  8 CPUs 118.5 min (job 1896290) and 107.2 min (job 1895826, pass-3 A2) → 16 CPUs 50.3 min
  (job 1896997, 1,508 % CPU, 0.44 GB RSS). `terminatornet.tsv` is byte-identical (`cmp`) to the
  8-CPU pass-3 A1 output (sha256 3018fa92…).
- Promotech chunk `windows_005.fa` (100,000 windows), rye-01: 3 CPUs / 10 GB 2.3 min wall,
  84 % CPU, 7.5 GiB peak RSS (job 1897274) vs 8 CPUs 2.9–3.7 min for the pass-3 A1 chunks on the
  same node (jobs 1896300, 1896301); `sequences_predictions.csv` byte-identical (`cmp`) to the
  pass-3 A1 chunk (job 1896305; only Promotech's own `sequences.log.txt` timing lines differ,
  which GIPPA does not read).
- Effective requests checked with a local-executor Nextflow run of base.config + farmshare.config
  (`--threads 14`): TERMINATORNET 14 CPUs/16 GB/12 h, PROMOTECH 3 CPUs/10 GB/12 h, OSTIR (control)
  8 CPUs/16 GB.
- `tests/unit/test_nf_modules_lint.py`, `tests/integration/test_nf_config.py`, the M7-15 retry
  policy and growth tests and `test_nf_conda_profile.py` pass unchanged.

## Consequences
- Expected critical-path cut on a 5-Mb genome: 55–70 min of 90–240 min (TerminatorNet), plus the
  Promotech queue wait when the cluster is busy. A 16-CPU request may wait longer than an 8-CPU
  one on a full cluster; the measured mean wait of 14–16-CPU GIPPA tasks was 6 min.
- R14 is unchanged: `queueSize` × the 16-CPU cap still bounds a run at 128 CPUs, and Promotech's
  11 chunks now ask for 33 CPUs instead of 88.
- If a Promotech chunk ever exceeds 10 GB it is killed (137) and retried at 20 GB.

## Amendment 1 (proposed 2026-10-27; M2-15b, docs/specs/M2-15b-terminatornet-sharding.md)
- `ANNOTATE:TERMINATORNET` is replaced by `TERMINATORNET_PLAN` → `TERMINATORNET_CHUNK` (one task
  per 250,000-window core of each contig strand) → `TERMINATORNET_MERGE`
  (`workflow/subworkflows/local/terminatornet_sharded.nf`, `workflow/bin/gippa_terminatornet_shard`)
  behind `params.terminatornet_sharded = true`. `false` runs the unsharded `TERMINATORNET` module
  unchanged, with the `cpus = 16`, `16.GB × attempt` request above; it stays the fallback and the
  acceptance reference.
- Chunk resources: `cpus = params.terminatornet_chunk_cpus` (default 4, `--terminatornet_chunk_cpus N`,
  clamped by `resourceLimits`), `memory = 4.GB × attempt`, and on FarmShare `time = 2.h × attempt`
  with the `max(cpus, ceil(memory / 4000 MB))` rule. PLAN and MERGE are `process_single`. Under
  `queueSize = 8` the step holds at most 8 × 4 = 32 CPUs at a time (16 before); the run's bound of
  128 CPUs is unchanged. A 137/140/143 retry re-runs one chunk, and `-resume` keeps finished chunks.
- Correction: "its output does not depend on `-n`" holds only for `-n ≥ 2`. With a sequence
  ≥ 501 nt, `scan_multiprocess` stable-sorts the calls by stop before `merge_overlapping`; `-n 1`
  runs `scan()`, which keeps window order, so its output can differ. The `-n ≥ 2` output (what the
  deployed module produces) is the reference; the sharded step reproduces it at any chunk CPU count.
- Resume: Nextflow does not hash `workflow/bin/` scripts, so after a change to
  `gippa_terminatornet_shard` a `-resume` reuses the old chunk pickles; re-run without `-resume`
  (or clear the TERMINATORNET_CHUNK work dirs).
- msmeg measurement (wall time from the first PLAN submit to the MERGE end, chunk count, max chunk
  wall time, max RSS and peak concurrent CPUs at `terminatornet_chunk_cpus = 4`, next to the failed
  unsharded attempt, exit 140 at 12 h): pending; recorded here and in PROGRESS.md once it exists.
