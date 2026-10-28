# M2-15b validation: sharded TerminatorNet on real genomes (2026-10-28)

This report checks the Acceptance section of `docs/specs/M2-15b-terminatornet-sharding.md` on real
genomes after the sharded scan merged as #455. It covers two things. First, the sharded
`terminatornet.tsv` must be byte-identical (`cmp`) to unsharded references on Nissle 1917,
K-12 MG1655, mhyo and gsulf. Second, a sharded msmeg scan is timed.

**Result: PASS.** All four sharded outputs are byte-identical to their unsharded references. msmeg
scans in **83 min 50 s** wall at 4-CPU chunks, 8 at a time (32 CPUs). The unsharded attempt failed
with exit 140 at 12 h.

sim100k is covered by `tests/integration/test_m2_15b_tnet_shard_real.py` (CORE = 10,000, 4-CPU
chunks, equal to the committed `-n 17` recording), which passed in job 1901547. It was not re-run.

## How it was run

Only the TerminatorNet step was run, not the pipeline. The driver follows the integration test.

- **Script.** `workflow/bin/gippa_terminatornet_shard` from main 5c6f2dab (worktree
  `local/worktrees/tnet-validation`, sha256 `d10d2a80…`). It runs with
  `$GIPPA_ROOT/micromamba/envs/gippa-terminatornet/bin/python`, with `TERMINATORNET_HOME` set to the
  env's `share/terminatornet` (pin `GIPPA_COMMIT` 72180309), `PYTHONNOUSERSITE=1` and no
  `PYTHONPATH`.
- **Working directory.** It holds `genome.fa` and a `model.pkl` symlink, as the module stages them.
- **Steps.** `plan` (default CORE = 250,000), then every `chunk <id> --cpus 4`, then `merge`.
- **Job layout.** Each sample is one Slurm job. Its chunks run `P` at a time with `xargs -P P`,
  and the job holds `4·P` CPUs. msmeg used P = 8 (32 CPUs), the same concurrency as the pipeline
  under `queueSize = 8`. The other samples used P = 4.
- **Timing.** Each subcommand is wrapped in `/usr/bin/time` and its start and end epochs go to
  `events.tsv`.
- **References.** The fresh unsharded references were recorded like
  `scripts/slurm/record_tools.sbatch` (`run_terminatornet`). The command was
  `TerminatorNet.py -f genome.fa -s -n 16 -o terminatornet.tsv` in `micromamba run -n gippa-terminatornet`.
- **Scripts.** The sbatch scripts are kept outside the repo, under
  `$GIPPA_ROOT/validation/m2-15b-tnet/bin/`: `shard.sbatch` takes `SAMPLE GENOME P` and
  `ref.sbatch` takes `SAMPLE GENOME`.
- **Outputs.** Sharded outputs, plans, pickles and logs are in
  `$GIPPA_ROOT/validation/m2-15b-tnet/<sample>/sharded/`.
- **Smoke test.** Before submitting, the harness was checked on a 4.7 kb excerpt of sim100k
  (CORE = 1,000, 10 chunks). Its output was identical to `-n 2`.

**Inputs.** Nissle and K-12 use the `genome.fa` of their pinned recording. For mhyo, gsulf and msmeg
the input is the pipeline's normalised genome, `local/more_datasets/<s>/out/work/preflight/genome.norm.fa`.
The #411 TERMINATORNET tasks staged exactly this file as `genome.fa` (their `.command.run`).

| sample | genome sha256 | records | cleaned length (nt) |
|---|---|---|---|
| nissle1917 | (recording symlink, `bakta_full/nissle1917/.../genome.norm.fa`) | 4 | 5,051,730 + 4,523 + 3,173 + 5,514 |
| k12_mg1655 | (recording symlink) | 1 (U00096.3) | 4,641,652 |
| mhyo | `923e34d48d692d70fe2f7a74916748d8c59aef32d9eabbe7de99dbe016cb929b` | 1 (NC_019552.1) | 836,897 |
| gsulf | `db5a9c0dc449f96ad7749236cbf5d339fd24ad54be5c60a8637db98127715534` | 1 (NC_002939.5) | 3,814,128 |
| msmeg | `456a563def9b963d7ac00df66bbef4c32f2bceb9896c1c636a35c3b0182ff6c4` | 1 (NC_008596.1) | 6,988,209 |

## 1. Byte identity against the pinned references (Nissle 1917, K-12 MG1655)

| sample | reference | sharded job | chunks | `cmp` | sha256 (both) | lines |
|---|---|---|---|---|---|---|
| nissle1917 | `$GIPPA_ROOT/recordings/terminatornet/nissle1917/terminatornet.tsv` (`-n 17`) | 1901829 | 48 | **identical** | `c319d52f8f852e95b157d7ec97c12fd6723413e923ae134d41e90ca6a4bf071e` | 8,644 |
| k12_mg1655 | `$GIPPA_ROOT/recordings/terminatornet/k12_mg1655/terminatornet.tsv` (`-n 17`) | 1901802 | 38 | **identical** | `673078256ba58086e91d2398d7fee9d5bfb49e3b3d539221191eb1923810cd4d` | 8,294 |

Nissle's 48 chunks are 2 × 21 for the chromosome plus 2 × 1 for each of the three plasmids. The
plasmid chunks are the 5 s minimum in the timing table below. Both sha256s match the prefixes in
the spec (`c319d52f…`, `67307825…`).

## 2. New unsharded references (mhyo, gsulf) and byte identity

**New references.** They are stored in `$GIPPA_ROOT/recordings/terminatornet/{mhyo,gsulf}/`. Each
directory has `genome.fa`, a real copy of the input rather than a symlink, so it survives work-dir
cleanup. Next to it are `model.pkl` (a symlink to the env), `command.txt`, `run.log`, `runtime.txt`,
`versions.txt`, `SHA256SUMS` and `.giap-complete`.

The `versions.txt` of both references:

```
terminatornet
d437c79205d1b313c5d060d0bff8c3fa676c6289de7c20011bbc598c8ba27c16  terminatornet.explicit.txt
90d0213b9cf22c8ed4d0603e5f0113ae7b7e0d64656027f4b1fb56d3cfb9cdfb  terminatornet.pip.txt
commit 72180309545988dba63b8bd5db5d0b09bf4dba70
```

The explicit lock and the commit equal the Nissle/K-12 recordings. The pip lock is the current
`workflow/envs/locks/terminatornet.pip.txt`, which the spec already noted differs from the recorded
`ab22b7d9…` only by packages TerminatorNet does not import.

| sample | reference job | `-n` | wall | node | reference sha256 | #411 cross-check | sharded job | chunks | `cmp` |
|---|---|---|---|---|---|---|---|---|---|
| mhyo | 1901796 | 16 | 395 s | rye-01 | `17fd23234143438bbd35d7704b66d60517d011309c0dba2ce93e108b2c9c0170` | identical to the `-n 3` work file `e9/564ba3…` (`17fd2323…`) | 1901798 | 8 | **identical** |
| gsulf | 1901795 | 16 | 3,182 s | rye-01 | `cfa108c9c90c3807a0919df1209a29ded0b181ba67095c28a0be1bf7e83f2d5f` | identical to the `-n 6` work file `dc/af3e1d…` (`cfa108c9…`) | 1901825 | 32 | **identical** |

The fresh `-n 16` references also equal the #411 pipeline outputs at `-n 3` and `-n 6`. This is
further evidence that the output does not depend on `-n` for `-n ≥ 2`. Line counts are 584 for mhyo
and 10,713 for gsulf.

gsulf's `run.log` holds only micromamba's "Could not set lock … Waiting for other mamba process"
warning, from the two reference jobs starting together. TerminatorNet itself printed nothing.

## 3. msmeg timing (sharded only, no reference)

msmeg ran as job **1901794** `gippa-tnetv-msmeg` on rye-01, with 32 CPUs and 40 GB allocated and
8 chunks at a time × 4 CPUs.

| metric | value |
|---|---|
| plan start → merge end | 00:58:17 → 02:22:07 PDT, **5,030 s (83 min 50 s)** |
| first chunk start → merge end | 5,027 s |
| chunks | **56** (2 strands × ⌈6,988,150 / 250,000⌉ = 2 × 28), all exit 0 |
| chunk wall: min / median / max | 454 s / 712 s / **960 s (16 min 0 s)** |
| sum of chunk wall | 39,491 s |
| plan / merge | 3.1 s / 1.3 s |
| CPUs | `--cpus 4` per chunk; **peak 8 concurrent chunks = 32 CPUs** (from the `events.tsv` intervals) |
| max RSS | 155 MB, the largest single process of any chunk (GNU time); 1.20 GB for the whole batch step with 8 chunks running (`sacct` MaxRSS) |
| total CPU (`sacct` TotalCPU) | 41 h 19 min, about 44 CPU-min per chunk |
| output | `terminatornet.tsv` 31,776 lines, sha256 `b44a5747bc09a262504fc29eb8c8d8145a138c39338f3a40e1fdf6fb911ee589` |
| unsharded comparison | the #411 TERMINATORNET task failed with exit 140 at the 12 h limit (PROGRESS.md) |

Notes:
- **Real-pipeline timing.** This is the step's own work at the pipeline's concurrency (8 slots ×
  4 CPUs). A Nextflow run adds Slurm queue wait for each chunk and the submission ramp:
  56 jobs at `submitRateLimit = 10/1min` take about 6 min, overlapped with execution.
- **Limits.** The longest chunk (16 min) is well inside the 2 h × attempt FarmShare time limit.
  The per-chunk RSS is far below the 4 GB × attempt request.
- **CPU cost.** msmeg's chunks cost about 2.2 times more CPU per Mb-strand than K-12: about
  177 CPU-min/Mb-strand against 81. The spec estimated about 76 CPU-min/Mb-strand and 19 CPU-min
  per chunk. The cost per Mb-strand rises with GC content: mhyo (GC 0.26) about 54, Nissle about
  73, K-12 about 81, gsulf (GC 0.61) about 111, msmeg (GC 0.67) about 177. Other jobs were running
  on rye-01 at the same time, so the wall times are indicative only. The spec's estimate of about
  5 min per chunk holds for E. coli (median 264–317 s) but not for msmeg (median 712 s).

## Timings of all sharded runs

All runs were on rye-01, CORE = 250,000, 4-CPU chunks.

| sample | job | P (CPUs) | chunks | plan → merge wall | chunk wall min / median / max | `sacct` TotalCPU |
|---|---|---|---|---|---|---|
| mhyo | 1901798 | 4 (16) | 8 | 441 s | 80 / 213 / 246 s | 1 h 30 min |
| k12_mg1655 | 1901802 | 4 (16) | 38 | 3,055 s | 154 / 317 / 381 s | 12 h 36 min |
| nissle1917 | 1901829 | 4 (16) | 48 | 2,948 s | 5 / 264 / 377 s | 12 h 21 min |
| gsulf | 1901825 | 4 (16) | 32 | 3,380 s | 109 / 433 / 551 s | 14 h 09 min |
| msmeg | 1901794 | 8 (32) | 56 | 5,030 s | 454 / 712 / 960 s | 41 h 19 min |

The unsharded gsulf reference used 13 h 18 min of CPU at `-n 16` (3,182 s wall). The sharded scan
used 14 h 09 min, so sharding adds little CPU overhead.

## Jobs and R14

Every job was named `gippa-tnetv-*`, ran with `--deadline=2026-10-28T14:00:00`, and used
partition `normal`, QOS `normal`.

| job | name | CPUs | mem | state | elapsed |
|---|---|---|---|---|---|
| 1901794 | gippa-tnetv-msmeg | 32 | 40G | COMPLETED 0:0 | 01:23:53 |
| 1901795 | gippa-tnetv-ref-gsulf | 16 | 32G | COMPLETED 0:0 | 00:53:05 |
| 1901796 | gippa-tnetv-ref-mhyo | 16 | 32G | COMPLETED 0:0 | 00:06:39 |
| 1901798 | gippa-tnetv-mhyo | 16 | 20G | COMPLETED 0:0 | 00:07:21 |
| 1901802 | gippa-tnetv-k12 | 16 | 20G | COMPLETED 0:0 | 00:50:56 |
| 1901823 | gippa-tnetv-gsulf | 16 | 20G | CANCELLED while pending (see below) | 0 |
| 1901825 | gippa-tnetv-gsulf | 16 | 20G | COMPLETED 0:0 | 00:56:21 |
| 1901829 | gippa-tnetv-nissle | 16 | 20G | COMPLETED 0:0 | 00:49:09 |

**CPU totals.** Before each sbatch, the `gippa-*`/`giap-*` running and pending totals (CPUs, jobs)
were 1/1, 49/3, 49/3, 98/7 and 66/6. These jobs held at most 64 CPUs at a time.

**Cancelled job.** 1901823 was submitted when the total was already 114 CPUs. The nightly jobs had
started in the meantime, so the sum would have been 130, over the 128 cap. It was cancelled before
it started and resubmitted as 1901825 once the total was 98.

## Not covered here

- **The Nextflow path.** The subworkflow (`TERMINATORNET_SHARDED` under `-profile farmshare,apptainer`)
  was not run end to end, and no Nextflow trace exists for msmeg. The measurement above is the
  step's script at the pipeline's concurrency. A full msmeg pipeline run would add queue waits and
  give the trace-based peak-concurrency figure the spec asks for.
- **The pinned reference sets.** The Nissle and K-12 references were not regenerated; they were used
  as pinned.
