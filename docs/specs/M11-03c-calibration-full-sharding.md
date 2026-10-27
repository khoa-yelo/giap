# M11-03c: make every full-tier calibration suite finish inside its time limit

Milestone M11 · **pass 2 (blocks M11-13)** · depends on: M11-03b (merged: `calibration_full_launch.sh`,
task records, `scripts/ci/calibration_rows.py`), M7-13 / ADR-0015 / ADR-0057 (replicate shards,
`--merge`, block 2) · estimate ≈ 640 code lines + ≈ 900 test lines in 5 PRs · **`stats` for PR 2
only** (it changes where block 2 runs, not what it computes; label `stats`, statistician review
before merge) · roles: implementer, test-author, lead (edits to existing launcher tests and
`docs/release.md`, see the PR split) · Khoa: none. Revision 2 folds in the #378 review round
(code review B1–B3, S1–S5; statistician B1, S1–S7) and the lead's Q1–Q8 decisions.

## Spec sections implemented

- §13.2 "Every statistical method passes its calibration suite (§8.4), and every attack null (§8.5)
  stays at level", read by `readiness.py` R3/R4 through the `calibration_full` and
  `calibration_attacks` rows. The rows cannot pass while a suite times out.
- §8.4 "The full suite runs nightly on FarmShare as a Slurm job array, one suite per task" and the
  "Nightly (full)" R column; §3.8 (long work goes to `sbatch`); §11.2.
- **R3:** this issue changes scheduling only. R, tolerances, α, bounds, seeds and the two-block rule
  (ADR-0015) stay as they are. Every replicate `i` of stream `s` still draws from
  `replicate_rng(s, i)` (or `records(i)` with the absolute index), whichever process runs it.
  **R14:** ≤ 128 CPUs and ≤ 20 queued or running `giap-*` jobs, checked before every `sbatch`; job
  names start with `giap-`.

## The problem (array 1894441, commit ef4a4cc77d7b, 2026-10-25 14:48)

`calibration_full_launch.sh` submitted `calibration_full.sbatch --array=0-14%15` (4 CPUs, 16 000 MB,
`--time=04:00:00`). Four tasks timed out at 4:00:09. The aggregator wrote both rows with exit 1
(`4 failed suites (placement no record, regions no record, determinants no record, novel_orfs no
record)`; `nightly.tsv` lines 37–38).

### Measured (sacct, the task JUnit, the harness records' `duration`, and probe jobs 1895351–1895355)

| Task | Suite | Elapsed | CPU-h used / 4 CPUs | MaxRSS | Longest tests (single process) | Single-process work, estimated | Basis |
|---|---|---|---|---|---|---|---|
| 0 | placement | **TIMEOUT** (41/58 tests) | 5.0 | 1.7 G | `attack5_reversed_stretches…shuffle_decoys_full`: **one check, 3.36 h** (R 300, harness, shardable). **28 shard-0-only tests** (conftest's predicates, below): the 17 under `test_placement_fdr_sim.py` / `test_readthrough_*` / `…readthrough_recoded`, each **≥ 69 min** for 13 of them, plus 11 more (`test_attack_denovo_precision` ×3, `…genome_decoys` ×2, the `…has_power…` tests, `…readthrough_decoys_{giapsim,stressed}`) | **≥ 13.4 h, likely 15–20 h**, of which ≥ 10 h shard-0-only | records; probe 1895351: 7 tests still running after 169 min under 8 workers |
| 1 | gene_fdr | 3:46:52 | 3.8 | 1.6 G | 30, 29, 23 min | 3.8 h | JUnit |
| 2 | detectability | 0:14:06 | 0.2 | 0.4 G | 4 min | 0.23 h | JUnit |
| 3 | regions | **TIMEOUT** (22/50) | 4.0 | 0.6 G | `a2_short…`, `a2_mid…`: 70 min each; `c1_level_full[rho0.0]` ≈ 69 min | **≥ 11 h, likely ≈ 12 h** | records; probe 1895352 |
| 4 | autocorrelation | 1:22:53 | 2.1 | 1.0 G | 16 min | 1.4 h | JUnit |
| 5 | islands | 1:47:23 | 1.8 | 0.5 G | 12 min | 1.8 h | JUnit |
| 6 | discordance | 2:25:22 | 2.5 | **13.8 G** | 30, 29, 29 min | 2.4 h | JUnit |
| 7 | missed_operons | 0:18:48 | 0.3 | 6.7 G | 2 min | 0.3 h | JUnit |
| 8 | ora | 3:17:52 | **7.8** | 1.2 G | 53, 34, 29 min (6 tests) | 3.3 h | JUnit |
| 9 | gsea | 3:18:40 | 3.7 | 0.7 G | 27, 26, 24 min | 3.3 h | JUnit |
| 10 | determinants | **TIMEOUT** (2/19) | 4.0 | 0.3 G | `test_null_level_under_attack[full-*-lightgbm]` ≈ 1–1.7 h each; all 19 are **direct** tests (own `for i in range(R)` loop, no harness call), so they run in shard 0 only | **≥ 8.5 h, likely ≈ 10 h** | the two finished tests (2.6 h); probe 1895353 |
| 11 | novel_orfs | **TIMEOUT** (3/56) | 4.6 | 1.0 G | `test_attack_gc_poor_intergenic_fdp_at_level_full` (`calibrate_fdr_counts`): **43.7 s per replicate**, a strict xfail, so it always runs both blocks: 2 × 300 × 43.7 s = **7.3 h**. It was at replicate 327 of 600 at the kill. **6 shard-0-only tests** (`test_planted_tir_calibration_300k_full` ×2, `test_sixframe_decoys_null_and_coin_full`, `test_planted_rbs_pass_detected_gene_quantile_full`, `test_attack_oracle_picked_class_fdp_dense`, the teeth test `…low_prior_attack_has_power…[full]`) ≈ 35 % of the suite (fast-tier audit) | **≥ 13 h, likely 13–18 h**, of which ≈ 4.5–6 h shard-0-only | records; probe 1895355 |
| 12 | denovo_precision | 0:13:47 | 0.2 | 0.7 G | 13 min | 0.22 h | JUnit |
| 13 | variants | 2:47:52 | 4.2 | 3.3 G | shards 0/2 and 1/2 ran concurrently, 1:07 each; `--merge` ran block 2 of the strict-xfail checks (g1 `calibrate_pooled` 80 min, a6b 12, a6c 8), **1:41 sequentially** | 2.2 h shards + 1.7 h merge | JUnit, records |
| 14 | compare | 0:26:57 | 0.4 | 0.4 G | 9 min | 0.44 h | JUnit |

**Total single-process work ≈ 75–90 h.** The array was allocated 144 CPU-h and used 44.8 (31 %).

Findings that shape the design:

1. **The harness runs single-process.** `_criteria.py` runs `pytest` without `-n`, and the sbatch
   adds nothing. CPU-h/elapsed is about 1.0 for 11 of the 15 tasks. Only `ora` is higher (2.4: its
   `function_gate_chain.workers()` prefetch uses `cpu_share()`). Three of every four allocated CPUs
   sat idle.
2. **A Slurm CPU is a hyperthread.** The nodes report `ThreadsPerCore=2` and
   `SelectTypeParameters=CR_CORE_MEMORY`, so 2 CPUs make 1 core and 4 CPUs make 2. In the probe,
   `xdist -n 8` on 8 CPUs (4 cores) slowed every test: `a2_mid_range_prophage…` took 70.8 min
   alone in 1894441 and 173.9 min in probe 1895352, i.e. **2.46× slower**, about 3.3× throughput
   for 8 CPUs. Two processes on 2 cores are fine: the variants shards (2 processes on 4 CPUs) ran
   50 replicates of g1 in 44.6 min, while block 2 alone ran 100 in 80.1 min. So more CPUs per task
   help only as **one process per core**, and only if each process keeps its threads to its own
   core (S1, below).
3. **xdist inside a task wastes work.** Per-process caches are built again in every worker.
   `@cache gate_record`/`gc_record` in `tests/helpers/m7_gate_novel_records.py` are shared by 30+
   novel_orfs tests; under 8 workers, none of the 42 heavy novel_orfs tests finished in 3 h. Tests
   also land on workers in an unpredictable way. This issue does not use xdist.
4. **Replicate shards alone cannot fix these suites.** Shard 0 runs every test that conftest keeps
   in shard 0 only: `calibration_unshardable`, `teeth_log`, and tests that reach neither a harness
   entry point nor `replicate_range` (`not _item_reaches_harness(item, SHARD_AWARE_NAMES)`). In the
   fast-tier shard audit (array 1886453, N = 3), shard 0's extra share of a suite's
   single-process time was 80 % for determinants, about 50 % for regions, islands, detectability
   and variants, and about 35 % for gene_fdr, gsea and novel_orfs. In the full tier, placement's
   shard 0 would hold 28 such tests (≥ 10 h), and novel_orfs' shard 0 would hold 6 (≈ 4.5–6 h).
   This design therefore moves **every** shard-0-only test to group lanes, using the same
   predicates conftest uses (rule "group selection").
5. **Block 2 in `--merge` runs in one process, and its big checks are counts checks.** novel_orfs
   `gc` (`calibrate_fdr_counts`) adds 3.6 h after the shards finish, variants g1
   (`calibrate_pooled`) 80 min. Both go through `_calibrate_counts`, whose confirm branch runs
   `range(R)` just as `calibrate`'s does.
6. **The lead's sketch does not fit the job cap.** With one array task per (suite, shard) and a
   merge job per sharded suite, every unit must stay ≤ 2 h 40 min (⅔ of 4 h). That needs about 35
   tasks, 8 merges and 1 aggregator, i.e. 44 jobs. R14 allows 20, and pending array tasks count
   (#349 B3; confirmed as lead decision Q1).
7. **A timed-out task left no record.** None of the four TIMEOUT tasks wrote `task-<i>.json`, and
   pytest wrote no JUnit, so the time spent left no timing evidence. The rows still failed closed
   ("no record"), which is correct. Bash runs a trap only after the current foreground command
   returns, so a trap cannot help while `uv run …` runs in the foreground. This design never runs
   a long command in the foreground (rule "signals").

## Definition

### Idea: one job per suite, one process per core ("lanes"), merge in the same job

Each of the 15 suites gets **one Slurm job** (`giap-cal-<suite>`) sized by a committed plan. Inside
the job, the suite runs as **L lanes**, one process per physical core (`cpus_per_lane` = 2 Slurm
CPUs: one core, two hyperthreads; ora and autocorrelation use 4). There are two lane kinds:

- **shard lane** `k` of `N`: `_criteria.py --shard k/N`, exactly as `calibration_sharded.sh` runs
  it today, except that the suite's *group tests* are deselected. Every harness check runs its
  replicate slice in every shard lane.
- **group lane** `g` of `G`: an **unsharded** `_criteria.py` run of a fixed subset of the suite's
  tests (block 2 inline, as in today's unsharded run). Each group test runs whole in exactly one
  group lane.

When the shard lanes are done, the job runs `--merge` (only when `N ≥ 1`). Block 2 of the pending
checks runs on `N` lanes (PR 2). The job then writes the suite's `task-<i>.json` itself, in the
schema `calibration_rows.py` reads today. Then the aggregator runs, unchanged.

Why this layout: it keeps **16 jobs** (15 suites + aggregator), as today, so R14's job arithmetic
is unchanged. It uses about 99 CPUs that are actually busy, instead of 60 of which 31 % were
busy. Every lane has its own core (finding 2). It needs no new job type and no cross-job
dependencies apart from the aggregator's. The lead accepted concurrent lanes in one job, one
process per physical core (Q2).

### The plan (committed scheduling data, outside `CALIB_PATHS`)

`scripts/slurm/calibration_plan.tsv`, one row per `replicates_full` key, in that order, with a header:

| Column | Meaning |
|---|---|
| `suite` | `replicates_full` key; the row index is the task index `i` (`task-<i>.json`) |
| `shards` | N ≥ 0 replicate-shard lanes |
| `groups` | G ≥ 0 group lanes; `shards + groups ≥ 1` |
| `group_select` | which tests are group tests: `*` = all (then `shards` must be 0); `-` = none (then `groups` must be 0, and the suite must have no shard-0-only test, see below); otherwise a comma-separated list of the token `direct` and/or **nodeid prefixes**. A test is a group test iff it is shard-0-only (when `direct` is given) or its nodeid starts with one of the prefixes |
| `cpus_per_lane` | even (whole cores) |
| `mem_per_lane_mb` | memory of one lane |
| `time` | the job's `--time` (`HH:MM:SS`, ≤ `08:00:00`) |
| `unknown_min` | duration assumed for a group test missing from the durations file |

**Group selection uses the harness's own rule.** PR 2 adds `shard0_only(item) -> bool` to
`tests/calibration/conftest.py`: `calibration_unshardable` marker, or the `teeth_log` fixture, or
`not _item_reaches_harness(item, SHARD_AWARE_NAMES)`. These are the predicates
`pytest_runtest_setup` already applies, and `pytest_runtest_setup` is refactored to call the helper,
so the two cannot drift. A conftest option `--calibration-lane-report=<path>` writes
`{nodeid: {"shard0_only": bool}}` for the collected items and changes nothing else. The job's
collection step uses it. **Whenever `shards > 0`, every shard-0-only test is a group test**, whatever
`group_select` says. `direct` states this in the plan, and `-` with a shard-0-only test present is
a `PlanError`. Shard 0 then runs only replicate slices, like every other shard lane. Hand prefixes
are an optional extra (for a heavy harness test that should run whole).

`scripts/slurm/calibration_durations.tsv` (`nodeid`, `seconds`, `source`): the single-process
duration of each full-tier test, from the newest full run's JUnit (`source` = job ID) or the probe
logs. Under lead decision Q4, the first run under this design is also the re-planning measurement;
the lead commits its durations TSV afterwards.

Both files live under `scripts/slurm/`, not in `CALIB_PATHS` (`src`, `tests/calibration`, `ci`,
`tests/sim`, `tests/helpers`); lead decision Q6. They decide scheduling only, never R, seeds or
which checks run. A test pins that every registered check runs in exactly one lane (in exactly one
group lane, or in every shard lane). **Assumption, stated for Q6:** a thread count can change BLAS
reductions at ULP level. So `cpus_per_lane` (hence every lane's thread limit, below) is fixed per
suite: 2, or 4 for ora and autocorrelation. A plan PR that changes a suite's `cpus_per_lane` must
say so in its description. Moving tests between lanes or changing `shards`/`groups` changes no
value (statistician review, #378).

**Group assignment (deterministic):** take the group tests from the collection. Sort them by
(duration descending, nodeid ascending) and put each on the group lane with the smallest load so
far, the lowest lane index on ties (LPT). The duration comes from `calibration_durations.tsv`,
else `unknown_min`. Shard lanes and the merge get `--deselect <nodeid>` for every group test. Group
lane `g` gets `--deselect` for every test not in `g`. Everything passes through `PYTEST_ADDOPTS`
(which the harness's pytest subprocesses inherit), each item quoted with `shlex.quote`.

**`--deselect` is a prefix match** (`_pytest/main.py`: `nodeid.startswith(deselect_prefixes)`), not
a set match. `lanes()` raises `PlanError` when a nodeid it deselects for a lane is a proper prefix
of a nodeid that lane keeps (e.g. `…::test_x` deselected, `…::test_x_full` kept). No such pair
exists among today's 436 full-tier nodeids (#378 review).

**Initial plan** (estimate per lane assumes one process per core; ⅔ rule: the job's estimate, shard
0's direct load included, is ≤ ⅔ of `time`; re-planned from the first run, Q4):

| Suite | shards | groups | group_select | lanes × CPUs | mem/lane | time | Estimate per job |
|---|---|---|---|---|---|---|---|
| placement | 2 | 8 | `direct` | 10 × 2 | 4000 | 06:00:00 | ≈ 2.5–3.5 h (attack5 → 1.7 h per shard lane; 28 group tests ≈ 15–20 h / 8, largest ≥ 1.2 h) |
| gene_fdr | 0 | 2 | `*` | 2 × 2 | 4000 | 06:00:00 (was 04:00:00, see note) | 1.9 h |
| detectability | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.25 h |
| regions | 0 | 5 | `*` | 5 × 2 | 4000 | 06:00:00 | ≈ 2.4 h (largest test 1.2 h) |
| autocorrelation | 0 | 1 | `*` | 1 × 4 | 4000 | 06:00:00 (was 04:00:00, see note) | 1.4 h |
| islands | 0 | 1 | `*` | 1 × 2 | 4000 | 06:00:00 (was 04:00:00, see note) | 1.8 h |
| discordance | 0 | 2 | `*` | 2 × 2 | 16000 | 06:00:00 (was 04:00:00, see note) | 1.2 h |
| missed_operons | 0 | 1 | `*` | 1 × 2 | 8000 | 02:00:00 | 0.3 h |
| ora | 0 | 2 | `*` | 2 × 4 | 4000 | 06:00:00 (was 04:00:00, see note) | 1.7 h (53-min test) |
| gsea | 0 | 2 | `*` | 2 × 2 | 4000 | 06:00:00 (was 04:00:00, see note) | 1.7 h |
| determinants | 0 | 5 | `*` | 5 × 2 | 4000 | 06:00:00 | ≈ 2.3 h (largest test ≈ 1.7 h) |
| novel_orfs | 6 | 3 | `direct` | 9 × 2 | 4000 | 06:00:00 | ≈ 2–2.5 h: group lanes 4.5–6 h / 3; shard lanes ≈ 1.1 h, then block 2 (`gc` 3.6 h and the other pending checks) on 6 lanes ≈ 0.7 h |
| denovo_precision | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.25 h |
| variants | 2 | 1 | `direct` | 3 × 2 | 4000 | 06:00:00 (was 04:00:00, see note) | ≈ 1.5 h: shard lanes ≈ 0.6 h, then block 2 on 2 lanes ≈ 0.9 h; group lane ≈ 1 h |
| compare | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.45 h |

Note (2026-10-27, lead): ora and variants measured 3 h 17 m and 2 h 47 m on rye-01 (full run 1894441_8, _13) and passed 3 h 50 m on wheat (1900073, 1900078, killed at the 4 h limit), so both now get the 6 h of the longest suite; the launch window (max time) is unchanged. The same run also killed autocorrelation at 4 h on wheat (1900069) and gene_fdr finished at 3 h 48 m on wheat (1900066), so every 4 h suite (gene_fdr, autocorrelation, islands, discordance, gsea) now gets 6 h as well.

Total 46 lanes, **98 CPUs** + 1 for the aggregator. The longest `time` is 6 h. The run's expected
wall time is about 3–3.5 h, set by placement and novel_orfs. A one-lane `*` row is exactly
today's unsharded task.

**Thread limits per lane (S1).** Lanes share the job's cpuset. Without limits, BLAS, OpenMP and
LightGBM size their pools from the whole allocation, and `src/giap/steps/determinants.py`,
`src/giap/determinants/sensitivity.py` size from `sched_getaffinity`. That would bring back
finding 2's slowdown. Each lane therefore:
- is pinned with `taskset -c <siblings>` to its own core's hyperthreads (`cpus_per_lane / 2`
  cores, sibling lists read from `/sys/devices/system/cpu/cpu<n>/topology/thread_siblings_list`
  within the job's `sched_getaffinity` set), so `sched_getaffinity` inside the lane is the lane's
  share;
- exports `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMBA_NUM_THREADS` and
  `SLURM_CPUS_PER_TASK` = `cpus_per_lane` (`cpu_share()` reads the last: determinants' `_threads`,
  ora's prefetch workers).

If the job's cpuset has fewer whole cores than lanes need, the job exits 2 with a record before any
pytest.

### `scripts/ci/calibration_lanes.py` (new; pure functions plus a thin CLI)

```text
read_plan(path, suites) -> list[PlanRow]                    # validates; raises PlanError
read_durations(path) -> dict[str, float]
lanes(row, report, durations) -> list[Lane]                  # report = {nodeid: shard0_only}
verify_partition(row, lanes, collected) -> list[str]         # problems; [] = ok
verify_counts(lanes, junit_counts) -> list[str]
fold_exit(lane_exits, merge_exit | None) -> int
record_parts(lanes, merge_parts) -> list[str]
durations_from_junit(paths) -> dict[str, float]
```

`Lane(kind, index, n, deselect, part, expected)`: `part` is the JUnit file name and `expected` the
nodeids it must run.

**Plan validation** (`PlanError`, exit 2 from the CLI):
- the suites are not exactly the `replicates_full` keys, in order;
- a count is negative or not an integer, or `shards + groups == 0`;
- `*` with `shards > 0`, `-` with `groups > 0`, or `direct` with `shards == 0` (use `*`);
- `cpus_per_lane` is odd, or `time` is unparsable or over 8 h;
- a prefix list is given with `groups == 0`.

Checks that need the collection:
- `-` while the suite holds a shard-0-only test;
- a nodeid matched by two lanes;
- the deselect prefix collision above;
- with `shards > 0`, no test left for the shard lanes (use `*`).

A group lane with no test is dropped; the record lists it with `tests: 0` and no part.

**In-job partition check (S2, statistician S7).** Before any lane starts, the job runs pytest
`--collect-only -q` once per lane, with that lane's `PYTEST_ADDOPTS` (about 3 s each). It then
calls `verify_partition`:
- every shard lane's set equals the planned shard set (so the shard lanes are identical);
- each group lane's set equals its `expected`;
- the group sets are disjoint from each other and from the shard set;
- the union equals the full collection.

Any problem ends the job with exit 2 and a record. After the lanes, `verify_counts` compares each
lane's JUnit testcase count with `len(expected)`, and a mismatch makes the job exit 2. This does
not apply to lanes that the USR1 path interrupted: those exit 143, and 143 wins in the fold.

**Exit folding** (as `calibration_full.sbatch` folds task 13 today): if any lane, the merge or a
check has an exit ∉ {0, 1}, the result is the first such exit, in order: lanes, merge, checks (a
crash, a timeout `143`, a merge error `2`, a partition or count problem `2`). `calibration_rows.py`
fails such suites closed (#349 B1). Otherwise the result is 1 if any lane or the merge returned 1,
else 0.

**Record** (`task-<i>.json`, written once by the job's EXIT/TERM/USR1 traps, same keys as
M11-03b):
- `job` = `$SLURM_JOB_ID` (a plain job has no `SLURM_ARRAY_JOB_ID`).
- `sharded` = `shards > 0`.
- `exit` = the fold; `args` = `[]`.
- `parts` = **the parts launched**, so a lane that died shows as `junit missing` in the rows:
  - `task-<i>.shard<k>.junit.xml` for every shard lane;
  - `task-<i>.group<g>.junit.xml` for every non-empty group lane;
  - the block-2 parts that `--merge` reports it launched (PR 2).

  All are plain names in `$GIAP_CALIBRATION_LOG_DIR`, as `_safe_part` requires.
- Extra keys, ignored by the rows script:
  - `lanes`: `[{"kind", "index", "exit", "tests", "junit_tests"}]`;
  - `plan`: the plan row as a dict;
  - `checks`: the partition and count problems.

### R4 under a sharded suite (S3; stricter than single-process, by design)

`calibration_rows.suite_result` (unchanged) fails `calibration_attacks` for a `sharded` record with
`exit == 1` and ≥ 1 attack testcase, because a merge exit 1 cannot be attributed to a test. In a
suite with `shards > 0`, a **non-attack** failure in a group lane also gives job exit 1. So
placement, novel_orfs or variants fail R4 as well as R3, whereas a single-process run would fail R3
only. This is conservative (no false pass), and `calibration_rows.py` stays byte-unchanged.
Criterion 4 pins it.

### `scripts/slurm/calibration_suite.sbatch` (new)

Takes `GIAP_SUITE`, `GIAP_WORKTREE` and `GIAP_CALIBRATION_LOG_DIR` from the environment. The
resources (`--cpus-per-task`, `--mem`, `--time`, `--job-name=giap-cal-<suite>`) come from the
launcher's `sbatch` options, so the script carries no `#SBATCH --array`.

1. Set up as `calibration_full.sbatch` does (`farmshare_env.sh`, its own `TMPDIR`, the worktree
   `PYTHONPATH`), but run every Python step with the **venv interpreter**
   `$UV_PROJECT_ENVIRONMENT/bin/python`, never `uv run`. The process tree is then batch shell →
   `_criteria.py` → pytest (→ block-2 pytests), with no uv in between for signals.
2. **Collect** with `--calibration-lane-report`, compute the lanes, run the per-lane
   `--collect-only` checks. A `PlanError`, a collection error or a partition problem ends the job
   with exit 2 and a record.
3. **Start every lane** in the background with its `taskset`, thread variables, `PYTEST_ADDOPTS`
   (`--junitxml=$GIAP_CALIBRATION_LOG_DIR/task-<i>.<shard<k>|group<g>>.junit.xml` plus the deselect
   list), and, for group lanes, `GIAP_CALIBRATION_LOG_DIR=<log>/lanes/<suite>-group<g>` (JUnit
   still goes to `<log>/`). The suite's `<suite>.jsonl` and `<suite>.shards.jsonl` then hold only
   shard records, which is what `merge_shards` and `shard_exit_verdict` expect. Diagnostics files
   (`*_gate_replicates.jsonl` and others) of group tests land under `lanes/<suite>-group<g>/`.
   Nothing reads them for verdicts; whoever rebuilds durations or evidence must look there too.
4. **Wait** for all lanes. If `shards > 0`, run `--merge` with the shard lanes' deselect list,
   `GIAP_CALIBRATION_CONFIRM_LANES=<shards>`, `GIAP_CALIBRATION_CONFIRM_JUNIT=<log>/task-<i>.merge`
   and the shard lanes' cpus (PR 2). Then check the counts and write the record.
5. **Signals and time limit.** Add `#SBATCH --signal=B:USR1@600`.
   - **Every step runs in the background** and is then waited for: collection, each lane and the
     merge, as `( trap - INT; exec taskset … "$py" … ) & pid=$!; wait "$pid"`, looping while
     `wait` returns > 128 because of a trapped signal. So USR1 is handled at once, even during the
     merge. `trap - INT` inside the subshell restores SIGINT, which bash otherwise sets to ignored
     for asynchronous commands without job control (tested in review: Python then never raises
     KeyboardInterrupt).
   - **On USR1** the job finds the **pytest descendants of the batch shell** by walking the
     `/proc/<pid>/stat` parent links from `$$` (or `pgrep -P`, recursively), including the block-2
     pytests that PR 2 starts. It sends them SIGINT; it never uses `pkill -f`, `-u $USER` or a
     process group shared with other jobs (a launch can put several suite jobs on one node). It
     signals pytest, not the `_criteria.py` wrapper, because `subprocess.run` kills its child
     0.25 s after its own KeyboardInterrupt.
   - The job then waits ≤ 120 s (`GIAP_CALIBRATION_GRACE`, for tests). pytest writes its JUnit with
     the tests finished so far. The job writes the record with `exit 143`, once (as in M11-03b).
     Any process still alive after the grace gets SIGKILL from Slurm; its JUnit is missing and the
     rows fail closed.

`calibration_full.sbatch` and `calibration_sharded.sh` stay byte-unchanged; they remain the tools
for hand runs and for re-running one suite. Retiring task 13's in-task `n=2` concurrent shards
once the lanes land goes to the backlog (lead decision Q7).

### `calibration_full_launch.sh` (changed)

1. Unchanged: argument refusal, `git fetch`, resolving the SHA, the detached worktree under flock,
   a fresh log directory, `GIAP_NOW` / `GIAP_NIGHTLY_AT`, and the fail-closed `squeue`.
2. Validate the plan of the **worktree** (`calibration_lanes.py check`); a `PlanError` gives exit 2
   and nothing submitted.
3. **R14 jobs:** refuse when `used_jobs + 16 > 20` (15 suite jobs + aggregator; same rule as today).
4. **R14 CPUs (replaces the throttle):** each suite is one job, so no array throttle exists. Refuse
   when `used_cpus + Σ_suites lanes × cpus_per_lane + 1 > 128`, counting the requested CPUs (Q2);
   the message names the numbers. With the initial plan this means `used_cpus ≤ 29`. On a busy
   queue the launch is refused, with no second wave (lead decision Q3).
5. **Nightly window:** refuse when `now + max_suites(time) + 20 min` is at or after the next nightly
   start (same python3 arithmetic as today). With the initial plan that is now + 6 h 20 min, so the
   last accepted launch is **19:09** for a 01:30 nightly (19:10 ends exactly at 01:30 and is
   refused).
6. Submit the 15 suite jobs in task order. If any `sbatch` fails, `scancel` **exactly** the jobs
   already submitted, by ID, and exit 1. Then submit the aggregator (`giap-calibration-rows`, 1 CPU,
   `--dependency=afterany:<j0>:<j1>:…:<j14>`, `--job <j0>`, `--expect 15`, otherwise as today). If
   that fails, `scancel` all 15 jobs and exit 1.
7. Print `jobs=<j0>,…,<j14> rows=<id> log=<dir> cpus=<Σ>`.

`calibration_rows.py` is byte-unchanged. The row's job field is `<j0>_calibration` and its text
says `array <j0>` (lead decision Q5: fine).

### PR 2 (stats): block 2 on the shard lanes, and `shard0_only`

`tests/calibration/_criteria.py`:
- **A separate variable.** `GIAP_CALIBRATION_CONFIRM_SLICE=k/L` is read **only** inside the confirm
  branches of `calibrate` **and** `_calibrate_counts` (the latter serves `calibrate_pooled` and
  `calibrate_fdr_counts`). There, `reps = shard_replicates(R, L, k)` of block 2 for a listed check
  (each check uses its own R), logged as a block-2 `shard` record. `GIAP_CALIBRATION_REPLICATES`
  and `GIAP_CALIBRATION_BLOCK` stay unset under CONFIRM, so `replicate_range`, `owns_direct_checks`
  and every conftest reader see exactly what the sequential confirm run sees today (e.g. `_gc`'s
  `first = gc_record(next(iter(C.replicate_range(R))))` vacuity guard still reads replicate 0).
  Without the slice variable, behaviour is byte-unchanged.
- **Empty slices** (R < L, e.g. R = 2, L = 3): the confirm branch logs nothing for that lane,
  mirroring "an empty shard logs nothing" (today it would `_append` a `start == stop` record that
  `_collect` rejects as malformed).
- **`_run_block2`** with `GIAP_CALIBRATION_CONFIRM_LANES=L`:
  - `L ≥ 2`: starts `L` pytest processes at once, each with `CONFIRM`, `CONFIRM_SLICE=k/L` and
    `--junitxml=<GIAP_CALIBRATION_CONFIRM_JUNIT>.lane<k>.junit.xml` **on the command line** (it wins
    over the merge's inherited `PYTEST_ADDOPTS` `--junitxml`, so the L processes never share a
    file). It waits for all of them, folds their exits with `fold_shard_exits`, then merges again.
  - Every run with `CONFIRM_LANES` set (`L ≥ 1`) writes `<log>/<suite>.confirm_parts.json`: the
    list of JUnit names it launched (empty when nothing was pending). The job lists exactly these in `parts`. A missing file after a
    merge that exited 0 or 1 is a check problem (exit 2).
  - A non-integer value, `L < 1`, or `CONFIRM_LANES` set without `CONFIRM_JUNIT` returns
    **`MERGE_ERROR`**; it never falls back silently to one lane. `L = 1`: one process (slice `0/1`,
    all of block 2) with its own lane-0 JUnit, listed in `confirm_parts.json`. Unset: unchanged
    (one process, today's JUnit, no `confirm_parts.json`).
  - Block-2 lanes never write `<suite>.shards.jsonl`, and they write no `deferred` outcome
    (conftest sets `shard = None` under CONFIRM).
- **Unchanged:** `check`, the criteria, R, α, the bounds, `_verdict`, `merge_shards`/`_collect`,
  `_merge_pooled`, the xfail logic (ADR-0057) and the seeds.

`tests/calibration/conftest.py`: `shard0_only(item)` and `--calibration-lane-report` (above);
`pytest_runtest_setup` calls the helper.

The **test-author audits** the 18 files calling `C.replicate_range` directly: under the new slice
variable nothing changes for them (it is not read there), and the audit confirms that no direct
caller relies on CONFIRM setting a slice, and that per-slice reach checks (`check_reach`) run only
where they do today.

Optional (statistician nit): concurrent `O_APPEND` writes to `<suite>.jsonl` on the shared
filesystem are not guaranteed atomic. A torn line raises a JSON error, giving `MERGE_ERROR` (fail
closed, already true for shard lanes). Per-lane record files would remove the risk; not required.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | The initial plan; `replicates_full` keys unchanged | `read_plan` returns 15 rows in key order; total lanes 46, CPUs 98; `max(time)` = 6 h |
| 2 | determinants row (`0 5 * 2 …`), collection = the 19 nodeids, durations: `[full-autocorr_pos-lightgbm]` 8366 s, `[full-detect-lightgbm]` 7784 s (probe job 1895353), the rest from the TSV | 5 group lanes; each nodeid in exactly one lane; the two largest go to lanes 0 and 1; each lane's `deselect` is the other 4 lanes' tests; no shard lane, no merge; record `sharded: false`, parts `task-10.group0..4.junit.xml` |
| 3 | placement row (`2 8 direct`), lane report marking the 28 shard-0-only tests | shard lanes 0/2 and 1/2 deselect exactly those 28 and keep the same 30 tests; the 28 are spread over 8 group lanes by LPT; `attack5_reversed_stretches…full` runs in both shard lanes (replicates 0–149, 150–299); the merge deselects the 28; record `sharded: true` |
| 4 | novel_orfs row with `group_select` `-` | `PlanError` at collection: 6 shard-0-only tests present (`test_planted_tir_calibration_300k_full[…]` ×2, …); exit 2, no lane starts |
| 5 | `fold_exit([0, 0, 1], 0)` / `([0, 143, 1], 0)` / `([0, 0], 2)` / `([0, 0], None)` | 1 / 143 / 2 / 0 |
| 6 | Launcher at 10:00, queue = `giap-nightly-launch` (1 CPU, pending) | 15 `sbatch` calls, `giap-cal-placement` with `--cpus-per-task=20 --mem=40000 --time=06:00:00`, …; aggregator `--dependency=afterany:<15 IDs>`; prints `cpus=98`; jobs 1 + 16 = 17 ≤ 20; CPUs 1 + 99 = 100 ≤ 128 |
| 7 | Launcher with 30 CPUs of other `giap-*` jobs queued | 30 + 98 + 1 = 129 > 128: exit 2, no `sbatch`, message `R14: … 30 CPUs queued + 99 > 128` |
| 8 | Launcher with `GIAP_NOW` = 19:09 / 19:10 local, nightly 01:30, max `time` 6 h | 19:09: end 01:29, submits. 19:10: end 01:30 ≥ 01:30, exit 2, nothing submitted. (Without PR 2, novel_orfs needs `08:00:00`; the last launch is then 17:09) |
| 9 | The 9th `sbatch` (suite `ora`) fails | `scancel` of exactly the 8 IDs already submitted; no aggregator; exit 1 |
| 10 | A lane and then the merge are still running at limit − 10 min (stub pytest sleeps, traps SIGINT, writes JUnit) | USR1 handled within a second (the step runs as `& wait`); SIGINT reaches only the job's own pytest descendants (a decoy `sleep` named pytest outside the tree survives); partial JUnit exists; `task-<i>.json` `exit: 143`, written once; the rows fail closed with `(<suite> exit 143)` |

## Edge cases

- **New test added after the plan:** with `*` it gets `unknown_min` and is placed by LPT; with
  `direct` it goes to a group lane if shard-0-only, else to the shard lanes; a matching prefix
  makes it a group test. No test is ever dropped (asserted, then verified in the job).
- **Fewer group tests than lanes:** the empty lanes are dropped (their CPUs sit idle; the record's
  `lanes` names them with `tests: 0`). **R < N** in a shard lane: an empty shard logs nothing
  (existing). **R < L** in a block-2 lane: logs nothing (PR 2).
- **Ties in LPT:** equal durations are ordered by nodeid, equal loads by the lowest lane index, so
  the assignment is a pure function of (plan, lane report, durations, collection).
- **Strict xfail in a group lane** (placement contaminants `flat_c200`): runs unsharded, block 2
  inline; pytest reports XFAIL in that lane's JUnit. In a shard lane, `--merge` judges it
  (ADR-0057), with block 2 on lanes giving the same `outcome` (criterion 7). A test deselected from
  every shard lane writes no `deferred` outcome, so `--merge` cannot mistake it for an XPASS.
- **A lane crashes or is OOM-killed:** its exit is ∉ {0, 1} and its JUnit is missing or short. The
  fold gives that exit or 2; `calibration_rows.py` fails the suite closed (`junit missing` / `exit …`).
- **A block-2 lane dies:** its replicates are missing, so `_collect` refuses the block and `--merge`
  returns `MERGE_ERROR`; its launched JUnit part is missing, so the rows say `junit missing`. Two
  independent fail-closed paths.
- **The whole job is SIGKILLed (node failure):** no record; the rows say `no record`, as today.
- **The merge times out:** the USR1 path covers it (the merge runs as `& wait`). Result: `exit 143`,
  fail closed, and a rerun is needed, never a pass.
- **Empty input:** an empty collection for a suite is a collection error, exit 2. An empty plan is a
  `PlanError`. An empty durations file means every test uses `unknown_min`.
- **Circular/linear replicons, minus strand:** not applicable; no genomic computation changes.

## Acceptance criteria

1. `read_plan` accepts the committed plan and rejects each plan-validation rule with a message
   naming the row and column (one parametrised test per rule).
2. `lanes()` property test over random collections (0–80 nodeids, including generated prefix pairs
   such as `t_x` / `t_x_full`), random shard-0-only flags, random durations and random valid rows.
   The test models pytest's `startswith` deselection, not set membership. It checks:
   - every nodeid is in exactly one group lane or in every shard lane, never both;
   - every shard-0-only test is a group test whenever `shards > 0`;
   - a prefix collision raises `PlanError`;
   - the output does not depend on input order;
   - examples 2–4 come out exactly.
3. `verify_partition` / `verify_counts` / `fold_exit`: example 5; the rule over all combinations of
   {0, 1, 2, 143} for ≤ 4 lanes plus the merge; a lane whose collected set or JUnit count differs
   from `expected` gives exit 2.
4. `calibration_suite.sbatch` with a stub venv python (pattern of `test_m11_03b_sbatch.py`):
   - examples 2 and 3 run every lane with its own `--junitxml`, deselect list, `taskset` and thread
     variables;
   - the per-lane `--collect-only` checks run before any lane;
   - the merge runs after every shard lane, and only when `shards > 0`;
   - there is one `task-<i>.json`, with `job` = `SLURM_JOB_ID` and `parts` = the parts launched
     (a stub lane that writes no JUnit gives `junit missing` in `calibration_rows.suite_result`);
   - **R4 strictness:** a sharded suite whose group lane fails one non-attack test fails both rows.
5. Example 10, with the stub trapping SIGINT and a shortened `GIAP_CALIBRATION_GRACE`. Also: a stub
   that ignores SIGINT still yields a record after the grace, and the decoy process outside the
   job's tree is not signalled.
6. Launcher, in the `slurm_sandbox` fakes:
   - examples 6–9;
   - `squeue` runs before the first `sbatch`;
   - a plan error refuses before any `sbatch`;
   - every job name starts with `giap-`;
   - the aggregator depends on exactly the 15 IDs;
   - harness arguments are still refused.
7. **PR 2 equivalence** (`test_criteria_shard_equivalence.py` pattern, on
   `_shard_fixture_suite.py` extended as needed). For `CONFIRM_LANES` = 1, 2 and 3 against the
   unsharded run, compare `values` byte for byte, plus `status`, `observed`, `bound`, `seed`, `ok`,
   `R`, `block`, `outcome` and the `--merge` exit code. Cover these checks:
   - a `calibrate` check that fails block 1;
   - a `calibrate_pooled` and a `calibrate_fdr_counts` check that fail block 1;
   - a check whose block 1 is `vacuous` (it must not reach block 2);
   - **both strict-xfail directions**: confirmed fail → `xfail`, exit 0; block 1 fails and block 2
     passes (`unconfirmed`) → `XPASS(strict)`, exit 1;
   - **R = 2 with L = 3** (empty slice, no malformed record).

   Also: bad `CONFIRM_LANES` values (`x`, `0`, `-1`, set without `CONFIRM_JUNIT`) give
   `MERGE_ERROR`; the L processes write L distinct JUnit files and `confirm_parts.json` lists
   them; `tolerances.yaml` is byte-unchanged.
8. **Group-lane equivalence (statistician S6):** on the fixture suite (or a small fast-tier suite),
   G = 1, 2 and 3 give the same per-check records (`values`, `status`), apart from order and
   `duration`, and the same JUnit outcomes. `shard0_only` agrees with `pytest_runtest_setup`'s
   skips on every collected item of the fixture suite.
9. `calibration_full.sbatch`, `calibration_sharded.sh`, `calibration_rows.py` and `readiness.py`
   are byte-unchanged (diff in review). `ruff`, `mypy` (new Python), `shellcheck` (if available)
   pass.
10. **The run itself (lead, after PR 5):** one launch on main finishes with every suite job
    `COMPLETED`, each within ⅔ of its `time`. Both rows are then written, and R3/R4 move from
    `fail` to whatever the suites say. This run is the re-planning measurement (Q4). The lead
    commits its durations TSV, and re-balances the plan in a data-only PR if a job exceeded ⅔.

## PR split (code lines per ADR-0007; tests and the durations TSV, ≈ 436 data rows, are outside the budget)

| PR | Content | Code | Tests | Roles | Depends |
|---|---|---|---|---|---|
| 1a | `scripts/ci/calibration_lanes.py` part 1: `read_plan`, validation, `fold_exit`, `record_parts`, `check` CLI; `scripts/slurm/calibration_plan.tsv` | ≈ 130 | ≈ 150 new (`tests/unit/test_m11_03c_plan.py`) | implementer, test-author | none |
| 1b | part 2: `lanes()` (LPT, deselect, prefix collisions, `direct`), `verify_partition`, `verify_counts`, `durations_from_junit`, `lanes`/`durations` CLI; `scripts/slurm/calibration_durations.tsv` (data) | ≈ 150 | ≈ 200 new (`tests/unit/test_m11_03c_lanes.py`) | implementer, test-author | 1a |
| 2 | `_criteria.py` confirm slice and lanes (both confirm branches), `confirm_parts.json`; conftest `shard0_only` and `--calibration-lane-report` (`stats` label) | ≈ 90 | ≈ 250 new (`tests/calibration/test_criteria_confirm_lanes.py`, the group-lane equivalence, the `replicate_range` audit notes) | implementer, test-author, statistician review | none (parallel with 1a/1b) |
| 3 | `scripts/slurm/calibration_suite.sbatch` (collection, partition check, lanes with taskset and thread limits, merge, counts, record, USR1) | ≈ 170 | ≈ 180 new (`tests/unit/test_m11_03c_suite_sbatch.py`) | implementer, test-author | 1b, 2 |
| 4 | `calibration_full_launch.sh` (per-suite jobs, CPU fit, window by `max(time)`, cancel on failure), changelog fragment | ≈ 100 | ≈ 120 new (`tests/unit/test_m11_03c_launcher.py`) **plus lead edits to the existing `tests/unit/test_m11_03b_launcher.py` (`0-14%15`, throttle, 4-CPU arithmetic) and `tests/unit/test_m11_03b_review_launcher.py:43` (`throttle=`)** and to **`docs/release.md:147`** (the old `array=… throttle=<T>` output line) | implementer, test-author (new file); **lead** edits existing tests and the runbook (role guard) | 3 |

`test_calibration_full_array.py`, `test_calibration_sharded_script.py`, `test_m11_03b_sbatch.py` and
`test_calibration_full_worktree.py` are not touched (their scripts are unchanged). Order: 1a → 1b,
2 in parallel, then 3 → 4. After PR 4 the lead launches the run (criterion 10), records the job
IDs in PROGRESS.md and commits the measured durations TSV.

## Lead decisions (2026-10-25, #378)

- **Q1:** no; pending array tasks count as jobs for R14 (#349 B3 stands). The one-job-per-suite lane
  design stays.
- **Q2:** yes; concurrent lanes in one job, one process per physical core; each lane requests 2 Slurm
  CPUs (ora and autocorrelation 4: two cores); the R14 sum counts the requested CPUs.
- **Q3:** refuse on a busy queue; no second wave.
- **Q4:** the first run under the new design is also the re-planning measurement; the lead commits
  its durations TSV afterwards.
- **Q5:** `<j0>_calibration` is fine.
- **Q6:** the plan and durations files stay outside `CALIB_PATHS`; they decide scheduling only; a
  test pins that every registered check runs in exactly one lane.
- **Q7:** retire task 13's in-task `n=2` once the lanes land; backlog.
- **Q8:** no separate investigation; `--signal=B:USR1@600` and the `& wait` steps replace the trap
  path, and PR 3's test covers a timeout writing `exit 143` plus partial JUnit.
- PR 2 gets the `stats` label and statistician review before merge.

## Notes for the backlog

- **Permanent block-2 cost of the accepted strict xfails.** ADR-0049 and ADR-0050 were
  `accepted-by-khoa` on 2026-10-23. ADR-0050 decision 1 keeps attack GC at q 0.01 (and A6c) as
  strict xfails. Under ADR-0015 and ADR-0057 such a check must run block 2 to be XFAIL. A block-2
  pass makes it `unconfirmed`, i.e. an XPASS(strict), which is how a fix shows up. So novel_orfs
  `gc`'s 2 × 300 × 43.7 s (7.3 h single-process) and variants g1's two blocks are **permanent**.
  The only allowed saving is running block 2 in parallel (PR 2). Skipping or stopping block 2 early
  would change the decision rule (R3).
- The `GC_FULL_XFAIL` reason text still reads "ADR-0050 (proposed, pending Khoa …)"; fix the text
  in a test PR (outside this issue).
- Retire `calibration_full.sbatch` task 13's in-task `n=2` (Q7).
- The 14 `slow`-but-not-`calibration_full` tests (M11-03b Q2) are still not scheduled.
