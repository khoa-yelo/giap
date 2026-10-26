# M11-03c: make every full-tier calibration suite finish inside its time limit

Milestone M11 · **pass 2 (blocks M11-13)** · depends on: M11-03b (merged: `calibration_full_launch.sh`,
task records, `scripts/ci/calibration_rows.py`), M7-13 / ADR-0015 / ADR-0057 (replicate shards,
`--merge`, block 2) · estimate ≈ 560 code lines + ≈ 620 test lines in 4 PRs · **`stats` for PR 2
only** (it changes where block 2 runs, not what it computes; statistician review) · roles:
implementer, test-author, lead (edits to existing launcher tests, see the PR split) · Khoa: none
blocking (Q1–Q8 are lead questions with reversible defaults).

## Spec sections implemented

- §13.2 "Every statistical method passes its calibration suite (§8.4), and every attack null (§8.5)
  stays at level", read by `readiness.py` R3/R4 through the `calibration_full` and
  `calibration_attacks` rows. The rows cannot pass while a suite times out.
- §8.4 "The full suite runs nightly on FarmShare as a Slurm job array, one suite per task" and the
  "Nightly (full)" R column; §3.8 (long work goes to `sbatch`); §11.2.
- **R3:** this issue changes scheduling only. R, tolerances, α, bounds, seeds and the two-block rule
  (ADR-0015) stay as they are. Every replicate `i` of stream `s` still draws from
  `replicate_rng(s, i)`, whichever process runs it. **R14:** ≤ 128 CPUs and ≤ 20 queued or running
  `giap-*` jobs, checked before every `sbatch`; job names start with `giap-`.

## The problem (array 1894441, commit ef4a4cc77d7b, 2026-10-25 14:48)

`calibration_full_launch.sh` submitted `calibration_full.sbatch --array=0-14%15` (4 CPUs, 16 000 MB,
`--time=04:00:00`). Four tasks timed out at 4:00:09. The aggregator wrote both rows with exit 1
(`4 failed suites (placement no record, regions no record, determinants no record, novel_orfs no
record)`; `nightly.tsv` lines 37–38).

### Measured (sacct, the task JUnit, the harness records' `duration`, and probe jobs 1895351–1895355)

| Task | Suite | Elapsed | CPU-h used / 4 CPUs | MaxRSS | Longest tests (single process) | Single-process work, estimated | Basis |
|---|---|---|---|---|---|---|---|
| 0 | placement | **TIMEOUT** (41/58 tests) | 5.0 | 1.7 G | `attack5_reversed_stretches…shuffle_decoys_full`: **one check, 3.36 h** (R 300, harness, shardable); 13 `calibration_unshardable` / direct tests (`test_placement_fdr_sim.py`, `test_readthrough_*`) each **≥ 69 min** | **≥ 13.4 h, likely 15–20 h** | records; probe 1895351: 7 tests still running after 169 min under 8 workers |
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
| 11 | novel_orfs | **TIMEOUT** (3/56) | 4.6 | 1.0 G | `test_attack_gc_poor_intergenic_fdp_at_level_full`: **43.7 s per replicate**, a strict xfail (ADR-0050), so it always runs both blocks: 2 × 300 × 43.7 s = **7.3 h**. It was at replicate 327 of 600 at the kill (`novel_orfs_gate_replicates.jsonl`) | **≥ 13 h, likely 13–18 h** | records; probe 1895355 |
| 12 | denovo_precision | 0:13:47 | 0.2 | 0.7 G | 13 min | 0.22 h | JUnit |
| 13 | variants | 2:47:52 | 4.2 | 3.3 G | shards 0/2 and 1/2 ran concurrently, 1:07 each; `--merge` ran block 2 of the strict-xfail checks (g1 80 min, a6b 12, a6c 8), **1:41 sequentially** | 2.2 h shards + 1.7 h merge | JUnit, records |
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
   help only as **one process per core**.
3. **xdist inside a task wastes work.** Per-process caches are built again in every worker.
   `@cache gate_record`/`gc_record` in `tests/helpers/m7_gate_novel_records.py` are shared by 30+
   novel_orfs tests; under 8 workers, none of the 42 heavy novel_orfs tests finished in 3 h. Tests
   also land on workers in an unpredictable way. This issue does not use xdist.
4. **Replicate shards alone cannot fix these suites.** Shard 0 runs every direct and every
   `calibration_unshardable` test. In the fast-tier shard audit (array 1886453, N = 3), shard 0's
   extra share of a suite's single-process time was 80 % for determinants, about 50 % for
   regions, islands, detectability and variants, and about 35 % for gene_fdr, gsea and novel_orfs.
   In the full tier, placement's shard 0 would hold ≥ 13 unshardable tests (≥ 15 h), and
   determinants' shard 0 would hold the whole suite.
5. **Block 2 in `--merge` runs in one process.** For novel_orfs that is ≥ 3.6 h of `gc` alone
   after the shards finish. For variants it was 1:41.
6. **The lead's sketch does not fit the job cap.** With one array task per (suite, shard) and a
   merge job per sharded suite, every unit must stay ≤ 2 h 40 min (⅔ of 4 h). That needs about 35
   tasks, 8 merges and 1 aggregator, i.e. 44 jobs. R14 allows 20, and pending array tasks count
   (#349 review B3). Even with 8 h tasks it is 29 jobs. Q1 asks whether a throttled array counts
   as one job; if it does, this sketch becomes possible.
7. **A timed-out task left no record.** The spec assumed SIGTERM at the limit would make the trap
   write `exit 143`. None of the four TIMEOUT tasks wrote `task-<i>.json`, and pytest wrote no
   JUnit, so the time spent left no timing evidence. The rows still failed closed ("no record"),
   which is correct. The root cause is not established (KillWait is 30 s, the trap waits for the
   foreground `uv`); Q8.

## Definition

### Idea: one job per suite, one process per core ("lanes"), merge in the same job

Each of the 15 suites gets **one Slurm job** (`giap-cal-<suite>`) sized by a committed plan. Inside
the job, the suite runs as **L lanes**, one process per physical core (`cpus_per_lane` = 2 Slurm
CPUs by default). There are two lane kinds:

- **shard lane** `k` of `N`: `_criteria.py --shard k/N`, exactly as `calibration_sharded.sh` runs
  it today, except that the suite's *group tests* are deselected (below). Every harness check
  runs its replicate slice in every shard lane.
- **group lane** `g` of `G`: an **unsharded** `_criteria.py` run of a fixed subset of the suite's
  tests (block 2 inline, as in today's unsharded run). Each group test runs whole in exactly one
  group lane.

When the shard lanes are done, the job runs `--merge` (only when `N ≥ 1`). With PR 2, block 2 of
the pending checks runs on the same `N` cores. The job then writes the suite's `task-<i>.json`
itself, in the same schema `calibration_rows.py` reads today. Then the aggregator runs, unchanged.

Why this layout: it keeps **16 jobs** (15 suites + aggregator), as today, so R14's job arithmetic
is unchanged. It uses about 89 CPUs that are actually busy, instead of 60 of which 31 % were
busy. Every lane has its own core (finding 2). It needs no new job type and no cross-job
dependencies apart from the aggregator's.

### The plan (committed scheduling data, outside `CALIB_PATHS`)

`scripts/slurm/calibration_plan.tsv`, one row per `replicates_full` key, in that order, with a header:

| Column | Meaning |
|---|---|
| `suite` | `replicates_full` key; the row index is the task index `i` (`task-<i>.json`) |
| `shards` | N ≥ 0 replicate-shard lanes |
| `groups` | G ≥ 0 group lanes; `shards + groups ≥ 1` |
| `group_select` | `*` = every test of the suite is a group test (then `shards` must be 0); `-` = none (then `groups` must be 0); else a comma-separated list of **nodeid prefixes**: a test is a group test iff its nodeid starts with one of them |
| `cpus_per_lane` | even (whole cores), default 2 |
| `mem_per_lane_mb` | memory of one lane |
| `time` | the job's `--time` (`HH:MM:SS`, ≤ `08:00:00`) |
| `unknown_min` | duration assumed for a group test missing from the durations file |

`scripts/slurm/calibration_durations.tsv` (`nodeid`, `seconds`, `source`): the single-process
duration of each full-tier test, taken from the newest full run's JUnit (`source` = job ID) or
probe logs. The planner tool rebuilds it from a log directory.

Both files live under `scripts/slurm/`, not in `CALIB_PATHS` (`src`, `tests/calibration`, `ci`,
`tests/sim`, `tests/helpers`). A plan change moves work between processes and changes no result,
so a code-equal calibration row stays code-equal (Q6).

**Group assignment (deterministic):** collect the suite's full-tier tests (`pytest --collect-only
-q` with `_criteria.pytest_args(suite, True)`) and take the group tests. Sort them by (duration
descending, nodeid ascending) and put each on the group lane with the smallest load so far, the
lowest lane index on ties (LPT). The duration comes from `calibration_durations.tsv`, else
`unknown_min`. Shard lanes and the merge get `--deselect <nodeid>` for every group test. Group lane
`g` gets `--deselect` for every test that is not in `g`. Everything passes through
`PYTEST_ADDOPTS`, which the harness's pytest subprocesses already inherit; nodeids are quoted for
`shlex`.

**Initial plan** (the estimate per lane assumes one process per core; ⅔ rule: estimate ≤ ⅔ of
`time`; revised after the first run, Q4):

| Suite | shards | groups | group_select | lanes × CPUs | mem/lane | time | Estimate per job |
|---|---|---|---|---|---|---|---|
| placement | 2 | 6 | `tests/calibration/test_placement_fdr_sim.py,tests/calibration/attacks/test_attack_readthrough_recoded.py,tests/calibration/test_readthrough_` | 8 × 2 | 4000 | 06:00:00 | ≈ 2.5–3.5 h (attack5 → 1.7 h per shard; ≥ 15 group tests ≈ 10–15 h / 6) |
| gene_fdr | 0 | 2 | `*` | 2 × 2 | 4000 | 04:00:00 | 1.9 h |
| detectability | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.25 h |
| regions | 0 | 5 | `*` | 5 × 2 | 4000 | 06:00:00 | ≈ 2.4 h (largest test 1.2 h) |
| autocorrelation | 0 | 1 | `*` | 1 × 4 | 4000 | 04:00:00 | 1.4 h |
| islands | 0 | 1 | `*` | 1 × 2 | 4000 | 04:00:00 | 1.8 h |
| discordance | 0 | 2 | `*` | 2 × 2 | 16000 | 04:00:00 | 1.2 h |
| missed_operons | 0 | 1 | `*` | 1 × 2 | 8000 | 02:00:00 | 0.3 h |
| ora | 0 | 2 | `*` | 2 × 4 | 4000 | 04:00:00 | 1.7 h (53-min test) |
| gsea | 0 | 2 | `*` | 2 × 2 | 4000 | 04:00:00 | 1.7 h |
| determinants | 0 | 5 | `*` | 5 × 2 | 4000 | 06:00:00 | ≈ 2.3 h (largest test ≈ 1.7 h) |
| novel_orfs | 6 | 0 | `-` | 6 × 2 | 4000 | 06:00:00 | ≈ 2.2 h shards + ≈ 0.7 h block 2 on 6 lanes (PR 2) |
| denovo_precision | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.25 h |
| variants | 3 | 0 | `-` | 3 × 2 | 4000 | 04:00:00 | 0.75 h shards + ≈ 0.6 h block 2 on 3 lanes |
| compare | 0 | 1 | `*` | 1 × 2 | 4000 | 02:00:00 | 0.45 h |

Total 41 lanes, **88 CPUs** + 1 for the aggregator. The longest `time` is 6 h. The run's expected
wall time is about 3–3.5 h, set by placement and novel_orfs. A one-lane `*` row is exactly
today's unsharded task. The variants row replaces the in-task `n=2` concurrent shards of
`calibration_full.sbatch` for launcher runs only.

Each lane exports `SLURM_CPUS_PER_TASK=<cpus_per_lane>` so that `cpu_share()` (determinants'
`_threads`, ora's prefetch workers) sees the lane's share of the job, not the whole job (5 lanes
× 8 LightGBM threads on 5 cores would oversubscribe the cores).

### `scripts/ci/calibration_lanes.py` (new; pure functions plus a thin CLI)

```text
read_plan(path, suites) -> list[PlanRow]                  # validates; raises PlanError
read_durations(path) -> dict[str, float]
lanes(row, nodeids, durations) -> list[Lane]               # Lane(kind, index, n, deselect, part)
fold_exit(lane_exits, merge_exit | None) -> int
record_parts(lanes, merged) -> list[str]
durations_from_junit(paths) -> dict[str, float]            # planner: rebuild the durations TSV
```

Plan validation (`PlanError`, exit 2 from the CLI): the suites are not exactly the `replicates_full`
keys in order; a count is negative or not an integer; `shards + groups == 0`; `*` with `shards > 0`;
`-` with `groups > 0`; odd `cpus_per_lane`; `time` unparsable or over 8 h; a prefix list given with
`groups == 0`. The checks that need the collection: a nodeid matched by two lanes; a test outside
every lane (impossible by construction; assert it); a group lane with no test (dropped, and the
record says so: no part for it); with `shards > 0`, no test left for the shard lanes (a `PlanError`:
use `*`).

**Exit folding** (the job's exit, as `calibration_full.sbatch` folds task 13 today): if any lane or
the merge has an exit ∉ {0, 1}, the result is the first such exit in lane order, merge last (a
crash, a timeout `143` or a merge error `2`; `calibration_rows.py` fails such suites closed, #349
B1). Otherwise it is 1 if any lane or the merge returned 1, else 0.

**Record** (`task-<i>.json`, written once by the job's EXIT/TERM/USR1 traps, same keys as
M11-03b): `sharded` = `shards > 0`; `parts` = `task-<i>.shard<k>.junit.xml` for every shard lane,
`task-<i>.group<g>.junit.xml` for every non-empty group lane, and `task-<i>.merge.junit.xml` when
the merge ran block 2 (as today); `exit` = the fold; `args` = `[]`. Extra keys are allowed and
ignored by the rows script: `lanes` (`[{"kind", "index", "exit", "tests"}]`) and `plan` (the plan
row as a dict).

### `scripts/slurm/calibration_suite.sbatch` (new)

Takes `GIAP_SUITE`, `GIAP_WORKTREE` and `GIAP_CALIBRATION_LOG_DIR` from the environment. The
resources (`--cpus-per-task`, `--mem`, `--time`, `--job-name=giap-cal-<suite>`) come from the
launcher's `sbatch` options, so the script carries no `#SBATCH --array`.

1. Same set-up as `calibration_full.sbatch`: `farmshare_env.sh`, its own `TMPDIR`, the worktree
   `PYTHONPATH`, `uv run --no-sync python` (never a plain `uv run`).
2. Collect, then compute the lanes with `calibration_lanes.py`. A `PlanError` or collection error
   ends the job with exit 2 and a record.
3. Start every lane in the background on its own core (`taskset` is not required; Slurm's
   cgroup holds the job). Each lane gets its own `PYTEST_ADDOPTS`: `--junitxml=…<part>`, the
   deselect list and `SLURM_CPUS_PER_TASK`. Group lanes get
   `GIAP_CALIBRATION_LOG_DIR=<log>/lanes/<suite>-group<g>`. The suite's `<suite>.jsonl` and
   `<suite>.shards.jsonl` then hold only shard records, which is what `merge_shards` and
   `shard_exit_verdict` expect. Shard lanes share `<log>`.
4. Wait for all lanes. If `shards > 0`, run `--merge` with the shard lanes' deselect list (and,
   after PR 2, `GIAP_CALIBRATION_CONFIRM_LANES=<shards>`).
5. **Time limit:** add `#SBATCH --signal=B:USR1@600`. On USR1 the job interrupts every running
   lane so that pytest writes its JUnit with the tests finished so far, waits ≤ 120 s, then writes
   the record with `exit 143`. Two traps to avoid, which the implementer must test with stubs:
   bash starts background commands with SIGINT ignored when job control is off (use `set -m` or
   an equivalent), and `subprocess.run` in `_criteria.py` kills its pytest child 0.25 s after its
   own KeyboardInterrupt (so signal the pytest process, not the wrapper). The record is still
   written once (as in M11-03b).

`calibration_full.sbatch` and `calibration_sharded.sh` stay byte-unchanged; they remain the tools
for hand runs and for re-running one suite.

### `calibration_full_launch.sh` (changed)

1. Unchanged: argument refusal, `git fetch`, resolving the SHA, the detached worktree under flock,
   a fresh log directory, `GIAP_NOW` / `GIAP_NIGHTLY_AT`, and the fail-closed `squeue`.
2. Validate the plan of the **worktree** (`calibration_lanes.py check`); a `PlanError` gives exit 2
   and nothing submitted.
3. **R14 jobs:** refuse when `used_jobs + 16 > 20` (15 suite jobs + aggregator; same rule as today).
4. **R14 CPUs (replaces the throttle):** each suite is one job, so no array throttle exists. Refuse
   when `used_cpus + Σ_suites lanes × cpus_per_lane + 1 > 128`; the message names the numbers.
   With the initial plan this means `used_cpus ≤ 39`. Q3 covers what to do on a busy queue.
5. **Nightly window:** refuse when `now + max_suites(time) + 20 min` is at or after the next nightly
   start (same python3 arithmetic as today). With the initial plan that is now + 6 h 20 min, so the
   last launch is at 19:10 for a 01:30 nightly.
6. Submit the 15 suite jobs in task order. If any `sbatch` fails, `scancel` **exactly** the jobs
   already submitted, by ID, and exit 1. Then submit the aggregator (`giap-calibration-rows`, 1 CPU,
   `--dependency=afterany:<j0>:<j1>:…:<j14>`, `--job <j0>`, `--expect 15`, otherwise as today). If
   that fails, `scancel` all 15 jobs and exit 1.
7. Print `jobs=<j0>,…,<j14> rows=<id> log=<dir> cpus=<Σ>`.

`calibration_rows.py` is byte-unchanged. The row's job field is `<j0>_calibration` and its text
says `array <j0>` (Q5).

### PR 2 (stats): block 2 on the shard lanes (`tests/calibration/_criteria.py`)

- `calibrate` under `GIAP_CALIBRATION_CONFIRM` with `GIAP_CALIBRATION_REPLICATES` set runs
  `replicate_range(R)` of block 2 for a listed check, as a block-2 `shard` record. Today it runs
  `range(R)`. Without the shard variable, behaviour is unchanged.
- `_run_block2` with `GIAP_CALIBRATION_CONFIRM_LANES=L` (`L ≥ 2`) starts `L` pytest processes at
  once with `CONFIRM` plus `SHARD_ENV=<shard_replicates(R, L, k)>@k/L` and `BLOCK_ENV=2`. It
  waits for all of them, folds their exits with `fold_shard_exits`, then merges again. Unset or
  `L = 1`: unchanged.
- Each block-2 process writes its own JUnit (`task-<i>.merge<k>.junit.xml`); the record lists the
  ones that exist.
- No change to `check`, the criteria, R, α, the bounds, `_verdict`, the xfail logic (ADR-0057) or
  the seeds.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | The initial plan; `replicates_full` keys unchanged | `read_plan` returns 15 rows in key order; total lanes 41, CPUs 88; `max(time)` = 6 h |
| 2 | determinants row (`0 5 * 2 …`), collection = the 19 nodeids, durations: `[full-autocorr_pos-lightgbm]` 8364 s, `[full-detect-lightgbm]` 7782 s, the rest from the TSV | 5 group lanes; each nodeid in exactly one lane; the first two (largest) go to lanes 0 and 1; each lane's `deselect` is the other 4 lanes' tests; no shard lane, no merge; record `sharded: false`, parts `task-10.group0..4.junit.xml` |
| 3 | placement row (`2 6 <prefixes>`) | shard lanes 0/2 and 1/2 deselect every test matching the prefixes (including both `test_attack_d3_recordings[k12]`, `[nissle]`); 6 group lanes split only those tests; `attack5_reversed_stretches…full` runs in both shard lanes (replicates 0–149, 150–299); merge deselects the group tests; record `sharded: true` |
| 4 | `fold_exit([0, 0, 1], 0)` / `([0, 143, 1], 0)` / `([0, 0], 2)` / `([0, 0], None)` | 1 / 143 / 2 / 0 |
| 5 | Launcher at 10:00, queue = `giap-nightly-launch` (1 CPU, pending) | 15 `sbatch` calls, `giap-cal-placement` with `--cpus-per-task=16 --mem=32000 --time=06:00:00`, …; aggregator `--dependency=afterany:<15 IDs>`; prints `cpus=88`; jobs 1 + 16 = 17 ≤ 20 |
| 6 | Launcher with 50 CPUs of other `giap-*` jobs queued | 50 + 88 + 1 = 139 > 128: exit 2, no `sbatch`, message `R14: … 50 CPUs queued + 89 > 128` |
| 7 | Launcher with `GIAP_NOW` = 18:00 / 19:15 local, nightly 01:30 | 18:00: end 00:20, submits. 19:15: end 01:35 ≥ 01:30, exit 2, nothing submitted |
| 8 | The 9th `sbatch` (suite `ora`) fails | `scancel` of exactly the 8 IDs already submitted; no aggregator; exit 1 |
| 9 | A lane is still running at limit − 10 min (stub pytest sleeps) | USR1 → the lanes are interrupted, the finished lanes' JUnit plus the interrupted lane's partial JUnit exist; `task-<i>.json` `exit: 143`, written once; the rows fail closed with `(<suite> exit 143)` |
| 10 | A plan prefix matches no collected test (a renamed file) / the same nodeid in two groups (a hand-edited lane list) | the first: those tests fall to the shard lanes (or, under `*`, nothing changes), with a warning printed; the second: `PlanError`, exit 2 before any pytest, record exit 2 |

## Edge cases

- **New test added after the plan:** with `*` it gets `unknown_min` and is placed by LPT; with a
  prefix list it goes to the shard lanes unless a prefix matches; with `-` it goes to the shard
  lanes. No test is ever dropped (asserted after the assignment).
- **Fewer group tests than lanes:** the empty lanes are dropped (their CPUs sit idle; the record's
  `lanes` names them with `tests: 0`). **R < N** for a check in shard mode: an empty shard logs
  nothing (existing harness behaviour).
- **Ties in LPT:** equal durations are ordered by nodeid, equal loads by the lowest lane index, so
  the assignment is a pure function of (plan, durations, collection).
- **Strict xfail in a group lane** (placement contaminants `flat_c200`): runs as unsharded, block 2
  inline, pytest reports XFAIL in that lane's JUnit. In a shard lane, `--merge` judges it
  (ADR-0057). A test deselected from every shard lane writes no `deferred` outcome, so `--merge`
  cannot mistake it for an XPASS. The test-author checks this on `_shard_fixture_suite.py`.
- **A lane crashes or is OOM-killed:** its exit is ∉ {0, 1} or its JUnit is missing. The fold
  gives that exit; `calibration_rows.py` fails the suite closed (`junit missing` / `exit …`).
- **The whole job is SIGKILLed (node failure):** no record; the rows say `no record`, as today.
- **The merge times out** (an unexpected block-1 failure of a large check, before PR 2): USR1 path,
  `exit 143`, fail closed. A release then needs a rerun, never a pass.
- **Empty input:** an empty collection for a suite is a collection error: exit 2. An empty plan or
  durations file is a `PlanError` (plan) or means every test uses `unknown_min` (durations).
- **Circular/linear replicons, minus strand:** not applicable; no genomic computation changes.

## Acceptance criteria

1. `read_plan` accepts the committed plan and rejects each rule under "Plan validation" with a
   message naming the row and column (one parametrised test per rule).
2. `lanes()` property test: over random collections (0–80 nodeids), random durations and random
   valid rows, every nodeid is in exactly one lane (or in every shard lane, never both); the
   `deselect` lists are complements; the output does not depend on the input order; it matches
   examples 2 and 3 exactly.
3. `fold_exit` matches example 4 and the rule over all combinations of {0, 1, 2, 143} for ≤ 4 lanes
   plus the merge.
4. `calibration_suite.sbatch` with the stub `uv` of `test_m11_03b_sbatch.py`: example 2 and
   example 3 run every lane with its own `--junitxml`, deselect list and `SLURM_CPUS_PER_TASK`; the
   merge runs after every shard lane and only when `shards > 0`; one `task-<i>.json` whose
   `parts`, `sharded` and `exit` follow the rules; `calibration_rows.suite_result` on that record
   and the stub JUnit gives the expected pass/fail.
5. Example 9 with a stub pytest that traps SIGINT and writes JUnit: the record says `exit 143` and
   lists the partial JUnit; a stub that ignores SIGINT still gets a record within the 120 s grace
   (test with a shortened grace through an environment variable).
6. Launcher, in the `slurm_sandbox` fakes: examples 5–8; `squeue` before the first `sbatch`; a plan
   error refuses before any `sbatch`; every job name starts with `giap-`; the aggregator depends on
   exactly the 15 IDs; harness arguments are still refused.
7. PR 2: the shard-equivalence test (`test_criteria_shard_equivalence.py` pattern, on
   `_shard_fixture_suite.py`) shows that `--merge` with `CONFIRM_LANES` = 1, 2 and 3 gives the same
   verdict records (status, observed, R, block) as the unsharded run, including a check that
   fails block 1 and a strict-xfail check. `tolerances.yaml` is byte-unchanged.
8. `calibration_full.sbatch`, `calibration_sharded.sh`, `calibration_rows.py` and `readiness.py`
   are byte-unchanged (diff in review). `ruff`, `mypy` (new Python), `shellcheck` (if available)
   pass.
9. **The run itself (lead, after PR 4):** one launch on main finishes with every suite job
   `COMPLETED`, each within ⅔ of its `time`. Both rows are then written, and R3/R4 move from
   `fail` to whatever the suites say. If a job exceeds ⅔, the planner rebuilds
   `calibration_durations.tsv` from that run's JUnit and a data-only plan PR re-balances the
   lanes before the next launch.

## PR split (≤ ≈ 400 code lines each)

| PR | Content | Code | Tests | Roles | Depends |
|---|---|---|---|---|---|
| 1 | `scripts/ci/calibration_lanes.py` (pure functions + `check`/`lanes`/`durations` CLI), `scripts/slurm/calibration_plan.tsv`, `scripts/slurm/calibration_durations.tsv` (from the 1894441 JUnit and the probe logs) | ≈ 230 | ≈ 250 new (`tests/unit/test_m11_03c_lanes.py`) | implementer, test-author | none |
| 2 | block 2 on lanes in `tests/calibration/_criteria.py` (`stats`) | ≈ 50 | ≈ 120 new (`tests/calibration/test_criteria_confirm_lanes.py`) | implementer, test-author, statistician review | none (parallel with 1) |
| 3 | `scripts/slurm/calibration_suite.sbatch` (lanes, merge, record, USR1) | ≈ 160 | ≈ 150 new (`tests/unit/test_m11_03c_suite_sbatch.py`, reusing the stub `uv`) | implementer, test-author | 1 |
| 4 | `calibration_full_launch.sh` (per-suite jobs, CPU fit, window by `max(time)`, cancel on failure), `docs/release.md` step 7 note, changelog fragment | ≈ 120 | ≈ 100 new (`tests/unit/test_m11_03c_launcher.py`) **plus edits to `tests/unit/test_m11_03b_launcher.py`** (`0-14%15`, throttle and the 4-CPU arithmetic are gone) | implementer, test-author (new file); **lead** edits the existing file (role guard) | 3 (2 soft: without it, set novel_orfs `time` to `08:00:00`) |

`test_calibration_full_array.py`, `test_calibration_sharded_script.py`, `test_m11_03b_sbatch.py` and
`test_calibration_full_worktree.py` are not touched (their scripts are unchanged). After PR 4 the
lead launches the run (criterion 9) and records the job IDs in PROGRESS.md.

## Questions for the lead (reversible defaults in use)

- **Q1. Does a throttled pending array count as one R14 job or as N?** Default: N, as in #349 B3,
  which is why the design uses one job per suite with lanes. If it counts as one, the lead's
  sketch (one array task per (suite, shard), `%T`, a merge job per suite) becomes possible:
  about 35 tasks at 2 CPUs. It is simpler to read but has more cross-job dependencies, and it does
  not fix findings 4–5 by itself.
- **Q2. Lanes run concurrently inside one job.** The brief asked for "not concurrent shards in one
  task". The measurement shows the cost comes from processes sharing a core, not from sharing a
  job: 2 processes on 2 cores ran at full speed, 8 on 4 cores at 1/2.46. Default: lanes, one per
  core.
- **Q3. Busy queue.** Today the throttle lets a launch start with T ≥ 4 when other jobs use CPUs.
  Default: refuse when the 89 CPUs do not fit (all or nothing; R14 always holds). Alternative: a
  second wave held by `--dependency=afterany` on the first wave's jobs (more code; the window grows
  to 2 × max(time)).
- **Q4. First-run budget.** The initial plan rests on estimates for placement, regions,
  determinants and novel_orfs: every probe job hit its 3 h limit, and the estimates have ±30 %
  uncertainty. Default: launch once after PR 4 (89 CPUs × ≤ 6 h), then re-plan from its JUnit
  (data-only PR). Alternative: a dedicated measurement run before PR 4.
- **Q5. Row label with 15 job IDs.** Default: `<j0>_calibration` and `array <j0>` (the IDs
  increase, so readiness's `newest` still works; `calibration_rows.py` stays unchanged).
  Alternative: change the summary to `jobs <j0>-<j14>` (touches the rows format and its tests).
- **Q6. The plan lives outside `CALIB_PATHS`.** Default: yes (`scripts/slurm/`), because it changes
  which process runs a replicate, not any result. A stricter reading would add the two TSVs to
  `CALIB_PATHS`, so every re-balance would invalidate the earlier rows.
- **Q7. `calibration_full.sbatch` task 13's in-task `n=2`.** Default: kept for hand runs (its tests
  pin it); the launcher no longer uses it. Retire it in a later clean-up?
- **Q8. Why the SIGTERM trap wrote nothing at TIMEOUT.** The USR1 path (PR 3) makes the record
  independent of the cause. Should someone still find the root cause? A probe: a 2-minute
  `--time` job with a sleeping stub.
- Not in scope, for the backlog: the full tier's block-2 cost for strict xfails (novel_orfs `gc`
  7.3 h single-process, variants g1) disappears once ADR-0049 and ADR-0050 are decided (Khoa). The
  14 `slow`-but-not-`calibration_full` tests (M11-03b Q2) are still not scheduled.
