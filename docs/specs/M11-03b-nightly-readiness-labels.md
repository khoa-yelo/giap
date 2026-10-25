# M11-03b: `calibration_full` and `calibration_attacks` rows for the readiness checker

Milestone M11 · **pass 1 (MVP follow-up of M11-03)** · depends on: M11-03 (merged:
`scripts/release/readiness.py`), M7-13 (the 15-suite `calibration_full.sbatch` array, task-13
sharding) · estimate ≈ 300 code lines + ≈ 300 test lines · not `stats` (no new statistics, no gate,
tolerance, α or R changes) · roles: implementer (`scripts/ci/calibration_rows.py`,
`scripts/slurm/calibration_full.sbatch`, `scripts/slurm/calibration_full_launch.sh`), test-author
(`tests/unit/test_m11_03b_*.py`), docs (the runbook note in `docs/release.md`, M11-05) · Khoa: none
blocking (Q1–Q6 are lead questions with reversible defaults).

## Spec sections implemented

- §13.2 "Every statistical method passes its calibration suite (§8.4), and every attack null (§8.5)
  stays at level": the evidence rows that M11-03's R3 and R4 read.
- §8.4 "The full suite runs nightly on FarmShare as a Slurm job array, one suite per task", the
  "Nightly (full)" R column. §8.5 "the method passes when it stays at level".
- §11.2 Calibration layer ("fast subset per PR, full nightly"). §12 "Every … FarmShare suite ends by
  posting a compact summary" (here: two `nightly.tsv` rows).
- R3: no threshold, tolerance, α or R changes; R14: ≤ 128 CPUs and ≤ 20 `giap-*` jobs (running +
  pending), checked before every `sbatch`.

## Scope (after the #344 review correction)

What `readiness.py` reads, checked against the merged code:

| Row | Label(s) it reads | Log | Writer today | This issue |
|---|---|---|---|---|
| R3 | `calibration_nightly` **and** `calibration_full` | `nightly.tsv` | `nightly.sbatch` task 2 writes `calibration_nightly: <pytest tail>` (`pyt calibration_nightly …`, line 178); rows observed 2026-10-20…24 | `calibration_nightly`: **unchanged**. `calibration_full`: **new** |
| R4 | `calibration_attacks` | `nightly.tsv` | none | **new** |
| R5 (TE8) | `k12_mg1655` | `weekly.tsv` | `weekly.sbatch` task 0 (`report "$rc" "$S: …"`, `S=k12_mg1655`) | unchanged (matches) |
| R6 | `integration` | `nightly.tsv` | `nightly.sbatch` task 0 (`pyt integration …`) | unchanged (matches) |

How readiness decides (unchanged, quoted from `readiness.py`):

- **Parse** (`parse_log`): tab-split into at most 5 fields; a line is malformed (skipped, noted) if
  it has < 5 fields, field 4 is not an integer, or field 3 is shorter than 7 characters.
- **Label** (`LogRow.label`): `summary.split(":", 1)[0].strip()`, compared with `==`.
- **Commit** (`code_equal`, `CALIB_PATHS = ("src", "tests/calibration", "ci")`): the commit must
  resolve in the release checkout, be HEAD or an ancestor, and `git diff --quiet <c> HEAD -- src
  tests/calibration ci` must hold.
- **Newest** (`newest`): max over (date string, list of integers in the job field, file order)
  among rows with that label and an accepted commit.
- **ok** (`LogRow.ok`): `exit == 0` and every match of `(\d+) (?:failed|errors?)\b` in the summary
  is 0.

So the readiness side needs nothing; the gap is two writers. The runbook's sentence "`nightly.sbatch`
labels its rows `task <n>`" (`docs/release.md` on `m11/05-runbook`, Known gaps) is inaccurate: only
the fallback rows of the EXIT/TERM traps (`task <n>: ended during …`, `task <n>: killed by SIGTERM
…`) and the MODE=record skip row use `task <n>:`. See Q5 for the one real consequence.

## Definition

### Which job writes which label

| Label | Written by | When | Covers |
|---|---|---|---|
| `calibration_nightly` | `scripts/slurm/nightly.sbatch` task 2 (existing, unchanged) | every night (scrontab) or `sbatch --array=2 …` | the whole fast tier, `tests -m "not slow and not calibration_full and not e2e"` with the coverage gate; this includes the 316 fast-tier attack tests (`calibration`, `calibration_nightly` markers) |
| `calibration_full` | **new** aggregator job `giap-calibration-rows`, submitted by **new** `scripts/slurm/calibration_full_launch.sh` with `--dependency=afterany:<array>` of `calibration_full.sbatch` | on demand (release step 7) and whenever the lead launches the full tier; not in the nightly scrontab (R14, below) | all 15 `replicates_full` suites (422 tests marked `calibration_full` today), at full R |
| `calibration_attacks` | the same aggregator job, same run | as above | the attack tests among them: testcases under `tests/calibration/attacks/` (230 tests in 14 suites today; `denovo_precision` has none) |

**Attack tier (decision, Q1):** `calibration_attacks` is the **full tier's** attack tests (the §8.4
"Nightly (full)" R). The fast-tier attack tests are not lost: they run in `nightly.sbatch` task 2,
so a failing fast-tier attack turns `calibration_nightly` red and fails R3. Every attack test in
both tiers therefore gates the release; R4's own row shows the full-R attack outcome.

**Why not inside the nightly array (R14):** the nightly holds 4 × 16 + 1 CPUs and 5 jobs
(+ the 1-CPU launcher). The full array is 15 × 4 CPUs. Together: 125–126 CPUs and 21–22 jobs, over
the 20-job cap, and with no room for anything else. The nightly launcher does not check the budget
today, so it can already overrun when other `giap-*` jobs run. This issue does not add load to it.

### `calibration_full.sbatch` (changed)

Per array task, in addition to what it does today:

1. **JUnit per pytest run.** Before each `harness` call it sets
   `PYTEST_ADDOPTS="${PYTEST_ADDOPTS:+$PYTEST_ADDOPTS }--junitxml=$GIAP_CALIBRATION_LOG_DIR/task-<i>.<part>.junit.xml"`
   with `<part>` = `run` (unsharded), `shard<k>` (task-13 shard k), `merge` (the `--merge` call;
   its block-2 pytest subprocess inherits the environment). Each concurrent shard gets its own file.
2. **Task record.** An EXIT trap (and a TERM trap: exit 143, record written once) writes
   `$GIAP_CALIBRATION_LOG_DIR/task-<i>.json`:

   ```json
   {"task": 3, "suite": "regions", "job": "1900000", "commit": "<40-hex HEAD of the checkout>",
    "dirty": false, "exit": 0, "args": [], "sharded": false,
    "parts": ["task-3.run.junit.xml"], "date": "2026-11-02"}
   ```

   `suite` is the `replicates_full` key of task i (same YAML lookup as `calibration_sharded.sh`);
   `dirty` is true when `git status --porcelain -- src tests/calibration tests/sim tests/helpers ci`
   is non-empty; `args` are the harness arguments (`"$@"`); `exit` is the task's own exit code as
   today. The record is written with one `write` to a temporary name and renamed (no partial JSON).
3. Nothing else changes: exit codes, task-13 sharding, `GIAP_WORKTREE`, `--array=0-14`,
   `--cpus-per-task=4`, `--time=04:00:00`. A hand run (`sbatch calibration_full.sbatch`) writes
   records but **no** `nightly.tsv` row; the aggregator can be run on its log directory by hand.

### `scripts/slurm/calibration_full_launch.sh [REF]` (new)

1. Resolve `REF` (default `origin/main`) after `git fetch -q origin` to a 40-hex SHA; refuse any
   further argument (harness options such as `--fast`, `--shard`, `--dry-run` would make the
   result something other than the full tier: exit 2, nothing submitted).
2. Check out SHA into a detached worktree `$GIAP_ROOT/ci/calibration-<sha12>` (flock on
   `$GIAP_ROOT/ci/.git.lock`, as the nightly does); fresh log directory
   `GIAP_CALIBRATION_LOG_DIR=$GIAP_ROOT/logs/calibration/full-<sha12>-<YYYYmmdd-HHMMSS>` (refuse if
   it exists).
3. **R14.** Sum running + pending `giap-*` jobs and CPUs of `$USER` with `squeue -h -r -u "$USER"
   -t RUNNING,PENDING -o '%C %j'` (`-r`: one line per array task); fail closed (exit 2) when squeue
   fails or prints a non-integer CPU count. Reserve the nightly footprint, 65 CPUs and 5 jobs, at
   all times (Q3). Throttle `T = min(15, ⌊(128 − 65 − used_cpus − 1) / 4⌋, 20 − 5 − used_jobs −
   1)` (the `− 1` is the aggregator job). Refuse (exit 2, message names the numbers) when `T < 4`.
   With nothing else queued, T = 14.
4. `sbatch --parsable --array=0-14%T --export=ALL,GIAP_WORKTREE=<wt>,GIAP_CALIBRATION_LOG_DIR=<dir>
   scripts/slurm/calibration_full.sbatch` (log dir via the environment, as `calibration_sharded.sh`
   does, so a comma is safe).
5. `sbatch --parsable --job-name=giap-calibration-rows --cpus-per-task=1 --mem=2G --time=00:20:00
   --dependency=afterany:<array> … python scripts/ci/calibration_rows.py --log-dir <dir> --job
   <array> --commit <sha> --expect 15 --tsv $GIAP_ROOT/logs/nightly.tsv`. If this submission fails,
   `scancel` the array just submitted (GIAP job by ID only) and exit 1: an array without its
   aggregator produces no evidence.
6. Print `array=<id> rows=<id> log=<dir> throttle=<T>`.

### `scripts/ci/calibration_rows.py` (new): the row writer

Pure part (no Slurm, no git, no clock; the test-author's target):

```text
suite_result(record: dict | None, junits: list[ET.Element], task: int, expect_commit: str)
    -> SuiteResult(task, suite, exit, sharded, all: Counts, attacks: Counts, problems: list[str])
full_rows(results: list[SuiteResult], job: str, commit: str, date: str) -> list[str]   # 2 lines
```

`Counts` = passed, failed, errors, skipped, xfailed, xpassed, from `<testcase>` children of the
JUnit files (xfail = `<skipped type="pytest.xfail">`; strict XPASS appears as a failure). A
testcase is an **attack** iff its `classname` starts with `tests.calibration.attacks.`.

A suite **fails** `calibration_full` when any of: record missing (problem `no record` plus the
sacct state if `--sacct-json` gives one); record unreadable; `record.commit != expect_commit`;
`dirty`; `args != []`; `exit != 0`; a JUnit file listed in `parts` missing or unparsable; any
failed or error testcase; zero testcases collected.

A suite **fails** `calibration_attacks` when any of: it fails `calibration_full` for a reason other
than a *non-attack* testcase failure, i.e. a missing/unreadable record or JUnit, commit mismatch,
dirty, args, or zero testcases (fail closed: the attack outcome is unknown); any attack testcase
failed or errored; it is unsharded and `exit` ∉ {0, 1} (crash, interrupt, collection error: JUnit
may be partial); it is unsharded, `exit == 1`, and no testcase at all failed (unattributable); it is
**sharded** (task 13), `exit != 0`, and it holds ≥ 1 attack testcase (shard tests skip without a
verdict; the merge exit cannot be attributed to a test). A suite with zero attack testcases and
none of the fail-closed problems is not an attack suite and is not counted in the attack row.

Row format (5 tab-separated fields, one `\n`; same layout as `nightly.sbatch` `tsv_line`):

```text
<date YYYY-MM-DD>\t<array>_calibration\t<40-hex commit>\t<exit>\t<summary>
calibration_full: <n> suites, <f> failed suites[ (<suite> <reason>, …)]; <N> tests: <a> failed, <b> errors, <c> passed, <d> skipped, <e> xfailed, <g> xpassed; array <array>
calibration_attacks: <m> suites, <f'> failed suites[ (…)]; <N'> attack tests: <a'> failed, <b'> errors, <c'> passed, <d'> skipped, <e'> xfailed, <g'> xpassed; array <array>
```

- `exit` is 0 when the row passes (no failed suite **and** `expect` suites present), else 1.
- A failing row always has `failed suites ≥ 1`, so `readiness.ok()` is false by exit **and** by
  count; a passing row has every `failed`/`errors` count 0. Reasons never contain
  `<digits> failed`/`<digits> error(s)` with a non-zero number other than the counts above (write
  `exit 1`, `no record`, `commit mismatch`, `dirty`, `args --fast`, `merge exit 2`, `TIMEOUT`).
- `date` is the aggregator's run date; `commit` is `--commit` (all records must equal it).
- Both rows are appended in **one** `write` call under `flock` on `<tsv>.lock`; a marker
  `<log-dir>/rows.done` (written after the append) makes a second run on the same directory exit 0
  without writing (idempotent).
- CLI exit: 0 rows written and both pass; 1 rows written, one fails; 2 usage error / unreadable log
  dir (still writes failing rows if the TSV is writable, so readiness sees a red run, not silence).

## Examples

Fixture: log dir with records `task-0.json … task-14.json`, JUnit files, commit `C` = 40-hex, job
`1900000`, date `2026-11-02`. Rows shown without the leading `2026-11-02\t1900000_calibration\tC\t`.

| # | Input | Expected rows (exit, summary) |
|---|---|---|
| 1 | all 15 records exit 0, commit C; 422 testcases: 400 passed, 20 skipped, 2 xfailed; 230 attack: 220 passed, 10 skipped | `0  calibration_full: 15 suites, 0 failed suites; 422 tests: 0 failed, 0 errors, 400 passed, 20 skipped, 2 xfailed, 0 xpassed; array 1900000` and `0  calibration_attacks: 14 suites, 0 failed suites; 230 attack tests: 0 failed, 0 errors, 220 passed, 10 skipped, 0 xfailed, 0 xpassed; array 1900000`; readiness R3 `pass` (with a passing `calibration_nightly` row), R4 `pass` |
| 2 | as 1, but `regions` (task 3) exit 1, one **non-attack** testcase failed | full `1  … 1 failed suites (regions exit 1); 422 tests: 1 failed, …`; attacks `0` (the failure is attributed outside attacks); R3 `fail`, R4 `pass` |
| 3 | as 1, but `islands` (task 5) one **attack** testcase failed, exit 1 | both rows exit 1, `1 failed suites (islands exit 1)`, attack `1 failed`; R3 `fail`, R4 `fail` |
| 4 | task 6 (`discordance`) record missing (SIGKILL / OOM), `--sacct-json` says `OUT_OF_MEMORY` | both rows `1  … 1 failed suites (discordance no record, OUT_OF_MEMORY)`; R3/R4 `fail` |
| 5 | task 13 (`variants`, sharded) exit 1 from `--merge` (confirmed fail), JUnit has 1 attack testcase skipped in every shard | both rows exit 1 naming `variants merge exit 1`; R4 `fail` (fail closed) |
| 6 | task 13 exit 1, its JUnit has no attack testcase | full exit 1; attacks exit 0 if all others pass |
| 7 | one record `args: ["--fast"]` (someone bypassed the launcher) | both rows exit 1, `(<suite> args --fast)`; never a full-tier pass from fast R |
| 8 | record commit ≠ `--commit` | both rows exit 1, `(<suite> commit mismatch)` |
| 9 | aggregator run twice on the same log dir | second run writes nothing, exit 0 (`rows.done`) |
| 10 | launcher with 70 CPUs / 6 jobs of other `giap-*` work queued | T = min(15, ⌊(128−65−70−1)/4⌋ < 0, …) → refuse, exit 2, no `sbatch` call |

## Edge cases

- **Array task that fails** (exit ≠ 0): its suite is named in `failed suites`; the other suites'
  counts are still reported. Both rows are written; readiness takes `exit 1`.
- **Task killed** by the time limit: the TERM trap writes the record with exit 143 → failed.
  SIGKILL/OOM/node failure: no record → `no record` (+ sacct state) → failed. Aggregator itself
  killed: no row → R3/R4 stay `pending` (never a pass).
- **Re-run** of the whole tier on the same commit: new array ID and log directory, two new rows;
  readiness's `newest` (date, then job ID integers, then file order) picks the re-run, so a later
  failure beats an earlier pass and vice versa. Re-running a **single** failed suite
  (`sbatch --array=3 …` or `calibration_sharded.sh`) produces no aggregate row by design (Q4).
- **Same date, two runs**: job ID decides (`1900500_calibration` > `1900000_calibration`).
- **Concurrent nightly**: different job IDs; `post_summary.py` only reads rows of its own array job
  (`jid != job` skipped), so the new rows never enter a nightly post. The nightly's own `>>` appends
  do not take the lock (unchanged); one-write appends keep lines whole on local disks (NFS caveat
  noted, pre-existing).
- **Dirty or moving checkout**: the launcher pins a detached worktree, so the main checkout's
  fast-forward mid-run cannot mix commits; a hand run from a dirty tree is refused by `dirty`.
- **Commit not in the release history** (a branch run): readiness's `code_equal` rejects it, so it
  never counts; no writer-side check needed.
- **R14**: launcher refuses when the throttle drops below 4; the array never exceeds T concurrent
  tasks; pending tasks count (squeue `-r`); squeue failure fails closed.
- **Empty input**: empty log dir → 15 `no record` suites → both rows exit 1. `--expect` ≠ number of
  `replicates_full` keys in the checkout → usage error, exit 2.
- **Ties at a bound / genomic cases** (strand, circular): not applicable; no statistic is computed.

## Acceptance criteria

1. The rows of examples 1–8, fed through the **merged** `readiness.parse_log` and `_calib`
   (R3 with a passing `calibration_nightly` row, R4), give exactly the R3/R4 statuses in the table
   (round-trip test; no change to `readiness.py`).
2. `suite_result`/`full_rows` are pure (no I/O beyond the parsed inputs): examples 1–8 are unit tests
   on fixture records and hand-written JUnit XML in `tmp_path`.
3. Every row `full_rows` emits has exactly 5 tab fields, a 40-hex commit, an integer exit, and a
   label equal to `calibration_full` or `calibration_attacks`; `exit == 0` iff the summary's
   `failed suites` is 0 and every `(\d+) (failed|errors?)` count is 0 (property-style test over
   random fixture outcomes).
4. Static checks of `calibration_full.sbatch`: `--array=0-14` and `--cpus-per-task=4` unchanged; it
   sets `PYTEST_ADDOPTS` with `--junitxml=` under `$GIAP_CALIBRATION_LOG_DIR`, distinct per shard;
   EXIT and TERM traps write `task-<i>.json`. Existing `tests/unit/test_calibration_full_*.py` still
   pass unchanged.
5. With the stub `uv` of `test_calibration_full_array.py` (writing a fake JUnit), a task exit 1 and a
   SIGTERM each leave one record with that exit (143 for TERM).
6. Static and sandbox checks of `calibration_full_launch.sh` (`tests/helpers/slurm_sandbox.py` fake
   `sbatch`/`squeue`/`scancel`): `squeue` runs before the first `sbatch`; example 10 refuses with no
   `sbatch` call; a squeue failure refuses; with an empty queue the array is `--array=0-14%14`, the
   aggregator has `--cpus-per-task=1`, `--dependency=afterany:<array>`, job names start `giap-`;
   extra arguments are refused; a failed aggregator submission `scancel`s exactly the array ID.
7. Example 9 (idempotent) and the `flock`ed single append are tested on a `tmp_path` TSV.
8. `nightly.sbatch`, `weekly.sbatch`, `post_summary.py` and `readiness.py` are byte-unchanged
   (a test or the diff in review).
9. `ruff check`, `ruff format --check`, `shellcheck` (if available) pass; type hints throughout.

Tests need no Slurm: everything above runs on fixtures, stubs and the existing sandbox.

## Estimate

Code ≈ 300: `calibration_rows.py` ≈ 170 (parse ≈ 50, rules ≈ 60, rows ≈ 30, CLI/lock ≈ 30);
`calibration_full.sbatch` +≈ 40; `calibration_full_launch.sh` ≈ 85; runbook note ≈ 5 lines.
Tests ≈ 300 (row writer + round-trip ≈ 170, sbatch ≈ 50, launcher sandbox ≈ 80). One PR, ≤ 400 code
lines.

## Questions for the lead (reversible defaults in use)

- **Q1. Which tier is `calibration_attacks`?** Default: the full tier's 230 attack tests; the 316
  fast-tier attack tests gate through `calibration_nightly` (R3). Alternative: the aggregator also
  runs the fast-tier attack subset itself so one row covers both (more CPU-hours, no gate change).
  A nightly-side `calibration_attacks` row would compete with this one under `newest` and is not
  proposed.
- **Q2. Calibration tests that no scheduled job runs (real gap in "every … suite").** 14 tests in
  `tests/calibration` are marked `slow` but not `calibration_full` (10 of them attacks: the R-engine
  GSEA/ORA attacks in `attacks/test_attack_gsea_r_engine.py`, `attacks/test_attack_function.py`,
  `attacks/test_attack_sim_function.py`; plus 4 R-engine tests in `test_gsea_r_engine_calibration.py`
  and `test_function_calibration.py`). The nightly excludes `slow`; the full array selects
  `-m calibration_full`. Default: out of scope here, listed in the backlog; R3/R4 pass without
  them. Options: mark them `calibration_full` with a suite (statistician), or a third row.
- **Q3. R14 reservation.** Default: the launcher always reserves the nightly's 65 CPUs / 5 jobs
  (T = 14 on an empty queue). Alternative: reserve only when the run could overlap 01:30–07:30.
  Separately, `nightly_launch.sh` checks no budget at all (pre-existing; not changed here).
- **Q4. Partial re-runs.** Default: only a full 15-suite run yields rows; a single-suite re-run
  (`calibration_sharded.sh`, `--array=<i>`) never does. Alternative: let the aggregator combine the
  newest record per suite across log dirs for one commit (more code, harder to audit).
- **Q5. Fallback rows of nightly task 2 are invisible to R3 (readiness reads them differently).**
  A night whose task 2 is killed writes `task 2: killed by SIGTERM …` (label `task 2`), so an
  *older* passing `calibration_nightly` row of a code-equal commit still decides R3 although the
  newest run died. Default: unchanged, per the #344 correction (the release evidence step submits a
  fresh `--array=2` on R, whose row is newest). Fix if wanted: label task 2's trap rows
  `calibration_nightly: task 2 killed …` (exit ≠ 0), keeping the MODE=record skip row as `task 2:`
  (it has exit 0 and would otherwise read as a pass).
- **Q6. `CALIB_PATHS` misses inputs of the calibration suites.** 64 files in `tests/calibration`
  import `giapsim` (`tests/sim/`) and many use `tests/helpers/`; a change there alters the nulls but
  not `CALIB_PATHS = ("src", "tests/calibration", "ci")`, so an older row stays code-equal. Default:
  readiness unchanged (outside this issue; the release freeze in runbook step 6 covers `src/`,
  `workflow/`, `tests/calibration/` only). Fix if wanted: add `tests/sim` and `tests/helpers` to
  `CALIB_PATHS` (stricter, not wider) and to the runbook freeze list.
- Also for the docs role: `docs/release.md` (M11-05) Known gaps should say "nothing writes
  `calibration_full` or `calibration_attacks`" until this lands, then point step 7 at
  `calibration_full_launch.sh` instead of `sbatch calibration_full.sbatch`.

## Lead decisions (2026-10-25)

- **Q1, Q3, Q4:** the defaults stand.
- **Q2:** the 14 `slow`-but-not-`calibration_full` calibration tests go to the backlog. They are a coverage gap to close before M11-13.
- **Q5, in scope (stricter):** nightly task 2's trap rows (killed or failed) are labelled `calibration_nightly` with their non-zero exit, so a killed run can no longer let an older pass decide R3. The MODE=record skip row stays `task 2` (exit 0, not a calibration result). This amends the acceptance criterion "nightly.sbatch byte-unchanged" to "unchanged apart from the task-2 trap label".
- **Q6, in scope (stricter):** `tests/sim` and `tests/helpers` join readiness.py's `CALIB_PATHS`. They also join the runbook's freeze list, through the M11-05 docs follow-up. This amends "readiness.py byte-unchanged" to "unchanged apart from CALIB_PATHS".
