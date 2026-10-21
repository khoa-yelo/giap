# CI-01: the `calibration_full` tier has no scheduled home — give it a sharded array job

Track: infrastructure (lead-raised 2026-10-20, not a milestone issue) · depends on: nothing

## The gap
`calibration_full` is **deselected everywhere** and therefore never runs except by hand:
- `.github/workflows/ci.yml:79` and `:96` — `-m "not … and not calibration_full and not …"`;
- `scripts/slurm/test_branch.sbatch:11,25` — the same exclusion in the default args;
- `scripts/slurm/nightly.sbatch` task 2 — `-m calibration_nightly` only;
- `scripts/slurm/weekly.sbatch` — end-to-end genome runs, no calibration.

161 of 557 calibration tests carry the mark today (2026-10-20), among them the two ADR-0019 estimate
checks that M3-25 just turned from strict xfail to passing, ADR-0020's three ~65-minute autocorrelation
attack arms, the M5-03 scan full arm, and the M5-02 engine full mode. M5-06's island attacks land next.
Nothing re-runs any of them when main changes, so a regression in a gated statistical claim would be
invisible until someone happened to run the tier by hand.

## Decisions (lead)
1. A new `scripts/slurm/calibration_full.sbatch`, a Slurm **array**, one shard per task, so no task
   approaches a wall-clock limit. Shard by test id with a stable assignment (hash of the node id modulo
   the shard count), not by file, because single tests are now the expensive unit (one arm ≈ 65 min).
2. Size the shards from measured durations, not guesses: the job records per-test durations
   (`--durations=0`) to `$GIAP_ROOT/logs/calibration_full.tsv`, and the shard count is a variable with a
   documented default chosen from the current total. R14: tasks × cpus-per-task ≤ 128 including anything
   else the account is running, so prefer more, smaller tasks over wide ones.
3. One row per shard in `$GIAP_ROOT/logs/calibration_full.tsv` (date, job_task, commit, exit, summary),
   the same shape as `nightly.tsv`, and a summary job with `afterany` that posts once every shard has a
   line — reuse `scripts/ci/post_summary.py` and the `nightly_launch.sh` pattern rather than inventing a
   second mechanism. `NIGHTLY_POST=dry` must suppress posting, as it does for the nightly.
4. **Do not touch the scrontab.** Scheduling needs Khoa's approval; this issue delivers a job that can be
   submitted by hand and is ready to schedule. Record in PROGRESS.md › Questions for Khoa: how often the
   full tier should run (weekly is the obvious candidate) and whether it may post to the nightly issue.
5. Fix the misleading backlog entry, which assumed the tier ran and was merely slow.

## Definition of done
The array job runs the whole `calibration_full` selection in shards, each well inside its time limit, on
current main; a dry-run summary is produced; durations are recorded; no scrontab change.
