# ADR-0052: PR CI under 20 minutes; calibration, integration and coverage move to the nightly
- Status: accepted-by-khoa (2026-10-23)
- Date: 2026-10-23
- Spec sections: §12 (the `ci.yml` and `nightly.sbatch` rows), §14 (feature loop)
- Needs Khoa: yes (Khoa asked for it: "we want 15 mins on all or 20 mins max … cut down on the
  tests, move to nightly", 2026-10-23)

## Context
- **PR CI is the bottleneck of the merge queue.** The single `calibration-and-coverage` job took
  85 min. Splitting it into 6 coverage shards (#273) still left 10–56 min shards. The slowest is
  one calibration test of about 25 min (`test_gene_fdr_one_replicate_through_real_mapping`, Slurm
  1879345 `--durations`).
- **The unsharded test matrix** (ubuntu 3.11, ubuntu 3.13, macOS 3.12) took 12–25 min per job.
- **Runner limits:** a free account runs about 20 jobs at once, 5 of them on macOS.
- **The Slurm fast suite already runs everything:** agents run the whole fast suite, calibration
  included, on Slurm (`scripts/slurm/test_branch.sbatch`) before a PR opens.

## Decision
1. **PR CI** (`ci.yml` on pull_request) runs:
   - lint;
   - guards;
   - nf-test in stub mode;
   - the unit, property, golden, contract and snapshot tests on ubuntu 3.13, in **4 parallel
     shards** (`scripts/ci/shard_tests.py`; `tests/calibration`, `tests/integration` and
     `tests/e2e` excluded).

   The target is every job ≤ 15 min, and the whole run ≤ 20 min.
2. **On pushes to main only:** Python 3.11 and macOS 3.12 (`test-compat`).
3. **Moved to the nightly:** the fast calibration subset and the coverage gate. Nightly task 2
   (`scripts/slurm/nightly.sbatch`) runs the whole fast suite, calibration_nightly included, with
   `--cov-fail-under` from `ci/gates.yaml` `coverage_min`. That value is unchanged (R3). A nightly
   failure becomes a top-priority issue (§12).
4. **Before a PR opens:** agents still run the full fast suite on Slurm, so calibration is checked
   before review, and the nightly is the backstop for anything merged since.

## Consequences
- **Later detection:** a calibration or coverage regression can merge and is found by the next
  nightly, at most about 24 h later, instead of blocking the PR. 3.11 and macOS regressions are
  found after the merge.
- **Gate values unchanged:** `ci/gates.yaml` and the tolerances are untouched; only where the
  checks run moves.
- **Budget:** the nightly's task 2 grows from the `calibration_nightly` tests to the whole fast
  suite: about 35–50 min on 16 CPUs, within its 6 h.
