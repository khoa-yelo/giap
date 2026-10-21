#!/bin/bash
#SBATCH --job-name=giap-calibration
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-12
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# The ADR-0015 calibration harness as a job array (GIAP_SPEC.md §8.4, §12; M4-05). Until CI-01 this
# file was scripts/slurm/calibration_full.sbatch; that name now belongs to the array that runs the
# whole `-m calibration_full` selection in shards, which is a different job. This one runs one
# *registered suite* per task through tests/calibration/_criteria.py (replicate blocks, --shard,
# --merge), i.e. only the tests that `calibrate()` owns.
# One suite per array task, task i = the i-th key
# of replicates_full in tests/calibration/tolerances.yaml (13 suites -> --array=0-12). Each test
# belongs to exactly one suite (tests/calibration/conftest.py SUITE_REGISTRY).
# Run by hand or from a weekly job; extra arguments go to the harness (tests/calibration/_criteria.py):
#   --fast, --dry-run
#   --shard K/N [--block 2]  run replicates block K of N (of confirmation block 2, ADR-0015); give
#                            every shard the same GIAP_CALIBRATION_LOG_DIR, e.g.
#                            GIAP_CALIBRATION_LOG_DIR=$GIAP_ROOT/logs/calibration/<tag> sbatch ...
#   --merge                  evaluate block 1 from the merged shard records on the full R and,
#                            for each failing check, run its block 2 (R..2R-1) here, sequentially
#                            (ADR-0015); same GIAP_CALIBRATION_LOG_DIR as the shards. Exit 0 pass
#                            or unconfirmed, 1 confirmed fail (both blocks), 2 merge error.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "${GIAP_WORKTREE:-$GIAP_REPO}"
export GIAP_CALIBRATION_LOG_DIR="${GIAP_CALIBRATION_LOG_DIR:-$GIAP_ROOT/logs/calibration/$SLURM_ARRAY_JOB_ID}"
uv run --no-sync python tests/calibration/_criteria.py --task-id "$SLURM_ARRAY_TASK_ID" "$@"
