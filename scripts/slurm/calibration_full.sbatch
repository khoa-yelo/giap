#!/bin/bash
#SBATCH --job-name=giap-calibration
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-14
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Full calibration suite (GIAP_SPEC.md §8.4, §12): one suite per array task, task i = the i-th key
# of replicates_full in tests/calibration/tolerances.yaml (15 suites -> --array=0-14). Each test
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
# Task 13 (variants, M7-13 lead decision 4) is sharded inside the task to fit --time: without
# --shard/--merge/--dry-run it runs --shard 1/2 and 2/2 concurrently, then --merge
# (block 2 of the failing checks). Other tasks run as before.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "${GIAP_WORKTREE:-$GIAP_REPO}"
export GIAP_CALIBRATION_LOG_DIR="${GIAP_CALIBRATION_LOG_DIR:-$GIAP_ROOT/logs/calibration/$SLURM_ARRAY_JOB_ID}"
harness() { uv run --no-sync python tests/calibration/_criteria.py --task-id "$SLURM_ARRAY_TASK_ID" "$@"; }
declare -A SHARDS=([13]=2)   # array task -> in-task shards
n="${SHARDS[$SLURM_ARRAY_TASK_ID]:-1}"
case " $* " in *" --shard "*|*" --merge "*|*" --dry-run "*) n=1 ;; esac
if [ "$n" -le 1 ]; then harness "$@"; exit $?; fi
pids=()
for k in $(seq 1 "$n"); do harness "$@" --shard "$k/$n" & pids+=("$!"); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=$?; done
echo "shards of task $SLURM_ARRAY_TASK_ID: pytest exit $rc; the verdict is the merge below"
harness "$@" --merge
