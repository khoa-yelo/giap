#!/bin/bash
#SBATCH --job-name=giap-calibration
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-14
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=04:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Full calibration suite (GIAP_SPEC.md §8.4, §12): one suite per array task, task i = the i-th key
# of replicates_full in tests/calibration/tolerances.yaml (15 suites -> --array=0-14). Each test
# belongs to exactly one suite (tests/calibration/conftest.py SUITE_REGISTRY).
# Run by hand or from a weekly job; extra arguments go to the harness (tests/calibration/_criteria.py):
#   --fast, --dry-run
#   --shard K/N [--block 2]  run replicates block K of N (of confirmation block 2, ADR-0015); give
#                            every shard the same GIAP_CALIBRATION_LOG_DIR, e.g.
#                            GIAP_CALIBRATION_LOG_DIR=$GIAP_ROOT/logs/calibration/<tag> sbatch ...
#   GIAP_WORKTREE=<dir>      run that checkout (PYTHONPATH set to its src, tests/sim, tests; the harness
#                            refuses a giap imported from elsewhere)
#   --merge                  evaluate block 1 from the merged shard records on the full R and,
#                            for each failing check, run its block 2 (R..2R-1) here, sequentially
#                            (ADR-0015); same GIAP_CALIBRATION_LOG_DIR as the shards. Exit 0 pass,
#                            unconfirmed or vacuous, 1 confirmed fail (both blocks), 2 merge error
#                            (including a non-zero pytest exit of the block-2 confirm run).
# Task 13 (variants, M7-13 lead decisions 4, 7) is sharded inside the task to fit --time: without
# --shard/--merge/--dry-run (either "--opt value" or "--opt=value" form) it runs --shard 0/2 and
# 1/2 concurrently (shard_replicates is 0-based), then --merge (block 2 of the failing checks).
# The task exits non-zero if any shard's pytest or the merge is non-zero. Other tasks run as before.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "${GIAP_WORKTREE:-$GIAP_REPO}"
# A worktree job must import the worktree's giap, not the shared venv's editable install (the main
# checkout's src; Slurm 1885838). _criteria.py prints giap.__file__ and exits 2 if it is not under
# $GIAP_WORKTREE/src.
if [ -n "${GIAP_WORKTREE:-}" ]; then
  export PYTHONPATH="$GIAP_WORKTREE/src:$GIAP_WORKTREE/tests/sim:$GIAP_WORKTREE/tests"
fi
export GIAP_CALIBRATION_LOG_DIR="${GIAP_CALIBRATION_LOG_DIR:-$GIAP_ROOT/logs/calibration/$SLURM_ARRAY_JOB_ID}"
harness() { uv run --no-sync python tests/calibration/_criteria.py --task-id "$SLURM_ARRAY_TASK_ID" "$@"; }
# array task -> in-task shards (a case, not declare -A: macOS /bin/bash is 3.2)
case "$SLURM_ARRAY_TASK_ID" in 13) n=2 ;; *) n=1 ;; esac
for a in "$@"; do
  case "$a" in --shard|--shard=*|--merge|--merge=*|--dry-run) n=1 ;; esac
done
if [ "$n" -le 1 ]; then rc=0; harness "$@" || rc=$?; exit "$rc"; fi
pids=()
for k in $(seq 0 $((n - 1))); do harness "$@" --shard "$k/$n" & pids+=("$!"); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=$?; done
echo "shards of task $SLURM_ARRAY_TASK_ID: pytest exit $rc"
mrc=0; harness "$@" --merge || mrc=$?
echo "merge of task $SLURM_ARRAY_TASK_ID: exit $mrc"
if [ "$rc" -ne 0 ]; then exit "$rc"; fi
exit "$mrc"
