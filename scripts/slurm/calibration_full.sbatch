#!/bin/bash
#SBATCH --job-name=giap-calfull
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# The `calibration_full` tier as a job array (CI-01; GIAP_SPEC.md §8.4, §12): task i runs shard i of
# the whole `-m calibration_full` selection. Nothing else runs this tier — CI, test_branch.sbatch and
# the nightly all deselect it — so it is submitted by hand (no scrontab entry; CI-01 decision 4):
#   sbatch scripts/slurm/calibration_full.sbatch              # the array alone, no summary job
#   bash   scripts/slurm/calibration_full_launch.sh           # the array + its afterany summary job
# Sharding (decision 1) is by test id, not by file: a single test is the expensive unit here. Each
# task collects the selection itself and keeps the ids with sha1(id) % CALFULL_SHARDS == task
# (scripts/ci/shard_tests.py) — a partition that every task computes identically, with no plan file.
# Shard count (decision 2, measured 2026-10-20 on commit 646da60, job 1861630: every one of the 161
# selected tests timed on its own, 8 at a time on one 16-CPU node): 158 tests under 20 s,
# test_related_dbstrain_stratum0_estimate_covers_false_placements 24 min,
# test_database_stressed_stratum0_placement_fdr_full 76 min, the third placement arm ~8 min; the
# tier totals ~2 h of CPU time and is dominated by those arms. 8 shards x 8 CPUs = 64 CPUs (+ 1 for
# the summary job) is the R14-safe size (tasks x cpus-per-task + 1 <= 128): the shard that draws the
# 76-minute arm sets the wall clock (~1.5 h), the others finish in minutes, and --time=04:00:00
# leaves room for the M5-06 island attacks and ADR-0020's hour-long arms as they land. Raise
# CALFULL_SHARDS (and --array to match) when the tier grows; more, smaller tasks beat wider ones.
# Every task appends one line to $GIAP_ROOT/logs/calibration_full.tsv (date, job_task, commit, exit,
# summary — the shape of nightly.tsv, decision 3), also when it fails early or is killed (SIGTERM at
# the time limit: EXIT and TERM traps), and the per-test durations (`--durations=0`) to
# $GIAP_ROOT/logs/calibration_full_durations.tsv (date, job_task, commit, test id, seconds, phase);
# the two files are kept apart so the duration rows cannot be mistaken for shard rows by
# scripts/ci/post_summary.py, which the task then calls to post the table once every shard has a
# line. NIGHTLY_POST=dry (the default here, decision 4: posting the full tier to the nightly issue
# is still Khoa's call) prints the comment instead of posting; NIGHTLY_POST=yes posts to the issue
# titled $CALFULL_TITLE.
# Each task checks out CALFULL_REF (default origin/main) in its own detached worktree
# $GIAP_ROOT/ci/calfull-<ref>-t<task>; scratch in $GIAP_ROOT/ci/runs/calfull_<job>_<task> (7 days).
# Overrides: CALFULL_REF, CALFULL_SHARDS (default = the array's task count), CALFULL_MARK (default
# calibration_full), CALFULL_PATHS (default tests/calibration), CALFULL_PYTEST_ARGS (extra pytest
# arguments, word-split), CALFULL_TITLE, CALFULL_ARRAY_JOB, CALFULL_EXPECT, NIGHTLY_POST.
# A cheap smoke run of the sharding, e.g. on four shards of the fast calibration tests:
#   sbatch --array=0-3 --cpus-per-task=2 --mem=8G --time=00:30:00 \
#     --export=ALL,CALFULL_SHARDS=4,CALFULL_MARK=calibration_smoke scripts/slurm/calibration_full.sbatch
# Never pushes and never commits. R14: 8 tasks x 8 CPUs = 64 CPUs + the 1-CPU summary job.
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # python for report() even before the checkout
REF="${CALFULL_REF:-origin/main}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
JOB="${CALFULL_ARRAY_JOB:-${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}}"
SHARDS="${CALFULL_SHARDS:-${SLURM_ARRAY_TASK_COUNT:-1}}"
EXPECT="${CALFULL_EXPECT:-$SHARDS}"
MARK="${CALFULL_MARK:-calibration_full}"
PATHS="${CALFULL_PATHS:-tests/calibration}"
TITLE="${CALFULL_TITLE:-Calibration full tier (FarmShare)}"
N="${SLURM_CPUS_PER_TASK:-8}"
TSV="$GIAP_ROOT/logs/calibration_full.tsv"
DURATIONS="$GIAP_ROOT/logs/calibration_full_durations.tsv"
CI="$GIAP_ROOT/ci"
W="${CALFULL_WORKTREE:-$CI/calfull-$(echo "${REF#origin/}" | tr '/' '-')-t$TASK}"
RUN="$CI/runs/calfull_${JOB}_$TASK"
commit=$(git -C "$GIAP_REPO" rev-parse -q --short --verify "$REF^0" 2>/dev/null || echo unknown)
stage="setup"; reported=""; numba_tmp=""
mkdir -p "$CI/runs" "$GIAP_ROOT/logs"

report() {  # report <exit> <summary>: this shard's calibration_full.tsv line, then the summary post
  reported=1
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "${JOB}_$TASK" "$commit" "$1" "$2" >> "$TSV"
  echo "CALFULL_EXIT=$1 $2"
  local dry=() ps="$W/scripts/ci/post_summary.py"
  [ -f "$ps" ] || ps="$GIAP_REPO/scripts/ci/post_summary.py"   # the checkout may have failed
  [ "${NIGHTLY_POST:-dry}" = dry ] && dry=(--dry-run)
  python "$ps" --job "$JOB" --tsv "$TSV" --expect "$EXPECT" --title "$TITLE" ${dry[@]+"${dry[@]}"}
}
on_exit() {  # a shard that ends without its line (early exit, error) still reports
  local rc=$?
  [ -n "$reported" ] || report "$(( rc ? rc : 1 ))" "shard $TASK/$SHARDS: ended during $stage (exit $rc)"
  [ -z "$numba_tmp" ] || rm -rf "$numba_tmp"   # only this task's cache, never the shared one
}
trap on_exit EXIT
# a row only if none was written yet (SIGTERM may land inside report itself)
trap '[ -n "$reported" ] || report 143 "shard $TASK/$SHARDS: killed by SIGTERM (time limit?) during $stage"; exit 143' TERM

case "$SHARDS" in ''|*[!0-9]*|0) echo "CALFULL_SHARDS must be a positive integer, not '$SHARDS'" >&2
  report 2 "shard $TASK: CALFULL_SHARDS='$SHARDS' is not a positive integer"; exit 2 ;;
esac
if [ "$TASK" -ge "$SHARDS" ]; then
  echo "array task $TASK outside 0..$((SHARDS - 1)): the array is wider than CALFULL_SHARDS" >&2
  report 2 "shard $TASK: outside 0..$((SHARDS - 1)) (CALFULL_SHARDS=$SHARDS)"; exit 2
fi
find "$CI/runs" -mindepth 1 -maxdepth 1 -mtime +7 -exec rm -rf {} + 2>/dev/null
rm -rf "$RUN"; mkdir -p "$RUN/tmp"
stage="checkout of $REF"
(
  flock 9   # the array tasks share one repository: fetch and worktree add one at a time
  git -C "$GIAP_REPO" fetch -q origin
  [ -d "$W" ] || git -C "$GIAP_REPO" worktree add -q --detach "$W" "$REF"
  git -C "$W" checkout -q -f --detach "$REF"
) 9>"$CI/.git.lock" || { echo "cannot check out $REF in $W" >&2; exit 2; }
cd "$W" || exit 2
numba_tmp="$(mktemp -d /tmp/giap-numba.XXXXXX)"
export TMPDIR="$RUN/tmp" NUMBA_CACHE_DIR="$numba_tmp"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export GIAP_CALIBRATION_LOG_DIR="${GIAP_CALIBRATION_LOG_DIR:-$GIAP_ROOT/logs/calibration/$JOB}"
commit=$(git -C "$W" rev-parse --short HEAD)
echo "calibration_full shard=$TASK/$SHARDS mark=$MARK ref=$REF commit=$commit node=$(hostname) start=$(date -Is)"

stage="collection of -m $MARK"
python -m pytest $PATHS -q -p no:cacheprovider --collect-only -m "$MARK" > "$RUN/collect.txt" 2>&1
crc=$?
if [ "$crc" -ne 0 ]; then
  tail -20 "$RUN/collect.txt"
  report "$crc" "shard $TASK/$SHARDS: collection of -m $MARK failed (exit $crc)"; exit "$crc"
fi
stage="sharding"
python "$W/scripts/ci/shard_tests.py" --shard "$TASK" --of "$SHARDS" "$RUN/collect.txt" \
  > "$RUN/shard.txt" || { report 2 "shard $TASK/$SHARDS: shard_tests.py failed"; exit 2; }
mapfile -t ids < "$RUN/shard.txt"
total=$(grep -c '::' "$RUN/collect.txt")
echo "shard $TASK/$SHARDS: ${#ids[@]} of $total collected tests"
if [ "${#ids[@]}" -eq 0 ]; then
  report 0 "shard $TASK/$SHARDS: 0 of $total tests, nothing to run"; exit 0
fi

stage="pytest ($MARK, ${#ids[@]} tests)"
# shellcheck disable=SC2086  # CALFULL_PYTEST_ARGS is meant to be word-split
python -m pytest "${ids[@]}" -q -p no:cacheprovider -n "$N" -m "$MARK" --durations=0 \
  ${CALFULL_PYTEST_ARGS:-} 2>&1 | tee "$RUN/pytest.txt"
rc=${PIPESTATUS[0]}
# decision 2: the per-test durations, so the next shard count is sized from measurement
awk -F'\n' -v d="$(date +%F)" -v j="${JOB}_$TASK" -v c="$commit" '
  /^[ ]*[0-9.]+s[ ]+(call|setup|teardown)[ ]+/ {
    line = $0
    sub(/^[ ]+/, "", line)
    secs = line + 0
    phase = line; sub(/^[0-9.]+s[ ]+/, "", phase); sub(/[ ].*$/, "", phase)
    id = line; sub(/^[0-9.]+s[ ]+[a-z]+[ ]+/, "", id)
    printf "%s\t%s\t%s\t%s\t%.2f\t%s\n", d, j, c, id, secs, phase
  }' "$RUN/pytest.txt" >> "$DURATIONS"
summary="shard $TASK/$SHARDS: ${#ids[@]} of $total tests, $(grep -E '[0-9]+ (passed|failed|error)' "$RUN/pytest.txt" | tail -1 | tr -d '=' | xargs)"
report "$rc" "$summary"
exit "$rc"
