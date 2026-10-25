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
# Evidence for scripts/ci/calibration_rows.py (M11-03b; launched by calibration_full_launch.sh):
# every pytest run writes $GIAP_CALIBRATION_LOG_DIR/task-<i>.<part>.junit.xml (part run, shard<k>
# or merge; through PYTEST_ADDOPTS, which the harness's pytest subprocesses inherit) and the EXIT
# and TERM traps write task-<i>.json once (task, suite, job, commit, dirty, exit, args, sharded,
# parts, date; SIGTERM: exit 143). The record uses the venv python or python3, never uv.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"
TASK="$SLURM_ARRAY_TASK_ID"
export GIAP_CALIBRATION_LOG_DIR="${GIAP_CALIBRATION_LOG_DIR:-$GIAP_ROOT/logs/calibration/$SLURM_ARRAY_JOB_ID}"
CHECKOUT="${GIAP_WORKTREE:-$GIAP_REPO}"
parts=(); sharded=false; recorded=""
RECORD_PY='import json, os, sys
out, task, job, commit, dirty, rc, sharded, tol, *rest = sys.argv[1:]
cut = rest.index("--")
suite = None
try:
    import yaml
    with open(tol) as fh:
        suites = list(yaml.safe_load(fh)["replicates_full"])
    suite = suites[int(task)] if 0 <= int(task) < len(suites) else None
except Exception:
    pass
import datetime
rec = {"task": int(task), "suite": suite, "job": job, "commit": commit or None,
       "dirty": dirty != "0", "exit": int(rc), "args": rest[cut + 1:], "sharded": sharded == "true",
       "parts": rest[:cut], "date": datetime.date.today().isoformat()}
tmp = f"{out}.tmp.{os.getpid()}"
with open(tmp, "w") as fh:
    fh.write(json.dumps(rec) + "\n")
os.replace(tmp, out)'
write_record() {  # write_record <exit>: task-<i>.json, once (one write to a temporary name, renamed)
  [ -z "$recorded" ] || return 0
  recorded=1
  local py="${UV_PROJECT_ENVIRONMENT:-}/bin/python" commit st dirty=1 have=() p
  [ -x "$py" ] || py=python3
  commit=$(git -C "$CHECKOUT" rev-parse HEAD 2>/dev/null) || commit=""
  if [ -n "$commit" ] && st=$(git -C "$CHECKOUT" status --porcelain -- src tests/calibration tests/sim tests/helpers ci 2>/dev/null); then
    [ -n "$st" ] || dirty=0
  fi
  for p in ${parts[@]+"${parts[@]}"}; do   # the merge writes JUnit only when it runs block 2
    if [ "$p" != "task-$TASK.merge.junit.xml" ] || [ -e "$GIAP_CALIBRATION_LOG_DIR/$p" ]; then have+=("$p"); fi
  done
  mkdir -p "$GIAP_CALIBRATION_LOG_DIR" &&
    "$py" -c "$RECORD_PY" "$GIAP_CALIBRATION_LOG_DIR/task-$TASK.json" "$TASK" "${SLURM_ARRAY_JOB_ID:-}" \
      "$commit" "$dirty" "$1" "$sharded" "$CHECKOUT/tests/calibration/tolerances.yaml" \
      ${have[@]+"${have[@]}"} -- ${ARGS[@]+"${ARGS[@]}"} ||
    echo "calibration_full: could not write the record of task $TASK" >&2
}
on_exit() {
  local rc=$?
  [ "${BASHPID:-$$}" = "$$" ] || exit "$rc"   # a forked shard subshell: no record of its own
  write_record "$rc"
  rm -rf "$TMPDIR"
  exit "$rc"
}
trap on_exit EXIT
trap '[ "${BASHPID:-$$}" = "$$" ] || exit 143; write_record 143; exit 143' TERM
ARGS=("$@")
cd "$CHECKOUT"
# A worktree job must import the worktree's giap, not the shared venv's editable install (the main
# checkout's src; Slurm 1885838). _criteria.py prints giap.__file__ and exits 2 if it is not under
# $GIAP_WORKTREE/src.
if [ -n "${GIAP_WORKTREE:-}" ]; then
  export PYTHONPATH="$GIAP_WORKTREE/src:$GIAP_WORKTREE/tests/sim:$GIAP_WORKTREE/tests"
fi
harness() {  # harness <part> <harness args...>: one pytest run with its own JUnit file
  local part="$1"; shift
  PYTEST_ADDOPTS="${PYTEST_ADDOPTS:+$PYTEST_ADDOPTS }--junitxml=$GIAP_CALIBRATION_LOG_DIR/task-$TASK.$part.junit.xml" \
    uv run --no-sync python tests/calibration/_criteria.py --task-id "$SLURM_ARRAY_TASK_ID" "$@"
}
# array task -> in-task shards (a case, not declare -A: macOS /bin/bash is 3.2)
case "$SLURM_ARRAY_TASK_ID" in 13) n=2 ;; *) n=1 ;; esac
for a in "$@"; do
  case "$a" in --shard|--shard=*|--merge|--merge=*|--dry-run) n=1 ;; esac
done
if [ "$n" -le 1 ]; then
  parts=("task-$TASK.run.junit.xml"); rc=0; harness run "$@" || rc=$?; exit "$rc"
fi
sharded=true
for k in $(seq 0 $((n - 1))); do parts+=("task-$TASK.shard$k.junit.xml"); done
parts+=("task-$TASK.merge.junit.xml")
pids=()
for k in $(seq 0 $((n - 1))); do harness "shard$k" "$@" --shard "$k/$n" & pids+=("$!"); done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=$?; done
echo "shards of task $SLURM_ARRAY_TASK_ID: pytest exit $rc"
mrc=0; harness merge "$@" --merge || mrc=$?
echo "merge of task $SLURM_ARRAY_TASK_ID: exit $mrc"
if [ "$rc" -ne 0 ]; then exit "$rc"; fi
exit "$mrc"
