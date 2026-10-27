#!/bin/bash
# Submit one calibration suite as N shard jobs plus a dependent --merge job (GIPPA_SPEC.md §8.4, §12;
# ADR-0015 two-block confirmation; ADR-0057 xfail-aware merge).
#
# Usage: scripts/slurm/calibration_sharded.sh SUITE N [harness args...]
#   SUITE  a key of replicates_full in tests/calibration/tolerances.yaml (its array task index is
#          looked up there, so calibration_full.sbatch runs exactly that suite)
#   N      number of shards, 1..10 (shard k runs replicates shard_replicates(R, N, k) of every check,
#          whatever its R; tests' direct parts run in shard 0 only)
#   extra  passed to tests/calibration/_criteria.py in every job, e.g. --fast. --shard, --merge and
#          --dry-run are refused here (this script sets them), and so is --block (a block-2 run never
#          owns the direct checks; --merge runs block 2 itself).
# Environment:
#   GIPPA_CALIBRATION_LOG_DIR  shared log directory of the shards and the merge (default
#                             $GIPPA_ROOT/logs/calibration/sharded-<suite>-<YYYYmmdd-HHMMSS>); it must
#                             not already hold any <suite>.* file (old records, outcomes or shard
#                             exits would be merged). It reaches the jobs through the environment
#                             (--export=ALL), so a comma in the path is safe.
#   GIPPA_SBATCH_EXTRA         extra sbatch options for every job (word-split), e.g. "--time=02:00:00".
# The shards are separate jobs gippa-cal-<suite>-s<k>; the merge (gippa-cal-<suite>-merge) waits for all of
# them with afterany, so it runs and reports MERGE_ERROR when a shard died. Exit of the merge job: 0 pass,
# unconfirmed, vacuous or expected failure (strict xfail), 1 confirmed fail or strict XPASS, 2 merge
# error. R14: N + 1 jobs of 4 CPUs (calibration_full.sbatch); before any sbatch the script sums the
# user's running + pending gippa-* jobs and CPUs (squeue) and refuses (exit 2) when the new jobs would
# pass 128 CPUs or 20 jobs; it fails closed (exit 2) when squeue fails.
# Prints "shards=<ids> merge=<id> log=<dir>".
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

usage() { echo "usage: $0 SUITE N [harness args...]  (N in 1..10)" >&2; exit 2; }
[ "$#" -ge 2 ] || usage
suite="$1"; n="$2"; shift 2
case "$n" in ''|*[!0-9]*) usage ;; esac
if [ "$n" -lt 1 ] || [ "$n" -gt 10 ]; then usage; fi
for a in "$@"; do
  case "$a" in
    --shard|--shard=*|--merge|--merge=*|--dry-run|--task-id|--task-id=*|--suite|--suite=*|--block|--block=*)
      echo "$0: $a is set by this script" >&2; exit 2 ;;
  esac
done

py="${GIPPA_PYTHON:-${UV_PROJECT_ENVIRONMENT:-}/bin/python}"
[ -x "$py" ] || py=python3
task="$("$py" - "$HERE/tests/calibration/tolerances.yaml" "$suite" <<'EOF'
import sys, yaml
suites = list(yaml.safe_load(open(sys.argv[1]))["replicates_full"])
print(suites.index(sys.argv[2]) if sys.argv[2] in suites else "")
EOF
)"
if [ -z "$task" ]; then echo "$0: unknown suite '$suite' (not in replicates_full)" >&2; exit 2; fi

log="${GIPPA_CALIBRATION_LOG_DIR:-${GIPPA_ROOT:?GIPPA_ROOT unset}/logs/calibration/sharded-$suite-$(date +%Y%m%d-%H%M%S)}"
for f in "$log/$suite".*; do
  if [ -e "$f" ] || [ -L "$f" ]; then
    echo "$0: $f exists; give a fresh GIPPA_CALIBRATION_LOG_DIR" >&2; exit 2
  fi
done

# R14 (GIPPA_SPEC.md §0.2): running + pending gippa-* jobs and CPUs of this user, plus the N + 1 new
# jobs of 4 CPUs each, must stay within 128 CPUs and 20 jobs. No squeue answer: refuse.
cpus_per_job=4
new_jobs=$((n + 1))
if ! queue="$(squeue -h -u "${USER:?USER unset}" -t RUNNING,PENDING -o '%C %j')"; then
  echo "$0: squeue failed; cannot check the R14 budget, nothing submitted" >&2; exit 2
fi
used_jobs=0; used_cpus=0
while read -r c name _; do
  [ -n "${c:-}" ] || continue
  case "$name" in giap-*|gippa-*) ;; *) continue ;; esac
  case "$c" in ''|*[!0-9]*) echo "$0: unexpected squeue CPU count '$c'; nothing submitted" >&2; exit 2 ;; esac
  used_jobs=$((used_jobs + 1)); used_cpus=$((used_cpus + c))
done <<< "$queue"
want_cpus=$((used_cpus + cpus_per_job * new_jobs)); want_jobs=$((used_jobs + new_jobs))
if [ "$want_cpus" -gt 128 ] || [ "$want_jobs" -gt 20 ]; then
  echo "$0: R14: $used_jobs gippa-* jobs / $used_cpus CPUs queued + $new_jobs jobs / $((cpus_per_job * new_jobs)) CPUs" \
       "= $want_jobs jobs / $want_cpus CPUs > 20 jobs / 128 CPUs; nothing submitted" >&2
  exit 2
fi
mkdir -p "$log"

sbatch_extra=()
if [ -n "${GIPPA_SBATCH_EXTRA:-}" ]; then
  # shellcheck disable=SC2206  # word splitting of the option string is intended
  sbatch_extra=(${GIPPA_SBATCH_EXTRA})
fi
job="$HERE/scripts/slurm/calibration_full.sbatch"
# the log dir travels in the environment (ALL): a comma in it would split an --export list
export GIPPA_CALIBRATION_LOG_DIR="$log"
export_vars="ALL,GIPPA_WORKTREE=$HERE"
submit() {  # submit NAME DEPENDENCY HARNESS-ARGS... -> prints the job id
  local name="$1" dep="$2"; shift 2
  local opts=(--parsable --array="$task" --job-name="$name" --export="$export_vars")
  if [ -n "$dep" ]; then opts+=(--dependency="$dep"); fi
  if [ "${#sbatch_extra[@]}" -gt 0 ]; then opts+=("${sbatch_extra[@]}"); fi
  sbatch "${opts[@]}" "$job" "$@" | cut -d';' -f1
}

ids=()
k=0
while [ "$k" -lt "$n" ]; do
  if ! id="$(submit "gippa-cal-$suite-s$k" "" --shard "$k/$n" "$@")" || [ -z "$id" ]; then
    echo "$0: sbatch failed for shard $k/$n; submitted so far: ${ids[*]:-none} (no merge job)" >&2
    exit 1
  fi
  ids+=("$id")
  k=$((k + 1))
done
dep="afterany$(printf ':%s' "${ids[@]}")"
if ! merge="$(submit "gippa-cal-$suite-merge" "$dep" --merge "$@")" || [ -z "$merge" ]; then
  echo "$0: sbatch failed for the merge job; shards ${ids[*]} run without a merge" >&2
  exit 1
fi
echo "shards=$(IFS=,; echo "${ids[*]}") merge=$merge log=$log"
