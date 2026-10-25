#!/bin/bash
# Run the full calibration tier and record it for the readiness checker (issue M11-03b; GIAP_SPEC.md
# §8.4, §13.2; scripts/release/readiness.py R3, R4).
#
# Usage: scripts/slurm/calibration_full_launch.sh [REF]     (REF default origin/main, after a fetch)
#   1. REF is resolved to a 40-hex SHA and checked out in the detached worktree
#      $GIAP_ROOT/ci/calibration-<sha12> (flock on $GIAP_ROOT/ci/.git.lock, as the nightly does);
#      a fresh log directory $GIAP_ROOT/logs/calibration/full-<sha12>-<YYYYmmdd-HHMMSS>.
#   2. R14 (spec rule 3, amended in the #349 review), before any worktree, log dir or sbatch: the
#      user's running + pending giap-* jobs and CPUs (squeue -r: one line per array task). Slurm
#      queues all 15 array tasks and the aggregator at once, so refuse when jobs + 16 > 20.
#      Throttle T = min(15, (128-cpus-1)/4) (4 CPUs a task, 1 for the aggregator); refuse when
#      T < 4. The nightly is protected by time, not by a reserve: refuse when the worst-case end
#      now + ceil(15/T) x 4 h + 20 min is at or after the next nightly start (GIAP_NIGHTLY_AT,
#      HH:MM local, default 01:30; a start equal to now counts). The clock is GIAP_NOW (epoch
#      seconds) else date +%s. A failing squeue, a non-integer CPU count, a bad GIAP_NOW or
#      GIAP_NIGHTLY_AT also refuse: exit 2, nothing submitted.
#   3. sbatch --array=0-14%T calibration_full.sbatch of the worktree (all 15 suites at full R), then
#      the 1-CPU giap-calibration-rows job (--dependency=afterany:<array>) running
#      scripts/ci/calibration_rows.py, which appends the calibration_full and calibration_attacks
#      rows to $GIAP_ROOT/logs/nightly.tsv. If that submission fails the array is scancelled (by
#      its ID) and the script exits 1: an array without its aggregator leaves no evidence.
# Any argument besides REF (harness options such as --fast, --shard, --dry-run) is refused (exit 2):
# the rows must describe the full tier. Prints "array=<id> rows=<id> log=<dir> throttle=<T>".
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
REPO="$GIAP_REPO"
SUITES=15   # replicates_full keys = calibration_full.sbatch --array=0-14 (the rows job checks it)
die() { echo "calibration_full_launch: $2" >&2; exit "$1"; }

[ "$#" -le 1 ] || die 2 "usage: $0 [REF]; harness arguments are refused (full tier only)"
ref="${1:-origin/main}"
case "$ref" in -*) die 2 "usage: $0 [REF]; '$ref' is not a ref (full tier only)" ;; esac
git -C "$REPO" fetch -q origin || die 2 "git fetch origin failed"
sha=$(git -C "$REPO" rev-parse -q --verify "$ref^{commit}") || die 2 "cannot resolve '$ref'"
short=${sha:0:12}

# R14 (GIAP_SPEC.md §0.2), fail closed
queue=$(squeue -h -r -u "${USER:?USER unset}" -t RUNNING,PENDING -o '%C %j') ||
  die 2 "squeue failed; cannot check the R14 budget, nothing submitted"
jobs=0; cpus=0
while read -r c name _; do
  [ -n "${c:-}" ] || continue
  case "$c" in *[!0-9]*) die 2 "R14: unexpected squeue CPU count '$c'; nothing submitted" ;; esac
  case "$name" in giap-*) jobs=$((jobs + 1)); cpus=$((cpus + c)) ;; esac
done <<< "$queue"
[ $((jobs + SUITES + 1)) -le 20 ] ||
  die 2 "R14: $jobs giap-* jobs / $cpus CPUs queued; this launch queues 16 jobs ($SUITES array tasks and the aggregator) and $jobs + 16 > 20 jobs; nothing submitted"
t=$SUITES
by_cpu=$(( (128 - cpus - 1) / 4 )); [ $((128 - cpus - 1)) -ge 0 ] || by_cpu=-1
[ "$by_cpu" -lt "$t" ] && t=$by_cpu
if [ "$t" -lt 4 ]; then
  die 2 "R14: $jobs giap-* jobs / $cpus CPUs queued leave a throttle of T=$t < 4 array tasks (4 CPUs each, 1 for the aggregator, of 128); nothing submitted"
fi
now="${GIAP_NOW:-$(date +%s)}"
case "$now" in ''|*[!0-9]*) die 2 "GIAP_NOW '$now' is not epoch seconds; nothing submitted" ;; esac
at="${GIAP_NIGHTLY_AT:-01:30}"
waves=$(( (SUITES + t - 1) / t ))
# python3 (not the venv) for local-time arithmetic; prints "<end> <start>", exit 1 = overlap
window=$(python3 -I -c '
import datetime as d, re, sys
now, secs, at = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
m = re.fullmatch(r"([01][0-9]|2[0-3]):([0-5][0-9])", at)
if not m:
    sys.exit(2)
t0 = d.datetime.fromtimestamp(now)
start = d.datetime.combine(t0.date(), d.time(int(m[1]), int(m[2])))
if start < t0:
    start = d.datetime.combine(t0.date() + d.timedelta(days=1), start.time())
end = d.datetime.fromtimestamp(now + secs)
print(end.strftime("%Y-%m-%d %H:%M"), start.strftime("%Y-%m-%d %H:%M"))
sys.exit(0 if now + secs < start.timestamp() else 1)
' "$now" $((waves * 4 * 3600 + 20 * 60)) "$at")
rc=$?
read -r end_d end_t start_d start_t _ <<< "$window"
case $rc in
  0) ;;
  1) die 2 "R14: $jobs giap-* jobs / $cpus CPUs queued, throttle T=$t: $waves wave(s) of 4 h + 20 min end $end_d $end_t, at or after the next nightly start $start_d $start_t (GIAP_NIGHTLY_AT $at); nothing submitted" ;;
  *) die 2 "cannot place the nightly window (GIAP_NIGHTLY_AT '$at' is not HH:MM, or python3 failed); nothing submitted" ;;
esac

wt="$GIAP_ROOT/ci/calibration-$short"
mkdir -p "$GIAP_ROOT/ci" "$GIAP_ROOT/logs/calibration" "$REPO/local/logs"
(
  flock 9
  [ -d "$wt" ] || git -C "$REPO" worktree add -q --detach "$wt" "$sha"
  git -C "$wt" checkout -q -f --detach "$sha"
) 9>"$GIAP_ROOT/ci/.git.lock" || die 2 "cannot check out $sha in $wt"
log="$GIAP_ROOT/logs/calibration/full-$short-$(date +%Y%m%d-%H%M%S)"
[ ! -e "$log" ] || die 2 "$log exists; nothing submitted"
mkdir -p "$log"

# the log dir travels in the environment (ALL): a comma in it would split an --export list
export GIAP_CALIBRATION_LOG_DIR="$log"
array=$(sbatch --parsable --array="0-$((SUITES - 1))%$t" --export="ALL,GIAP_WORKTREE=$wt" \
  "$wt/scripts/slurm/calibration_full.sbatch") || die 1 "the calibration array was not submitted"
array=${array%%;*}
py="${UV_PROJECT_ENVIRONMENT:-}/bin/python"; [ -x "$py" ] || py=python3
cmd=$(printf '%q ' "$py" "$wt/scripts/ci/calibration_rows.py" --log-dir "$log" --job "$array" \
  --commit "$sha" --expect "$SUITES" --tsv "$GIAP_ROOT/logs/nightly.tsv")
if ! rows=$(sbatch --parsable --job-name=giap-calibration-rows --cpus-per-task=1 --mem=2G \
  --time=00:20:00 --dependency="afterany:$array" --output="$REPO/local/logs/%x.%j.out" \
  --wrap="$cmd"); then
  scancel "$array"
  die 1 "array $array submitted, but its giap-calibration-rows job was refused; array cancelled"
fi
echo "array=$array rows=${rows%%;*} log=$log throttle=$t"
