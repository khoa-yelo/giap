#!/bin/bash
# Run the full calibration tier and record it for the readiness checker (issues M11-03b, M11-03c;
# GIPPA_SPEC.md §8.4, §13.2; scripts/release/readiness.py R3, R4; spec
# docs/specs/M11-03c-calibration-full-sharding.md, "calibration_full_launch.sh (changed)").
#
# Usage: scripts/slurm/calibration_full_launch.sh [REF]     (REF default origin/main, after a fetch)
#   1. REF is resolved to a 40-hex SHA. R14 jobs, before any worktree, log dir or sbatch: the
#      user's running + pending gippa-* jobs and CPUs (squeue -r: one line per array task); a launch
#      queues 16 jobs (15 suite jobs and the aggregator), so refuse when jobs + 16 > 20.
#   2. The SHA is checked out in the detached worktree $GIPPA_ROOT/ci/calibration-<sha12> (flock on
#      $GIPPA_ROOT/ci/.git.lock, as the nightly does) and the worktree's plan
#      (scripts/slurm/calibration_plan.tsv) is validated with its calibration_lanes.py check (venv
#      python: it needs yaml). Each suite's job asks for lanes x cpus_per_lane CPUs (lanes = shards
#      + groups), lanes x mem_per_lane_mb MB and the plan time.
#   3. R14 CPUs: refuse when cpus + sum(suite CPUs) + 2 > 144, the hard cap (no throttle, no second
#      wave; the 1-CPU aggregator holds a whole core: 2), and warn on stderr ("R14: ... soft cap
#      128"), then submit, when that sum is over 128: 128 is a soft cap (Khoa 2026-10-29: "It's okay
#      if it exceeds 128 a little bit. It's not a hard threshold."). A suite's CPUs are what Slurm allocates (and squeue shows): max(lanes x cpus_per_lane,
#      ceil(lanes x mem_per_lane_mb / GIPPA_MAX_MEM_PER_CPU_MB)), the partition's MaxMemPerCPU,
#      default 4000 (a 2 x 16000 MB discordance job asks for 4 CPUs and gets 8). The
#      nightly is protected by time: refuse when now + max(time) + 20 min is at or after the next
#      nightly start (GIPPA_NIGHTLY_AT, HH:MM local, default 01:30; a start equal to now counts). The
#      clock is GIPPA_NOW (epoch seconds) else date +%s. A failing squeue, a non-integer CPU count, a
#      bad plan, GIPPA_NOW or GIPPA_NIGHTLY_AT also refuse: exit 2, nothing submitted.
#   4. A fresh log directory $GIPPA_ROOT/logs/calibration/full-<sha12>-<YYYYmmdd-HHMMSS>; one
#      calibration_suite.sbatch job per plan row, in plan order (gippa-cal-<suite>, GIPPA_SUITE,
#      GIPPA_WORKTREE, GIPPA_CALIBRATION_LOG_DIR), then the 1-CPU gippa-calibration-rows job
#      (--dependency=afterany:<j0>:...:<j14>) running scripts/ci/calibration_rows.py --job <j0>,
#      which appends the calibration_full and calibration_attacks rows to
#      $GIPPA_ROOT/logs/nightly.tsv. If any submission fails, the jobs already submitted are
#      scancelled (by ID) and the script exits 1: suites without their aggregator leave no evidence.
# Any argument besides REF (harness options such as --fast, --shard, --dry-run) is refused (exit 2):
# the rows must describe the full tier. Prints "jobs=<j0>,...,<j14> rows=<id> log=<dir> cpus=<sum>"
# (<sum> = the allocated CPUs of the suite jobs, as summed for R14).
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
REPO="$GIPPA_REPO"
SUITES=15   # replicates_full keys = plan rows (calibration_lanes.py check; the rows job checks it)
die() { echo "calibration_full_launch: $2" >&2; exit "$1"; }

[ "$#" -le 1 ] || die 2 "usage: $0 [REF]; harness arguments are refused (full tier only)"
ref="${1:-origin/main}"
case "$ref" in -*) die 2 "usage: $0 [REF]; '$ref' is not a ref (full tier only)" ;; esac
git -C "$REPO" fetch -q origin || die 2 "git fetch origin failed"
sha=$(git -C "$REPO" rev-parse -q --verify "$ref^{commit}") || die 2 "cannot resolve '$ref'"
short=${sha:0:12}

# R14 (GIPPA_SPEC.md §0.2), fail closed
queue=$(squeue -h -r -u "${USER:?USER unset}" -t RUNNING,PENDING -o '%C %j') ||
  die 2 "squeue failed; cannot check the R14 budget, nothing submitted"
jobs=0; cpus=0
while read -r c name _; do
  [ -n "${c:-}" ] || continue
  case "$c" in *[!0-9]*) die 2 "R14: unexpected squeue CPU count '$c'; nothing submitted" ;; esac
  # whole cores (CR_CORE_MEMORY, 2 threads per core): a pending 1-CPU job will hold 2 (#445 S1)
  case "$name" in giap-*|gippa-*) jobs=$((jobs + 1)); cpus=$((cpus + (10#$c + 1) / 2 * 2)) ;; esac
done <<< "$queue"
[ $((jobs + SUITES + 1)) -le 20 ] ||
  die 2 "R14: $jobs gippa-* jobs / $cpus CPUs queued; this launch queues 16 jobs ($SUITES suite jobs and the aggregator) and $jobs + 16 > 20 jobs; nothing submitted"
now="${GIPPA_NOW:-$(date +%s)}"
case "$now" in ''|*[!0-9]*) die 2 "GIPPA_NOW '$now' is not epoch seconds; nothing submitted" ;; esac
at="${GIPPA_NIGHTLY_AT:-01:30}"
mpc="${GIPPA_MAX_MEM_PER_CPU_MB:-4000}"
case "$mpc" in ''|*[!0-9]*|0|0*) die 2 "GIPPA_MAX_MEM_PER_CPU_MB '$mpc' is not a positive integer; nothing submitted" ;; esac

wt="$GIPPA_ROOT/ci/calibration-$short"
mkdir -p "$GIPPA_ROOT/ci" "$GIPPA_ROOT/logs/calibration" "$REPO/local/logs"
(
  flock 9
  [ -d "$wt" ] || git -C "$REPO" worktree add -q --detach "$wt" "$sha"
  git -C "$wt" checkout -q -f --detach "$sha"
) 9>"$GIPPA_ROOT/ci/.git.lock" || die 2 "cannot check out $sha in $wt"

# the worktree's plan, validated before it sizes anything
plan="$wt/scripts/slurm/calibration_plan.tsv"
py="${UV_PROJECT_ENVIRONMENT:-}/bin/python"; [ -x "$py" ] || py=python3
"$py" "$wt/scripts/ci/calibration_lanes.py" check --plan "$plan" ||
  die 2 "the plan $plan of $sha is invalid; nothing submitted"
names=(); sizes=(); total=0; max_s=0
# `|| [ -n "$suite" ]` keeps a last row without a final newline (check accepts it; #399 S1)
suite=""
while IFS=$'\t' read -r suite shards groups _ cpl mpl time _ || [ -n "$suite" ]; do
  lanes=$((10#$shards + 10#$groups))
  IFS=: read -r hh mm ss <<< "$time"
  secs=$((10#$hh * 3600 + 10#$mm * 60 + 10#$ss))
  [ "$secs" -gt "$max_s" ] && max_s=$secs
  names+=("$suite")
  req=$((lanes * 10#$cpl)); mem=$((lanes * 10#$mpl))
  alloc=$(( (mem + 10#$mpc - 1) / 10#$mpc )); [ "$alloc" -lt "$req" ] && alloc=$req
  alloc=$(( (alloc + 1) / 2 * 2 ))  # whole cores, as Slurm allocates them
  sizes+=("$req $mem $time")
  total=$((total + alloc))  # R14 counts what Slurm allocates, not what is asked (MaxMemPerCPU)
done < <(tail -n +2 "$plan")
[ "${#names[@]}" -eq "$SUITES" ] || die 2 "the plan $plan has ${#names[@]} suites, not $SUITES; nothing submitted"

# R14 CPUs: every suite job and the aggregator are queued at once
# (the 1-CPU aggregator holds a whole core: 2 CPUs)
# 128 is a soft cap (warn, then submit), 144 the hard one (refuse); Khoa 2026-10-29
sum=$((cpus + total + 2))
[ "$sum" -le 144 ] ||
  die 2 "R14: $jobs gippa-* jobs / $cpus CPUs queued; this launch requests $((total + 2)) CPUs ($total allocated for $SUITES suite jobs, MaxMemPerCPU $mpc MB, and 2 for the aggregator) and $cpus + $((total + 2)) = $sum > 144 CPUs; nothing submitted"
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
' "$now" $((max_s + 20 * 60)) "$at")
rc=$?
read -r end_d end_t start_d start_t _ <<< "$window"
case $rc in
  0) ;;
  1) die 2 "R14: the longest suite ($((max_s / 60)) min) + 20 min ends $end_d $end_t, at or after the next nightly start $start_d $start_t (GIPPA_NIGHTLY_AT $at); nothing submitted" ;;
  *) die 2 "cannot place the nightly window (GIPPA_NIGHTLY_AT '$at' is not HH:MM, or python3 failed); nothing submitted" ;;
esac
# every refusal is behind us: warn once over the soft cap, then submit
[ "$sum" -le 128 ] ||
  echo "R14: $jobs gippa-* jobs / $cpus CPUs queued + $((total + 2)) for this launch = $sum CPUs, over the soft cap 128 (hard cap 144); submitting" >&2

log="$GIPPA_ROOT/logs/calibration/full-$short-$(date +%Y%m%d-%H%M%S)"
[ ! -e "$log" ] || die 2 "$log exists; nothing submitted"
mkdir -p "$log"

# the log dir travels in the environment (ALL): a comma in it would split an --export list
export GIPPA_CALIBRATION_LOG_DIR="$log"
ids=()
cancel() {
  if [ "${#ids[@]}" -gt 0 ] && ! scancel "${ids[@]}"; then  # orphans would hold R14 (#399 S2)
    die 1 "$1; scancel FAILED for the job(s) already submitted (${ids[*]}): cancel them by hand"
  fi
  die 1 "$1; cancelled the ${#ids[@]} job(s) already submitted (${ids[*]:-none})"
}
for k in "${!names[@]}"; do
  read -r c mem time <<< "${sizes[$k]}"
  id=$(sbatch --parsable --job-name="gippa-cal-${names[$k]}" --cpus-per-task="$c" --mem="${mem}M" \
    --time="$time" --export="ALL,GIPPA_SUITE=${names[$k]},GIPPA_WORKTREE=$wt" \
    "$wt/scripts/slurm/calibration_suite.sbatch") || cancel "suite ${names[$k]} was not submitted"
  ids+=("${id%%;*}")
done
cmd=$(printf '%q ' "$py" "$wt/scripts/ci/calibration_rows.py" --log-dir "$log" --job "${ids[0]}" \
  --commit "$sha" --expect "$SUITES" --tsv "$GIPPA_ROOT/logs/nightly.tsv")
dep=$(IFS=:; echo "afterany:${ids[*]}")
rows=$(sbatch --parsable --job-name=gippa-calibration-rows --cpus-per-task=1 --mem=2G \
  --time=00:20:00 --dependency="$dep" --output="$REPO/local/logs/%x.%j.out" \
  --wrap="$cmd") || cancel "the gippa-calibration-rows job was refused"
echo "jobs=$(IFS=,; echo "${ids[*]}") rows=${rows%%;*} log=$log cpus=$total"
