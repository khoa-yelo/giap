#!/bin/bash
#SBATCH --job-name=giap-nightly
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Nightly job array on the latest origin/main (GIAP_SPEC.md §12; M2-21 lead decisions 3-5):
#   0  tests/integration + tests/contract with --require-capabilities <ci/gates.yaml farmshare_nightly>
#      (a missing DB fails, never skips), then nf-test test workflow/tests --profile $NF_PROFILE
#      (default test,stub: the suites assert fake-tool argv logs and stub outputs; under
#      farmshare,apptainer the pinned images replace the fakes and 9 of 30 tests fail, job 1852060)
#   1  benchmarks pinned to the wheat node class (CPU_SKU:6148): a task that lands elsewhere
#      resubmits itself with --constraint and exits; > 15 % slower than benchmarks/baseline.json fails
#   2  the calibration_nightly tests (the pre-M2-21 nightly)
#   3  real `giap run` on sim100k (giapsim SIM100K + its proteins): Bakta with the light DB, and
#      with the recorded full-DB JSON (annotation.precomputed.bakta_json; its GBFF feeds
#      antiSMASH); measured runtimes via scripts/measure_runtimes.py
# MODE=test (default) | record: record runs the scripts/slurm/record_tools.sbatch logic in task 0
# (sim100k into a fresh directory, compared with $GIAP_ROOT/recordings: identical / changed); the
# other tasks only log "skipped". Submit: sbatch --array=0 --export=ALL,MODE=record nightly.sbatch
# Each task checks out NIGHTLY_REF (default origin/main; decision 5: main only) in its own detached
# worktree $GIAP_ROOT/ci/<ref>-t<task>; work dirs in $GIAP_ROOT/ci/runs/<job>_<task> (kept 7 days).
# Every task writes one line to $GIAP_ROOT/logs/nightly.tsv (date, job_task, commit, exit,
# summary), also when it fails early or is killed (SIGTERM at the time limit; decision 14: EXIT and
# TERM traps), then calls scripts/ci/post_summary.py, which posts the table to the "Nightly
# (FarmShare)" issue once every task has a line (NIGHTLY_POST=dry: print it instead). The 1-CPU
# summary job of scripts/slurm/nightly_launch.sh (afterany) posts whatever is there when a task
# could not report at all (SIGKILL, OOM). Never pushes and never commits.
# Overrides: NIGHTLY_REF, NIGHTLY_PROFILE (default farmshare_local,apptainer: every Nextflow task
# stays inside this allocation; a list with the Slurm-executor profile `farmshare` is refused with
# exit 2, decisions 12 and 22), NIGHTLY_NF_PROFILE, NIGHTLY_POST, NIGHTLY_RECORD_SAMPLES.
# Runs nightly from scripts/slurm/nightly_launch.sh (scrontab: scripts/slurm/nightly.scrontab) or by
# hand: sbatch scripts/slurm/nightly.sbatch (then nothing posts for a task that is SIGKILLed).
# R14: 4 tasks x 16 CPUs = 64 CPUs + the 1-CPU summary job; no Nextflow child jobs (farmshare_local).
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"   # python for report() even before the checkout
MODE="${MODE:-test}"
REF="${NIGHTLY_REF:-origin/main}"
PROFILE="${NIGHTLY_PROFILE:-farmshare_local,apptainer}"
NF_PROFILE="${NIGHTLY_NF_PROFILE:-test,stub}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
JOB="${NIGHTLY_ARRAY_JOB:-${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}}"
EXPECT="${NIGHTLY_EXPECT:-${SLURM_ARRAY_TASK_COUNT:-1}}"
N="${SLURM_CPUS_PER_TASK:-16}"
TSV="$GIAP_ROOT/logs/nightly.tsv"
CI="$GIAP_ROOT/ci"
W="${NIGHTLY_WORKTREE:-$GIAP_ROOT/ci/$(echo "${REF#origin/}" | tr '/' '-')-t$TASK}"
RUN="$CI/runs/${JOB}_$TASK"
commit=$(git -C "$GIAP_REPO" rev-parse -q --short --verify "$REF^0" 2>/dev/null || echo unknown)
stage="setup"; reported=""; numba_tmp=""
mkdir -p "$CI/runs" "$GIAP_ROOT/logs"

tsv_line() {  # tsv_line <job_task> <exit> <summary>
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "$1" "$commit" "$2" "$3" >> "$TSV"
}
report() {  # report <exit> <sub-task: summary>: this task's nightly.tsv line, then the summary post
  reported=1
  tsv_line "${JOB}_$TASK" "$1" "$2"
  echo "NIGHTLY_EXIT=$1 $2"
  local dry=() ps="$W/scripts/ci/post_summary.py"
  [ -f "$ps" ] || ps="$GIAP_REPO/scripts/ci/post_summary.py"   # the checkout may have failed
  [ "${NIGHTLY_POST:-yes}" = dry ] && dry=(--dry-run)
  python "$ps" --job "$JOB" --tsv "$TSV" --expect "$EXPECT" ${dry[@]+"${dry[@]}"}
}
on_exit() {  # decision 14: a task that ends without its line (early exit, error) still reports
  local rc=$?
  [ "$BASHPID" = "$$" ] || return 0   # a forked child hit by the group SIGTERM before its exec: not this task's exit
  [ -n "$reported" ] || report "$(( rc ? rc : 1 ))" "task $TASK: ended during $stage without a result (exit $rc)"
  [ -z "$numba_tmp" ] || rm -rf "$numba_tmp"   # only this task's cache, never the shared one
}
trap on_exit EXIT
# decision 23: a row only if none was written yet (SIGTERM may land inside report itself)
trap '[ "$BASHPID" = "$$" ] || exit 143; [ -n "$reported" ] || report 143 "task $TASK: killed by SIGTERM (time limit?) during $stage"; exit 143' TERM

case "$MODE" in test|record) ;; *) echo "MODE must be test or record, not '$MODE'" >&2; exit 2 ;; esac
case ",$PROFILE," in
  *,farmshare,*)   # decision 22: Slurm-executor child jobs would break the R14 CPU budget
    echo "NIGHTLY_PROFILE=$PROFILE: the farmshare (Slurm executor) profile is refused (R14); use farmshare_local" >&2
    report 2 "task $TASK: NIGHTLY_PROFILE=$PROFILE refused (R14: use farmshare_local)"; exit 2 ;;
esac
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
export TMPDIR="$RUN/tmp" NXF_ANSI_LOG=false NUMBA_CACHE_DIR="$numba_tmp"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
commit=$(git -C "$W" rev-parse --short HEAD)
stage="task $TASK"
echo "nightly task=$TASK mode=$MODE ref=$REF commit=$commit node=$(hostname) profile=$PROFILE start=$(date -Is)"

pyt() {  # pyt <label> <pytest args...>: sets rc and summary
  python -m pytest -p no:cacheprovider "${@:2}" 2>&1 | tee "$RUN/$1.txt"
  rc=${PIPESTATUS[0]}
  summary="$1: $(grep -E '[0-9]+ (passed|failed|error)' "$RUN/$1.txt" | tail -1 | tr -d '=' | xargs)"
}

task_integration() {
  local caps r0 s0
  caps=$(python -c 'import yaml; print(",".join(yaml.safe_load(open("ci/gates.yaml"))["required_capabilities"]["farmshare_nightly"]))')
  pyt integration tests/integration tests/contract -q -n 8 --require-capabilities "$caps"
  r0=$rc; s0=$summary
  nf-test test workflow/tests --profile "$NF_PROFILE" 2>&1 | tee "$RUN/nf-test.txt"
  local r1=${PIPESTATUS[0]}
  rc=$(( r0 > r1 ? r0 : r1 ))
  summary="$s0; nf-test $NF_PROFILE: $(sed 's/\x1b\[[0-9;]*m//g' "$RUN/nf-test.txt" | grep -E '^(SUCCESS|FAILURE):' | tail -1)"
}

task_record() {
  local out="$RUN/recordings" same=0 changed=() nfail
  RECORD_OUT="$out" RECORD_SAMPLES="${NIGHTLY_RECORD_SAMPLES:-sim100k}" \
    bash scripts/slurm/record_tools.sbatch 2>&1 | tee "$RUN/record_tools.txt"
  nfail=$(grep -oP '^failures: \K[0-9]+' "$RUN/record_tools.txt" || echo 1)
  for d in "$out"/*/*/; do   # byte-identical with the accepted recording? (logs and timings aside)
    [ -f "$d/SHA256SUMS" ] || continue
    local ts; ts=$(basename "$(dirname "$d")")/$(basename "$d")
    if [ "$(grep -vE '\./(run\.log|runtime\.txt)$' "$d/SHA256SUMS")" = \
         "$(grep -vE '\./(run\.log|runtime\.txt)$' "$GIAP_ROOT/recordings/$ts/SHA256SUMS" 2>/dev/null)" ]; then
      same=$((same + 1))
    else changed+=("$ts"); fi
  done
  rc=$(( nfail > 0 ? 1 : 0 ))
  summary="record: failures $nfail; identical $same; changed ${#changed[@]}${changed[0]+ (${changed[*]})}"
}

extend_summary() {  # extend_summary <job>: decision 19, the night's summary also waits for <job>
  local sid
  sid=$(squeue -u "${USER:-$(id -un)}" -h -t PENDING -n giap-nightly-summary -o '%i %k' 2>/dev/null \
    | awk -v c="nightly:$JOB" '$2 == c {print $1; exit}')
  [ -n "$sid" ] && scontrol update JobId="$sid" Dependency="afterany:$JOB:$1"
}

task_benchmarks() {
  case "$(hostname -s)" in
    wheat-*) ;;
    *) if [ -z "${NIGHTLY_BENCH_PINNED:-}" ]; then   # array tasks cannot carry their own constraint
         local j   # decision 15: the same commit; the pinned job writes this task's line
         if ! j=$(sbatch --parsable --array="$TASK" --constraint=CPU_SKU:6148 \
           --export="ALL,MODE=$MODE,NIGHTLY_BENCH_PINNED=1,NIGHTLY_ARRAY_JOB=$JOB,NIGHTLY_EXPECT=$EXPECT,NIGHTLY_REF=$commit,NIGHTLY_WORKTREE=$W" \
           scripts/slurm/nightly.sbatch); then
           rc=1; summary="benchmarks: not a wheat node ($(hostname -s)) and the resubmission to wheat failed"
           return
         fi
         j=${j%%;*}
         echo "not a wheat node ($(hostname -s)): benchmarks resubmitted as job $j"
         extend_summary "$j" || tsv_line "${JOB}_note" 0 \
           "benchmarks resubmitted as job $j; no pending giap-nightly-summary job of $JOB to extend: it may post before them"
         reported=1; exit 0
       fi
       rc=1; summary="benchmarks: not on a wheat node ($(hostname -s))"; return ;;
  esac
  pyt benchmarks benchmarks --benchmark-only -q -p no:xdist --benchmark-json "$RUN/benchmarks.json"
  summary="$summary (wheat $(hostname -s))"
}

task_calibration() {
  pyt calibration_nightly tests/calibration -q -n "$N" -m calibration_nightly --durations=15
}

task_giap_run() {
  local sim="$RUN/sim" full="$GIAP_ROOT/recordings/bakta_full/sim100k/bakta/sim100k.json" r out
  python -c 'import dataclasses, pathlib, sys, giapsim; from giapsim.proteomics import ProteomicsConfig
giapsim.simulate(dataclasses.replace(giapsim.SIM100K, proteomics=ProteomicsConfig()), giapsim.SIM100K_SEED, pathlib.Path(sys.argv[1]))' "$sim" || { rc=2; summary="giap run: giapsim failed"; return; }
  local common=(--genome "$sim/genome.fasta" --proteins "$sim/proteins.faa" --threads "$N" --memory 56.GB --profile "$PROFILE")
  rc=0; summary="giap run sim100k:"
  for variant in bakta_light bakta_full_json; do
    out="$RUN/$variant"
    if [ "$variant" = bakta_light ]; then
      giap run "${common[@]}" --outdir "$out" --set annotation.bakta.db=light
    else
      giap run "${common[@]}" --outdir "$out" --set "annotation.precomputed.bakta_json=$full"
    fi
    r=$?
    [ "$r" -gt "$rc" ] && rc=$r
    summary="$summary $variant exit $r ($(awk -F'\t' 'NR>1 {n[$6]++} END {for (s in n) printf "%s %d ", s, n[s]}' "$out/logs/trace.tsv" 2>/dev/null | xargs));"
    python scripts/measure_runtimes.py "$out/logs/trace.tsv" --sacct "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}${SLURM_ARRAY_JOB_ID:+_$SLURM_ARRAY_TASK_ID}" --genome "sim100k ($variant)" \
      | tee "$RUN/runtimes_$variant.md"
  done
}

if [ "$MODE" = record ] && [ "$TASK" != 0 ]; then
  report 0 "task $TASK: skipped (MODE=record)"; exit 0
fi
rc=2; summary="task $TASK: no such sub-task"
case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) if [ "$MODE" = record ]; then task_record; else task_integration; fi ;;
  1) task_benchmarks ;;
  2) task_calibration ;;
  3) task_giap_run ;;
esac
report "$rc" "$summary"
exit "$rc"
