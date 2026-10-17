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
#      (a missing DB fails, never skips), then nf-test test workflow/tests --profile $PROFILE
#   1  benchmarks pinned to the wheat node class (CPU_SKU:6148): a task that lands elsewhere
#      resubmits itself with --constraint and exits; > 15 % slower than benchmarks/baseline.json fails
#   2  the calibration_nightly tests (the pre-M2-21 nightly)
#   3  real `giap run` on sim100k (giapsim SIM100K + its proteins): Bakta with the light DB, and
#      with the recorded full-DB JSON (--annotation); measured runtimes via scripts/measure_runtimes.py
# MODE=test (default) | record: record runs the scripts/slurm/record_tools.sbatch logic in task 0
# (sim100k into a fresh directory, compared with $GIAP_ROOT/recordings: identical / changed); the
# other tasks only log "skipped". Submit: sbatch --array=0 --export=ALL,MODE=record nightly.sbatch
# Each task checks out NIGHTLY_REF (default origin/main; decision 5: main only) in its own detached
# worktree $GIAP_ROOT/ci/<ref>-t<task>; work dirs in $GIAP_ROOT/ci/runs/<job>_<task> (kept 7 days).
# One line per sub-task goes to $GIAP_ROOT/logs/nightly.tsv (date, job_task, commit, exit, summary);
# the last task to finish posts the table with scripts/ci/post_summary.py to the "Nightly
# (FarmShare)" issue (NIGHTLY_POST=dry: print it instead). Never pushes and never commits.
# Overrides: NIGHTLY_REF, NIGHTLY_PROFILE (default farmshare,apptainer; farmshare_local,apptainer
# keeps every Nextflow task inside this allocation), NIGHTLY_POST, NIGHTLY_RECORD_SAMPLES.
# Runs from scrontab (scripts/slurm/nightly.scrontab) or by hand: sbatch scripts/slurm/nightly.sbatch
# R14: 4 tasks x 16 CPUs = 64 CPUs; Nextflow queueSize 8 x 16 CPUs under -profile farmshare.
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
MODE="${MODE:-test}"
case "$MODE" in test|record) ;; *) echo "MODE must be test or record, not '$MODE'" >&2; exit 2 ;; esac
REF="${NIGHTLY_REF:-origin/main}"
PROFILE="${NIGHTLY_PROFILE:-farmshare,apptainer}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
JOB="${NIGHTLY_ARRAY_JOB:-${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}}"
EXPECT="${NIGHTLY_EXPECT:-${SLURM_ARRAY_TASK_COUNT:-1}}"
N="${SLURM_CPUS_PER_TASK:-16}"
TSV="$GIAP_ROOT/logs/nightly.tsv"
CI="$GIAP_ROOT/ci"
W="$GIAP_ROOT/ci/$(echo "${REF#origin/}" | tr '/' '-')-t$TASK"
RUN="$CI/runs/${JOB}_$TASK"
mkdir -p "$CI/runs" "$GIAP_ROOT/logs"
find "$CI/runs" -mindepth 1 -maxdepth 1 -mtime +7 -exec rm -rf {} + 2>/dev/null
rm -rf "$RUN"; mkdir -p "$RUN/tmp"
(
  flock 9   # the array tasks share one repository: fetch and worktree add one at a time
  git -C "$GIAP_REPO" fetch -q origin
  [ -d "$W" ] || git -C "$GIAP_REPO" worktree add -q --detach "$W" "$REF"
  git -C "$W" checkout -q -f --detach "$REF"
) 9>"$CI/.git.lock" || { echo "cannot check out $REF in $W" >&2; exit 2; }
cd "$W" || exit 2
export TMPDIR="$RUN/tmp" NXF_ANSI_LOG=false NUMBA_CACHE_DIR="$(mktemp -d /tmp/giap-numba.XXXXXX)"
trap 'rm -rf "$NUMBA_CACHE_DIR"' EXIT
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
commit=$(git -C "$W" rev-parse --short HEAD)
echo "nightly task=$TASK mode=$MODE ref=$REF commit=$commit node=$(hostname) profile=$PROFILE start=$(date -Is)"

report() {  # report <exit> <sub-task: summary>: one nightly.tsv line, then the summary post
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "${JOB}_$TASK" "$commit" "$1" "$2" >> "$TSV"
  echo "NIGHTLY_EXIT=$1 $2"
  local dry=(); [ "${NIGHTLY_POST:-yes}" = dry ] && dry=(--dry-run)
  python scripts/ci/post_summary.py --job "$JOB" --tsv "$TSV" --expect "$EXPECT" "${dry[@]}"
}
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
  nf-test test workflow/tests --profile "$PROFILE" 2>&1 | tee "$RUN/nf-test.txt"
  local r1=${PIPESTATUS[0]}
  rc=$(( r0 > r1 ? r0 : r1 ))
  summary="$s0; nf-test: $(grep -E '^(SUCCESS|FAILURE):' "$RUN/nf-test.txt" | tail -1)"
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
  summary="record: failures $nfail; identical $same; changed ${#changed[@]}${changed:+ (${changed[*]})}"
}

task_benchmarks() {
  case "$(hostname -s)" in
    wheat-*) ;;
    *) if [ -z "${NIGHTLY_BENCH_PINNED:-}" ]; then   # array tasks cannot carry their own constraint
         local j; j=$(sbatch --parsable --array="$TASK" --constraint=CPU_SKU:6148 \
           --export="ALL,MODE=$MODE,NIGHTLY_BENCH_PINNED=1,NIGHTLY_ARRAY_JOB=$JOB,NIGHTLY_EXPECT=$EXPECT" \
           scripts/slurm/nightly.sbatch)
         echo "not a wheat node ($(hostname -s)): benchmarks resubmitted as job $j"; exit 0
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
      giap run "${common[@]}" --outdir "$out" --annotation "$full"
    fi
    r=$?
    [ "$r" -gt "$rc" ] && rc=$r
    summary="$summary $variant exit $r ($(awk -F'\t' 'NR>1 {n[$6]++} END {for (s in n) printf "%s %d ", s, n[s]}' "$out/logs/trace.tsv" 2>/dev/null | xargs));"
    python scripts/measure_runtimes.py "$out/logs/trace.tsv" --sacct "${SLURM_JOB_ID:-0}" --genome "sim100k ($variant)" \
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
