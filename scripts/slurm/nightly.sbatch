#!/bin/bash
#SBATCH --job-name=gippa-nightly
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --signal=B:TERM@600
#SBATCH --array=0-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=06:00:00
#SBATCH --exclude=rye-02
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# --exclude=rye-02 is a stopgap: on that node Python/JVM start-up stalls 50-250x (nightly 1870584,
# 1896018); a slow-node guard replaces it.
# Nightly job array on the latest origin/main (GIPPA_SPEC.md §12; M2-21 lead decisions 3-5):
#   0  tests/integration + tests/contract with --require-capabilities <ci/gates.yaml farmshare_nightly>
#      (a missing DB fails, never skips), then nf-test test workflow/tests --profile $NF_PROFILE
#      (default test,stub: the suites assert fake-tool argv logs and stub outputs; under
#      farmshare,apptainer the pinned images replace the fakes and 9 of 30 tests fail, job 1852060)
#   1  benchmarks pinned to the wheat node class (CPU_SKU:6148): a task that lands elsewhere
#      resubmits itself with --constraint and exits; > 15 % slower than benchmarks/baseline.json fails
#   2  the whole fast suite, i.e. the fast calibration tier including calibration_nightly, with the
#      coverage gate (ci/gates.yaml coverage_min): moved here from PR CI by ADR-0052
#   3  real `gippa run` on sim100k (gippasim SIM100K + its proteins): Bakta with the light DB, and
#      with the recorded full-DB JSON (annotation.precomputed.bakta_json; its GBFF feeds
#      antiSMASH); measured runtimes via scripts/measure_runtimes.py
# MODE=test (default) | record: record runs the scripts/slurm/record_tools.sbatch logic in task 0
# (sim100k into a fresh directory, compared with $GIPPA_ROOT/recordings: identical / changed); the
# other tasks only log "skipped". Submit: sbatch --array=0 --export=ALL,MODE=record nightly.sbatch
# Each task checks out NIGHTLY_REF (default origin/main; decision 5: main only) in its own detached
# worktree $GIPPA_ROOT/ci/<ref>-t<task>; work dirs in $GIPPA_ROOT/ci/runs/<job>_<task> (kept 7 days).
# Every task writes one line to $GIPPA_ROOT/logs/nightly.tsv (date, job_task, commit, exit,
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
# Slow-node guard: a task whose `python -c 'import gippa'` takes longer than NIGHTLY_SLOW_NODE_SECS
# (default 60) resubmits itself with --exclude=<NIGHTLY_EXCLUDE, the #SBATCH excludes, this host>
# and writes no row; at most NIGHTLY_SLOW_MAX_HOPS (default 2) hops; NIGHTLY_SLOW_NODE_GUARD=0 off.
# --signal=B:TERM@600 and `cmd & wait $!` for the long steps: the TERM trap writes the row in time.
# R14: 4 tasks x 16 CPUs = 64 CPUs + the 1-CPU summary job; no Nextflow child jobs (farmshare_local).
# Task 2 shards: NIGHTLY_CAL_SHARDS=N (default 1; nightly_launch.sh sets 2) splits the fast suite
# into N file shards (scripts/ci/shard_tests.py): shard 1 (the array task) submits shards 2..N as
# 16-CPU jobs on the same commit (own worktree <W>-s<k>, run dir <job>_2-s<k>, summary extended),
# each records its result and coverage data in $CI/runs/<job>_2.shards, and the last to finish
# writes the single calibration_nightly row (scripts/ci/nightly_shards.py: counts summed, largest
# exit) after `coverage combine` + the coverage_min gate. R14: + 16 CPUs per extra shard (81 at
# N=2). The 2026-10-24/27 unsharded runs took 5 h 43 / > 5 h 50 of the 6 h limit on rye-01.
# Output is -v, so a killed shard's row counts the FAILED lines so far ($RUN/failed_so_far.txt).
# N > 4 is refused (R14: 64 + 16 (N-1) + 1 <= 128). A shard that never records (SIGKILL, OOM,
# cancelled) is recorded as failed by nightly_summary.sh, which then writes the failing row.
# Shards 2..N run the tested commit's nightly.sbatch (submitted from $W), shard 1 the main
# checkout's spooled copy, as the benchmark resubmission: the protocol is the same unless the main
# checkout is stale.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
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
SHARD="${NIGHTLY_CAL_SHARD:-1}"; SHARDS_N="${NIGHTLY_CAL_SHARDS:-1}"
SHARDS="$GIPPA_ROOT/ci/runs/${JOB}_2.shards"
sharded=""   # set once this task-2 run is one of several shards (shard 1: after submitting the rest)
TSV="$GIPPA_ROOT/logs/nightly.tsv"
CI="$GIPPA_ROOT/ci"
W="${NIGHTLY_WORKTREE:-$GIPPA_ROOT/ci/$(echo "${REF#origin/}" | tr '/' '-')-t$TASK}"
RUN="$CI/runs/${JOB}_$TASK"
if [ "$TASK" = 2 ] && [ "$SHARD" != 1 ]; then RUN="$RUN-s$SHARD"; sharded=1; fi
commit=$(git -C "$GIPPA_REPO" rev-parse -q --short --verify "$REF^0" 2>/dev/null || echo unknown)
stage="setup"; reported=""; numba_tmp=""; tmp_dir=""
# M11-03b Q5: task 2's early-exit and killed rows are failing calibration_nightly rows, so readiness
# R3 cannot fall back to an older pass; the MODE=record skip row keeps "task 2" (exit 0)
trap_label="task $TASK:"
if [ "$TASK" = 2 ] && [ "$MODE" != record ]; then trap_label="calibration_nightly: task $TASK"; fi
mkdir -p "$CI/runs" "$GIPPA_ROOT/logs"

tsv_line() {  # tsv_line <job_task> <exit> <summary>
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "$1" "$commit" "$2" "$3" >> "$TSV"
}
report() {  # report <exit> <sub-task: summary>: this task's nightly.tsv line, then the summary post
  reported=1
  tsv_line "${JOB}_$TASK" "$1" "$2"
  echo "NIGHTLY_EXIT=$1 $2"
  local dry=() ps="$W/scripts/ci/post_summary.py"
  [ -f "$ps" ] || ps="$GIPPA_REPO/scripts/ci/post_summary.py"   # the checkout may have failed
  [ "${NIGHTLY_POST:-yes}" = dry ] && dry=(--dry-run)
  python "$ps" --job "$JOB" --tsv "$TSV" --expect "$EXPECT" ${dry[@]+"${dry[@]}"}
}
finish() {  # finish <exit> <summary>: the task's row, or (a task-2 shard) its shard record
  [ -n "$sharded" ] || { report "$@"; return; }
  local ns="$W/scripts/ci/nightly_shards.py" merged mrc msum
  [ -f "$ns" ] || ns="$GIPPA_REPO/scripts/ci/nightly_shards.py"
  echo "shard $SHARD/$SHARDS_N: exit $1 $2"
  mkdir -p "$SHARDS"
  # reported only once recorded (or the row written): SIGTERM during combine still gets a row
  merged=$( { flock 7 &&
    python "$ns" record --dir "$SHARDS" --shard "$SHARD" --exit "$1" --summary "$2" ||
    { echo "shard $SHARD/$SHARDS_N: record failed (nightly_summary.sh records it)" >&2; exit 4; }
    python "$ns" merge --dir "$SHARDS" --of "$SHARDS_N"; } 7>"$SHARDS/.lock" )
  case $? in 0) ;; 3) reported=1; return 0 ;; *) return 0 ;; esac
  mrc=${merged%%$'\t'*}; msum=${merged#*$'\t'}
  coverage_gate; report "$(( mrc > cov_rc ? mrc : cov_rc ))" "$msum$cov_note"
}
coverage_gate() {  # the last shard: combine the shards' coverage data, then ci/gates.yaml coverage_min
  local min k files=() rc_file="$SHARDS/coveragerc" data="$SHARDS/combined.coverage"
  local txt="$SHARDS/coverage_report.txt"   # neither may start with coverage.: report combines <data>.*
  cov_rc=0; cov_note=""
  min=$(awk '/^coverage_min:/ {print $2*100}' ci/gates.yaml 2>/dev/null)
  for k in $(seq 1 "$SHARDS_N"); do [ -s "$SHARDS/coverage.s$k" ] && files+=("$SHARDS/coverage.s$k"); done
  if [ -z "$min" ] || [ "${#files[@]}" != "$SHARDS_N" ]; then
    cov_rc=1; cov_note="; coverage not measured (${#files[@]} of $SHARDS_N shard data files)"; return
  fi
  # ignore_errors: a module a test wrote under its tmp dir (test_steps_review's zz_review_toy.py)
  # is gone by report time, and "No source for code" failed the gate (run 1905610)
  printf '[paths]\ngippa =\n    %s/src/gippa\n    */src/gippa\n[report]\nignore_errors = True\n' \
    "$PWD" > "$rc_file"
  rm -f "$data"
  python -m coverage combine --keep --rcfile="$rc_file" --data-file="$data" "${files[@]}" &&
    python -m coverage report --keep-combined --rcfile="$rc_file" --data-file="$data" \
      --skip-covered --show-missing --fail-under="$min" > "$txt"
  cov_rc=$?
  cov_note="; coverage $(awk '$1 == "TOTAL" {print $NF}' "$txt" 2>/dev/null) (min $min%$( [ "$cov_rc" = 0 ] || echo ", FAIL"))"
}
failed_so_far() {  # on SIGTERM: the FAILED/ERROR lines pytest -v printed so far
  local f="$RUN/calibration_nightly.txt" n
  [ "$TASK" = 2 ] && [ -f "$f" ] || return 0
  grep -oE '(FAILED|ERROR) tests/[^ ]+' "$f" | sort -u > "$RUN/failed_so_far.txt"
  n=$(wc -l < "$RUN/failed_so_far.txt")
  echo "; $n failed so far ($RUN/failed_so_far.txt)"
}
on_exit() {  # decision 14: a task that ends without its line (early exit, error) still reports
  local rc=$?
  [ "${BASHPID:-$$}" = "$$" ] || return 0   # a forked child hit by the group SIGTERM before its exec (bash < 4 has no BASHPID: old behaviour)
  [ -n "$reported" ] || finish "$(( rc ? rc : 1 ))" "$trap_label ended during $stage without a result (exit $rc)"
  [ -z "$numba_tmp" ] || rm -rf "$numba_tmp"   # only this task's cache, never the shared one
  [ -z "$tmp_dir" ] || rm -rf "$tmp_dir"       # this task's node-local TMPDIR
}
trap on_exit EXIT
# decision 23: a row only if none was written yet (SIGTERM may land inside report itself)
trap '[ "${BASHPID:-$$}" = "$$" ] || exit 143; [ -n "$reported" ] || finish 143 "$trap_label killed by SIGTERM (time limit?) during $stage$(failed_so_far)"; exit 143' TERM

case "$MODE" in test|record) ;; *) echo "MODE must be test or record, not '$MODE'" >&2; exit 2 ;; esac
case ",$PROFILE," in
  *,farmshare,*)   # decision 22: Slurm-executor child jobs would break the R14 CPU budget
    echo "NIGHTLY_PROFILE=$PROFILE: the farmshare (Slurm executor) profile is refused (R14); use farmshare_local" >&2
    finish 2 "$trap_label NIGHTLY_PROFILE=$PROFILE refused (R14: use farmshare_local)"; exit 2 ;;
esac
if [ "$TASK" = 2 ] && [ "$MODE" != record ]; then
  case "$SHARDS_N" in
    1|2|3|4) ;;
    *) echo "NIGHTLY_CAL_SHARDS=$SHARDS_N: refused (R14: 1-4 shards of 16 CPUs)" >&2
       finish 2 "$trap_label NIGHTLY_CAL_SHARDS=$SHARDS_N refused (R14: at most 4)"; exit 2 ;;
  esac
fi
find "$CI/runs" -mindepth 1 -maxdepth 1 -mtime +7 -exec rm -rf {} + 2>/dev/null
find "$CI" -maxdepth 1 -name 'summary-deps-*' -mtime +7 -delete 2>/dev/null
rm -rf "$RUN"; mkdir -p "$RUN"
stage="checkout of $REF"
(
  flock 9   # the array tasks share one repository: fetch and worktree add one at a time
  git -C "$GIPPA_REPO" fetch -q origin
  target=$(git -C "$GIPPA_REPO" rev-parse -q --verify "$REF^0") || exit 1
  if [ -d "$W" ]; then
    # a killed task can leave index.lock behind; under this flock no other git runs in $W
    gd=$(git -C "$W" rev-parse --absolute-git-dir 2>/dev/null) && rm -f "$gd/index.lock"
    # skip checkout -f only for a clean tree at the target: the row must name the tested tree
    if { [ "$(git -C "$W" rev-parse -q --verify HEAD 2>/dev/null)" != "$target" ] ||
         ! git -C "$W" diff --quiet HEAD --; } &&
       ! git -C "$W" checkout -q -f --detach "$target"; then
      echo "checkout in $W failed: removing and re-adding the worktree" >&2
      git -C "$GIPPA_REPO" worktree remove --force "$W" || rm -rf "$W"
      git -C "$GIPPA_REPO" worktree prune
    fi
  fi
  [ -d "$W" ] || git -C "$GIPPA_REPO" worktree add -q --detach "$W" "$target"
) 9>"$CI/.git.lock" || { echo "cannot check out $REF in $W" >&2; exit 2; }
cd "$W" || exit 2
numba_tmp="$(mktemp -d /tmp/gippa-numba.XXXXXX)"
tmp_dir="$(mktemp -d /tmp/gippa.XXXXXX)"   # node-local, as test_branch.sbatch (not NFS under $RUN)
export TMPDIR="$tmp_dir" NXF_ANSI_LOG=false NUMBA_CACHE_DIR="$numba_tmp"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export GIPPA_PROBE_TIMEOUT="${GIPPA_PROBE_TIMEOUT:-180}"   # tests/conftest.py's tool probes (cold NFS)
commit=$(git -C "$W" rev-parse --short HEAD)
stage="task $TASK"
echo "nightly task=$TASK mode=$MODE ref=$REF commit=$commit node=$(hostname) profile=$PROFILE start=$(date -Is)"

pyt() {  # pyt <label> <pytest args...>: sets rc and summary
  # in the background, so the TERM trap runs at once (bash defers traps during a foreground wait)
  { python -m pytest -p no:cacheprovider "${@:2}" 2>&1 | tee "$RUN/$1.txt"; exit "${PIPESTATUS[0]}"; } &
  wait $!; rc=$?
  summary="$1: $(grep -E '[0-9]+ (passed|failed|error)' "$RUN/$1.txt" | tail -1 | tr -d '=' | xargs)"
}

task_integration() {
  local caps r0 s0 s
  # Warm the gippa-wf JVM first, as test_branch.sbatch (M9-13 decision 5): a cold `nextflow
  # -version` from NFS outlasts tests/conftest.py's 30 s probe and every needs_nextflow test skips.
  # Not `nf-test version`: decision 23 pins that nf-test never starts before pytest ends.
  s=$SECONDS; nextflow -version >/dev/null 2>&1; echo "warm-up: nextflow -version rc=$? $((SECONDS - s))s"
  caps=$(python -c 'import yaml; print(",".join(yaml.safe_load(open("ci/gates.yaml"))["required_capabilities"]["farmshare_nightly"]))')
  pyt integration tests/integration tests/contract -q -n 8 --require-capabilities "$caps"
  r0=$rc; s0=$summary
  # docker_smoke.nf.test and docker_db.nf.test run real images (docker-smoke.yml), as in ci.yml
  { nf-test test $(find workflow/tests -name '*.nf.test' ! -name docker_smoke.nf.test \
    ! -name docker_db.nf.test | sort) --profile "$NF_PROFILE" 2>&1 | tee "$RUN/nf-test.txt"
    exit "${PIPESTATUS[0]}"; } &
  local r1; wait $!; r1=$?
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
         "$(grep -vE '\./(run\.log|runtime\.txt)$' "$GIPPA_ROOT/recordings/$ts/SHA256SUMS" 2>/dev/null)" ]; then
      same=$((same + 1))
    else changed+=("$ts"); fi
  done
  rc=$(( nfail > 0 ? 1 : 0 ))
  summary="record: failures $nfail; identical $same; changed ${#changed[@]}${changed[0]+ (${changed[*]})}"
}

extend_summary() {  # extend_summary <job>: decision 19, the night's summary also waits for <job>
  local sid
  sid=$(squeue -u "${USER:-$(id -un)}" -h -t PENDING -n gippa-nightly-summary -o '%i %k' 2>/dev/null \
    | awk -v c="nightly:$JOB" '$2 == c {print $1; exit}')
  [ -n "$sid" ] || return 1
  # scontrol replaces the whole dependency list: keep every hop of the night in one file (flock)
  (
    flock 8
    echo "$1" >> "$CI/summary-deps-$JOB"
    scontrol update JobId="$sid" Dependency="afterany:$JOB:$(paste -sd: "$CI/summary-deps-$JOB")"
  ) 8>"$CI/summary-deps-$JOB.lock"
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
           "benchmarks resubmitted as job $j; no pending gippa-nightly-summary job of $JOB to extend: it may post before them"
         reported=1; exit 0
       fi
       rc=1; summary="benchmarks: not on a wheat node ($(hostname -s))"; return ;;
  esac
  pyt benchmarks benchmarks --benchmark-only -q -p no:xdist --benchmark-json "$RUN/benchmarks.json"
  summary="$summary (wheat $(hostname -s))"
}

# ADR-0052: the fast calibration tier (with calibration_nightly) and the coverage gate
task_calibration() {
  local min k j f sel=(tests); min=$(awk '/^coverage_min:/ {print $2*100}' ci/gates.yaml)
  if [ "$SHARDS_N" -gt 1 ]; then
    if [ "$SHARD" = 1 ]; then   # once per night: clear the shard records, submit shards 2..N
      rm -rf "$SHARDS"; mkdir -p "$SHARDS"; sharded=1
      printf '%s\t%s\n' "$SHARDS_N" "$commit" > "$SHARDS/plan"   # nightly_summary.sh: N and commit
      echo "1 ${SLURM_JOB_ID:-}" > "$SHARDS/jobs"
      for k in $(seq 2 "$SHARDS_N"); do
        if j=$(sbatch --parsable --array="$TASK" \
          --export="ALL,MODE=$MODE,NIGHTLY_ARRAY_JOB=$JOB,NIGHTLY_EXPECT=$EXPECT,NIGHTLY_CAL_SHARD=$k,NIGHTLY_CAL_SHARDS=$SHARDS_N,NIGHTLY_SLOW_HOPS=0,NIGHTLY_REF=$commit,NIGHTLY_WORKTREE=$W-s$k" \
          scripts/slurm/nightly.sbatch); then
          j=${j%%;*}; echo "task 2 shard $k/$SHARDS_N submitted as job $j"; echo "$k $j" >> "$SHARDS/jobs"
          extend_summary "$j" || tsv_line "${JOB}_note" 0 \
            "task 2 shard $k/$SHARDS_N is job $j; no pending gippa-nightly-summary job of $JOB to extend: it may post before it"
        else   # the night's row still comes: a refused shard is a failed shard
          python scripts/ci/nightly_shards.py record --dir "$SHARDS" --shard "$k" --exit 2 \
            --summary "calibration_nightly: shard submission refused by sbatch"
        fi
      done
    fi
    mkdir -p "$SHARDS"; stage="task 2 shard $SHARD/$SHARDS_N"
    sel=(); while IFS= read -r f; do sel+=("$f"); done < <(python scripts/ci/shard_tests.py "$SHARD" "$SHARDS_N" "$PWD")
    [ "${#sel[@]}" -gt 0 ] || { rc=2; summary="calibration_nightly: no test files in shard $SHARD/$SHARDS_N"; return; }
    export COVERAGE_FILE="$SHARDS/coverage.s$SHARD"
    pyt calibration_nightly "${sel[@]}" -v -n "$N" \
      -m "not slow and not calibration_full and not e2e" --cov=gippa --cov-report= \
      --durations=30 --junitxml="$RUN/junit.xml"
    return
  fi
  pyt calibration_nightly "${sel[@]}" -v -n "$N" -m "not slow and not calibration_full and not e2e" \
    --cov=gippa --cov-report=term-missing:skip-covered --cov-fail-under="$min" --durations=30 \
    --junitxml="$RUN/junit.xml"
}

task_gippa_run() {
  local sim="$RUN/sim" full="$GIPPA_ROOT/recordings/bakta_full/sim100k/bakta/sim100k.json" r out
  python -c 'import dataclasses, pathlib, sys, gippasim; from gippasim.proteomics import ProteomicsConfig
gippasim.simulate(dataclasses.replace(gippasim.SIM100K, proteomics=ProteomicsConfig()), gippasim.SIM100K_SEED, pathlib.Path(sys.argv[1]))' "$sim" || { rc=2; summary="gippa run: gippasim failed"; return; }
  local common=(--genome "$sim/genome.fasta" --proteins "$sim/proteins.faa" --threads "$N" --memory 56.GB --profile "$PROFILE")
  rc=0; summary="gippa run sim100k:"
  for variant in bakta_light bakta_full_json; do
    out="$RUN/$variant"
    if [ "$variant" = bakta_light ]; then
      gippa run "${common[@]}" --outdir "$out" --set annotation.bakta.db=light &
    else
      gippa run "${common[@]}" --outdir "$out" --set "annotation.precomputed.bakta_json=$full" &
    fi
    wait $!; r=$?
    [ "$r" -gt "$rc" ] && rc=$r
    summary="$summary $variant exit $r ($(awk -F'\t' 'NR>1 {n[$6]++} END {for (s in n) printf "%s %d ", s, n[s]}' "$out/logs/trace.tsv" 2>/dev/null | xargs));"
    python scripts/measure_runtimes.py "$out/logs/trace.tsv" --sacct "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}${SLURM_ARRAY_JOB_ID:+_$SLURM_ARRAY_TASK_ID}" --genome "sim100k ($variant)" \
      | tee "$RUN/runtimes_$variant.md"
  done
}

if [ "$MODE" = record ] && [ "$TASK" != 0 ]; then
  report 0 "task $TASK: skipped (MODE=record)"; exit 0
fi
slow_node_guard() {  # resubmits this task off a node where Python start-up stalls (rye-02, 10-26)
  local secs="${NIGHTLY_SLOW_NODE_SECS:-60}" hops="${NIGHTLY_SLOW_HOPS:-0}" t0 took host excl j
  [ "${NIGHTLY_SLOW_NODE_GUARD:-1}" = 0 ] && return 0
  t0=$(date +%s%N)
  timeout "$((secs + 60))" python -c 'import gippa' >/dev/null 2>&1   # timed whether or not it imports
  took=$(( ($(date +%s%N) - t0) / 1000000 ))
  echo "slow-node probe: import gippa took ${took} ms on $(hostname -s) (limit ${secs} s, hop $hops)"
  [ "$took" -le $((secs * 1000)) ] && return 0
  if [ "$hops" -ge "${NIGHTLY_SLOW_MAX_HOPS:-2}" ]; then
    slow_note=" [slow node $(hostname -s): import gippa ${took} ms; hop limit reached]"; return 0
  fi
  host=$(hostname -s)
  # --exclude on the command line replaces the #SBATCH one: carry those hosts over
  excl=$(printf '%s\n' ${NIGHTLY_EXCLUDE:+${NIGHTLY_EXCLUDE//,/ }} \
           $(grep -oP '^#SBATCH\s+(--exclude[= ]|-x\s*)\K\S+' scripts/slurm/nightly.sbatch | tr ',' ' ') \
           "$host" | awk 'NF && !seen[$0]++' | paste -sd,)
  # NIGHTLY_EXCLUDE rides in the environment (--export=ALL): sbatch splits --export at commas
  local pin=()   # a pinned benchmark task keeps its wheat constraint across the hop
  [ -n "${NIGHTLY_BENCH_PINNED:-}" ] && pin=(--constraint=CPU_SKU:6148)
  if ! j=$(NIGHTLY_EXCLUDE="$excl" sbatch --parsable --array="$TASK" --exclude="$excl" "${pin[@]}" \
    --export="ALL,MODE=$MODE,NIGHTLY_ARRAY_JOB=$JOB,NIGHTLY_EXPECT=$EXPECT,NIGHTLY_SLOW_HOPS=$((hops + 1)),NIGHTLY_REF=$commit,NIGHTLY_WORKTREE=$W" \
    scripts/slurm/nightly.sbatch); then
    slow_note=" [slow node $host: import gippa ${took} ms; the resubmission failed]"; return 0
  fi
  j=${j%%;*}
  echo "slow node $host (import gippa ${took} ms > ${secs} s): task $TASK resubmitted as job $j, excluding $excl"
  extend_summary "$j" || tsv_line "${JOB}_note" 0 \
    "task $TASK resubmitted off slow node $host as job $j; no pending gippa-nightly-summary job of $JOB to extend: it may post before it"
  reported=1; exit 0
}

slow_note=""
stage="slow-node probe"; slow_node_guard; stage="task $TASK"
rc=2; summary="task $TASK: no such sub-task"
case "${SLURM_ARRAY_TASK_ID:-0}" in
  0) if [ "$MODE" = record ]; then task_record; else task_integration; fi ;;
  1) task_benchmarks ;;
  2) task_calibration ;;
  3) task_gippa_run ;;
esac
finish "$rc" "$summary$slow_note"
exit "$rc"
