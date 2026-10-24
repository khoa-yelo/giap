#!/bin/bash
#SBATCH --job-name=giap-bench
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-2%1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --constraint=CPU_SKU:6148
#SBATCH --time=22:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# The standard-run benchmark (GIAP_SPEC.md §13.2, §15.2 step 9; docs/specs/M10-12-benchmark-standard.md).
# Array task r is repeat r+1: the cells k12|sim5m x standard|core run one after another in a
# rotated order, each as its own srun step (16 CPUs, 64,000 MB = MaxMemPerCPU x 16, so Slurm
# allocates exactly 16 CPUs) under `timeout`, from a fresh launch directory, with a fresh numba
# cache, on BENCH_COMMIT checked out detached in $GIAP_ROOT/ci/bench-<sha7>. Every run writes
# RUN_DIR=$BENCH_DIR/<bench_id>/r<k>[.<attempt>]/<genome>_<variant>/{bench_run.json,runtimes.md,
# cotenancy.tsv}, recorded by $GIAP_REPO/scripts/dev/bench_record.py; wall_s is the step's sacct
# Elapsed. scripts/dev/benchmark_summary.py summarises.
# Variables: BENCH_COMMIT (required, 40 hex, on origin), BENCH_ID (<submit date>_<sha7>; required for a replacement),
# BENCH_DIR ($GIAP_ROOT/bench), BENCH_GENOMES, BENCH_VARIANTS, BENCH_ONLY=<genome>/<variant> with
# BENCH_ATTEMPT=b|c (a replacement), BENCH_RUN_TIMEOUT (5h; the cells x it fit in --time), BENCH_PROFILE (farmshare_local,
# apptainer; the Slurm-executor profile `farmshare` is refused, R14), BENCH_DRY_RUN=1 (preconditions
# and the planned lines only). Exit 0 all runs ok, 1 any failed/timeout/invalid, 2 a
# precondition, 143 on SIGTERM (the current run is written as timeout). Never pushes or commits.
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
K12_GENOME_REL=tests/data/genomes/k12_mg1655.fna.gz
K12_GENOME_SHA256=4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c
K12_JSON_REL=tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz
K12_JSON_SHA256=c08c6466edf99d16cbe79679a8971b6e06e081cfc6adb9cc6f1ea9d51d3f083d
MV="$GIAP_ROOT/mvp04/20261023"           # MVP-04 K-12 peptides and config.yaml
PROFILE="${BENCH_PROFILE:-farmshare_local,apptainer}"
TIMEOUT="${BENCH_RUN_TIMEOUT:-5h}"
GENOMES="${BENCH_GENOMES-k12 sim5m}"
VARIANTS="${BENCH_VARIANTS-standard core}"
ATT="${BENCH_ATTEMPT:-a}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"; REPEAT=$((TASK + 1))
JOBRAW="${SLURM_JOB_ID:-0}"
JOBID="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_$TASK}"; JOBID="${JOBID:-$JOBRAW}"
NODE="${SLURMD_NODENAME:-$(hostname -s)}"
NODETMP="${TMPDIR:-/tmp}"
ME="${USER:-$(id -un)}"
cur=""; spid=""; numba_tmp=""; SACCT=""; SKU=""; nsteps=0; worst=0

die() { echo "benchmark_standard: $*" >&2; exit 2; }
stop_sampler() { [ -z "$spid" ] || { kill "$spid" 2>/dev/null; wait "$spid" 2>/dev/null; }; spid=""; }
on_exit() {
  [ "${BASHPID:-$$}" = "$$" ] || return 0
  stop_sampler
  [ -z "$numba_tmp" ] || rm -rf "$numba_tmp"
}
on_term() {
  [ "${BASHPID:-$$}" = "$$" ] || exit 143   # a forked child hit by the group SIGTERM
  stop_sampler
  if [ -n "$cur" ]; then
    t1=$(date +%s.%N); SACCT=""
    cd "$W" || exit 143
    set -- $cur
    record "$1" "$2" "$3" 143 timeout >/dev/null
  fi
  echo "benchmark_standard: SIGTERM during ${cur:-set-up}; not-started runs are absent" >&2
  exit 143
}
trap on_exit EXIT
trap on_term TERM

case ",$(printf '%s' "$PROFILE" | tr -d ' \t')," in
  *,farmshare,*) die "BENCH_PROFILE=$PROFILE: the farmshare (Slurm executor) profile is refused (R14); use farmshare_local" ;;
esac
printf '%s\n' "${BENCH_COMMIT:-}" | grep -qE '^[0-9a-f]{40}$' || die "BENCH_COMMIT must be a full 40-hex SHA (got '${BENCH_COMMIT:-}')"
[ -n "$(echo $GENOMES)" ] || die "BENCH_GENOMES is empty"
[ -n "$(echo $VARIANTS)" ] || die "BENCH_VARIANTS is empty"
for g in $GENOMES; do case "$g" in k12|sim5m) ;; *) die "unknown genome '$g' (k12, sim5m)" ;; esac; done
for v in $VARIANTS; do case "$v" in standard|core) ;; *) die "unknown variant '$v' (standard, core)" ;; esac; done
case "$ATT" in a|b|c) ;; *) die "BENCH_ATTEMPT must be a, b or c (R9: 3 attempts)" ;; esac
[ -n "${BENCH_ID:-}" ] || { [ -z "${BENCH_ONLY:-}" ] && [ "$ATT" = a ]; } || die "a replacement (BENCH_ONLY or BENCH_ATTEMPT b|c) needs BENCH_ID"
case "${BENCH_ONLY:-k12/core}" in k12/standard|k12/core|sim5m/standard|sim5m/core) ;; *) die "BENCH_ONLY='$BENCH_ONLY' is not <genome>/<variant>" ;; esac
SHA7="${BENCH_COMMIT:0:7}"
if [ -z "${BENCH_ID:-}" ] && [ -n "${SLURM_ARRAY_JOB_ID:-}" ]; then
  sub=$(scontrol show job "$SLURM_ARRAY_JOB_ID" 2>/dev/null | sed -n 's/.*SubmitTime=\([0-9]\{4\}-[0-9][0-9]-[0-9][0-9]\).*/\1/p' | head -1)
fi
BENCH_ID="${BENCH_ID:-${sub:-$(date +%F)}_$SHA7}"
REP_NAME="r$REPEAT"; [ "$ATT" = a ] || REP_NAME="r$REPEAT.$ATT"
REP="${BENCH_DIR:-$GIAP_ROOT/bench}/$BENCH_ID/$REP_NAME"
IN="${BENCH_DIR:-$GIAP_ROOT/bench}/$BENCH_ID/inputs"   # the K-12 JSON (the same bytes every task); SIM5M per repeat
W="$GIAP_ROOT/ci/bench-$SHA7"

ALL=(k12/standard k12/core sim5m/standard sim5m/core)   # rotated by the task: rows of the spec's table
CELLS=(); ORDER=(); i=0
while [ "$i" -lt 4 ]; do
  c="${ALL[$(((i + TASK) % 4))]}"
  if [ -n "${BENCH_ONLY:-}" ]; then
    [ "$c" = "$BENCH_ONLY" ] && { CELLS+=("$c"); ORDER+=("$i"); }
  else
    case " $(echo $GENOMES) " in *" ${c%/*} "*)
      case " $(echo $VARIANTS) " in *" ${c#*/} "*) CELLS+=("$c"); ORDER+=("$i") ;; esac ;;
    esac
  fi
  i=$((i + 1))
done
[ "${#CELLS[@]}" -gt 0 ] || die "no cell selected"
tsec=$(printf '%s\n' "$TIMEOUT" | awk '/^[0-9]+[smhd]?$/ {n = $0 + 0; u = substr($0, length($0)); print n * (u == "d" ? 86400 : u == "h" ? 3600 : u == "m" ? 60 : 1)}')
[ -n "$tsec" ] || die "BENCH_RUN_TIMEOUT='$TIMEOUT' is not N[smhd]"
[ $((tsec * ${#CELLS[@]})) -le 79200 ] || die "${#CELLS[@]} runs x BENCH_RUN_TIMEOUT=$TIMEOUT overrun --time=22:00:00"

mkdir -p "$GIAP_ROOT/ci"
(
  flock 9
  git -C "$GIAP_REPO" fetch -q origin || exit 1
  [ -n "$(git -C "$GIAP_REPO" branch -r --contains "$BENCH_COMMIT" 2>/dev/null)" ] || exit 1
  [ -d "$W" ] || git -C "$GIAP_REPO" worktree add -q --detach "$W" "$BENCH_COMMIT" || exit 1
  git -C "$W" checkout -q -f --detach "$BENCH_COMMIT"
) 9>"$GIAP_ROOT/ci/.git.lock" || die "BENCH_COMMIT=$BENCH_COMMIT is not on origin or cannot be checked out in $W"
cd "$W" || die "cannot enter $W"
export PYTHONPATH="$W/src:$W/tests/sim" NXF_ANSI_LOG=false
giap setup --check --db-dir "$GIAP_DB_DIR" --bakta full --genomad --antismash --promotech-models \
  --go --kegg --cog --contaminants >&2 || die "giap setup --check --db-dir $GIAP_DB_DIR failed"
case " ${CELLS[*]} " in *" k12/"*)
  [ -f "$MV/peptides.tsv" ] && [ -f "$MV/config.yaml" ] || die "K-12 peptides.tsv/config.yaml missing under $MV (mvp04_inputs.sbatch)"
  [ "$(sha256sum < "$K12_GENOME_REL" | cut -d' ' -f1)" = "$K12_GENOME_SHA256" ] || die "$K12_GENOME_REL differs from the pinned sha256"
  [ "$(sha256sum < "$K12_JSON_REL" | cut -d' ' -f1)" = "$K12_JSON_SHA256" ] || die "$K12_JSON_REL differs from the pinned sha256"
  ;;
esac
for c in "${CELLS[@]}"; do
  [ ! -e "$REP/${c%/*}_${c#*/}" ] || die "$REP/${c%/*}_${c#*/} already exists (a new BENCH_ID or BENCH_ATTEMPT)"
done

plan() {  # plan <genome> <variant>: the run's command line in CMD
  local gen pep ann extra
  if [ "$1" = k12 ]; then
    gen="$W/$K12_GENOME_REL"; pep="$MV/peptides.tsv"; ann="$IN/k12_mg1655.json"; extra=(--config "$MV/config.yaml")
  else
    gen="$IN/sim5m_$REP_NAME/genome.fasta"; pep="$IN/sim5m_$REP_NAME/peptides.tsv"; ann="$IN/sim5m_$REP_NAME/bakta_like.json"; extra=()
  fi
  if [ "$2" = core ]; then extra+=(--annotation "$ann" --modules core); else extra+=(--set annotation.bakta.db=full); fi
  CMD=(srun --ntasks=1 --cpus-per-task=16 --mem=64000M --export=ALL
    timeout --signal=TERM --kill-after=5m "$TIMEOUT"
    giap run --genome "$gen" --peptides "$pep" "${extra[@]}" --db-dir "$GIAP_DB_DIR" --threads 16
    --memory 56.GB --profile "$PROFILE" --seed 20261014 --outdir "$REP/$1_$2/out")
}

echo "benchmark bench_id=$BENCH_ID repeat=$REP_NAME commit=$BENCH_COMMIT node=$NODE job=$JOBID profile=$PROFILE start=$(date -Is)"
if [ "${BENCH_DRY_RUN:-}" = 1 ]; then
  for c in "${CELLS[@]}"; do plan "${c%/*}" "${c#*/}"; printf 'PLAN %s:' "$c"; printf ' %q' "${CMD[@]}"; echo; done
  exit 0
fi

snap() {  # the image and tool-environment caches, one entry per line
  local d
  for d in "${NXF_APPTAINER_CACHEDIR:-}" "${NXF_CONDA_CACHEDIR:-}"; do
    [ -n "$d" ] && [ -d "$d" ] && ls -1A "$d" | while IFS= read -r f; do printf '%s/%s\n' "$d" "$f"; done
  done
  return 0
}
sample() {  # one cotenancy.tsv row: epoch, node CPUAlloc, CPULoad, AllocMem, our other giap-* jobs
  local n q
  n=$(scontrol show node "$NODE" 2>/dev/null | awk '{for (i = 1; i <= NF; i++) {split($i, a, "="); v[a[1]] = a[2]}}
    END {printf "%s\t%s\t%s", v["CPUAlloc"], v["CPULoad"], v["AllocMem"]}')
  q=$(squeue -u "$ME" -t R -h -o '%A|%j|%C|%N' 2>/dev/null | awk -F'|' -v me="$JOBRAW" '$1 != me && $2 ~ /^giap-/' | paste -sd';' -)
  printf '%s\t%s\t%s\n' "$(date +%s)" "$n" "$q" >> "$1"
}
sacct_step() {  # SACCT = the step's complete row; 3 tries over 2 min (accounting can lag or be down)
  local i
  for i in 1 2 3; do
    SACCT=$(sacct -j "$JOBRAW.$1" -P -n -o JobIDRaw,ElapsedRaw,MaxRSS,TotalCPU,AllocCPUS,State,ExitCode 2>/dev/null \
      | grep "^$JOBRAW\.$1|" | awk -F'|' 'NF >= 7 && $2 ~ /^[0-9]+$/ && $5 ~ /^[0-9]+$/ && $6 != ""' | head -1)
    [ -z "$SACCT" ] || break
    [ "$i" = 3 ] || sleep 60
  done
}
record() {  # record <genome> <variant> <order> <exit> [<forced status>]: RUN_DIR/bench_run.json; prints the status
  local rd="$REP/$1_$2" st rc
  st=$(python "$GIAP_REPO/scripts/dev/bench_record.py" "$rd" "$1" "$2" "$3" "$step" "$4" "$t0" "$t1" "${5:-}" \
    "$SACCT" "$BENCH_ID" "$BENCH_COMMIT" "$REPEAT" "$ATT" "$JOBID" "$TASK" "$NODE" "$SKU" "${CMD[@]}"); rc=$?
  if [ "$rc" != 0 ] || [ -z "$st" ]; then  # the minimal record: the summary names this repeat failed
    st=failed
    printf '{"schema": "giap-bench-run/1", "bench_id": "%s", "commit": "%s", "genome": "%s", "variant": "%s", "repeat": %s, "attempt": "%s", "job_id": "%s", "node": "%s", "exit_code": %s, "status": "failed", "invalid_reasons": ["record_error:bench_record.py exit %s"], "warnings": []}\n' \
      "$BENCH_ID" "$BENCH_COMMIT" "$1" "$2" "$REPEAT" "$ATT" "$JOBID" "$NODE" "$4" "$rc" > "$rd/bench_run.json"
  fi
  echo "$st"
}
run_cell() {  # run_cell <genome> <variant> <order index>
  local rd="$REP/$1_$2" rc status
  plan "$1" "$2"; mkdir -p "$rd/launch"
  numba_tmp=$(mktemp -d "$NODETMP/giap-bench-numba.XXXXXX")
  snap > "$rd/caches.before"
  sample "$rd/cotenancy.tsv"
  ( while :; do sleep 300; sample "$rd/cotenancy.tsv"; done ) & spid=$!
  t0=$(date +%s.%N); step=$nsteps; nsteps=$((nsteps + 1)); cur="$1 $2 $3"
  cd "$rd/launch" || exit 1
  NUMBA_CACHE_DIR="$numba_tmp" "${CMD[@]}"
  rc=$?; t1=$(date +%s.%N)
  cd "$W" || exit 1
  stop_sampler; sample "$rd/cotenancy.tsv"; snap > "$rd/caches.after"
  sacct_step "$step"
  status=$(record "$1" "$2" "$3" "$rc"); cur=""
  echo "RUN $1/$2 exit=$rc status=${status:-?} step=$JOBRAW.$step"
  [ "$status" = ok ] || worst=1
  python "$GIAP_REPO/scripts/measure_runtimes.py" "$rd/out/logs/trace.tsv" --sacct "$JOBRAW" --genome "$1 ($2)" > "$rd/runtimes.md" 2>/dev/null || true
  rm -rf "$rd/out/work" "$numba_tmp"; numba_tmp=""
}

SKU=$(scontrol show node "$NODE" 2>/dev/null | sed -n 's/.*ActiveFeatures=\([^ ]*\).*/\1/p' | head -1 | tr ',' '\n' | grep '^CPU_SKU:' | head -1)
mkdir -p "$IN" || die "cannot create $IN"
case " ${CELLS[*]} " in *" k12/core "*)  # shared by the repeats: skip the same bytes, else temp file + mv
  want=$(gzip -dc "$K12_JSON_REL" | sha256sum | cut -d' ' -f1)
  if [ "$(sha256sum < "$IN/k12_mg1655.json" 2>/dev/null | cut -d' ' -f1)" != "$want" ]; then
    gzip -dc "$K12_JSON_REL" > "$IN/.k12json.$$" && mv -f "$IN/.k12json.$$" "$IN/k12_mg1655.json" \
      || { rm -f "$IN/.k12json.$$"; die "cannot decompress $K12_JSON_REL"; }
  fi ;;
esac
case " ${CELLS[*]} " in *" sim5m/"*)
  python -c 'import sys; from giapsim import SIM5M, SIM5M_SEED, simulate; simulate(SIM5M, SIM5M_SEED, sys.argv[1])' \
    "$IN/sim5m_$REP_NAME" || die "giapsim SIM5M generation failed" ;;
esac
i=0
while [ "$i" -lt "${#CELLS[@]}" ]; do
  c="${CELLS[$i]}"
  run_cell "${c%/*}" "${c#*/}" "${ORDER[$i]}"
  i=$((i + 1))
done
echo "BENCH_EXIT=$worst bench_id=$BENCH_ID repeat=$REP_NAME end=$(date -Is)"
exit "$worst"
