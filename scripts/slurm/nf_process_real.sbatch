#!/bin/bash
#SBATCH --job-name=gippa-nf-process-real
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-8%4
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Each Nextflow tool process run ALONE with its real tool, one array task per process, so one failure
# never blocks the others and a rerun takes minutes. Task i runs workflow/tests/real/<tool>.nf (one
# process of workflow/modules/local/<tool>/main.nf) with -profile farmshare_local,apptainer (pinned
# images; conda-only tools in the prebuilt envs, prebuilt_envs.config) and compares its outputs with
# the recording (scripts/dev/nf_real_compare.py documents the normalised comparisons).
#   0 BAKTA  1 GENOMAD  2 ANTISMASH  3 TERMINATORNET  4 TRANSTERMHP  5 OSTIR  6 PROMOTECH
#   7 PROMOTER_CALCULATOR  8 UNIOP       (BACTERMFINDER: no env built, not covered)
# Inputs: $GIPPA_ROOT/recordings/bakta_full/<sample> (preflight genome + replicons, Bakta JSON/GBFF);
# genes.coords, ostir_input.csv and TU regions rebuilt by scripts/dev/record_tools_inputs.py as
# record_tools.sbatch did. PROMOTECH windows (promoter_windows step, M2-13) and the UNIOP GFF
# (uniop_input step, M2-12) are rebuilt by its promotech-windows / uniop subcommands, which run the
# merged gippa steps of this checkout, and compared with $GIPPA_ROOT/recordings/{promotech,uniop}/<sample>;
# sim100k keeps the tests/data/recorded fixtures (their inputs predate the bakta_full sim100k recording).
# NF_REAL_RECORD=1 (PROMOTECH, UNIOP): with no recording yet, the run itself is saved as the recording
# (inputs, outputs, command.txt, versions.txt, runtime.txt, SHA256SUMS, .giap-complete); a second
# run then checks it reproduces. Existing recordings are never touched. PROMOTECH compares only
# promotech_NNN/sequences_predictions.csv (sequences.log.txt holds elapsed times).
# Usage (from the repo or worktree root; the code under test is that checkout):
#   sbatch scripts/slurm/nf_process_real.sbatch [sample]          # default sim100k; all 9 processes
#   sbatch --array=4 --export=ALL,NF_REAL_JOB=<array job> scripts/slurm/nf_process_real.sbatch
#       reruns one process (here TRANSTERMHP) into that job's directory, replacing its verdict
#   sbatch --array=6,8 --export=ALL,NF_REAL_RECORD=1 scripts/slurm/nf_process_real.sbatch k12_mg1655
#       records PROMOTECH (E. coli: 10 chunks of 100k windows, ~4 min) and UNIOP (~10 s)
# Each task writes "PROCESS <name>: ok|mismatch|failed (<reason>)" to its log and to
# $GIPPA_ROOT/scratch/nf-process-real/<array job>/verdicts/; summary.txt there lists all of them.
# Work dirs live under that directory (never node-local /tmp); $GIPPA_DB_DIR is only read.
# R14: 64G is billed as 18 CPUs (MaxMemPerCPU), so %4 keeps the array at 72 CPUs; check squeue first.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
REPO="${NF_REAL_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"   # sbatch runs a spooled copy: use the submit dir
[ -d "$REPO/workflow/tests/real" ] || { echo "no workflow/tests/real in $REPO: submit from the repo root" >&2; exit 2; }
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH" NXF_ANSI_LOG=false
S="${1:-sim100k}"
PROCS=(BAKTA GENOMAD ANTISMASH TERMINATORNET TRANSTERMHP OSTIR PROMOTECH PROMOTER_CALCULATOR UNIOP)
P="${PROCS[${SLURM_ARRAY_TASK_ID:-0}]}"; t=$(printf '%s' "$P" | tr '[:upper:]' '[:lower:]')
REC="$GIPPA_ROOT/recordings"; PRE="$REC/bakta_full/$S/preflight/work/preflight"; B="$REC/bakta_full/$S/bakta"
FIX="$REPO/tests/data/recorded"
J="$GIPPA_ROOT/scratch/nf-process-real/${NF_REAL_JOB:-${SLURM_ARRAY_JOB_ID:-manual}}"
T="$J/$t"; rm -rf "$T"; mkdir -p "$T/in" "$T/want" "$J/verdicts" "$T/tmp" && cd "$T" || exit 1
export TMPDIR="$T/tmp"
INPUTS="python $REPO/scripts/dev/record_tools_inputs.py"
GIPPA_PY="env PYTHONPATH=$REPO/src $INPUTS"   # subcommands that run this checkout's gippa steps
echo "process=$P sample=$S repo=$REPO ($(git -C "$REPO" rev-parse --short HEAD)) work=$T"

verdict() {
  local line="PROCESS $P: $1"
  echo "$line" | tee "$J/verdicts/$(printf %02d "${SLURM_ARRAY_TASK_ID:-0}")_$t.txt"
  cat "$J"/verdicts/*.txt > "$J/summary.txt"
  echo "--- summary so far ($J/summary.txt)"; cat "$J/summary.txt"
  exit 0
}
[ -f "$REC/bakta_full/$S/.giap-complete" ] || verdict "failed (no complete Bakta recording for $S)"

want=("$T/want")   # the recorded outputs this run must reproduce (nf_real_compare.py layout)
case "$P" in
  BAKTA) args=(--genome "$PRE/genome.norm.fa" --replicons "$PRE/replicons.tsv" --bakta_db "$GIPPA_DB_DIR/bakta/db"); want=("$B") ;;
  GENOMAD) args=(--genome "$PRE/genome.norm.fa" --db_dir "$GIPPA_DB_DIR"); want=("$REC/genomad/$S") ;;
  ANTISMASH) args=(--gbff "$REC/antismash/$S/bakta.gbff" --db_dir "$GIPPA_DB_DIR"); want=("$REC/antismash/$S") ;;
  TERMINATORNET) args=(--genome "$PRE/genome.norm.fa"); ln -s "$REC/terminatornet/$S/terminatornet.tsv" want/ ;;
  TRANSTERMHP) $INPUTS coords "$B/$S.json" "$T/in/genes.coords" || verdict "failed (record_tools_inputs.py coords)"
    args=(--genome "$PRE/genome.norm.fa" --coords "$T/in/genes.coords"); ln -s "$REC/transtermhp/$S/transtermhp.tt" want/ ;;
  OSTIR) $INPUTS ostir "$B/$S.json" "$T/in/ostir_input.csv" || verdict "failed (record_tools_inputs.py ostir)"
    args=(--csv "$T/in/ostir_input.csv"); ln -s "$REC/ostir/$S/ostir.csv" want/ ;;
  PROMOTECH) if [ "$S" = sim100k ]; then
      d=$(ls -d "$FIX"/promotech/*/*/ | head -1)
      args=(--windows "$d/windows_000.fa" --db_dir "$GIPPA_DB_DIR")
      mkdir want/promotech_000; ln -s "$d/sequences_predictions.csv" want/promotech_000/
    else
      $GIPPA_PY promotech-windows "$B/$S.json" "$PRE/genome.norm.fa" "$T/in" || verdict "failed (record_tools_inputs.py promotech-windows)"
      args=(--windows "$T/in/windows/windows_*.fa" --db_dir "$GIPPA_DB_DIR")
      for x in "$REC/promotech/$S"/promotech_*/sequences_predictions.csv; do   # not the timed log
        [ -f "$x" ] && c=$(basename "$(dirname "$x")") && mkdir "want/$c" && ln -s "$x" "want/$c/"
      done
      recdir="$REC/promotech/$S"
    fi ;;
  PROMOTER_CALCULATOR) $INPUTS regions "$B/$S.json" "$T/in/regions" || verdict "failed (record_tools_inputs.py regions)"
    args=(--regions "$T/in/regions"); want=("$REC/promoter_calculator/$S/csv") ;;
  UNIOP) if [ "$S" = sim100k ]; then
      d=$(ls -d "$FIX"/uniop/*/sim100k)
      args=(--gff "$d/uniop_cds.gff"); ln -s "$d/uniop.pred" want/
    else
      $GIPPA_PY uniop "$B/$S.json" "$PRE/genome.norm.fa" "$T/in" || verdict "failed (record_tools_inputs.py uniop)"
      args=(--gff "$T/in/uniop_cds.gff")
      [ -f "$REC/uniop/$S/uniop.pred" ] && ln -s "$REC/uniop/$S/uniop.pred" want/
      recdir="$REC/uniop/$S"
    fi ;;
esac
recording=0
if [ -n "${recdir:-}" ] && [ "${NF_REAL_RECORD:-0}" = 1 ] && [ ! -e "$recdir" ]; then recording=1
elif [ -n "${recdir:-}" ] && [ ! -f "$recdir/.giap-complete" ]; then verdict "failed (no complete recording at $recdir)"
else [ -n "$(ls -A "${want[0]}")" ] || verdict "failed (no recording at ${want[0]})"; fi

SECONDS=0
nextflow -log "$T/nextflow.log" run "$REPO/workflow/tests/real/$t.nf" -c "$REPO/workflow/tests/real/real.config" \
  -profile farmshare_local,apptainer -work-dir "$T/work" --sample "$S" "${args[@]}"
rc=$?
echo "nextflow exit=$rc wall=${SECONDS}s"
if [ "$rc" -ne 0 ]; then
  err=$(grep -m1 -A1 -E "ERROR ~|Caused by:" "$T/nextflow.log" | tr '\n' ' ' | cut -c1-200)
  verdict "failed (nextflow exit $rc: ${err:-see $T/nextflow.log})"
fi
if [ "$recording" = 1 ]; then   # save this run as the recording (inputs, outputs, provenance)
  rm -rf "$recdir.partial"; mkdir -p "$recdir.partial" && cp -r "$T/in/." "$T/out/." "$recdir.partial/" && rm -rf "$recdir.partial/parse_bakta"
  { echo "# inputs: record_tools_inputs.py ${t/promotech/promotech-windows} $S.json genome.norm.fa . (gippa $(git -C "$REPO" rev-parse --short HEAD))"
    echo "# nextflow run workflow/tests/real/$t.nf -profile farmshare_local,apptainer; per task (.command.sh):"
    for c in $(grep -l . "$T"/work/*/*/.command.sh); do sed -n '/^[^#]/p' "$c" | grep -m1 -E "promotech.py|UniOP"; done | sort -u
  } > "$recdir.partial/command.txt"
  { cat "$(ls "$T"/work/*/*/versions.yml | head -1)"; envl="$REPO/workflow/envs/locks"; (cd "$envl" && sha256sum "$t".explicit.txt "$t".pip.txt)
    [ "$t" = promotech ] && echo "RF-HOT.zip $(cat "$GIPPA_DB_DIR/promotech/RF-HOT.zip.sha256")"; } > "$recdir.partial/versions.txt"
  echo "wall_s=$SECONDS cpus=${SLURM_CPUS_PER_TASK:-?} node=$(hostname) job=${SLURM_ARRAY_JOB_ID:-?}_${SLURM_ARRAY_TASK_ID:-?}" > "$recdir.partial/runtime.txt"
  (cd "$recdir.partial" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
  touch "$recdir.partial/.giap-complete" && mv "$recdir.partial" "$recdir" || verdict "failed (could not save the recording $recdir)"
  verdict "recorded ($recdir, ${SECONDS}s)"
fi
v=$(python "$REPO/scripts/dev/nf_real_compare.py" "$P" "$T/out" "${want[0]}")
verdict "$v (${SECONDS}s)"
