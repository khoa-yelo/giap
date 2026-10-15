#!/bin/bash
#SBATCH --job-name=giap-nf-process-real
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
# Inputs: $GIAP_ROOT/recordings/bakta_full/<sample> (preflight genome + replicons, Bakta JSON/GBFF);
# genes.coords, ostir_input.csv and TU regions rebuilt by scripts/dev/record_tools_inputs.py as
# record_tools.sbatch did; PROMOTECH windows and the UNIOP GFF from tests/data/recorded (sim100k only).
# Usage (from the repo or worktree root; the code under test is that checkout):
#   sbatch scripts/slurm/nf_process_real.sbatch [sample]          # default sim100k; all 9 processes
#   sbatch --array=4 --export=ALL,NF_REAL_JOB=<array job> scripts/slurm/nf_process_real.sbatch
#       reruns one process (here TRANSTERMHP) into that job's directory, replacing its verdict
# Each task writes "PROCESS <name>: ok|mismatch|failed (<reason>)" to its log and to
# $GIAP_ROOT/scratch/nf-process-real/<array job>/verdicts/; summary.txt there lists all of them.
# Work dirs live under that directory (never node-local /tmp); $GIAP_DB_DIR is only read.
# R14: 64G is billed as 18 CPUs (MaxMemPerCPU), so %4 keeps the array at 72 CPUs; check squeue first.
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
REPO="${NF_REAL_REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"   # sbatch runs a spooled copy: use the submit dir
[ -d "$REPO/workflow/tests/real" ] || { echo "no workflow/tests/real in $REPO: submit from the repo root" >&2; exit 2; }
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH" NXF_ANSI_LOG=false
S="${1:-sim100k}"
PROCS=(BAKTA GENOMAD ANTISMASH TERMINATORNET TRANSTERMHP OSTIR PROMOTECH PROMOTER_CALCULATOR UNIOP)
P="${PROCS[${SLURM_ARRAY_TASK_ID:-0}]}"; t="${P,,}"
REC="$GIAP_ROOT/recordings"; PRE="$REC/bakta_full/$S/preflight/work/preflight"; B="$REC/bakta_full/$S/bakta"
FIX="$REPO/tests/data/recorded"
J="$GIAP_ROOT/scratch/nf-process-real/${NF_REAL_JOB:-${SLURM_ARRAY_JOB_ID:-manual}}"
T="$J/$t"; rm -rf "$T"; mkdir -p "$T/in" "$T/want" "$J/verdicts" "$T/tmp" && cd "$T" || exit 1
export TMPDIR="$T/tmp"
INPUTS="python $REPO/scripts/dev/record_tools_inputs.py"
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
  BAKTA) args=(--genome "$PRE/genome.norm.fa" --replicons "$PRE/replicons.tsv" --bakta_db "$GIAP_DB_DIR/bakta/db"); want=("$B") ;;
  GENOMAD) args=(--genome "$PRE/genome.norm.fa" --db_dir "$GIAP_DB_DIR"); want=("$REC/genomad/$S") ;;
  ANTISMASH) args=(--gbff "$REC/antismash/$S/bakta.gbff" --db_dir "$GIAP_DB_DIR"); want=("$REC/antismash/$S") ;;
  TERMINATORNET) args=(--genome "$PRE/genome.norm.fa"); ln -s "$REC/terminatornet/$S/terminatornet.tsv" want/ ;;
  TRANSTERMHP) $INPUTS coords "$B/$S.json" "$T/in/genes.coords" || verdict "failed (record_tools_inputs.py coords)"
    args=(--genome "$PRE/genome.norm.fa" --coords "$T/in/genes.coords"); ln -s "$REC/transtermhp/$S/transtermhp.tt" want/ ;;
  OSTIR) $INPUTS ostir "$B/$S.json" "$T/in/ostir_input.csv" || verdict "failed (record_tools_inputs.py ostir)"
    args=(--csv "$T/in/ostir_input.csv"); ln -s "$REC/ostir/$S/ostir.csv" want/ ;;
  PROMOTECH) [ "$S" = sim100k ] || verdict "failed (recorded windows exist for sim100k only)"
    d=$(ls -d "$FIX"/promotech/*/*/ | head -1)
    args=(--windows "$d/windows_000.fa" --db_dir "$GIAP_DB_DIR")
    mkdir want/promotech_000; ln -s "$d/sequences_predictions.csv" want/promotech_000/ ;;
  PROMOTER_CALCULATOR) $INPUTS regions "$B/$S.json" "$T/in/regions" || verdict "failed (record_tools_inputs.py regions)"
    args=(--regions "$T/in/regions"); want=("$REC/promoter_calculator/$S/csv") ;;
  UNIOP) [ "$S" = sim100k ] || verdict "failed (recorded UniOP input exists for sim100k only)"
    d=$(ls -d "$FIX"/uniop/*/sim100k)
    args=(--gff "$d/uniop_cds.gff"); ln -s "$d/uniop.pred" want/ ;;
esac
[ -e "${want[0]}" ] || verdict "failed (no recording at ${want[0]})"

SECONDS=0
nextflow -log "$T/nextflow.log" run "$REPO/workflow/tests/real/$t.nf" -c "$REPO/workflow/tests/real/real.config" \
  -profile farmshare_local,apptainer -work-dir "$T/work" --sample "$S" "${args[@]}"
rc=$?
echo "nextflow exit=$rc wall=${SECONDS}s"
if [ "$rc" -ne 0 ]; then
  err=$(grep -m1 -A1 -E "ERROR ~|Caused by:" "$T/nextflow.log" | tr '\n' ' ' | cut -c1-200)
  verdict "failed (nextflow exit $rc: ${err:-see $T/nextflow.log})"
fi
v=$(python "$REPO/scripts/dev/nf_real_compare.py" "$P" "$T/out" "${want[0]}")
verdict "$v (${SECONDS}s)"
