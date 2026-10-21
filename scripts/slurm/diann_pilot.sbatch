#!/bin/bash
#SBATCH --job-name=giap-diann
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# DIA-NN for the pilot arms A1 and A2 (GIAP_SPEC.md §15.2 step 5; M3-24 decisions 1-6).
#
#   sbatch -J giap-diann-A1 scripts/slurm/diann_pilot.sbatch A1
#   sbatch -J giap-diann-A2 scripts/slurm/diann_pilot.sbatch A2
#
# The search settings live in validation/pilot/diann_<arm>.cfg, not here: this script only resolves
# the FASTA paths, the raw runs and the output names. DIA-NN 2.3.1 Academia runs from the Apptainer
# image Khoa supplied (decision 1); the image is never copied into the repository (decision 3) and
# is taken from $DIANN_SIF, defaulting to the path below.
#
# Two calls, in order (decision 4): DIA-NN 2.x disables the FASTA digest when --fasta-search,
# --predictor and raw files are combined in one call, so step 1 predicts the in-silico library from
# the FASTA alone and step 2 searches the runs against it. Step 1 is skipped when its library is
# already there, so a re-queued job resumes.
#
# Environment: DIANN_SIF (image), RAW_DIR (the *.d runs of PXD078459, diaPASEF, §15.1),
# OUT_DIR (parent of the per-arm output directory). Everything is written under $GIAP_ROOT.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
module load apptainer/1.4.4 >/dev/null 2>&1 || true

usage() {
  echo "usage: sbatch [-J giap-diann-A1] scripts/slurm/diann_pilot.sbatch A1|A2" >&2
  exit 2
}
[ $# -eq 1 ] || usage
ARM="$1"
case "$ARM" in
  A1 | A2) ;;
  *) echo "unknown arm '$ARM': the pilot arms searched with DIA-NN are A1 and A2 (§15.1)" >&2; usage ;;
esac

SIF="${DIANN_SIF:-/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif}"
DIANN="/diann-2.3.1/diann-linux"   # the binary inside the image (decision 1)
if [ ! -f "$SIF" ]; then
  echo "DIA-NN image not found: $SIF" >&2
  echo "set DIANN_SIF to the 2.3.1 Academia image (validation/pilot/README.md); it is never copied into the repository" >&2
  exit 2
fi
command -v apptainer >/dev/null 2>&1 || { echo "apptainer is not on PATH (module load apptainer/1.4.4)" >&2; exit 2; }

CFG="${GIAP_REPO}/validation/pilot/diann_${ARM}.cfg"
[ -f "$CFG" ] || { echo "missing settings file: $CFG" >&2; exit 2; }
RAW_DIR="${RAW_DIR:-$GIAP_DATA_DIR/pxd078459}"
OUT="${OUT_DIR:-$GIAP_ROOT/pilot}/$ARM"
FASTA_DIR="$OUT/fasta"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
PREDICTED_LIB="$OUT/${ARM}.predicted.speclib"   # DIA-NN names it after --out-lib
mkdir -p "$OUT" "$FASTA_DIR"

# The timsTOF runs of PXD078459 (§15.1); checked before any DIA-NN call.
RAW_FILES=()
for d in "$RAW_DIR"/*.d; do
  [ -e "$d" ] || continue
  RAW_FILES+=(--f "$d")
done
if [ "${#RAW_FILES[@]}" -eq 0 ]; then
  echo "no *.d runs in $RAW_DIR (run scripts/slurm/fetch_pxd078459.sbatch first)" >&2
  exit 2
fi

expand() {  # $VAR / ${VAR} of the GIAP environment in a cfg value
  local s="$1" name val
  for name in GIAP_REPO GIAP_ROOT GIAP_DATA_DIR GIAP_DB_DIR GIAP_TOOLS_DIR; do
    val="${!name:-}"
    s="${s//\$\{$name\}/$val}"
    s="${s//\$$name/$val}"
  done
  printf '%s' "$s"
}

CFG_ARGS=()
load_step() {  # flags of [step1]/[step2] of $CFG, in file order, into CFG_ARGS
  local want="$1" section="" line tok
  CFG_ARGS=()
  set -f  # K*,R* and UniMod:35,… must not be globbed
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [ -n "$line" ] || continue
    case "$line" in
      "["*"]")
        section="${line//[![:alnum:]]/}"
        section="${section,,}"
        continue
        ;;
    esac
    [ "$section" = "$want" ] || continue
    for tok in $line; do CFG_ARGS+=("$tok"); done
  done < "$CFG"
  set +f
  [ "${#CFG_ARGS[@]}" -gt 0 ] || { echo "$CFG: section [$want] is empty" >&2; exit 2; }
}

STEP_ARGS=()
resolve_fastas() {  # expand every --fasta of CFG_ARGS; unpack a gzipped one outside the repository
  local i=0 tok val base
  STEP_ARGS=()
  while [ "$i" -lt "${#CFG_ARGS[@]}" ]; do
    tok="${CFG_ARGS[$i]}"
    if [ "$tok" != "--fasta" ]; then
      STEP_ARGS+=("$tok")
      i=$((i + 1))
      continue
    fi
    val="$(expand "${CFG_ARGS[$((i + 1))]}")"
    case "$val" in
      /*) ;;
      *) echo "$CFG: --fasta $val is not an absolute path" >&2; exit 2 ;;
    esac
    if [ "${val%.gz}" != "$val" ]; then
      [ -f "$val" ] || { echo "missing FASTA: $val" >&2; exit 2; }
      base="$(basename "${val%.gz}")"
      [ -s "$FASTA_DIR/$base" ] || gunzip -c "$val" >"$FASTA_DIR/$base"
      val="$FASTA_DIR/$base"
    fi
    [ -f "$val" ] || { echo "missing FASTA: $val (giap setup / fetch_resources.sbatch)" >&2; exit 2; }
    STEP_ARGS+=(--fasta "$val")
    i=$((i + 2))
  done
}

diann() { apptainer exec -B /home/classes -B /scratch -B /farmshare "$SIF" "$DIANN" "$@"; }

echo "== arm $ARM; image $SIF; settings $CFG; output $OUT"
diann --version >"$OUT/diann_version.txt" 2>&1 || echo "warning: the image did not answer --version" >&2
head -n 2 "$OUT/diann_version.txt" || true

# ---- step 1: predict the in-silico spectral library from the FASTA files (no raw files)
if [ -f "$PREDICTED_LIB" ]; then
  echo "step 1: the predicted library is already there, skipping: $PREDICTED_LIB"
else
  echo "step 1: predicting the spectral library"
  load_step step1
  resolve_fastas
  diann "${STEP_ARGS[@]}" --threads "$THREADS" --out-lib "$OUT/${ARM}.parquet"
  [ -f "$PREDICTED_LIB" ] || { echo "step 1 wrote no $PREDICTED_LIB" >&2; exit 1; }
fi

# ---- step 2: search the runs against the predicted library
echo "step 2: searching ${#RAW_FILES[@]} run(s) from $RAW_DIR"
load_step step2
resolve_fastas
diann "${STEP_ARGS[@]}" "${RAW_FILES[@]}" --threads "$THREADS" --lib "$PREDICTED_LIB" \
  --out "$OUT/${ARM}_report.parquet" --out-lib "$OUT/${ARM}_empirical.parquet"

echo "== $ARM done: $OUT/${ARM}_report.parquet"
[ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null || true
