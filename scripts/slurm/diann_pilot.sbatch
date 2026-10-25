#!/bin/bash
#SBATCH --job-name=giap-diann
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# DIA-NN for the pilot arms A1 and A2 (GIAP_SPEC.md §15.2 step 5; M3-24 decisions 1-6; M10-02).
#
#   sbatch -J giap-diann-A1 scripts/slurm/diann_pilot.sbatch A1
#   sbatch -J giap-diann-A2 scripts/slurm/diann_pilot.sbatch A2
#
# The search settings live in validation/pilot/diann_<arm>.cfg, not here: this script only resolves
# the FASTA paths, the raw runs and the output names. DIA-NN 2.3.1 Academia runs only from the
# verified read-only copy that scripts/slurm/diann_image.sbatch installs (M10-02, index lead
# decision 2): $DIANN_SIF, defaulting to $GIAP_ROOT/images/diann_2.3.1.sif. Its sha256 (pin:
# validation/pilot/diann_image.sha256), its mode and the acquisition mode are checked before any
# container call.
#
# Two calls, in order (decision 4): DIA-NN 2.x disables the FASTA digest when --fasta-search,
# --predictor and raw files are combined in one call, so step 1 predicts the in-silico library from
# the FASTA alone and step 2 searches the runs against it. Step 1 is skipped only when step1.done
# records the sha256 of the current library and the sha256s of the current cfg, FASTA files and
# image, so a re-queued job resumes; step 2 always starts fresh.
#
# --mem=64000M is 16 x MaxMemPerCPU (4000 MB on partition normal): any more and Slurm raises the
# CPU count (--mem=64G gave 17 CPUs per task, 18 allocated; PR #293 review S1).
#
# Environment: DIANN_SIF (image), RAW_DIR (the *.d runs of PXD078459, §15.1), ACQUISITION_JSON
# (M10-01's run list and mode), OUT_DIR (parent of the per-arm output directory), GIAP_PYTHON.
# Everything is written under $GIAP_ROOT; provenance.json is written on every exit.
# Exit 0 report written; 1 DIA-NN failed or wrote no report; 2 a precondition; 3 the acquisition
# mode forbids these settings (§17, ADR-0051); 143/130/129 killed by TERM/INT/HUP (scancel,
# TIMEOUT), recorded as status failed with FAILED.
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

SIF="${DIANN_SIF:-$GIAP_ROOT/images/diann_2.3.1.sif}"
DIANN="/diann-2.3.1/diann-linux"   # the binary inside the image (decision 1)
PIN_FILE="$GIAP_REPO/validation/pilot/diann_image.sha256"
CFG="${GIAP_REPO}/validation/pilot/diann_${ARM}.cfg"
ACQ_JSON="${ACQUISITION_JSON:-$GIAP_REPO/validation/pilot/inputs/acquisition.json}"
PY="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
REPORT_PY="$GIAP_REPO/validation/pilot/diann_report.py"
RAW_DIR="${RAW_DIR:-$GIAP_DATA_DIR/pxd078459}"
OUT="${OUT_DIR:-$GIAP_ROOT/pilot}/$ARM"
FASTA_DIR="$OUT/fasta"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
PREDICTED_LIB="$OUT/${ARM}.predicted.speclib"   # DIA-NN names it after --out-lib
REPORT="$OUT/${ARM}_report.parquet"
MARKER="$OUT/step1.done"
NL='
'
mkdir -p "$OUT" "$FASTA_DIR"
rm -f "${OUT:?}/diann_version.txt"

# ---- provenance on every exit after the argument check (spec rule 7)
START_EPOCH="$(date +%s)"
# STEP<N>_RC stays empty until DIA-NN has returned: a job killed during a call never records 0
STEP1=not_reached STEP1_RC="" STEP2=not_reached STEP2_RC=""
INPUT_FASTAS="" APPTAINER_VERSION="" UNPACK_TMP=""
write_failed() {  # step, exit code[, reason]
  printf 'step=%s exit_code=%s%s utc=%s\n' "$1" "$2" "${3:+ reason=$3}" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >"$OUT/FAILED"
}
on_exit() {
  local rc=$? status=failed
  trap - TERM INT HUP
  [ -z "$UNPACK_TMP" ] || rm -f "$UNPACK_TMP" "$UNPACK_TMP.sha"
  case "$rc" in 0) status=ok ;; 2) status=precondition ;; 3) status=acquisition_mode ;; esac
  # a step still running when the job ended (a signal): it takes the job's exit code, and FAILED
  if [ "$rc" -ne 0 ]; then
    if [ "$STEP2" = ran ] && [ -z "$STEP2_RC" ]; then
      STEP2_RC=$rc
      write_failed step2 "$rc" interrupted
    elif [ "$STEP1" = ran ] && [ -z "$STEP1_RC" ]; then
      STEP1_RC=$rc
      write_failed step1 "$rc" interrupted
    fi
  fi
  "$PY" -I "$REPORT_PY" provenance --arm "$ARM" --out "$OUT" --status "$status" \
    --exit-code "$rc" --step1 "$STEP1" --step1-rc "$STEP1_RC" --step2 "$STEP2" \
    --step2-rc "$STEP2_RC" --image "$SIF" --pin-file "$PIN_FILE" --cfg "$CFG" \
    --repo "$GIAP_REPO" --acquisition "$ACQ_JSON" --raw-dir "$RAW_DIR" \
    --fastas "$INPUT_FASTAS" --apptainer "$APPTAINER_VERSION" --node "$(hostname)" \
    --start "$START_EPOCH" || echo "warning: provenance.json was not written" >&2
  exit "$rc"
}
trap on_exit EXIT
# Slurm's scancel and TIMEOUT send TERM: without these traps the EXIT trap would see $? = 0
trap 'exit 143' TERM
trap 'exit 130' INT
trap 'exit 129' HUP
fail() {  # step, exit code[, reason]: record FAILED and exit 1; FAILED never says exit code 0
  local code="$2"
  [ "$code" -ne 0 ] || code=1
  write_failed "$1" "$code" "${3:-}"
  exit 1
}

# ---- rule 1: the verified read-only copy only, checked before any container call
if [ ! -f "$SIF" ]; then
  echo "DIA-NN image copy not found: $SIF (DIANN_SIF); submit scripts/slurm/diann_image.sbatch first" >&2
  exit 2
fi
if [ ! -f "$PIN_FILE" ] || [ "$(grep -c '' "$PIN_FILE")" -ne 1 ] \
  || ! grep -Eqx '[0-9a-f]{64}  diann_2\.3\.1\.sif' "$PIN_FILE"; then
  echo "missing or malformed pin file $PIN_FILE: want one line '<sha256>  diann_2.3.1.sif'" >&2
  exit 2
fi
PIN="$(cut -d' ' -f1 "$PIN_FILE")"
SIF_SHA="$(sha256sum "$SIF" | cut -d' ' -f1)"
if [ "$SIF_SHA" != "$PIN" ]; then
  echo "$SIF hashes to $SIF_SHA, not the pin $PIN ($PIN_FILE); inspect and delete the copy, then re-run scripts/slurm/diann_image.sbatch" >&2
  exit 2
fi
SIF_MODE="$(stat -c %a "$SIF")"
if [ $((0$SIF_MODE & 0222)) -ne 0 ]; then
  echo "$SIF has mode $SIF_MODE: the copy must be read-only, mode 0444" >&2
  exit 2
fi
[ -f "$CFG" ] || { echo "missing settings file: $CFG" >&2; exit 2; }

# ---- rule 2: the acquisition mode (exit 3); rule 3: the runs, in the file's order
if [ ! -f "$PY" ] || [ ! -x "$PY" ]; then
  echo "GIAP_PYTHON=$PY is not an executable file (default \$UV_PROJECT_ENVIRONMENT/bin/python)" >&2
  exit 2
fi
RUN_LIST="" rc=0
RUN_LIST="$("$PY" -I "$REPORT_PY" runs "$ACQ_JSON")" || rc=$?
case "$rc" in
  0) ;;
  3) exit 3 ;;
  2) echo "cannot read the runs of $ACQ_JSON (GIAP_PYTHON=$PY)" >&2; exit 2 ;;
  *) echo "GIAP_PYTHON=$PY failed with exit $rc reading the runs of $ACQ_JSON" >&2; exit 2 ;;
esac
RAW_FILES=()
RAW_N=0
while IFS= read -r run; do
  [ -n "$run" ] || continue
  for f in analysis.tdf analysis.tdf_bin; do
    [ -f "$RAW_DIR/$run.d/$f" ] || { echo "run $run: no $RAW_DIR/$run.d/$f" >&2; exit 2; }
  done
  RAW_FILES+=(--f "$RAW_DIR/$run.d")
  RAW_N=$((RAW_N + 1))
done <<EOF
$RUN_LIST
EOF
if [ "$RAW_N" -eq 0 ]; then
  echo "no runs listed in $ACQ_JSON for $RAW_DIR (fetch_pxd078459.sbatch and M10-01 first)" >&2
  exit 2
fi
for d in "$RAW_DIR"/*.d; do
  [ -e "$d" ] || continue
  case "$NL$RUN_LIST$NL" in
    *"$NL$(basename "$d" .d)$NL"*) ;;
    *) echo "warning: $d is not in $ACQ_JSON; ignored" >&2 ;;
  esac
done
command -v apptainer >/dev/null 2>&1 || { echo "apptainer is not on PATH (module load apptainer/1.4.4)" >&2; exit 2; }

# bash 3.2 (macOS CI) runs these functions too: no ${v,,}, no ${!v}, no empty-array expansion.
expand() {  # $VAR / ${VAR} of the GIAP environment in a cfg value
  local s="$1" name val
  for name in GIAP_REPO GIAP_ROOT GIAP_DATA_DIR GIAP_DB_DIR GIAP_TOOLS_DIR; do
    eval "val=\${$name:-}"
    s="${s//\$\{$name\}/$val}"
    s="${s//\$$name/$val}"
  done
  printf '%s' "$s"
}

CFG_ARGS=()
CFG_N=0
load_step() {  # flags of [step1]/[step2] of $CFG, in file order, into CFG_ARGS
  local want="$1" section="" line tok
  CFG_ARGS=()
  CFG_N=0
  set -f  # K*,R* and UniMod:35,… must not be globbed
  while IFS= read -r line || [ -n "$line" ]; do
    line="${line%%#*}"
    case "$line" in
      *\[*\]*)
        section=$(printf '%s' "$line" | tr -cd '[:alnum:]' | tr '[:upper:]' '[:lower:]')
        continue
        ;;
    esac
    [ "$section" = "$want" ] || continue
    for tok in $line; do
      CFG_ARGS+=("$tok")
      CFG_N=$((CFG_N + 1))
    done
  done < "$CFG"
  set +f
  [ "$CFG_N" -gt 0 ] || { echo "$CFG: section [$want] is empty" >&2; exit 2; }
}

sha() { sha256sum "$1" | cut -d' ' -f1; }
unpack() {  # .gz, target: reuse the target only if <target>.src.sha256 is the .gz's sha256
  local gz="$1" dst="$2" gz_sha
  gz_sha="$(sha "$gz")"
  if [ -f "$dst" ] && [ -f "$dst.src.sha256" ] \
    && [ "$(cut -d' ' -f1 "$dst.src.sha256")" = "$gz_sha" ]; then
    return 0
  fi
  echo "unpacking $gz to $dst"
  rm -f "$dst" "$dst.src.sha256"
  UNPACK_TMP="$(dirname "$dst")/.$(basename "$dst").tmp.$$"
  if ! gunzip -c "$gz" >"$UNPACK_TMP"; then
    rm -f "$UNPACK_TMP"
    echo "cannot unpack $gz" >&2
    exit 2
  fi
  printf '%s  %s\n' "$gz_sha" "$(basename "$gz")" >"$UNPACK_TMP.sha"
  mv "$UNPACK_TMP" "$dst"
  mv "$UNPACK_TMP.sha" "$dst.src.sha256"
  UNPACK_TMP=""
}

STEP_ARGS=()
resolve_fastas() {  # expand every --fasta of CFG_ARGS; unpack a gzipped one outside the repository
  local i=0 tok val base
  STEP_ARGS=()
  while [ "$i" -lt "$CFG_N" ]; do
    tok="${CFG_ARGS[$i]}"
    if [ "$tok" != "--fasta" ]; then
      STEP_ARGS+=("$tok")
      i=$((i + 1))
      continue
    fi
    if [ $((i + 1)) -ge "$CFG_N" ]; then
      echo "$CFG: a --fasta at the end of a section has no value" >&2
      usage
    fi
    val="$(expand "${CFG_ARGS[$((i + 1))]}")"
    case "$val" in
      /*) ;;
      *) echo "$CFG: --fasta $val is not an absolute path" >&2; exit 2 ;;
    esac
    if [ "${val%.gz}" != "$val" ]; then
      [ -f "$val" ] || { echo "missing FASTA: $val" >&2; exit 2; }
      base="$(basename "${val%.gz}")"
      unpack "$val" "$FASTA_DIR/$base"
      INPUT_FASTAS="$INPUT_FASTAS$val$NL"
      val="$FASTA_DIR/$base"
    fi
    [ -f "$val" ] || { echo "missing FASTA: $val (giap setup / fetch_resources.sbatch)" >&2; exit 2; }
    STEP_ARGS+=(--fasta "$val")
    INPUT_FASTAS="$INPUT_FASTAS$val$NL"
    i=$((i + 2))
  done
}

diann() { apptainer exec -B /home/classes -B /scratch -B /farmshare "$SIF" "$DIANN" "$@"; }
run_step() {  # N, then the DIA-NN arguments: record the argv one token per line, then call it
  local n="$1"
  shift
  printf '%s\n' apptainer exec -B /home/classes -B /scratch -B /farmshare "$SIF" "$DIANN" "$@" \
    >"$OUT/command_step$n.txt"
  diann "$@"
}
step1_inputs() {  # the inputs a step-1 library depends on, one 'key sha256 [path]' per line
  local f
  printf 'cfg_sha256 %s\n' "$(sha "$CFG")"
  printf 'image_sha256 %s\n' "$SIF_SHA"
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    printf 'fasta_sha256 %s %s\n' "$(sha "$f")" "$f"
  done <<EOF
$INPUT_FASTAS
EOF
}

echo "== arm $ARM; image $SIF; settings $CFG; output $OUT"
APPTAINER_VERSION="$(apptainer --version 2>/dev/null || true)"
diann --version >"$OUT/diann_version.txt" 2>&1 || true
grep -m 2 . "$OUT/diann_version.txt" || true
case "$(grep -m 1 . "$OUT/diann_version.txt" || true)" in  # first non-empty line
  *"DIA-NN 2.3.1 Academia"*) ;;
  *) echo "the image's banner is not 'DIA-NN 2.3.1 Academia' (ADR-0051)" >&2; exit 2 ;;
esac
rm -f "${OUT:?}/FAILED" "${OUT:?}/outputs.sha256"

# ---- step 1: predict the in-silico spectral library from the FASTA files (no raw files)
# step1.done: '<sha256>  <arm>.predicted.speclib', job_id=, commit=, then the step1_inputs lines
load_step step1
resolve_fastas
STEP1_INPUTS="$(step1_inputs)"
if [ -f "$MARKER" ] && [ -f "$PREDICTED_LIB" ] \
  && [ "$(head -n 1 "$MARKER" | cut -d' ' -f1)" = "$(sha "$PREDICTED_LIB")" ] \
  && [ "$(sed -e 1d -e '/^job_id=/d' -e '/^commit=/d' "$MARKER")" = "$STEP1_INPUTS" ]; then
  echo "step 1: step1.done matches the predicted library and its inputs, skipping: $PREDICTED_LIB"
  sed -n -e '2,3s/^/  built by /p' "$MARKER"
  STEP1=skipped
else
  echo "step 1: predicting the spectral library"
  rm -f "${MARKER:?}" "${OUT:?}/${ARM:?}.parquet" "${PREDICTED_LIB:?}" "${OUT:?}/${ARM:?}.log.txt"
  STEP1=ran STEP1_RC=""
  rc=0
  run_step 1 "${STEP_ARGS[@]}" --threads "$THREADS" --out-lib "$OUT/${ARM}.parquet" || rc=$?
  STEP1_RC=$rc
  [ "$STEP1_RC" -eq 0 ] || fail step1 "$STEP1_RC"
  [ -s "$PREDICTED_LIB" ] || { echo "step 1 wrote no $PREDICTED_LIB" >&2; fail step1 1 no_output; }
  COMMIT="$(git -C "$GIAP_REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
  {
    printf '%s  %s\n' "$(sha "$PREDICTED_LIB")" "$(basename "$PREDICTED_LIB")"
    printf 'job_id=%s\ncommit=%s\n' "${SLURM_JOB_ID:-}" "$COMMIT"
    printf '%s\n' "$STEP1_INPUTS"
  } >"$MARKER.tmp"
  mv "$MARKER.tmp" "$MARKER"
fi

# ---- step 2: search the runs against the predicted library, always from scratch
echo "step 2: searching $RAW_N run(s) from $RAW_DIR"
rm -f "${OUT:?}/${ARM:?}_report"* "${OUT:?}/${ARM:?}_empirical"*
load_step step2
resolve_fastas
STEP2=ran STEP2_RC=""
# --temp: DIA-NN writes <run>.d.quant next to the raw runs by default, so two arms searching the
# same runs at once read each other's (job 1879914: "library size mismatch"); each arm keeps its own
QUANT="$OUT/quant"
rm -rf "${QUANT:?}"
mkdir -p "$QUANT"
rc=0
run_step 2 "${STEP_ARGS[@]}" "${RAW_FILES[@]}" --threads "$THREADS" --lib "$PREDICTED_LIB" \
  --temp "$QUANT" --out "$REPORT" --out-lib "$OUT/${ARM}_empirical.parquet" || rc=$?
STEP2_RC=$rc
[ "$STEP2_RC" -eq 0 ] || fail step2 "$STEP2_RC"
[ -s "$REPORT" ] || { echo "step 2 wrote no report $REPORT" >&2; fail step2 1 no_output; }

# ---- the files DIA-NN wrote, sorted by name (fasta/ and this job's own records excluded)
(cd "$OUT" && find . -path ./fasta -prune -o -type f ! -name provenance.json \
  ! -name 'outputs.sha256*' ! -name FAILED ! -name step1.done ! -name diann_version.txt \
  ! -name 'command_step*.txt' ! -name check.json ! -name '.*' -print \
  | sed 's|^\./||' | LC_ALL=C sort | while IFS= read -r f; do sha256sum "$f"; done) \
  >"$OUT/outputs.sha256.tmp"
mv "$OUT/outputs.sha256.tmp" "$OUT/outputs.sha256"

echo "== $ARM done: $REPORT"
[ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null || true
