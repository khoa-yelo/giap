#!/bin/bash
#SBATCH --job-name=giap-denovo
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000M
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Pilot arm A3, de novo sequencing (GIAP_SPEC.md §15.1, §15.2 step 6; M3-24).
#
#   sbatch scripts/slurm/denovo_pilot.sbatch              # Cascadia 0.0.7 on WT_1, WT_2, WT_3
#   A3_RUN=WT_1 sbatch --time=04:00:00 scripts/slurm/denovo_pilot.sbatch cascadia   # one run
#   sbatch scripts/slurm/denovo_pilot.sbatch cascadia WT_2                      # the same
#   sbatch scripts/slurm/denovo_pilot.sbatch casanovo     # the documented fallback
#
# Main branch (M10-08, deliverable 8): Cascadia 0.0.7 (DIA de novo) on one GPU (1 x L40S on the
# `gpu` partition), one call per run, on $GIAP_DATA_DIR/pxd078459/mzml/<RUN>.mzML (from
# mzml_pilot.sbatch). The weights and the Unimod table must match the sha256 committed in
# validation/pilot/a3_assets.tsv (exit 2 otherwise, R6). Cascadia runs through
# scripts/pilot/cascadia_offline.py (no Unimod download) in a per-job mktemp -d
# ${TMPDIR:-/tmp}/giap-a3.XXXXXX/<RUN> (its scratch goes to the cwd; removed on exit), with its
# default score threshold 0.8, and writes $GIAP_ROOT/pilot/A3_cascadia/<RUN>/<RUN>.ssl. Exit 1
# when CUDA is not available, when its log shows a unimod.org URL (the launcher did not take
# effect), when Cascadia fails, or when it writes no <RUN>.ssl. 64000M = 16 x MaxMemPerCPU 4000M.
#
# Fallback branch (§15.1): if DIA de novo fails on timsTOF data, run Casanovo 5.2.1 on three UTI89
# DDA runs of PXD053727 (genome GCA_000013265.1) instead, and say so in validation/pilot/REPORT.md.
# UTI89 is a different strain, so the A3 comparison then measures de novo against A2's peptides
# rather than against the same spectra; §15.5 metric 6 is reported with that caveat.
#
# Environment: DENOVO_ENV (micromamba env name), DENOVO_MODEL (model weights), A3_RUN (Cascadia
# runs), MZML_DIR (Cascadia input), RAW_DIR (Casanovo input runs), OUT_DIR (parent of the output directory); outputs go under $GIAP_ROOT.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

TOOL="${1:-cascadia}"
case "$TOOL" in
  cascadia)
    ENV_NAME="${DENOVO_ENV:-giap-cascadia}"
    DEFAULT_RAW="$GIAP_DATA_DIR/pxd078459"
    DEFAULT_MODEL="$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt"
    PATTERN="*.d"
    ;;
  casanovo)
    ENV_NAME="${DENOVO_ENV:-giap-casanovo}"
    DEFAULT_RAW="$GIAP_DATA_DIR/pxd053727"
    DEFAULT_MODEL="$GIAP_DB_DIR/casanovo/casanovo_massivekb_5.2.1.ckpt"
    PATTERN="*.mzML"
    ;;
  *)
    echo "usage: sbatch scripts/slurm/denovo_pilot.sbatch [cascadia|casanovo]" >&2
    exit 2
    ;;
esac

RAW_DIR="${RAW_DIR:-$DEFAULT_RAW}"
MODEL="${DENOVO_MODEL:-$DEFAULT_MODEL}"
OUT="${OUT_DIR:-$GIAP_ROOT/pilot}/A3_$TOOL"

if [ ! -d "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME" ]; then
  echo "the $TOOL environment is missing: $MAMBA_ROOT_PREFIX/envs/$ENV_NAME" >&2
  echo "add workflow/envs/${TOOL}.yml and run: sbatch scripts/slurm/build_envs.sbatch $TOOL" >&2
  exit 2
fi
if [ ! -f "$MODEL" ]; then
  echo "the $TOOL model weights are missing: $MODEL (set DENOVO_MODEL)" >&2
  exit 2
fi

if [ "$TOOL" = "cascadia" ]; then
  ASSETS="$GIAP_REPO/validation/pilot/a3_assets.tsv"
  UNIMOD="$GIAP_DB_DIR/cascadia/unimod_tables.xml"
  MZML_DIR="${MZML_DIR:-$GIAP_DATA_DIR/pxd078459/mzml}"
  field() {  # field ASSET COLUMN: one cell of a3_assets.tsv
    awk -F'\t' -v a="$1" -v c="$2" \
      'NR == 1 { for (i = 1; i <= NF; i++) if ($i == c) k = i; next } $1 == a { print $k; exit }' \
      "$ASSETS"
  }
  sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum -- "$1"; else shasum -a 256 -- "$1"; fi \
      | cut -d' ' -f1
  }
  verify() {  # verify ASSET FILE: the committed sha256, never an unverified file (R6)
    [ -f "$2" ] || { echo "$1 is missing: $2 (fetch_a3_assets.sbatch)" >&2; exit 2; }
    want="$(field "$1" sha256)"
    case "$want" in
      *[!0-9a-f]* | "") want="" ;;
    esac
    if [ "${#want}" -ne 64 ]; then
      echo "$1: no sha256 in $ASSETS: not running on an unverified file" >&2
      exit 2
    fi
    got="$(sha256_of "$2")"
    [ "$got" = "$want" ] || { echo "$1: sha256 $got != $want ($ASSETS)" >&2; exit 2; }
    echo "$1 sha256 $got ok"
  }
  verify cascadia_weights "$MODEL"
  verify unimod_table "$UNIMOD"
  RUN_LIST="${A3_RUN:-${2:-WT_1 WT_2 WT_3}}"
  for r in $RUN_LIST; do
    [ -f "$MZML_DIR/$r.mzML" ] || { echo "no mzML for $r: $MZML_DIR/$r.mzML" >&2; exit 2; }
  done
  echo "== A3 de novo with Cascadia 0.0.7 on: $RUN_LIST; score threshold 0.8 (CLI default)"
  micromamba run -n "$ENV_NAME" cascadia --help 2>&1 | head -n 40 || true
  nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true
  cuda="$(micromamba run -n "$ENV_NAME" python -c 'import torch; print(torch.cuda.is_available())')" \
    || cuda="error"
  echo "torch.cuda.is_available(): $cuda"
  [ "$cuda" = "True" ] || { echo "CUDA is not available: Cascadia would run on the CPU" >&2; exit 1; }
  SCRATCH_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/giap-a3.XXXXXX")"
  trap 'rm -rf "$SCRATCH_ROOT"' EXIT
  for r in $RUN_LIST; do
    mkdir -p "$OUT/$r"
    SCRATCH="$SCRATCH_ROOT/$r"
    mkdir "$SCRATCH"
    log="$OUT/$r/cascadia.log"
    rc=0
    (cd "$SCRATCH" && micromamba run -n "$ENV_NAME" python -I \
      "$GIAP_REPO/scripts/pilot/cascadia_offline.py" "$UNIMOD" sequence "$MZML_DIR/$r.mzML" \
      "$MODEL" -o "$OUT/$r/$r" -b 32) 2>&1 | tee "$log" || rc=$?
    if grep -q 'unimod\.org' "$log"; then
      echo "the job log shows a unimod.org URL: the offline launcher did not take effect" >&2
      exit 1
    fi
    [ "$rc" -eq 0 ] || { echo "Cascadia failed on $r (exit $rc): R9" >&2; exit 1; }
    [ -f "$OUT/$r/$r.ssl" ] || { echo "Cascadia exited 0 but wrote no $OUT/$r/$r.ssl: R9" >&2; exit 1; }
    ls -l "$OUT/$r/$r.ssl"
  done
  echo "== A3 done: $OUT (next: a3_feasibility.py evaluate, then giap convert --from cascadia)"
  [ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null || true
  exit 0
fi

RUNS=()
RUN_N=0   # a counter, so bash 3.2 (macOS CI) never expands an empty array under `set -u`
for f in "$RAW_DIR"/$PATTERN; do
  [ -e "$f" ] || continue
  RUNS+=("$f")
  RUN_N=$((RUN_N + 1))
done
if [ "$RUN_N" -eq 0 ]; then
  echo "no $PATTERN runs in $RAW_DIR" >&2
  exit 2
fi

mkdir -p "$OUT"
echo "== A3 de novo with $TOOL on $RUN_N run(s) from $RAW_DIR; output $OUT"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

micromamba run -n "$ENV_NAME" casanovo sequence \
  --model "$MODEL" --output "$OUT/a3_casanovo.mztab" "${RUNS[@]}"

echo "== A3 done: $OUT (next: giap convert --from cascadia, §15.2 step 7)"
[ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null || true
