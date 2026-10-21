#!/bin/bash
#SBATCH --job-name=giap-denovo
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Pilot arm A3, de novo sequencing (GIAP_SPEC.md §15.1, §15.2 step 6; M3-24).
#
#   sbatch scripts/slurm/denovo_pilot.sbatch              # Cascadia 0.0.7 on the PXD078459 DIA runs
#   sbatch scripts/slurm/denovo_pilot.sbatch cascadia     # the same, explicitly
#   sbatch scripts/slurm/denovo_pilot.sbatch casanovo     # the documented fallback
#
# Main branch: Cascadia 0.0.7 (DIA de novo) on one GPU (1 x L40S on the `gpu` partition) over the
# same three wild-type runs of PXD078459 that arms A1 and A2 search.
#
# Fallback branch (§15.1): if DIA de novo fails on timsTOF data, run Casanovo 5.2.1 on three UTI89
# DDA runs of PXD053727 (genome GCA_000013265.1) instead, and say so in validation/pilot/REPORT.md.
# UTI89 is a different strain, so the A3 comparison then measures de novo against A2's peptides
# rather than against the same spectra; §15.5 metric 6 is reported with that caveat.
#
# Prepared only: M3 writes this script, M10 runs it. Neither tool is installed yet, so the script
# stops with exit 2 and a clear message until the environment and the model weights are there
# (add workflow/envs/cascadia.yml, then `sbatch scripts/slurm/build_envs.sbatch cascadia`).
#
# Environment: DENOVO_ENV (micromamba env name), DENOVO_MODEL (model weights),
# RAW_DIR (input runs), OUT_DIR (parent of the output directory); outputs go under $GIAP_ROOT.
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
THREADS="${SLURM_CPUS_PER_TASK:-8}"

if [ ! -d "$MAMBA_ROOT_PREFIX/envs/$ENV_NAME" ]; then
  echo "the $TOOL environment is missing: $MAMBA_ROOT_PREFIX/envs/$ENV_NAME" >&2
  echo "add workflow/envs/${TOOL}.yml and run: sbatch scripts/slurm/build_envs.sbatch $TOOL" >&2
  exit 2
fi
if [ ! -f "$MODEL" ]; then
  echo "the $TOOL model weights are missing: $MODEL (set DENOVO_MODEL)" >&2
  exit 2
fi

RUNS=()
for f in "$RAW_DIR"/$PATTERN; do
  [ -e "$f" ] || continue
  RUNS+=("$f")
done
if [ "${#RUNS[@]}" -eq 0 ]; then
  echo "no $PATTERN runs in $RAW_DIR" >&2
  exit 2
fi

mkdir -p "$OUT"
echo "== A3 de novo with $TOOL on ${#RUNS[@]} run(s) from $RAW_DIR; output $OUT"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader || true

case "$TOOL" in
  cascadia)
    micromamba run -n "$ENV_NAME" cascadia sequence \
      --model "$MODEL" --num_workers "$THREADS" --output "$OUT/a3_cascadia" "${RUNS[@]}"
    ;;
  casanovo)
    micromamba run -n "$ENV_NAME" casanovo sequence \
      --model "$MODEL" --output "$OUT/a3_casanovo.mztab" "${RUNS[@]}"
    ;;
esac

echo "== A3 done: $OUT (next: giap convert --from cascadia, §15.2 step 7)"
[ -z "${SLURM_JOB_ID:-}" ] || seff "$SLURM_JOB_ID" 2>/dev/null || true
