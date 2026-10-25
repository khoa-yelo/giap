#!/bin/bash
#SBATCH --job-name=giap-cascadia-cpu
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000M
#SBATCH --time=03:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Cascadia 0.0.7 on CPU as prep -> Slurm array -> merge (lead, 2026-10-24; Khoa: "split the job
# such that the task can be finished within 2 to 3 hours max"). Submit with
# scripts/pilot/cascadia_array_submit.sh RUN N_SHARDS. STEP=prep|predict|merge (set by the
# submitter); A3_RUN names the run (never $RUN: FarmShare ~/.bashrc exports it). The weights and
# Unimod table must match validation/pilot/a3_assets.tsv (exit 2 otherwise, R6). Work files go to
# $GIAP_ROOT/pilot/A3_cascadia_cpu/<RUN>/shards; the merged table is <RUN>/<RUN>.ssl.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
R="${A3_RUN:?A3_RUN}"; STEP="${STEP:?STEP}"
ASSETS="$GIAP_REPO/validation/pilot/a3_assets.tsv"
MODEL="$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt"; UNIMOD="$GIAP_DB_DIR/cascadia/unimod_tables.xml"
OUT="$GIAP_ROOT/pilot/A3_cascadia_cpu/$R"; WORK="$OUT/shards"; mkdir -p "$OUT"
for pair in "cascadia_weights $MODEL" "unimod_table $UNIMOD"; do
  set -- $pair
  want=$(awk -F'\t' -v a="$1" 'NR==1{for(i=1;i<=NF;i++)if($i=="sha256")k=i;next}$1==a{print $k;exit}' "$ASSETS")
  got=$(sha256sum -- "$2" | cut -d' ' -f1)
  [ "${#want}" -eq 64 ] && [ "$got" = "$want" ] || { echo "$1 sha256 $got != '$want'" >&2; exit 2; }
done
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}" MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export HDF5_USE_FILE_LOCKING=FALSE
PY=(micromamba run -n giap-cascadia python -I "$GIAP_REPO/scripts/pilot/cascadia_shards.py" "$UNIMOD")
case "$STEP" in
  prep) "${PY[@]}" prep "$GIAP_DATA_DIR/pxd078459/mzml/$R.mzML" "$MODEL" "$WORK" ;;
  predict) "${PY[@]}" predict "$MODEL" "$WORK" "$SLURM_ARRAY_TASK_ID" "$N_SHARDS" ${LIMIT:+--limit "$LIMIT"} ;;
  merge) "${PY[@]}" merge "$WORK" "$OUT/$R"; ls -l "$OUT/$R.ssl" ;;
  *) echo "STEP must be prep, predict or merge" >&2; exit 2 ;;
esac
