#!/bin/bash
# Launcher for the `calibration_full` tier (CI-01 decision 3), the nightly_launch.sh pattern:
#   1. sbatch --array=0-<CALFULL_SHARDS-1> scripts/slurm/calibration_full.sbatch   (the shards)
#   2. scripts/slurm/nightly_summary.sh <array> as a 1-CPU giap-calfull-summary job with
#      --dependency=afterany:<array>, pointed at calibration_full.tsv: post_summary.py --final
#      posts the table even when a shard was killed before it could write its line (its sacct
#      state fills the row). The same summary script serves both jobs; its post_summary.py
#      options (--tsv, --expect, --title) are overridden by the ones appended here.
# No scrontab entry (decision 4): run this by hand. CALFULL_SHARDS (default 8, see the header of
# calibration_full.sbatch) sets both the array range and --expect, so the summary waits for exactly
# the shards that were submitted. NIGHTLY_POST=dry (the default) prints the comment instead of
# posting it; NIGHTLY_POST=yes posts to the issue titled $CALFULL_TITLE. Other CALFULL_* variables
# reach the array (--export=ALL). A refused array submission exits 1 and submits no summary job.
# R14: CALFULL_SHARDS x --cpus-per-task of the array + 1 CPU for the summary job (8 x 8 + 1 = 65).
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
REPO="$GIAP_REPO"
SHARDS="${CALFULL_SHARDS:-8}"
TITLE="${CALFULL_TITLE:-Calibration full tier (FarmShare)}"
TSV="$GIAP_ROOT/logs/calibration_full.tsv"
case "$SHARDS" in ''|*[!0-9]*|0)
  echo "calibration_full_launch: CALFULL_SHARDS='$SHARDS' is not a positive integer" >&2; exit 2 ;;
esac
mkdir -p "$REPO/local/logs" "$GIAP_ROOT/logs"
array=$(sbatch --parsable --array="0-$((SHARDS - 1))" --export="ALL,CALFULL_SHARDS=$SHARDS" \
  "$REPO/scripts/slurm/calibration_full.sbatch") || {
  echo "calibration_full_launch: the array was not submitted" >&2; exit 1; }
array=${array%%;*}
dry=(); [ "${NIGHTLY_POST:-dry}" = dry ] && dry=(--dry-run)
summary=$(sbatch --parsable --job-name=giap-calfull-summary --cpus-per-task=1 \
  --dependency="afterany:$array" --comment="calfull:$array" --output="$REPO/local/logs/%x.%j.out" \
  "$REPO/scripts/slurm/nightly_summary.sh" "$array" \
  --tsv "$TSV" --expect "$SHARDS" --title "$TITLE" ${dry[@]+"${dry[@]}"}) || {
  echo "calibration_full_launch: array $array submitted, but its summary job was refused" >&2
  exit 1; }
echo "calibration_full_launch: array ${array} ($SHARDS shards), summary job ${summary%%;*} (afterany:$array)"
