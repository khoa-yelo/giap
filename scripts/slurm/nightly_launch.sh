#!/bin/bash
# Nightly launcher (M2-21 lead decisions 14, 18, 19), run by scrontab (scripts/slurm/nightly.scrontab)
# in a 1-CPU job, because scrontab may not accept --array:
#   1. sbatch --array=0-3 scripts/slurm/nightly.sbatch               (the four sub-tasks, 16 CPUs each)
#   2. scripts/slurm/nightly_summary.sh <array> as a 1-CPU giap-nightly-summary job,
#      --dependency=afterany:<array>: post_summary.py --final posts the night even when a task was
#      killed before it could write its nightly.tsv line (missing tasks get their sacct state). Its
#      --comment nightly:<array> lets a resubmitted benchmark task find it and extend its dependency.
# A refused array submission exits 1 and submits no summary job. NIGHTLY_POST=dry: the summary job
# prints the comment instead of posting. Other NIGHTLY_* variables reach the array (--export=ALL).
# R14: 4 x 16 + 1 CPUs (+ this launcher's 1 CPU, which ends before the array starts).
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
REPO="$GIAP_REPO"
mkdir -p "$REPO/local/logs" "$GIAP_ROOT/logs"
array=$(sbatch --parsable --array=0-3 "$REPO/scripts/slurm/nightly.sbatch") || {
  echo "nightly_launch: the nightly array was not submitted" >&2; exit 1; }
array=${array%%;*}
dry=(); [ "${NIGHTLY_POST:-yes}" = dry ] && dry=(--dry-run)
summary=$(sbatch --parsable --job-name=giap-nightly-summary --cpus-per-task=1 \
  --dependency="afterany:$array" --comment="nightly:$array" --output="$REPO/local/logs/%x.%j.out" \
  "$REPO/scripts/slurm/nightly_summary.sh" "$array" "${dry[@]}") || {
  echo "nightly_launch: array $array submitted, but its summary job was refused" >&2; exit 1; }
echo "nightly_launch: array ${array}, summary job ${summary%%;*} (afterany:$array)"
