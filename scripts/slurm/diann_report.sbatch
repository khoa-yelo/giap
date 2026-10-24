#!/bin/bash
#SBATCH --job-name=giap-diann-report
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000M
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Check one arm's DIA-NN report, then write its stratified sample (docs/specs/M10-02-diann-searches.md).
#
#   sbatch -J giap-diann-report-A1 --dependency=afterok:<search job> scripts/slurm/diann_report.sbatch A1
#
# check.json goes next to the report; the sample (report.sample.parquet, schema.tsv, summary.json)
# goes to $GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/ and is never committed (Khoa, 2026-10-23).
# Exit 0 ok; 1 the report fails a check; 2 usage or no report.
# --mem=16000M is 4 x MaxMemPerCPU (4000 MB): 16G would make Slurm allocate a fifth CPU (R14).
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

[ $# -eq 1 ] || { echo "usage: sbatch scripts/slurm/diann_report.sbatch A1|A2" >&2; exit 2; }
ARM="$1"
case "$ARM" in
  A1 | A2) ;;
  *) echo "unknown arm '$ARM' (A1 or A2)" >&2; exit 2 ;;
esac
PY="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
REPORT_PY="$GIAP_REPO/validation/pilot/diann_report.py"
ACQ_JSON="${ACQUISITION_JSON:-$GIAP_REPO/validation/pilot/inputs/acquisition.json}"
OUT="${OUT_DIR:-$GIAP_ROOT/pilot}/$ARM"
SAMPLE_DIR="$GIAP_ROOT/pilot/recorded/diann/2.3.1/$ARM"
REPORT="$OUT/${ARM}_report.parquet"

[ -f "$REPORT" ] || { echo "no report $REPORT: run scripts/slurm/diann_pilot.sbatch $ARM first" >&2; exit 2; }
"$PY" -I "$REPORT_PY" check "$REPORT" --acquisition "$ACQ_JSON" --sums "$OUT/outputs.sha256" \
  --out "$OUT/check.json"
"$PY" -I "$REPORT_PY" sample "$REPORT" --out "$SAMPLE_DIR"
echo "== $ARM: $OUT/check.json ok; sample under $SAMPLE_DIR"
