#!/bin/bash
#SBATCH --job-name=giap-pilot-compare
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# giap compare of the pass-1 pilot arms A1 and A2, twice, then rule set C (M10-04).
#
#   sbatch -J giap-pilot-compare scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 CMP
#
# CMP must not exist (no merging into stale output). The second compare goes to $TMPDIR/cmp2 and
# must be byte-identical (rule C6). The evidence goes to $(dirname CMP)/evidence/compare.json.
# Exit codes: 2 bad arguments; the first giap compare's code when it fails (1: CompareInputError);
# else run_evidence.py's (0 pass, 1 a rule of set C failed).
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
set -f

usage() {
  echo "usage: sbatch scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 CMP" >&2
  exit 2
}
[ $# -eq 3 ] || usage
OUT_A1="$1"; OUT_A2="$2"; CMP="$3"
[ -d "$OUT_A1" ] || { echo "not a run directory: $OUT_A1" >&2; usage; }
[ -d "$OUT_A2" ] || { echo "not a run directory: $OUT_A2" >&2; usage; }
if [ -e "$CMP" ]; then
  echo "CMP exists: $CMP (use a new name or remove it; no merging into stale output)" >&2
  exit 2
fi

GIAP_PYTHON="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
export PYTHONPATH="$GIAP_REPO/src"
TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/giap.XXXXXX")"; export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
EVIDENCE="$(dirname "$CMP")/evidence/compare.json"

giap_compare() {  # OUTDIR
  "$GIAP_PYTHON" -c 'from giap.cli import cli; cli()' compare "$OUT_A1" "$OUT_A2" \
    --outdir "$1" --labels A1,A2
}

rc=0
giap_compare "$CMP" || rc=$?
if [ "$rc" -ne 0 ]; then
  echo "giap compare exited $rc" >&2
  exit "$rc"
fi
rc=0
giap_compare "$TMPDIR/cmp2" || rc=$?
if [ "$rc" -ne 0 ]; then
  echo "the second giap compare exited $rc" >&2
  exit "$rc"
fi

rc=0
"$GIAP_PYTHON" "$GIAP_REPO/validation/pilot/run_evidence.py" compare "$CMP" \
  --second "$TMPDIR/cmp2" --labels A1,A2 --out "$EVIDENCE" || rc=$?
echo "evidence: $EVIDENCE (exit $rc)"
exit "$rc"
