#!/bin/bash
#SBATCH --job-name=gippa-pilot-compare
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# gippa compare of the pilot arms, twice, then rule set C (M10-04; three arms: M10-10).
#
#   sbatch -J gippa-pilot-compare scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 CMP
#   sbatch -J gippa-pilot-compare3 scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 OUT_A3 CMP
#
# CMP must not exist (no merging into stale output). The second compare goes to $TMPDIR/cmp2 and
# must be byte-identical (rule C6). The evidence goes to $(dirname CMP)/evidence/compare.json
# (labels A1,A2) or compare3.json (labels A1,A2,A3). With three arms, pass 1's
# $(dirname OUT_A1)/cmp_a1_a2, when it exists, is the reference of rule C7 (a warning only).
# Exit codes: 2 bad arguments; the first gippa compare's code when it fails (1: CompareInputError);
# else run_evidence.py's (0 pass, 1 a rule of set C failed).
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
set -f

usage() {
  echo "usage: sbatch scripts/slurm/pilot_compare.sbatch OUT_A1 OUT_A2 [OUT_A3] CMP" >&2
  exit 2
}
case $# in
  3) RUNS=("$1" "$2"); CMP="$3"; LABELS=A1,A2; NAME=compare.json ;;
  4) RUNS=("$1" "$2" "$3"); CMP="$4"; LABELS=A1,A2,A3; NAME=compare3.json ;;
  *) usage ;;
esac
OUT_A1="$1"
for d in "${RUNS[@]}"; do
  [ -d "$d" ] || { echo "not a run directory: $d" >&2; usage; }
done
if [ -e "$CMP" ]; then
  echo "CMP exists: $CMP (use a new name or remove it; no merging into stale output)" >&2
  exit 2
fi

GIPPA_PYTHON="${GIPPA_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
export PYTHONPATH="$GIPPA_REPO/src"
TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/gippa.XXXXXX")"; export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
EVIDENCE="$(dirname "$CMP")/evidence/$NAME"
PAIR_ARGS=()
PAIR_CMP="$(dirname "$OUT_A1")/cmp_a1_a2"
if [ $# -eq 4 ] && [ -d "$PAIR_CMP" ]; then
  PAIR_ARGS=(--pair-cmp "$PAIR_CMP")
fi

gippa_compare() {  # OUTDIR
  "$GIPPA_PYTHON" -c 'from gippa.cli import cli; cli()' compare "${RUNS[@]}" \
    --outdir "$1" --labels "$LABELS"
}

rc=0
gippa_compare "$CMP" || rc=$?
if [ "$rc" -ne 0 ]; then
  echo "gippa compare exited $rc" >&2
  exit "$rc"
fi
rc=0
gippa_compare "$TMPDIR/cmp2" || rc=$?
if [ "$rc" -ne 0 ]; then
  echo "the second gippa compare exited $rc" >&2
  exit "$rc"
fi

rc=0
# ${PAIR_ARGS[@]+...}: bash 3.2 (macOS CI) never expands an empty array under `set -u`
"$GIPPA_PYTHON" "$GIPPA_REPO/validation/pilot/run_evidence.py" compare "$CMP" \
  --second "$TMPDIR/cmp2" --labels "$LABELS" --genes "$OUT_A1/tables/genes.tsv" \
  ${PAIR_ARGS[@]+"${PAIR_ARGS[@]}"} --out "$EVIDENCE" || rc=$?
echo "evidence: $EVIDENCE (exit $rc)"
exit "$rc"
