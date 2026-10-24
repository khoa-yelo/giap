#!/bin/bash
#SBATCH --job-name=giap-mzml-pilot
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000M
#SBATCH --time=02:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# .d -> mzML for ONE PXD078459 run, then the structural mzML check (M10-08 deliverable 7).
#
#   sbatch scripts/slurm/mzml_pilot.sbatch WT_1
#
# MZML_ROUTE=timsconvert (default, route B2: TIMSCONVERT 2.0.0 + pyTDFSDK's TDF-SDK, ADR-0054)
# or pwiz (route B1: ProteoWizard msconvert under Wine; needs PWIZ_SIF, an image pulled by digest,
# and Khoa's acceptance of the five non-Bruker EULAs, which he has not given). Reads
# $GIAP_DATA_DIR/pxd078459/<RUN>.d, writes $GIAP_DATA_DIR/pxd078459/mzml/<RUN>.mzML and
# <RUN>.mzml_check.json (validation/pilot/a3_feasibility.py mzml-check). No network (R6).
# Exit codes: 0 converted and the check passed; 1 conversion or check failed; 2 a precondition
# (no run, unknown run or route, less free space than 4 x the .d size).
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"

RUN="${1:-}"
ROUTE="${MZML_ROUTE:-timsconvert}"
DATA="$GIAP_DATA_DIR/pxd078459"
SRC="$DATA/$RUN.d"
OUT="$DATA/mzml"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"

[ -n "$RUN" ] || { echo "usage: mzml_pilot.sbatch RUN   (e.g. WT_1)" >&2; exit 2; }
case "$ROUTE" in
  timsconvert) ;;
  pwiz) [ -f "${PWIZ_SIF:-}" ] || { echo "MZML_ROUTE=pwiz needs PWIZ_SIF (an image)" >&2; exit 2; } ;;
  *) echo "unknown MZML_ROUTE: $ROUTE (timsconvert | pwiz)" >&2; exit 2 ;;
esac
[ -d "$SRC" ] || { echo "no such run: $SRC" >&2; exit 2; }
ACQ="$GIAP_REPO/validation/pilot/inputs/acquisition.json"
"$PY" -I -c 'import json, sys
sys.exit(sys.argv[2] not in {r["run"] for r in json.load(open(sys.argv[1]))["runs"]})' "$ACQ" "$RUN" \
  || { echo "run $RUN not in $ACQ: not converting" >&2; exit 2; }

# Disk: free space under $GIAP_DATA_DIR >= 4 x the .d size (apparent bytes).
size=$("$PY" -I -c 'import os, sys
print(sum(os.path.getsize(os.path.join(d, f)) for d, _, fs in os.walk(sys.argv[1]) for f in fs))' \
  "$SRC")
avail=$(( $(df -Pk "$DATA" | awk 'NR == 2 { print $4 }') * 1024 ))
echo "== $RUN.d: $size bytes; free under $DATA: $avail bytes"
if [ "$avail" -lt $(( 4 * size )) ]; then
  echo "less than 4 x $size bytes free under $DATA: not converting" >&2
  exit 2
fi

mkdir -p "$OUT"
rm -f "$OUT/$RUN.mzML" "$OUT/$RUN.mzml_check.json"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/giap-mzml.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
cd "$WORK"
rc=0
if [ "$ROUTE" = "timsconvert" ]; then
  micromamba run -n giap-timsconvert timsconvert --help 2>&1 | head -n 40 || true
  micromamba run -n giap-timsconvert timsconvert --input "$SRC" --outdir "$OUT" \
    --mode centroid --compression zlib --exclude_mobility || rc=$?
else
  pwiz() { apptainer exec --cleanenv --writable-tmpfs -B "$GIAP_DATA_DIR" "$PWIZ_SIF" "$@"; }
  pwiz wine msconvert --help 2>&1 | head -n 40 || true
  pwiz wine msconvert "$SRC" --mzML --zlib --64 --combineIonMobilitySpectra \
    --filter "peakPicking vendor msLevel=1-" -o "$OUT" --outfile "$RUN.mzML" || rc=$?
fi
if [ "$rc" -ne 0 ] || [ ! -s "$OUT/$RUN.mzML" ]; then
  echo "conversion failed (route $ROUTE, exit $rc): R9" >&2
  exit 1
fi
ls -l "$OUT/$RUN.mzML"

PYTHONPATH="$GIAP_REPO/src" "$PY" "$GIAP_REPO/validation/pilot/a3_feasibility.py" mzml-check \
  "$OUT/$RUN.mzML" --acquisition "$ACQ" \
  --run "$RUN" --out "$OUT/$RUN.mzml_check.json" || rc=$?
cat "$OUT/$RUN.mzml_check.json" 2>/dev/null || true
[ "$rc" -eq 0 ] || echo "mzml-check failed for $RUN (route $ROUTE): a route failure, not a verdict" >&2
exit "$rc"
