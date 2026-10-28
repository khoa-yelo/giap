#!/bin/bash
#SBATCH --job-name=gippa-cascadia-convert
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000M
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Per-run L6 check and gippa convert --from cascadia of the three A3 Cascadia runs (M10-10a,
# docs/specs/M10-10-a3-run.md, "cascadia_convert.sbatch (step 3)"). No arguments.
#
#   sbatch -J gippa-cascadia-convert --dependency=afterok:$d1:$d2:$d3 scripts/slurm/cascadia_convert.sbatch
#
# Reads $A3/<RUN>/<RUN>.ssl (A3=${A3_OUT_DIR:-$GIPPA_ROOT/pilot/A3}; runs from acquisition.json);
# writes $A3/ssl.sha256, $A3/<RUN>/feasibility.json and $GIPPA_ROOT/pilot/inputs/a3.tsv.gz,
# a3.convert.log, a3.convert.json. Nothing is written under $GIPPA_REPO, nothing is committed.
# Exit 0 converted and checked; 1 L6 evaluation error, conversion, determinism or summary failure
# (a3.FAILED); 2 a precondition; 3 the gate (M10-08 verdict fail) or L6 failing on >= 2 runs.
set -euo pipefail
set -f
die2() { echo "cascadia_convert: $*" >&2; exit 2; }
[ $# -eq 0 ] || die2 "usage: sbatch scripts/slurm/cascadia_convert.sbatch (no arguments)"
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT
GIPPA_JOB_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
export GIPPA_JOB_START

A3="${A3_OUT_DIR:-$GIPPA_ROOT/pilot/A3}"
OUT="$GIPPA_ROOT/pilot/inputs"
GIPPA_PYTHON="${GIPPA_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
PILOT="$GIPPA_REPO/validation/pilot"
ACQ="$PILOT/inputs/acquisition.json"

# ---- 1. the entry gate (exit 3 passed through; nothing written)
rc=0
"$GIPPA_PYTHON" -I "$PILOT/cascadia_convert.py" gate || rc=$?
[ "$rc" -ne 3 ] || exit 3
[ "$rc" -eq 0 ] || die2 "entry gate exited $rc (M10-08 not done)"

# ---- 2. the runs, in acquisition.json order; each .ssl exists and is non-empty
RUNS="$("$GIPPA_PYTHON" -I -c 'import json, sys
print(" ".join(r["run"] for r in json.load(open(sys.argv[1]))["runs"]))' "$ACQ")" \
  || die2 "cannot read the runs of $ACQ"
[ -n "$RUNS" ] || die2 "$ACQ lists no run"
SSLS=""
for r in $RUNS; do
  [ -s "$A3/$r/$r.ssl" ] || die2 "run $r: $A3/$r/$r.ssl is missing or empty"
  SSLS="$SSLS $A3/$r/$r.ssl"
done
R_REAL="$(cd "$GIPPA_REPO" && pwd -P)" || die2 "cannot resolve GIPPA_REPO=$GIPPA_REPO"
mkdir -p "$OUT"
O_REAL="$(cd "$OUT" && pwd -P)" || die2 "cannot resolve $OUT"
case "$O_REAL/" in
  "$R_REAL"/*) git -C "$R_REAL" check-ignore -q "$O_REAL" || die2 "$OUT is inside the checkout" ;;
esac

# ---- from here on every exit-1 failure leaves a3.FAILED and no temporaries
JOB="${SLURM_JOB_ID:-$$}"
NAMES="a3.tsv.gz a3.convert.log a3.convert.json"
T_PEP="$OUT/.a3.tsv.gz.tmp.$JOB.gz"
T_LOG="$OUT/.a3.convert.log.tmp.$JOB"
T_JSON="$OUT/.a3.convert.json.tmp.$JOB"
cleanup() { rm -f "$T_PEP" "$T_LOG" "$T_JSON"; rm -rf "$TMPDIR"; }
trap cleanup EXIT
fail() {
  trap - ERR
  echo "cascadia_convert: $*" >&2
  printf 'exit 1\n%s\n%s\n' "$*" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/a3.FAILED"
  exit 1
}
trap 'fail "command failed at line $LINENO"' ERR
rm -f "$OUT/a3.FAILED"
mkdir "$TMPDIR/previous"  # older outputs go aside (discarded): no stale table after a failure
for name in $NAMES; do
  if [ -e "$OUT/$name" ]; then mv -f "$OUT/$name" "$TMPDIR/previous/$name"; fi
done
for s in $SSLS; do sha256sum "$s"; done | sort -k2 > "$TMPDIR/ssl.sha256"
mv -f "$TMPDIR/ssl.sha256" "$A3/ssl.sha256"

# ---- 3. per-run L6 check (M10-08's evaluate); >= 2 failing runs stop the job (exit 3)
nfail=0
for r in $RUNS; do
  rc=0
  "$GIPPA_PYTHON" -I "$PILOT/a3_feasibility.py" evaluate "$A3/$r/$r.ssl" --run "$r" \
    --out "$A3/$r/feasibility.json" || rc=$?
  case "$rc" in
    0) echo "cascadia_convert: L6 passes on $r" >&2 ;;
    3) nfail=$((nfail + 1)); echo "cascadia_convert: L6 fails on $r (reported)" >&2 ;;
    *) fail "L6 evaluation of $r exited $rc" ;;
  esac
done
n_runs=$(echo $RUNS | wc -w | tr -d ' ')
if [ "$nfail" -ge 2 ]; then
  echo "cascadia_convert: L6 fails on $nfail of $n_runs runs: lead decision before A3 runs" >&2
  exit 3
fi

# ---- 4-5. convert, then the same conversion again outside $OUT (determinism)
convert() {  # convert <peptides> <log>
  # shellcheck disable=SC2086
  PYTHONPATH="$GIPPA_REPO/src" "$GIPPA_PYTHON" -P -c 'from gippa.cli import cli; cli()' \
    convert --from cascadia $SSLS --out "$1" --source denovo 2> "$2"
}
rc=0
convert "$T_PEP" "$T_LOG" || rc=$?
cat "$T_LOG" >&2
[ "$rc" -eq 0 ] || fail "gippa convert exited $rc"
rc=0
convert "$TMPDIR/b.tsv.gz" "$TMPDIR/b.log" || rc=$?
[ "$rc" -eq 0 ] || fail "second gippa convert exited $rc"
cmp "$T_PEP" "$TMPDIR/b.tsv.gz" || fail "the peptide table differs between two conversions"

# ---- 6. summary
rc=0
"$GIPPA_PYTHON" -I "$PILOT/cascadia_convert.py" summary --ssl-sums "$A3/ssl.sha256" \
  --peptides "$T_PEP" --log "$T_LOG" --acquisition "$ACQ" --out "$T_JSON" || rc=$?
[ "$rc" -eq 0 ] || fail "summary rules failed (exit $rc)"

# ---- 7. move into place, only now (convert.json last); a failed move unpublishes the rest
published=""
for name in $NAMES; do
  case "$name" in
    a3.tsv.gz) src="$T_PEP" ;; a3.convert.log) src="$T_LOG" ;; *) src="$T_JSON" ;;
  esac
  if ! mv -f "$src" "$OUT/$name"; then
    for p in $published; do rm -f "$OUT/$p"; done
    fail "could not move $name into $OUT"
  fi
  published="$published $name"
done
echo "cascadia_convert: $n_runs runs converted and checked into $OUT"
