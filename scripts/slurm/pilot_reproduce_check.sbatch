#!/bin/bash
#SBATCH --job-name=giap-pilot-repro-check
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# S-check of the pilot runbook at the freeze commit F (docs/specs/M10-14-pilot-reproducibility.md,
# "S-check"): 1. reproduce_check.py static; 2. reproduce_check.py verify into
# $GIAP_ROOT/pilot/reproduce/verify-<F7>.json; 3. the re-conversion of each DIA-NN report at F into
# $TMPDIR, compared with the tables under $PILOT/inputs (identical | changed | absent).
#
#   sbatch scripts/slurm/pilot_reproduce_check.sbatch
#
# RUNBOOK (default $GIAP_REPO/validation/pilot/reproduce.md); PILOT (default $GIAP_ROOT/pilot).
# Nothing is written under $GIAP_REPO. Exit 0 checked (a changed table is a warning, example 8);
# 1 static, verify or a re-conversion failed (named in the JSON's failures[]); 2 a precondition.
set -euo pipefail
set -f
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"
export TMPDIR
trap 'rm -rf "$TMPDIR"' EXIT

die2() { echo "pilot_reproduce_check: $*" >&2; exit 2; }

GIAP_PYTHON="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
RUNBOOK="${RUNBOOK:-$GIAP_REPO/validation/pilot/reproduce.md}"
CHECK="$GIAP_REPO/validation/pilot/reproduce_check.py"
PILOT="${PILOT:-$GIAP_ROOT/pilot}"
[ -f "$RUNBOOK" ] || die2 "no runbook at $RUNBOOK"
[ -x "$GIAP_PYTHON" ] || die2 "GIAP_PYTHON=$GIAP_PYTHON is not executable"

HEAD="$(git -C "$GIAP_REPO" rev-parse HEAD)" || die2 "cannot read HEAD of $GIAP_REPO"
F="$("$GIAP_PYTHON" -I "$CHECK" protected "$RUNBOOK" --repo "$GIAP_REPO" | head -n 1)"
if [ "$F" = "pending" ]; then
  echo "pilot_reproduce_check: freeze_commit=pending; checking HEAD $HEAD" >&2
  F="$HEAD"
fi
[ "$F" = "$HEAD" ] || die2 "the runbook names F=$F but $GIAP_REPO is at $HEAD"
[ -z "$(git -C "$GIAP_REPO" status --porcelain --untracked-files=no)" ] \
  || die2 "$GIAP_REPO has uncommitted changes"
F7="$(printf '%s' "$F" | cut -c1-7)"
OUTDIR="$GIAP_ROOT/pilot/reproduce"
OUTJ="$OUTDIR/verify-$F7.json"
mkdir -p "$OUTDIR"
echo "== S-check at $F ($RUNBOOK)"

# ---- 1. static, 2. verify (rc: 0 or 1 for the job; static's goes to the JSON as static_ok)
rc=0
static_ok=1
"$GIAP_PYTHON" -I "$CHECK" static "$RUNBOOK" --repo "$GIAP_REPO" || { static_ok=0; rc=1; }
v=0
"$GIAP_PYTHON" -I "$CHECK" verify "$RUNBOOK" --repo "$GIAP_REPO" --out "$OUTJ" || v=$?
[ "$v" -le 1 ] || die2 "verify exited $v"
[ "$v" -eq 0 ] || rc=1

# ---- 3. re-conversion at F (before S15 deletes anything)
RESULTS="$TMPDIR/reconvert.tsv"
: > "$RESULTS"
for ARM in A1 A2; do
  case "$ARM" in
    A1) STEM=a1; SEARCH_DB=other ;;
    A2) STEM=a2; SEARCH_DB=same ;;
  esac
  REPORT="$PILOT/$ARM/${ARM}_report.parquet"
  if [ ! -f "$REPORT" ]; then
    printf '%s\tabsent\n' "$STEM" >> "$RESULTS"
    continue
  fi
  if ! PYTHONPATH="$GIAP_REPO/src" "$GIAP_PYTHON" -P -c 'from giap.cli import cli; cli()' \
    convert --from diann "$REPORT" --out "$TMPDIR/$STEM.tsv.gz" \
    --abundance-out "$TMPDIR/${STEM}_abundance.tsv.gz" --search-db "$SEARCH_DB" \
    2> "$TMPDIR/$STEM.log"; then
    cat "$TMPDIR/$STEM.log" >&2
    echo "pilot_reproduce_check: giap convert $ARM failed" >&2
    printf '%s\tfailed\n' "$STEM" >> "$RESULTS"
    rc=1
    continue
  fi
  result=identical
  for name in "$STEM.tsv.gz" "${STEM}_abundance.tsv.gz"; do
    new="$(sha256sum "$TMPDIR/$name" | awk '{print $1}')"
    old="$(sha256sum "$PILOT/inputs/$name" 2>/dev/null | awk '{print $1}')" || old=""
    [ "$new" = "$old" ] || result=changed
  done
  [ "$result" = identical ] || echo "WARNING: $ARM re-converts to a changed table (example 8)" >&2
  printf '%s\t%s\n' "$STEM" "$result" >> "$RESULTS"
done
"$GIAP_PYTHON" -I -c '
import json, sys
doc = json.load(open(sys.argv[1]))
doc["reconvert"] = dict(ln.rstrip("\n").split("\t") for ln in open(sys.argv[2]) if ln.strip())
doc["failures"] += [f"reconvert:{k}" for k, v in doc["reconvert"].items() if v == "failed"]
doc["static_ok"] = sys.argv[3] == "1"
json.dump(doc, open(sys.argv[1], "w"), indent=2, sort_keys=True)
print("reconvert:", doc["reconvert"])' "$OUTJ" "$RESULTS" "$static_ok" || die2 "cannot update $OUTJ"

echo "== $OUTJ"
[ "$rc" -eq 0 ] || { echo "pilot_reproduce_check: a check failed (see $OUTJ)" >&2; exit 1; }
echo "pilot_reproduce_check: ok at $F"
