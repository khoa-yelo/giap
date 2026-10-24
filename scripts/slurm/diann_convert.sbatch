#!/bin/bash
#SBATCH --job-name=giap-diann-convert
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# giap convert --from diann on one pilot arm's DIA-NN 2.3.1 report (M10-03,
# docs/specs/M10-03-diann-convert.md, "The conversion job").
#
#   sbatch -J giap-diann-convert-A1 scripts/slurm/diann_convert.sbatch A1
#
# Reads $PILOT/<arm>/ (PILOT=${OUT_DIR:-$GIAP_ROOT/pilot}); writes $PILOT/inputs/<a>.tsv.gz,
# <a>_abundance.tsv.gz, <a>.convert.log and <a>.convert.json. Nothing is written under $GIAP_REPO
# and nothing here is committed (Khoa 2026-10-23). Exit 0 converted and checked; 1 conversion,
# determinism, summary or contract failure (<a>.FAILED holds the exit code); 2 a precondition.
set -euo pipefail
set -f
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
GIAP_JOB_START="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
export GIAP_JOB_START

die2() { echo "diann_convert: $*" >&2; exit 2; }

[ $# -eq 1 ] || die2 "usage: sbatch scripts/slurm/diann_convert.sbatch A1|A2"
ARM="$1"
case "$ARM" in
  A1) STEM=a1; SEARCH_DB=other ;;
  A2) STEM=a2; SEARCH_DB=same ;;
  *) die2 "unknown arm '$ARM' (A1 or A2)" ;;
esac
PILOT="${OUT_DIR:-$GIAP_ROOT/pilot}"
GIAP_PYTHON="${GIAP_PYTHON:-$UV_PROJECT_ENVIRONMENT/bin/python}"
REPORT="$PILOT/$ARM/${ARM}_report.parquet"
CHECK="$PILOT/$ARM/check.json"
SUMS="$PILOT/$ARM/outputs.sha256"

# ---- preconditions (exit 2, nothing written)
[ -s "$REPORT" ] || die2 "report missing or empty: $REPORT"
[ -f "$CHECK" ] || die2 "missing $CHECK (M10-02 check not run)"
"$GIAP_PYTHON" -I -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1])).get("ok") is True else 1)' "$CHECK" \
  || die2 "$CHECK is not \"ok\": true"
[ -f "$SUMS" ] || die2 "missing $SUMS"
WANT="$(awk -v n="${ARM}_report.parquet" '{f=$2; sub(/^\*/, "", f); sub(/.*\//, "", f)} f == n {print $1; exit}' "$SUMS")"
[ -n "$WANT" ] || die2 "$SUMS has no line for ${ARM}_report.parquet"
GOT="$(sha256sum "$REPORT" | awk '{print $1}')"
[ "$GOT" = "$WANT" ] || die2 "sha256 of $REPORT is $GOT, $SUMS lists $WANT"
[ ! -e "$PILOT/$ARM/FAILED" ] || die2 "$PILOT/$ARM/FAILED exists (the M10-02 search failed)"

# ---- outputs under temporary names
OUT="$PILOT/inputs"
mkdir -p "$OUT"
chmod 0755 "$OUT"
rm -f "$OUT/$STEM.FAILED"
JOB="${SLURM_JOB_ID:-$$}"
NAMES="$STEM.tsv.gz ${STEM}_abundance.tsv.gz $STEM.convert.log $STEM.convert.json"
tmp() {  # .<name>.tmp.<job>, keeping .gz so the CLI still gzips
  case "$1" in *.gz) echo "$OUT/.$1.tmp.$JOB.gz" ;; *) echo "$OUT/.$1.tmp.$JOB" ;; esac
}
T_PEP="$(tmp "$STEM.tsv.gz")"
T_AB="$(tmp "${STEM}_abundance.tsv.gz")"
T_LOG="$(tmp "$STEM.convert.log")"
T_JSON="$(tmp "$STEM.convert.json")"
cleanup() {
  rm -f "$T_PEP" "$T_AB" "$T_LOG" "$T_JSON"
  rm -rf "$TMPDIR"
}
trap cleanup EXIT
fail() {
  echo "diann_convert: $*" >&2
  printf 'exit 1\n%s\n%s\n' "$*" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$OUT/$STEM.FAILED"
  exit 1
}

convert() {  # convert <peptides> <abundance> <log>
  PYTHONPATH="$GIAP_REPO/src" "$GIAP_PYTHON" -c 'from giap.cli import cli; cli()' \
    convert --from diann "$REPORT" --out "$1" --abundance-out "$2" --search-db "$SEARCH_DB" \
    2> "$3"
}

rc=0
convert "$T_PEP" "$T_AB" "$T_LOG" || rc=$?
cat "$T_LOG" >&2
[ "$rc" -eq 0 ] || fail "giap convert exited $rc"

# ---- determinism: the same conversion again, outside $OUT
rc=0
convert "$TMPDIR/b.tsv.gz" "$TMPDIR/b_abundance.tsv.gz" "$TMPDIR/b.log" || rc=$?
[ "$rc" -eq 0 ] || fail "second giap convert exited $rc"
cmp "$T_PEP" "$TMPDIR/b.tsv.gz" || fail "the peptide table differs between two conversions"
cmp "$T_AB" "$TMPDIR/b_abundance.tsv.gz" || fail "the abundance table differs between two conversions"

# ---- summary
rc=0
"$GIAP_PYTHON" -I "$GIAP_REPO/validation/pilot/diann_convert.py" summary --arm "$ARM" \
  --report "$REPORT" --peptides "$T_PEP" --abundance "$T_AB" --log "$T_LOG" --out "$T_JSON" || rc=$?
[ "$rc" -eq 0 ] || fail "summary rules failed (exit $rc)"

# ---- full-report contract checks
rc=0
(cd "$GIAP_REPO" && PYTHONPATH="$GIAP_REPO/src:$GIAP_REPO/tests/sim:$GIAP_REPO/tests" \
  "$GIAP_PYTHON" -m pytest -m slow -p no:cacheprovider \
  tests/contract/test_m10_03_diann_contract.py -k full_report) || rc=$?
[ "$rc" -eq 0 ] || fail "contract checks failed (pytest exit $rc)"

# ---- move into place, only now
for name in $NAMES; do
  mv -f "$(tmp "$name")" "$OUT/$name"
done
echo "diann_convert: $ARM converted and checked into $OUT"
