#!/bin/bash
#SBATCH --job-name=giap-weekly
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Weekly end-to-end runs on origin/main (GIAP_SPEC.md §12; M2-21), one array task per genome:
#   0 = k12_mg1655 (U00096.3), 1 = nissle1917 (CP171242-CP171245), from $GIAP_ROOT/data/genomes
#   (every replicon marked [topology=circular]: all are closed).
# Per genome: (a) the farmshare_weekly capability gate (pytest --require-capabilities, the Bakta
# full-DB contract tests; a missing DB fails); (b) real `giap run` with live full-DB Bakta and every
# annotation tool, including the antiSMASH live run (its trace row must be COMPLETED); (c) `giap
# run` with the recorded full-DB JSON ($GIAP_ROOT/recordings/bakta_full/<genome>, --annotation);
# (d) K-12 only: TE6 (§11.6) on the live Bakta JSON: an oriC call whose nearest gene is mnmG or
# mioC with mnmG within 2 kb. Measured runtimes: scripts/measure_runtimes.py on each trace.
# Evidence stand-in (M2 annotation track only; MAPPING is M3): `giap run` requires --proteins or
# --peptides, so the evidence is a deterministic 60 % subset of the recorded Bakta proteins
# (records 0,1,2 of every 5 of bakta/<genome>.faa), as in the lead's hand run of 2026-10-16.
# MODE=test (default) | record: record also keeps the live Bakta outputs (JSON, GBFF, FAA, TSV with
# SHA256SUMS) under $GIAP_ROOT/recordings/weekly/<genome>/<job>/; both modes report whether the
# live JSON's features equal the recorded one's.
# Overrides: WEEKLY_REF (default origin/main), WEEKLY_PROFILE (default farmshare,apptainer).
# One line per task goes to $GIAP_ROOT/logs/weekly.tsv. Never pushes and never commits.
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
MODE="${MODE:-test}"
case "$MODE" in test|record) ;; *) echo "MODE must be test or record, not '$MODE'" >&2; exit 2 ;; esac
REF="${WEEKLY_REF:-origin/main}"
PROFILE="${WEEKLY_PROFILE:-farmshare,apptainer}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
N="${SLURM_CPUS_PER_TASK:-16}"
SAMPLES=(k12_mg1655 nissle1917)
S="${SAMPLES[$TASK]}"
CI="$GIAP_ROOT/ci"
W="$GIAP_ROOT/ci/weekly-$(echo "${REF#origin/}" | tr '/' '-')-t$TASK"
RUN="$CI/runs/weekly_${JOB}_$TASK"
REC="$GIAP_ROOT/recordings/bakta_full/$S"
mkdir -p "$CI/runs" "$GIAP_ROOT/logs"
rm -rf "$RUN"; mkdir -p "$RUN/tmp"
(
  flock 9
  git -C "$GIAP_REPO" fetch -q origin
  [ -d "$W" ] || git -C "$GIAP_REPO" worktree add -q --detach "$W" "$REF"
  git -C "$W" checkout -q -f --detach "$REF"
) 9>"$CI/.git.lock" || { echo "cannot check out $REF in $W" >&2; exit 2; }
cd "$W" || exit 2
export TMPDIR="$RUN/tmp" NXF_ANSI_LOG=false NUMBA_CACHE_DIR="$(mktemp -d /tmp/giap-numba.XXXXXX)"
trap 'rm -rf "$NUMBA_CACHE_DIR"' EXIT
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
commit=$(git -C "$W" rev-parse --short HEAD)
echo "weekly genome=$S mode=$MODE ref=$REF commit=$commit node=$(hostname) profile=$PROFILE start=$(date -Is)"
rc=0; summary="$S:"
step() { [ "$1" -gt "$rc" ] && rc=$1; summary="$summary $2;"; echo "STEP exit=$1 $2"; }

# (a) capability gate
caps=$(python -c 'import yaml; print(",".join(yaml.safe_load(open("ci/gates.yaml"))["required_capabilities"]["farmshare_weekly"]))')
python -m pytest -p no:cacheprovider -q tests/contract/test_recording_bakta_full.py tests/contract/test_te6_oric.py \
  --require-capabilities "$caps" 2>&1 | tee "$RUN/gate.txt"
step "${PIPESTATUS[0]}" "gate $(grep -E '[0-9]+ (passed|failed|error)' "$RUN/gate.txt" | tail -1 | tr -d '=' | xargs)"

# inputs: genome (circular replicons) and the 60 % protein stand-in
zcat "$GIAP_ROOT"/data/genomes/"$S"/*.fna.gz | sed -E 's/^>(\S+)/>\1 [topology=circular]/' > "$RUN/genome.fa"
awk '/^>/ {n++} (n - 1) % 5 < 3' "$REC/bakta/$S.faa" > "$RUN/proteins.faa"
common=(--genome "$RUN/genome.fa" --proteins "$RUN/proteins.faa" --threads "$N" --memory 56.GB --profile "$PROFILE")

# (b) live full-DB Bakta + annotation track (antiSMASH live); (c) recorded full-DB JSON
for variant in live recorded_json; do
  out="$RUN/$variant"
  if [ "$variant" = live ]; then
    giap run "${common[@]}" --outdir "$out" --set annotation.bakta.db=full
  else
    giap run "${common[@]}" --outdir "$out" --annotation "$REC/bakta/$S.json"
  fi
  r=$?
  step "$r" "giap run $variant exit $r"
  python scripts/measure_runtimes.py "$out/logs/trace.tsv" --sacct "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}${SLURM_ARRAY_JOB_ID:+_$SLURM_ARRAY_TASK_ID}" --genome "$S ($variant)" \
    | tee "$RUN/runtimes_$variant.md"
done
grep -qP '\tANNOTATE:ANTISMASH \([^)]*\)\t[^\t]*\tCOMPLETED\t' "$RUN/live/logs/trace.tsv" 2>/dev/null
step $? "antismash live"

# the live Bakta outputs: compare with the recording; MODE=record keeps them
bdir=$(awk -F'\t' '$4 ~ /:BAKTA \(/ && $6 == "COMPLETED" {print $NF}' "$RUN/live/logs/trace.tsv" 2>/dev/null | tail -1)/bakta
live_json=$(ls "$bdir"/*.json 2>/dev/null | head -1)
if [ -z "$live_json" ]; then
  step 1 "no live Bakta JSON"
else
  same=$(python -c 'import json, sys
f = lambda p: [(x.get("type"), x.get("contig"), x.get("start"), x.get("stop"), x.get("strand"), x.get("product")) for x in json.load(open(p))["features"]]
print("identical" if f(sys.argv[1]) == f(sys.argv[2]) else "changed")' "$live_json" "$REC/bakta/$S.json")
  step 0 "live Bakta features vs recording: $same"
  if [ "$MODE" = record ]; then
    keep="$GIAP_ROOT/recordings/weekly/$S/$JOB"; mkdir -p "$keep"
    p=$(basename "$live_json" .json)
    for x in json gbff faa tsv; do cp "$bdir/$p.$x" "$keep/$S.$x"; done
    (cd "$keep" && sha256sum "$S".* > SHA256SUMS)
    step 0 "kept in $keep"
  fi
  # (d) TE6 on the live K-12 JSON
  if [ "$S" = k12_mg1655 ]; then
    python - "$live_json" "$RUN/genome.fa" <<'PY'
import sys
from pathlib import Path

from giap.annotate.bakta import parse_bakta_json
from giap.annotate.bakta_layers import nearest_genes
from giap.io.fasta import read_genome_fasta

seqs = {r.id: r.sequence for r in read_genome_fasta(Path(sys.argv[2]))}
ann = parse_bakta_json(Path(sys.argv[1]), genome=seqs)
genes = list(ann.model.genes)
name = {g.locus_tag: g.gene for g in genes}
ok = []
for o in (f for f in ann.model.features if f.type == "oric"):
    near = nearest_genes(o, genes, ann.model.replicons[o.interval.contig], n=20)
    mnmg = [d for loc, d in near if name[loc] == "mnmG"]
    ok.append(name[near[0][0]] in ("mnmG", "mioC") and bool(mnmg) and mnmg[0] <= 2_000)
print(f"TE6: {len(ok)} oriC calls, at the true origin: {sum(ok)}")
sys.exit(0 if any(ok) else 1)
PY
    step $? "TE6 live"
  fi
fi
printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "${JOB}_$TASK" "$commit" "$rc" "$summary" >> "$GIAP_ROOT/logs/weekly.tsv"
echo "WEEKLY_EXIT=$rc $summary"
exit "$rc"
