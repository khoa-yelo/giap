#!/bin/bash
#SBATCH --job-name=gippa-weekly
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --array=0-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%A_%a.out
# Weekly end-to-end runs on origin/main (GIPPA_SPEC.md §12; M2-21), one array task per genome:
#   0 = k12_mg1655 (U00096.3), 1 = nissle1917 (CP171242-CP171245), from $GIPPA_ROOT/data/genomes
#   (every replicon marked [topology=circular]: all are closed); 2 = k12_distance_model (14s2).
# R14: 3 x 16 = 48 CPUs; with the Sunday nightly (65) that is 113 of 128 during the overlap.
# Per genome: (a) the farmshare_weekly capability gate (pytest --require-capabilities, the Bakta
# full-DB contract tests; a missing DB fails); (b) real `gippa run` with live full-DB Bakta and every
# annotation tool, including the antiSMASH live run (its trace row must be COMPLETED); (c) `gippa
# run` with the recorded full-DB JSON (annotation.precomputed.bakta_json =
# $GIPPA_ROOT/recordings/bakta_full/<genome>/bakta/<genome>.json);
# (d) K-12 only: TE6 (§11.6) on the live Bakta JSON: an oriC call whose nearest gene is mnmG or
# mioC with mnmG within 2 kb; (e) K-12 only: TE8 (§11.6, M3-22): validation/fetch_uniprot_k12.py,
# then validation/proxies/te_proxies.py te8 -> $GIPPA_ROOT/ci/validation/<date>_te8_<job>.json;
# (f) Nissle only: TE2 (§11.6, M7-14): validation/proxies/te_proxies.py te2 on the committed Nissle
# genome and recorded full-DB JSON -> $GIPPA_ROOT/ci/validation/<date>_te2_<job>.json.
# (g) task 2 only (M11-14 14s2, AC 10): K-12 `gippa run` with the recorded full-DB JSON and
# `annotation.uniop.enabled=false` (the distance-model operon fallback), row label
# k12_distance_model (never k12_mg1655, so readiness keeps reading TE8 from task 0). The row has
# the run's exit, operon_source from its run_manifest.json and the TU pair precision / recall vs
# RegulonDB (scripts/tu_agreement.py) next to the UniOP K-12 arm's (task 0's recorded_json of this
# job, else the newest earlier one; "not available" if none). Reported, not gated: task 2 exits 0.
# Measured runtimes: scripts/measure_runtimes.py on each trace.
# Evidence stand-in (M2 annotation track only; MAPPING is M3): `gippa run` requires --proteins or
# --peptides, so the evidence is a deterministic 60 % subset of the recorded Bakta proteins
# (records 0,1,2 of every 5 of bakta/<genome>.faa), as in the lead's hand run of 2026-10-16.
# MODE=test (default) | record: record also keeps the live Bakta outputs (JSON, GBFF, FAA, TSV with
# SHA256SUMS) under $GIPPA_ROOT/recordings/weekly/<genome>/<job>/; both modes report whether the
# live JSON's features equal the recorded one's.
# Overrides: WEEKLY_REF (default origin/main), WEEKLY_PROFILE (default farmshare_local,apptainer:
# every Nextflow task inside this 16-CPU allocation; a list with the Slurm-executor profile
# `farmshare` is refused with exit 2, R14 and lead decisions 12, 22).
# Every task writes exactly one line to $GIPPA_ROOT/logs/weekly.tsv and prints WEEKLY_EXIT=<rc>, also
# on an early exit or SIGTERM (exit 143; EXIT/TERM traps, decision 27). Never pushes or commits.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
MODE="${MODE:-test}"
REF="${WEEKLY_REF:-origin/main}"
PROFILE="${WEEKLY_PROFILE:-farmshare_local,apptainer}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"
N="${SLURM_CPUS_PER_TASK:-16}"
SAMPLES=(k12_mg1655 nissle1917 k12_mg1655)
LABELS=(k12_mg1655 nissle1917 k12_distance_model)
S="${SAMPLES[$TASK]:-none}"
L="${LABELS[$TASK]:-none}"
CI="$GIPPA_ROOT/ci"
W="$GIPPA_ROOT/ci/weekly-$(echo "${REF#origin/}" | tr '/' '-')-t$TASK"
RUN="$CI/runs/weekly_${JOB}_$TASK"
REC="$GIPPA_ROOT/recordings/bakta_full/$S"
commit=$(git -C "$GIPPA_REPO" rev-parse -q --short --verify "$REF^0" 2>/dev/null || echo unknown)
rc=0; summary="$L:"; stage="setup"; reported=""; numba_tmp=""
mkdir -p "$CI/runs" "$GIPPA_ROOT/logs"

report() {  # report <exit> <summary>: this task's one weekly.tsv line
  reported=1
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "${JOB}_$TASK" "$commit" "$1" "$2" >> "$GIPPA_ROOT/logs/weekly.tsv"
  echo "WEEKLY_EXIT=$1 $2"
}
on_exit() {
  local r=$?
  [ "${BASHPID:-$$}" = "$$" ] || return 0   # a forked child hit by the group SIGTERM before its exec (bash < 4 has no BASHPID: old behaviour)
  [ -n "$reported" ] || report "$(( r ? r : 1 ))" "$summary ended during $stage (exit $r)"
  [ -z "$numba_tmp" ] || rm -rf "$numba_tmp"   # only this task's cache, never the shared one
}
trap on_exit EXIT
trap '[ "${BASHPID:-$$}" = "$$" ] || exit 143; [ -n "$reported" ] || report 143 "$summary killed by SIGTERM (time limit?) during $stage"; exit 143' TERM

case "$MODE" in test|record) ;; *) echo "MODE must be test or record, not '$MODE'" >&2; exit 2 ;; esac
case ",$PROFILE," in
  *,farmshare,*)   # decision 22: Slurm-executor child jobs would break the R14 CPU budget
    echo "WEEKLY_PROFILE=$PROFILE: the farmshare (Slurm executor) profile is refused (R14); use farmshare_local" >&2
    report 2 "$L: WEEKLY_PROFILE=$PROFILE refused (R14: use farmshare_local)"; exit 2 ;;
esac
rm -rf "$RUN"; mkdir -p "$RUN/tmp"
stage="checkout of $REF"
(
  flock 9
  git -C "$GIPPA_REPO" fetch -q origin
  [ -d "$W" ] || git -C "$GIPPA_REPO" worktree add -q --detach "$W" "$REF"
  git -C "$W" checkout -q -f --detach "$REF"
) 9>"$CI/.git.lock" || { echo "cannot check out $REF in $W" >&2; exit 2; }
cd "$W" || exit 2
numba_tmp="$(mktemp -d /tmp/gippa-numba.XXXXXX)"
export TMPDIR="$RUN/tmp" NXF_ANSI_LOG=false NUMBA_CACHE_DIR="$numba_tmp"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
commit=$(git -C "$W" rev-parse --short HEAD)
stage="the run"
echo "weekly genome=$S arm=$L mode=$MODE ref=$REF commit=$commit node=$(hostname) profile=$PROFILE start=$(date -Is)"
step() { [ "$1" -gt "$rc" ] && rc=$1; summary="$summary $2;"; echo "STEP exit=$1 $2"; }

# inputs: genome (circular replicons) and the 60 % protein stand-in
zcat "$GIPPA_ROOT"/data/genomes/"$S"/*.fna.gz | sed -E 's/^>(\S+)/>\1 [topology=circular]/' > "$RUN/genome.fa"
awk '/^>/ {n++} (n - 1) % 5 < 3' "$REC/bakta/$S.faa" > "$RUN/proteins.faa"
common=(--genome "$RUN/genome.fa" --proteins "$RUN/proteins.faa" --threads "$N" --memory 56.GB --profile "$PROFILE")

# (g) task 2: the distance-model fallback arm, reported, not gated
if [ "$TASK" = 2 ]; then
  stage="the distance-model fallback run"
  out="$RUN/distance_model"
  gippa run "${common[@]}" --outdir "$out" --set "annotation.precomputed.bakta_json=$REC/bakta/$S.json" \
    --set annotation.uniop.enabled=false
  r=$?
  src=$(python -c 'import json, sys
try: print(json.load(open(sys.argv[1])).get("operon_source") or "missing")
except (OSError, ValueError, AttributeError): print("missing")' "$out/run_manifest.json")
  uni="$CI/runs/weekly_${JOB}_0/recorded_json"
  if [ ! -f "$uni/tables/operons.tsv" ]; then
    last=$(ls -td "$CI"/runs/weekly_*_0/recorded_json/tables/operons.tsv 2>/dev/null | head -1)
    [ -z "$last" ] || uni="${last%/tables/operons.tsv}"
  fi
  # a crash reports "not available", never its last stderr line (stderr stays in the job log)
  agree=$(python scripts/tu_agreement.py --regulondb-dir "$GIPPA_DATA_DIR/regulondb/14.5.0" --arm "fallback=$out" \
    --arm "uniop ($(basename "$(dirname "$uni")"))=$uni" | tail -1 | tr '\t' ' '; exit "${PIPESTATUS[0]}") || agree=""
  report "$r" "$summary gippa run distance_model exit $r; operon_source=$src; ${agree:-TU pairs vs regulondb: not available; uniop: not available};"
  exit 0
fi

# (a) capability gate
caps=$(python -c 'import yaml; print(",".join(yaml.safe_load(open("ci/gates.yaml"))["required_capabilities"]["farmshare_weekly"]))')
python -m pytest -p no:cacheprovider -q tests/contract/test_recording_bakta_full.py tests/contract/test_te6_oric.py \
  --require-capabilities "$caps" 2>&1 | tee "$RUN/gate.txt"
step "${PIPESTATUS[0]}" "gate $(grep -E '[0-9]+ (passed|failed|error)' "$RUN/gate.txt" | tail -1 | tr -d '=' | xargs)"

# (b) live full-DB Bakta + annotation track (antiSMASH live); (c) recorded full-DB JSON
for variant in live recorded_json; do
  out="$RUN/$variant"
  if [ "$variant" = live ]; then
    gippa run "${common[@]}" --outdir "$out" --set annotation.bakta.db=full
  else
    gippa run "${common[@]}" --outdir "$out" --set "annotation.precomputed.bakta_json=$REC/bakta/$S.json"
  fi
  r=$?
  step "$r" "gippa run $variant exit $r"
  python scripts/measure_runtimes.py "$out/logs/trace.tsv" --sacct "${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}${SLURM_ARRAY_JOB_ID:+_$SLURM_ARRAY_TASK_ID}" --genome "$S ($variant)" \
    | tee "$RUN/runtimes_$variant.md"
done
grep -qP '\tANNOTATE:ANTISMASH \([^)]*\)\t[^\t]*\tCOMPLETED\t' "$RUN/live/logs/trace.tsv" 2>/dev/null
step $? "antismash live"

# the live Bakta outputs: compare with the recording; MODE=record keeps them
# the workdir column by name: trace fields end in hostname,container (workflow/nextflow.config)
bdir=$(awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "workdir") w = i; next}
  w && $4 ~ /:BAKTA \(/ && $6 == "COMPLETED" {print $w}' "$RUN/live/logs/trace.tsv" 2>/dev/null | tail -1)/bakta
live_json=$(ls "$bdir"/*.json 2>/dev/null | head -1)
if [ -z "$live_json" ]; then
  step 1 "no live Bakta JSON"
else
  same=$(python -c 'import json, sys
f = lambda p: [(x.get("type"), x.get("contig"), x.get("start"), x.get("stop"), x.get("strand"), x.get("product")) for x in json.load(open(p))["features"]]
print("identical" if f(sys.argv[1]) == f(sys.argv[2]) else "changed")' "$live_json" "$REC/bakta/$S.json")
  step 0 "live Bakta features vs recording: $same"
  if [ "$MODE" = record ]; then
    keep="$GIPPA_ROOT/recordings/weekly/$S/$JOB"; mkdir -p "$keep"
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

from gippa.annotate.bakta import parse_bakta_json
from gippa.annotate.bakta_layers import nearest_genes
from gippa.io.fasta import read_genome_fasta

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
# (e) K-12 task only: TE8 (§11.6, M3-22) on the live UniProt K-12 proteome (rest.uniprot.org) with
# the committed genomes and recorded full-DB Bakta JSONs; result JSON kept as gate-M3 evidence
if [ "$S" = k12_mg1655 ]; then
  stage="TE8"
  te8="$CI/validation/$(date +%F)_te8_${JOB}.json"; mkdir -p "$CI/validation"
  D=tests/data; B=$D/recorded/bakta/1.12.2-full
  python validation/fetch_uniprot_k12.py \
    && python validation/proxies/te_proxies.py te8 "$D/genomes/k12_mg1655.fna.gz" \
         "$B/k12/k12_mg1655.json.gz" "$D/genomes/nissle1917.fna.gz" "$B/nissle/nissle1917.json.gz" \
         "$GIPPA_DATA_DIR/uniprot/UP000000625.fasta" --out "$te8" 2>&1 | tee "$RUN/te8.log" | tail -1
  step "${PIPESTATUS[0]}" "TE8 $te8"
fi
# (f) Nissle task only: TE2 (§11.6, M7-14) on the committed genome and the recorded full-DB Bakta
# JSON (no download); a violation turns the weekly line red. Result JSON kept as gate-M7 evidence
if [ "$S" = nissle1917 ]; then
  stage="TE2"
  te2="$CI/validation/$(date +%F)_te2_${JOB}.json"; mkdir -p "$CI/validation"
  D=tests/data; B=$D/recorded/bakta/1.12.2-full
  python validation/proxies/te_proxies.py te2 "$D/genomes/nissle1917.fna.gz" \
    "$B/nissle/nissle1917.json.gz" --out "$te2" 2>&1 | tee "$RUN/te2.log" | tail -1
  step "${PIPESTATUS[0]}" "TE2 $te2"
fi
report "$rc" "$summary"
exit "$rc"
