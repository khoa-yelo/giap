#!/bin/bash
#SBATCH --job-name=gippa-record-tools
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Real runs of the annotation tools of GIPPA_SPEC.md §7.1.1 (M2-13..M2-18) on the three samples
# recorded by record_bakta_full.sbatch (sim100k, k12_mg1655, nissle1917), for the contract tests.
# Commands follow the tool table and the workflow/modules/local/<tool>/main.nf scripts; tools with
# a pinned image run in it (workflow/conf/containers.config), the others in their micromamba env.
# Inputs the GIPPA steps will build (genes.coords, ostir_input.csv, TU regions) come from
# scripts/dev/record_tools_inputs.py (spec definitions; TUs from the distance-model fallback).
# Output: $GIPPA_ROOT/recordings/<tool>/<sample>/ with command.txt, versions.txt, runtime.txt,
# SHA256SUMS and .giap-complete (pairs already complete are skipped). A failing pair is logged as
# "FAILED <tool> <sample>: <reason>" and the job continues.
# Usage: sbatch record_tools.sbatch [tool ...]   (default: all tools below)
# Not here: Bakta, UniOP (recorded elsewhere), codon-bias (in-process library, no tool output),
# BacTermFinder (opt-in, no env). Promotech runs genome-wide on sim100k only: the targeted windows
# writer (M2-13) does not exist yet, and genome-wide mode is ~5.8 GB per strand on E. coli.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
module load apptainer 2>/dev/null || true
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd)"
# sbatch runs a spooled copy of this script: fall back to the record-tools worktree, then the repo.
[ -f "$REPO/scripts/dev/record_tools_inputs.py" ] || REPO="$GIPPA_REPO/local/worktrees/record-tools"
[ -f "$REPO/scripts/dev/record_tools_inputs.py" ] || REPO="$GIPPA_REPO"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
INPUTS="$PY $REPO/scripts/dev/record_tools_inputs.py"
BAKTA="$GIPPA_ROOT/recordings/bakta_full"
OUT="${RECORD_OUT:-$GIPPA_ROOT/recordings}"   # RECORD_OUT, RECORD_SAMPLES: dry checks
ENVS="$MAMBA_ROOT_PREFIX/envs"
N="${SLURM_CPUS_PER_TASK:-16}"
SAMPLES=(${RECORD_SAMPLES:-sim100k k12_mg1655 nissle1917})
TOOLS=("$@"); [ "${#TOOLS[@]}" -gt 0 ] || TOOLS=(terminatornet transtermhp ostir promoter_calculator promotech genomad antismash)
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
echo "repo=$REPO cpus=$N tools=${TOOLS[*]} samples=${SAMPLES[*]}"

image() {  # pinned image path for a containers.config process name
  local ref; ref=$(grep -oP "withName: '$1' \\{ container = '(docker://)?\\K[^']+" "$REPO/workflow/conf/containers.config")
  echo "$GIPPA_ROOT/containers/$(echo "$ref" | tr '/:' '--').img"
}
digest() { grep -oP "withName: '$1' \\{ container = '[^']*@sha256:\\K[0-9a-f]{64}" "$REPO/workflow/conf/containers.config"; }
envlock() {  # env lock sha256s and pinned source commit
  local t="$1"
  (cd "$REPO/workflow/envs/locks" && sha256sum "$t".explicit.txt "$t".pip.txt)
  [ -f "$ENVS/gippa-$t/share/$t/GIPPA_COMMIT" ] && echo "commit $(cat "$ENVS/gippa-$t/share/$t/GIPPA_COMMIT")"
  return 0
}
inenv() { local e="$1"; shift; micromamba run -n "gippa-$e" "$@"; }
cmdline() { printf '%q ' "$@"; echo; }

# ---- one function per tool: runs in $d (cwd), sample inputs in $w (preflight) and $b (bakta) ----
run_terminatornet() {
  local H="$ENVS/gippa-terminatornet/share/terminatornet"
  ln -sf "$w/genome.norm.fa" genome.fa
  ln -sf "$H/model.pkl" model.pkl   # TerminatorNet loads model.pkl from the cwd
  local c=(python "$H/TerminatorNet.py" -f genome.fa -s -n "$N" -o terminatornet.tsv)
  { echo "# cwd holds model.pkl -> \$TERMINATORNET_HOME/model.pkl"; cmdline "${c[@]}"; } > command.txt
  inenv terminatornet "${c[@]}"
  { echo "terminatornet"; envlock terminatornet; } > versions.txt
}
run_transtermhp() {
  local img; img=$(image TRANSTERMHP)
  ln -sf "$w/genome.norm.fa" genome.fa
  $INPUTS coords "$b/$s.json" genes.coords
  local c=(transterm -p /usr/local/data/expterm.dat genome.fa genes.coords)
  { echo "# genes.coords: record_tools_inputs.py coords $s.json"; cmdline "${c[@]}"; } > command.txt
  apptainer exec --bind "$d,$w" "$img" "${c[@]}" > transtermhp.tt
  { echo "transtermhp 2.09"; echo "image sha256:$(digest TRANSTERMHP)"; } > versions.txt
}
run_ostir() {
  $INPUTS ostir "$b/$s.json" ostir_input.csv
  local c=(ostir -i ostir_input.csv -o ostir.csv -j "$N")
  { echo "# ostir_input.csv: record_tools_inputs.py ostir $s.json"; cmdline "${c[@]}"; } > command.txt
  inenv ostir "${c[@]}"
  { echo "ostir $(inenv ostir python -c 'import importlib.metadata as m; print(m.version("ostir"))')"; envlock ostir; } > versions.txt
}
run_promoter_calculator() {
  $INPUTS regions "$b/$s.json" regions
  mkdir -p csv
  { echo "# regions/: record_tools_inputs.py regions $s.json (TU upstream -300..+50, distance-model TUs)"
    echo "# per region, $N in parallel:"; cmdline promoter-calculator -i 'regions/<TU>.fa' -o 'csv/<TU>.csv' -j 1; } > command.txt
  ls regions | sed 's/\.fa$//' |
    xargs -P "$N" -I{} micromamba run -n gippa-promoter_calculator promoter-calculator -i regions/{}.fa -o csv/{}.csv -j 1 > pc.log 2>&1
  local nr nc; nr=$(ls regions | wc -l); nc=$(ls csv | wc -l)
  echo "regions=$nr csv=$nc"; [ "$nr" -eq "$nc" ] || { echo "only $nc of $nr regions produced a CSV" >&2; return 1; }
  { echo "promotercalculator $(inenv promoter_calculator python -c 'import importlib.metadata as m; print(m.version("promotercalculator"))')"; envlock promoter_calculator; } > versions.txt
}
run_promotech() {
  [ "$s" = sim100k ] || { echo "SKIPPED promotech $s: genome-wide mode on sim100k only (targeted windows need M2-13)"; return 97; }
  local H="$ENVS/gippa-promotech/share/promotech"
  # Promotech loads models/ relative to the cwd: link the env's models plus RF-HOT.model from the DB.
  mkdir -p models; for f in "$H"/models/*; do ln -sf "$f" models/; done
  ln -sf "$GIPPA_DB_DIR/promotech/RF-HOT.model" models/RF-HOT.model
  ln -sf "$w/genome.norm.fa" genome.fa
  local c1=(python "$H/promotech.py" -pg -m RF-HOT -f genome.fa -o parsed/)
  local c2=(python "$H/promotech.py" -g -m RF-HOT -i parsed/ -o out/)
  { echo "# genome-wide (opt-in mode); cwd holds models/ -> env models + \$GIPPA_DB_DIR/promotech/RF-HOT.model"
    cmdline "${c1[@]}"; cmdline "${c2[@]}"; } > command.txt
  OMP_NUM_THREADS="$N" inenv promotech "${c1[@]}"
  mkdir -p out; OMP_NUM_THREADS="$N" inenv promotech "${c2[@]}"
  rm -rf parsed models   # large intermediates; out/genome_predictions.csv is the recording
  { echo "promotech RF-HOT"; envlock promotech; echo "RF-HOT.zip $(cat "$GIPPA_DB_DIR/promotech/RF-HOT.zip.sha256")"; } > versions.txt
}
run_genomad() {
  local img; img=$(image GENOMAD)
  ln -sf "$w/genome.norm.fa" genome.fa
  local c=(genomad end-to-end --cleanup --splits 8 --threads "$N" genome.fa genomad "$GIPPA_DB_DIR/genomad_db")
  cmdline "${c[@]}" > command.txt
  apptainer exec --bind "$d,$w,$GIPPA_DB_DIR" "$img" "${c[@]}"
  { apptainer exec "$img" genomad --version 2>&1 | tail -1; echo "image sha256:$(digest GENOMAD)"
    echo "genomad_db $(cat "$GIPPA_DB_DIR/genomad_db/version.txt")"; } > versions.txt
}
run_antismash() {
  local img; img=$(image ANTISMASH)
  ln -sf "$b/$s.gbff" bakta.gbff
  local c=(antismash --genefinding-tool none --minimal --cpus "$N" --databases "$GIPPA_DB_DIR/antismash" --output-dir antismash bakta.gbff)
  cmdline "${c[@]}" > command.txt
  apptainer exec --bind "$d,$b,$GIPPA_DB_DIR" "$img" "${c[@]}"
  { apptainer exec "$img" antismash --version 2>&1 | tail -1; echo "image sha256:$(digest ANTISMASH)"; } > versions.txt
}

nfail=0
for t in "${TOOLS[@]}"; do
  for s in "${SAMPLES[@]}"; do
    d="$OUT/$t/$s"; w="$BAKTA/$s/preflight/work/preflight"; b="$BAKTA/$s/bakta"
    [ -f "$d/.giap-complete" ] && { echo "$t $s: done"; continue; }
    [ -f "$BAKTA/$s/.giap-complete" ] || { echo "FAILED $t $s: no complete Bakta recording"; nfail=$((nfail+1)); continue; }
    rm -rf "$d"; mkdir -p "$d"
    t0=$(date +%s)
    ( set -eo pipefail; cd "$d"; "run_$t" ) 2>&1 | tee "$d/run.log"
    rc=${PIPESTATUS[0]}
    if [ "$rc" -eq 97 ] && [ "$t" = promotech ]; then rm -rf "$d"; continue; fi
    if [ "$rc" -ne 0 ]; then
      echo "FAILED $t $s: exit $rc; last log line: $(tail -n 1 "$d/run.log" 2>/dev/null)"; nfail=$((nfail+1)); continue
    fi
    echo "wall_s=$(( $(date +%s) - t0 )) cpus=$N node=$(hostname)" > "$d/runtime.txt"
    (cd "$d" && find . -type f ! -name SHA256SUMS ! -name .giap-complete -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
    touch "$d/.giap-complete"; echo "$t $s: $(cat "$d/runtime.txt")"
  done
done
du -sh "$OUT"/*/ 2>/dev/null
echo "failures: $nfail"
