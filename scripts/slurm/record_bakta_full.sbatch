#!/bin/bash
#SBATCH --job-name=giap-record-bakta-full
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Bakta 1.12.2 with the full DB v6.0 on K-12 MG1655, Nissle 1917 and the 100-kb simulator genome
# (SIM100K): the recorded-fixture runs of spec §3.7 / G17. Inputs are the committed, circular-tagged
# genomes tests/data/genomes/<sample>.fna.gz (M2-10 lead decision 1). Each sample is normalised with
# `giap preflight` (replicons.tsv -> Bakta --replicons), then annotated with the exact command of
# `giap.annotate.bakta_cmd.bakta_command` in the pinned BAKTA image. Raw outputs, the command,
# versions and sha256 sums go to $OUT/<sample>/; scripts/import_recording.py copies them into
# tests/data/recorded/bakta/1.12.2-<db>/.
# Env: DBTYPE=full|light (default full; light uses $GIAP_DB_DIR/bakta/db-light),
#      SAMPLES="sim100k k12_mg1655 nissle1917", OUT (default $GIAP_ROOT/recordings/bakta_<db>_circ).
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
cd "$GIAP_REPO"
export PYTHONPATH="$GIAP_REPO/src:$GIAP_REPO/tests/sim"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
DBTYPE="${DBTYPE:-full}"
DB="$GIAP_DB_DIR/bakta/db"; [ "$DBTYPE" = light ] && DB="$GIAP_DB_DIR/bakta/db-light"
OUT="${OUT:-$GIAP_ROOT/recordings/bakta_${DBTYPE}_circ}"
DIGEST=$(grep -oP "withName: 'BAKTA' \\{ container = 'docker://oschwengers/bakta@sha256:\\K[0-9a-f]{64}" \
  workflow/conf/containers.config)
IMG="$GIAP_ROOT/containers/oschwengers-bakta@sha256-$DIGEST.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 1; }
grep -q "\"type\": \"$DBTYPE\"" "$DB/version.json" || { echo "$DBTYPE Bakta DB missing at $DB" >&2; exit 1; }
mkdir -p "$OUT"

for s in ${SAMPLES:-sim100k k12_mg1655 nissle1917}; do
  GENOME="$GIAP_REPO/tests/data/genomes/$s.fna.gz"
  [ -f "$GENOME" ] || { echo "missing $GENOME (scripts/fetch_test_genomes.py)" >&2; exit 1; }
  d="$OUT/$s"; [ -f "$d/.giap-complete" ] && { echo "$s: done"; continue; }
  rm -rf "$d/preflight" "$d/bakta"; mkdir -p "$d"
  # preflight requires evidence; a one-protein placeholder satisfies `inputs.evidence_present` and does not
  # affect genome normalisation (genome.norm.fa, replicons.tsv), which is all Bakta reads.
  printf '>placeholder_evidence\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n' > "$d/placeholder.faa"
  "$UV_PROJECT_ENVIRONMENT/bin/giap" preflight --genome "$GENOME" --proteins "$d/placeholder.faa" \
    --outdir "$d/preflight"
  w="$d/preflight/work/preflight"
  mapfile -t CMD < <("$PY" - "$w/genome.norm.fa" "$d/bakta" "$s" "$DB" "$w/replicons.tsv" <<'EOF'
import sys
from pathlib import Path
from giap.annotate.bakta_cmd import bakta_command
g, o, p, db, rep = sys.argv[1:]
for a in bakta_command(genome=Path(g), outdir=Path(o), prefix=p, db=Path(db), threads=16,
                       all_circular=False, complete=False, replicons_tsv=Path(rep), regions_gff3=None,
                       extra_args=[]):
    print(a)
EOF
)
  printf '%q ' "${CMD[@]}" > "$d/command.txt"; echo >> "$d/command.txt"
  t0=$(date +%s)
  apptainer exec --bind "$GIAP_DB_DIR,$d" "$IMG" "${CMD[@]}"
  echo "wall_s=$(( $(date +%s) - t0 )) cpus=16 node=$(hostname -s) job=${SLURM_JOB_ID:-none}" > "$d/runtime.txt"
  { apptainer exec "$IMG" bakta --version; echo "image sha256:$DIGEST"; cat "$DB/version.json"; } > "$d/versions.txt" 2>&1
  (cd "$d" && sha256sum "$GENOME" > input.sha256; find bakta preflight/work -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
  touch "$d/.giap-complete"; echo "$s: $(cat "$d/runtime.txt")"
done
du -sh "$OUT"/*
