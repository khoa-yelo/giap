#!/bin/bash
#SBATCH --job-name=giap-record-bakta-full
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Bakta 1.12.2 with the full DB v6.0 on K-12 MG1655, Nissle 1917 and the 100-kb simulator genome
# (SIM100K): the recorded-fixture runs of spec §3.7 / G17. Each sample is normalised with `giap
# preflight`, then annotated with the exact command of `giap.annotate.bakta_cmd.bakta_command` in the
# pinned BAKTA image. Raw outputs, the command, versions and sha256 sums go to
# $GIAP_ROOT/recordings/bakta_full/<sample>/ for M2-10 to import into tests/data/recorded/.
# Needs the genomes from fetch_resources.sbatch (submit with --dependency=afterok:<id>).
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
cd "$GIAP_REPO"
export PYTHONPATH="$GIAP_REPO/src:$GIAP_REPO/tests/sim"
PY="$UV_PROJECT_ENVIRONMENT/bin/python"
DATA="${GIAP_DATA_DIR:-$GIAP_ROOT/data}"
OUT="$GIAP_ROOT/recordings/bakta_full"
DIGEST=$(grep -oP "withName: 'BAKTA' \\{ container = 'docker://oschwengers/bakta@sha256:\\K[0-9a-f]{64}" \
  workflow/conf/containers.config)
IMG="$GIAP_ROOT/containers/oschwengers-bakta@sha256-$DIGEST.img"
[ -f "$IMG" ] || { echo "missing $IMG" >&2; exit 1; }
[ -f "$GIAP_DB_DIR/bakta/db/version.json" ] || { echo "full Bakta DB missing" >&2; exit 1; }
mkdir -p "$OUT"

# Inputs: two real genomes and SIM100K (written by the simulator with its fixed seed).
declare -A GENOME
GENOME[k12_mg1655]=$(ls "$DATA"/genomes/k12_mg1655/*_genomic.fna.gz)
GENOME[nissle1917]=$(ls "$DATA"/genomes/nissle1917/*_genomic.fna.gz)
if [ ! -f "$OUT/sim100k/input/genome.fasta" ]; then
  "$PY" -c "from giapsim import SIM100K, SIM100K_SEED, simulate; simulate(SIM100K, SIM100K_SEED, '$OUT/sim100k/input')"
fi
GENOME[sim100k]="$OUT/sim100k/input/genome.fasta"

for s in sim100k k12_mg1655 nissle1917; do
  d="$OUT/$s"; [ -f "$d/.giap-complete" ] && { echo "$s: done"; continue; }
  rm -rf "$d/preflight" "$d/bakta"; mkdir -p "$d"
  # preflight requires evidence; a one-protein placeholder satisfies `inputs.evidence_present` and does not
  # affect genome normalisation (genome.norm.fa, replicons.tsv), which is all Bakta reads.
  printf '>placeholder_evidence\nMSTNPKPQRKTKRNTNRRPQDVKFPGG\n' > "$d/placeholder.faa"
  "$UV_PROJECT_ENVIRONMENT/bin/giap" preflight --genome "${GENOME[$s]}" --proteins "$d/placeholder.faa" \
    --outdir "$d/preflight"
  w="$d/preflight/work/preflight"
  mapfile -t CMD < <("$PY" - "$w/genome.norm.fa" "$d/bakta" "$s" "$GIAP_DB_DIR/bakta/db" "$w/replicons.tsv" <<'EOF'
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
  echo "wall_s=$(( $(date +%s) - t0 ))" > "$d/runtime.txt"
  { apptainer exec "$IMG" bakta --version; echo "image sha256:$DIGEST"; cat "$GIAP_DB_DIR/bakta/db/version.json"; } > "$d/versions.txt" 2>&1
  (cd "$d" && find bakta preflight/work -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
  touch "$d/.giap-complete"; echo "$s: $(cat "$d/runtime.txt")"
done
du -sh "$OUT"/*
