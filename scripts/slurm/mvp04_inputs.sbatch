#!/bin/bash
#SBATCH --job-name=gippa-mvp04-inputs
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# MVP-04 input step (docs/specs/MVP-04-ecoli-evidence-run.md "The run", rule 2): gippasim
# (scenarios.k12_evidence) on the committed K-12 genome and recorded full-DB Bakta JSON, with
# RegulonDB 14.5.0 TUs when present (gap rule otherwise). Two calls with the same seed: the inputs
# in $GIPPA_ROOT/mvp04/<seed>/ and the TE7 copy in $GIPPA_ROOT/mvp04/<seed>_te7/ (gate 9).
# Usage: sbatch scripts/slurm/mvp04_inputs.sbatch [SEED]   (default K12_EVIDENCE_SEED = 20261023)
# MVP04_OUT_ROOT replaces $GIPPA_ROOT/mvp04 (smoke tests into a scratch directory).
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIPPA_REPO"
export PYTHONPATH="$GIPPA_REPO/src:$GIPPA_REPO/tests/sim"
SEED="${1:-${SEED:-20261023}}"   # K12_EVIDENCE_SEED (gippasim.scenarios)
INPUTS="$GIPPA_ROOT/mvp04/$SEED"
if [[ -n "${MVP04_OUT_ROOT:-}" ]]; then INPUTS="$MVP04_OUT_ROOT/$SEED"; fi
TU_SET="$GIPPA_DATA_DIR/regulondb/14.5.0/TUSet.txt"
TU_ARGS=()
if [[ -f "$TU_SET" ]]; then TU_ARGS=(--tu-set "$TU_SET"); else echo "WARNING: $TU_SET absent: gap rule TUs"; fi
for OUT in "$INPUTS" "${INPUTS}_te7"; do
  rm -rf "${OUT:?}"
  "$UV_PROJECT_ENVIRONMENT/bin/python" -m gippasim.real_genome \
    --genome tests/data/genomes/k12_mg1655.fna.gz \
    --annotation tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz \
    --seed "$SEED" --outdir "$OUT" "${TU_ARGS[@]}"
done
echo "inputs: $INPUTS  te7: ${INPUTS}_te7"
