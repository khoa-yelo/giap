#!/bin/bash
#SBATCH --job-name=gippa-entrapment
#SBATCH --partition=normal
#SBATCH --qos=dev
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# §15.5 metric 8: A1's free entrapment (docs/specs/M10-07-free-entrapment.md, criterion 12).
#
#   sbatch scripts/slurm/entrapment.sbatch
#
# FarmShare's `dev` is a QOS on the `normal` partition (GIPPA_SPEC.md §3: 8 h, 8 CPU, 32 GB, one job).
# Inputs default to the spec's table and can be overridden from the environment: PEP_A1, K12_FASTA,
# CONTAMINANTS, GENOME, STRAIN_DB, RUNS_MANIFEST, METRICS_DIR, and the optional OUT_A1 (default: the
# M10-04 A1 directory, resolved through ATTEMPTS.tsv by entrapment.py), DIANN_CMD and PROVENANCE.
# An optional input set explicitly but missing exits 2 before anything runs; a default one that is
# absent is skipped and the skipped check is logged. Everything is written under $METRICS_DIR
# (default $GIPPA_ROOT/pilot/metrics): outside the checkout or in a git-ignored part of it
# (`git check-ignore`; farmshare_env.sh puts GIPPA_ROOT at $GIPPA_REPO/local). The expected run is
# minutes, dominated by the tolerant search of the K-12-only peptides of >= 10 aa.
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba" NUMBA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
cd "$GIPPA_REPO"
export PYTHONPATH="$PWD/src:$PWD"

PEP_A1="${PEP_A1:-$GIPPA_ROOT/pilot/inputs/a1.tsv.gz}"
K12_FASTA="${K12_FASTA:-$GIPPA_DATA_DIR/proteomes/UP000000625.fasta}"
CONTAMINANTS="${CONTAMINANTS:-$GIPPA_DB_DIR/contaminants/contaminants.fasta}"
GENOME="${GENOME:-tests/data/genomes/nissle1917.fna.gz}"
STRAIN_DB="${STRAIN_DB:-validation/pilot/db_a2/strain_db.faa.gz}"
RUNS_MANIFEST="${RUNS_MANIFEST:-validation/pilot/inputs/pxd078459_files.tsv}"
METRICS_DIR="${METRICS_DIR:-$GIPPA_ROOT/pilot/metrics}"
for v in OUT_A1 DIANN_CMD PROVENANCE; do  # S5: an explicit path must exist
  if [ -n "${!v+set}" ] && [ ! -e "${!v}" ]; then
    echo "precondition: $v=${!v} does not exist" >&2; exit 2
  fi
done
OUT_A1="${OUT_A1:-$GIPPA_ROOT/runs/pilot/pass1/out_a1}"
DIANN_CMD="${DIANN_CMD:-$GIPPA_ROOT/pilot/A1/command_step1.txt}"
PROVENANCE="${PROVENANCE:-$GIPPA_ROOT/pilot/A1/provenance.json}"
case "$(realpath -m "$METRICS_DIR")/" in
  "$(realpath "$GIPPA_REPO")"/*)
    git -C "$GIPPA_REPO" check-ignore -q "$METRICS_DIR" || {
      echo "precondition: METRICS_DIR=$METRICS_DIR is inside the checkout and not git-ignored" >&2
      exit 2
    } ;;
esac

opt=()
if [ -d "$OUT_A1" ]; then opt+=(--out-a1 "$OUT_A1")
else echo "cross-checks 1 and 2 skipped: no A1 run directory $OUT_A1"; fi
if [ -f "$DIANN_CMD" ]; then opt+=(--diann-cmd "$DIANN_CMD")
else echo "cross-check 4 (DIA-NN settings) skipped: no $DIANN_CMD"; fi
if [ -f "$PROVENANCE" ]; then opt+=(--provenance "$PROVENANCE")
else echo "provenance check skipped: no $PROVENANCE"; fi

echo "node=$(hostname) head=$(git rev-parse --short HEAD) out=$METRICS_DIR optional=${opt[*]:-none}"
rc=0
"$UV_PROJECT_ENVIRONMENT/bin/python" validation/pilot/entrapment.py compute \
  --pep-a1 "$PEP_A1" --k12-fasta "$K12_FASTA" --contaminants "$CONTAMINANTS" \
  --genome "$GENOME" --strain-db "$STRAIN_DB" --runs-manifest "$RUNS_MANIFEST" \
  --out "$METRICS_DIR" "${opt[@]}" || rc=$?
echo "ENTRAPMENT_EXIT=$rc"
exit $rc
