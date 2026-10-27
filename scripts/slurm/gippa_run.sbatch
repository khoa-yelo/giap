#!/bin/bash
#SBATCH --job-name=gippa-run
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Nextflow head job (spec §3.8): Nextflow submits every task as its own gippa-* Slurm job.
# Usage: sbatch scripts/slurm/gippa_run.sbatch <gippa run arguments, e.g. --peptides P --genome G --outdir OUT>
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIPPA_REPO"
# pin the code to this checkout too: the shared venv's editable gippa is the main checkout
export PYTHONPATH="$GIPPA_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
uv run --no-sync gippa run --profile farmshare,apptainer "$@"
