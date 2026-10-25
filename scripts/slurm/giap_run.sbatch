#!/bin/bash
#SBATCH --job-name=giap-run
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Nextflow head job (spec §3.8): Nextflow submits every task as its own giap-* Slurm job.
# Usage: sbatch scripts/slurm/giap_run.sbatch <giap run arguments, e.g. --peptides P --genome G --outdir OUT>
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIAP_REPO"
# pin the code to this checkout too: the shared venv's editable giap is the main checkout
export PYTHONPATH="$GIAP_REPO/src${PYTHONPATH:+:$PYTHONPATH}"
uv run --no-sync giap run --profile farmshare,apptainer "$@"
