#!/bin/bash
#SBATCH --job-name=gippa-setup
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# `gippa setup` into $GIPPA_DB_DIR (spec §7.9.1, §3.9; M2-08). Adopts what is already there (prefetch),
# downloads the rest from the official sources only (R6), writes $GIPPA_DB_DIR/db_manifest.json.
# Extra arguments replace the defaults, e.g. `sbatch scripts/slurm/gippa_setup.sbatch --check`.
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIPPA_REPO"
if [ $# -eq 0 ]; then set -- --bakta full --genomad --promotech-models --go --kegg --cog --contaminants; fi
uv run --no-sync gippa setup --db-dir "$GIPPA_DB_DIR" "$@"
