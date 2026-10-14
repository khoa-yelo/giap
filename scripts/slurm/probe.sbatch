#!/bin/bash
#SBATCH --job-name=giap-probe
#SBATCH --partition=normal
#SBATCH --qos=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Proves the batch path (spec §13.1 step 4): environment, tools and network from a compute node.
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIAP_REPO"
uv run --no-sync python scripts/probe_env.py --out "$GIAP_ROOT/logs/probe.compute.$SLURM_JOB_ID.json"
