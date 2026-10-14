#!/bin/bash
#SBATCH --job-name=giap-build-envs
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Builds the tool envs (workflow/envs/*.yml) and pulls every pinned image by digest into
# $NXF_APPTAINER_CACHEDIR (spec §3.8, §3.9, M2-01).  Usage: sbatch build_envs.sbatch [tool ...]
set -euo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIAP_REPO"
if [ "$#" -eq 0 ]; then set -- all; fi
SECONDS=0
bash scripts/build_envs.sh "$@"

module load apptainer 2>/dev/null || true
mkdir -p "$NXF_APPTAINER_CACHEDIR"
grep -o "docker://[^'\"]*@sha256:[0-9a-f]\{64\}" workflow/conf/containers.config | sort -u |
while read -r ref; do
  # Nextflow's cache naming: '/', ':' and '@' become '-'.
  sif="$NXF_APPTAINER_CACHEDIR/$(echo "${ref#docker://}" | tr '/:@' '---').img"
  if [ -s "$sif" ]; then echo "present: $sif"; else apptainer pull "$sif" "$ref"; fi
done

echo "wall time: ${SECONDS}s"
du -sh "$MAMBA_ROOT_PREFIX"/envs/giap-* "$NXF_APPTAINER_CACHEDIR" 2>/dev/null || true
