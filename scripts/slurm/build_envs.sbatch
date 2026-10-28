#!/bin/bash
#SBATCH --job-name=gippa-build-envs
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Builds the tool envs (workflow/envs/*.yml) and pulls every pinned image by digest into
# $NXF_APPTAINER_CACHEDIR (spec §3.8, §3.9, M2-01).  Usage: sbatch build_envs.sbatch [tool ...]
set -euo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
cd "$GIPPA_REPO"
if [ "$#" -eq 0 ]; then set -- all; fi
SECONDS=0
status=0
bash scripts/build_envs.sh "$@" || status=$?

module load apptainer 2>/dev/null || true
mkdir -p "$NXF_APPTAINER_CACHEDIR"
# Pins are bare refs (ADR-0059); a legacy docker:// prefix is accepted. `|| true`: no pins is no error.
refs=$(sed -nE "s|.*container = '(docker://)?([^'\"]*@sha256:[0-9a-f]{64})'.*|\2|p" \
  workflow/conf/containers.config | sort -u) || true
for ref in $refs; do
  # Nextflow 26.10 SingularityCache.simpleName: ':' and '/' become '-', '@' kept.
  sif="$NXF_APPTAINER_CACHEDIR/$(echo "$ref" | tr '/:' '--').img"
  if [ -s "$sif" ]; then echo "present: $sif"; else apptainer pull "$sif" "docker://$ref"; fi
done

echo "wall time: ${SECONDS}s"
du -sh "$MAMBA_ROOT_PREFIX"/envs/gippa-* "$NXF_APPTAINER_CACHEDIR" 2>/dev/null || true
exit "$status"
