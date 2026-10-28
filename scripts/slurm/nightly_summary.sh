#!/bin/bash
#SBATCH --job-name=gippa-nightly-summary
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# The night's summary job (M2-21 lead decisions 14, 21), submitted by scripts/slurm/nightly_launch.sh
# with --dependency=afterany:<array>: posts the array's nightly.tsv rows with
# scripts/ci/post_summary.py --final, so a task killed before it could report (SIGKILL, OOM) still
# shows up (its sacct state); a task-2 shard without a record fails the calibration_nightly row
# (scripts/ci/nightly_shards.py close). Usage: nightly_summary.sh ARRAY_JOB [post_summary.py options]
# Slurm runs a spooled copy of this file, so the repository comes from GIPPA_REPO (exported by the
# launcher), never from this script's own path.
set -uo pipefail
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
job="${1:?usage: nightly_summary.sh ARRAY_JOB [post_summary.py options]}"
shift
# sharded task 2: shards that left no record fail the night's calibration_nightly row
python "$GIPPA_REPO/scripts/ci/nightly_shards.py" close --dir "$GIPPA_ROOT/ci/runs/${job}_2.shards" \
  --job "$job" --tsv "$GIPPA_ROOT/logs/nightly.tsv" || echo "nightly_shards.py close failed" >&2
exec python "$GIPPA_REPO/scripts/ci/post_summary.py" --job "$job" --tsv "$GIPPA_ROOT/logs/nightly.tsv" \
  --expect "${NIGHTLY_EXPECT:-4}" --final "$@"
