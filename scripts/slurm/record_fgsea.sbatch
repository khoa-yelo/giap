#!/bin/bash
#SBATCH --job-name=giap-record-fgsea
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# M8-17a §6 (GIAP_SPEC.md §11.4): real fgsea 1.36.2 recordings of workflow/bin/giap_fgsea.R in the
# giap-fgsea env (sbatch scripts/slurm/build_envs.sbatch fgsea first), for
# tests/data/recorded/fgsea/1.36.2/: sim_main, sim_retest, ties, unbalanced, tiny_p, log2err_na
# (scripts/dev/record_fgsea.py; every input synthetic, every file < 1 MB).
# Output: $GIAP_ROOT/recordings/fgsea/1.36.2/<name>/; copy into the worktree after review.
# Usage: sbatch record_fgsea.sbatch <worktree> [name ...]
set -euo pipefail
WT="${1:?worktree dir}"; shift || true
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
D="$GIAP_ROOT/recordings/fgsea/1.36.2"
mkdir -p "$D"
echo "worktree=$WT head=$(git -C "$WT" rev-parse --short HEAD) node=$(hostname)"
SECONDS=0
PYTHONPATH="$WT/src:$WT/tests/sim:$WT/tests" "$UV_PROJECT_ENVIRONMENT/bin/python" \
  "$WT/scripts/dev/record_fgsea.py" "$D" "$@"
echo "wall time: ${SECONDS}s"
du -sh "$D"/*
