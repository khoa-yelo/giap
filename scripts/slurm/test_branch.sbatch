#!/bin/bash
#SBATCH --job-name=giap-test
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Run pytest for one worktree on a batch node instead of the shared 8-CPU session allocation.
# Usage: sbatch scripts/slurm/test_branch.sbatch <worktree-dir> [pytest args...]
#   default pytest args: tests -q -n 16 --cov=giap -m "not slow and not calibration_full and not calibration_nightly and not e2e"
# The last log line is "PYTEST_EXIT=<code>"; poll with `squeue -j <id>` and read the log. Read-only on the
# worktree apart from pytest caches/coverage files. R14: 16 CPUs; keep at most 4 giap-test jobs at a time.
set -uo pipefail
WT="${1:?worktree dir}"; shift || true
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba"   # node-local; avoids cache races with the session node
cd "$WT"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
echo "worktree=$WT branch=$(git rev-parse --abbrev-ref HEAD) head=$(git rev-parse --short HEAD) node=$(hostname)"
if [ "$#" -eq 0 ]; then
  set -- tests -q -n "${SLURM_CPUS_PER_TASK:-16}" --cov=giap -m "not slow and not calibration_full and not calibration_nightly and not e2e"
fi
python -m pytest -p no:cacheprovider "$@"
rc=$?
echo "PYTEST_EXIT=$rc"
exit $rc
