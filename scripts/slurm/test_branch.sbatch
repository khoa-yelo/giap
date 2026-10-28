#!/bin/bash
#SBATCH --job-name=gippa-test
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Run pytest for one worktree on a batch node instead of the shared 8-CPU session allocation.
# Usage: sbatch scripts/slurm/test_branch.sbatch <worktree-dir> [pytest args...]
#   default pytest args: tests -q -n 16 --durations=25 -m "not slow and not calibration_full and not calibration_nightly and not e2e"
#   No coverage by default (CI measures it, and it slows every test); GIPPA_TEST_COV=1 adds --cov=gippa.
# The last log line is "PYTEST_EXIT=<code>"; poll with `squeue -j <id>` and read the log. Read-only on the
# worktree apart from pytest caches/coverage files. R14: 16 CPUs; keep at most 4 gippa-test jobs at a time.
set -uo pipefail
WT="${1:?worktree dir}"; shift || true
source "${GIPPA_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
export TMPDIR="$(mktemp -d /tmp/gippa.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba"   # node-local; avoids cache races with the session node
cd "$WT"
export PYTHONPATH="$PWD/src:$PWD/tests/sim"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
echo "worktree=$WT branch=$(git rev-parse --abbrev-ref HEAD) head=$(git rev-parse --short HEAD) node=$(hostname)"
# Warm the gippa-wf JVM before pytest (M9-13 decision 5): on a node that has not run it, the first
# `nextflow -version` reads Java and the Nextflow jars cold from NFS and took 29 s on rye-02
# (job 1869186; 1.4 s once warm), so tests/conftest.py's 30 s capability probe timed out and every
# needs_nextflow test skipped (jobs 1850601, 1860124, 1869158; with -n 16 the workers cold-start together).
for probe in "nextflow -version" "nf-test version"; do
  s=$SECONDS; $probe >/dev/null 2>&1; echo "warm-up: $probe rc=$? $((SECONDS - s))s"
done
if [ "$#" -eq 0 ]; then
  set -- tests -q -n "${SLURM_CPUS_PER_TASK:-16}" --durations=25 -m "not slow and not calibration_full and not calibration_nightly and not e2e"
  [ "${GIPPA_TEST_COV:-0}" = 1 ] && set -- "$@" --cov=gippa
fi
python -m pytest -p no:cacheprovider "$@"
rc=$?
echo "PYTEST_EXIT=$rc"
exit $rc
