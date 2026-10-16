#!/bin/bash
#SBATCH --job-name=giap-nightly
#SBATCH --partition=normal
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home/classes/bios/270/khoa/giap/local/logs/%x.%j.out
# Nightly calibration on the latest origin/main (GIAP_SPEC.md §12): the fast-mode calibration tests
# marked calibration_nightly (over 60 s each; GitHub CI skips them, see .github/workflows/ci.yml).
# Runs from scrontab (scripts/slurm/nightly.scrontab) or by hand: sbatch scripts/slurm/nightly.sbatch
# Writes NIGHTLY_EXIT=<pytest exit> and the commit to its log and a one-line summary to
# $GIAP_ROOT/logs/nightly.tsv (date, job, commit, exit, passed/failed counts).
set -uo pipefail
source "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}/scripts/farmshare_env.sh"
unset PYTHONHOME
W="$GIAP_ROOT/worktrees/nightly"            # dedicated detached checkout; never a feature worktree
git -C "$GIAP_REPO" fetch -q origin
[ -d "$W" ] || git -C "$GIAP_REPO" worktree add -q --detach "$W" origin/main
git -C "$W" checkout -q --detach origin/main
cd "$W" || exit 2
export TMPDIR="$(mktemp -d /tmp/giap.XXXXXX)"; trap 'rm -rf "$TMPDIR"' EXIT
export NUMBA_CACHE_DIR="$TMPDIR/numba" PYTHONPATH="$PWD/src:$PWD/tests/sim"
export PATH="$UV_PROJECT_ENVIRONMENT/bin:$PATH"
commit=$(git rev-parse --short HEAD)
echo "nightly commit=$commit node=$(hostname) start=$(date -Is)"
pytest tests/calibration -q -p no:cacheprovider -n "${SLURM_CPUS_PER_TASK:-16}" -m calibration_nightly --durations=15
rc=$?
echo "NIGHTLY_EXIT=$rc"
summary=$(grep -E '[0-9]+ (passed|failed)' "/home/classes/bios/270/khoa/giap/local/logs/giap-nightly.${SLURM_JOB_ID:-x}.out" 2>/dev/null | tail -1 | tr -d '=' | xargs)
printf '%s\t%s\t%s\t%s\t%s\n' "$(date +%F)" "${SLURM_JOB_ID:-manual}" "$commit" "$rc" "$summary" >> "$GIAP_ROOT/logs/nightly.tsv"
exit "$rc"
