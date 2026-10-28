#!/bin/bash
# replant.sh WORKTREE FORKPOINT PR — rebase a branch's own commits onto origin/main, verify, retarget the PR, push.
set -euo pipefail
wt=$1; fork=$2; pr=$3
cd "$wt"
source /home/classes/bios/270/khoa/giap/scripts/farmshare_env.sh
export GIPPA_REPO=/home/classes/bios/270/khoa/giap PYTHONPATH=$PWD/src:$PWD/tests/sim
git fetch -q origin
branch=$(git branch --show-current)
git rebase --onto origin/main "$fork" "$branch"
uv run --no-sync pytest tests -q -n 4 --ignore=tests/calibration -p no:cacheprovider > $GIPPA_LOG_DIR/replant_$pr.log 2>&1 \
  || { tail -20 $GIPPA_LOG_DIR/replant_$pr.log; echo "TESTS FAILED"; exit 1; }
tail -1 $GIPPA_LOG_DIR/replant_$pr.log
uv run --no-sync ruff check -q .
uv run --no-sync ruff format --check -q .
uv run --no-sync mypy src/gippa | tail -1
python3 scripts/ci/check_roles.py origin/main
python3 scripts/ci/check_protected.py origin/main
gh pr edit "$pr" --base main > /dev/null   # retarget first, so CI diffs against main
git push -q --force-with-lease
echo "replanted $branch (#$pr) onto main: $(git diff --shortstat origin/main...HEAD)"
