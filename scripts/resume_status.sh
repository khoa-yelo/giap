#!/bin/bash
# What was in flight when the last session stopped (GIPPA_SPEC.md R8, §14.10).
# Run first in every new session: lists each worktree's issue, the last feature-loop step that was
# committed, uncommitted work left behind, the next step, open PRs, GIPPA Slurm jobs and time left.
set -uo pipefail
source "$(dirname "$0")/farmshare_env.sh"
cd "$GIPPA_REPO" || exit 1
git fetch -q origin --prune 2>/dev/null || echo "(git fetch failed; showing local state)"

next_step() {  # role of the last commit -> next step of the feature loop (§14.3)
  case "$1" in
    spec-writer) echo "test-author: write failing tests" ;;
    test-author) echo "implementer: make the tests pass (or lead: align tests if a contradiction was reported)" ;;
    implementer) echo "verify suite, push, open/update PR, code-reviewer (or merge if review already done)" ;;
    lead)        echo "read the last commit subject; usually push/PR or next loop step" ;;
    *)           echo "inspect" ;;
  esac
}

echo "== GIPPA resume status ($(date '+%F %T')) =="
if [ -n "${SLURM_JOB_ID:-}" ]; then
  echo "allocation: job $SLURM_JOB_ID, time left $(squeue -h -j "$SLURM_JOB_ID" -o %L 2>/dev/null)"
else
  echo "allocation: NONE — start one before working (R14)"
fi
echo "main: $(git log --oneline -1 origin/main 2>/dev/null)"
echo

echo "== Worktrees (local/worktrees/*) =="
git worktree list --porcelain | awk '/^worktree /{print $2}' | grep '/local/worktrees/' | while read -r wt; do
  branch=$(git -C "$wt" branch --show-current)
  subject=$(git -C "$wt" log -1 --format=%s)
  role=$(git -C "$wt" log -1 --format=%B | sed -n 's/^Agent-Role: *//p' | head -1)
  ahead=$(git -C "$wt" rev-list --count origin/main..HEAD 2>/dev/null || echo "?")
  dirty=$(git -C "$wt" status --porcelain | grep -c -v '^??' || true)
  untracked=$(git -C "$wt" status --porcelain | grep -c '^??' || true)
  pushed="pushed"
  if ! git -C "$wt" rev-parse -q --verify "origin/$branch" >/dev/null; then pushed="NOT PUSHED"
  elif [ "$(git -C "$wt" rev-parse HEAD)" != "$(git -C "$wt" rev-parse "origin/$branch")" ]; then pushed="local commits not pushed"; fi
  pr=$(gh pr list --head "$branch" --state open --json number --jq '.[0].number' 2>/dev/null)
  echo "- $branch  ($wt)"
  echo "    last commit [$role]: $subject"
  echo "    $ahead commits ahead of main; $pushed; PR: ${pr:-none}"
  if [ "$dirty" -gt 0 ] || [ "$untracked" -gt 0 ]; then
    echo "    UNCOMMITTED: $dirty modified, $untracked untracked (an agent stopped mid-step; inspect before continuing):"
    git -C "$wt" status --porcelain | head -8 | sed 's/^/      /'
  fi
  echo "    next: $(next_step "$role")"
done
echo

echo "== Open PRs =="
gh pr list --state open --json number,title,headRefName --jq '.[] | "- #\(.number) \(.headRefName): \(.title)"' 2>/dev/null || echo "(gh unavailable)"
echo

echo "== GIPPA Slurm jobs =="
squeue --me -h -o '%i %j %T %M %l %R' 2>/dev/null | grep -E ' (giap|gippa)-' || echo "(none running)"
for log in $(ls -t "$GIPPA_LOG_DIR"/gippa-*.out "$GIPPA_LOG_DIR"/giap-*.out 2>/dev/null | head -3); do
  echo "  recent log $(basename "$log"): $(tail -1 "$log")"
done
echo
echo "Next: read PROGRESS.md (Next actions, Questions for Khoa), then continue each worktree's next step."
