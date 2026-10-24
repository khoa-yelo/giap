#!/bin/bash
# waitmerge.sh PR TITLE — wait until all CI checks on a PR are finished; squash-merge only if all passed.
# Decides from ONE snapshot per poll, and requires the full set of checks to be present.
set -uo pipefail
pr=$1; title=$2; expected=${3:-8}
cd "${GIAP_REPO:-/home/classes/bios/270/khoa/giap}"
source scripts/farmshare_env.sh
# A draft skips the coverage shards (ci.yml), and skipped checks count as passed below: never merge one.
if [ "$(gh pr view "$pr" --json isDraft -q .isDraft 2>/dev/null)" = true ]; then
  echo "NOT MERGED: #$pr is a draft (gh pr ready $pr first)"; exit 1
fi
sleep 30
while true; do
  # A PR that conflicts with its base gets no CI run at all: stop instead of waiting forever.
  if [ "$(gh pr view "$pr" --json mergeable -q .mergeable 2>/dev/null)" = CONFLICTING ]; then
    echo "NOT MERGED: Pull Request has merge conflicts (no CI run); merge the base branch first"; exit 1
  fi
  snap=$(gh pr checks "$pr" 2>/dev/null | cut -f1,2)
  n=$(printf '%s\n' "$snap" | grep -c .)
  if [ "$n" -ge "$expected" ] && ! printf '%s\n' "$snap" | cut -f2 | grep -qE 'pending|queued|in_progress'; then
    break
  fi
  sleep 20
done
printf '%s\n' "$snap"
if printf '%s\n' "$snap" | cut -f2 | grep -qvE '^(pass|skipping)$'; then
  echo "NOT MERGED: a check did not pass"; exit 1
fi
gh pr merge "$pr" --squash --subject "$title (#$pr)" && git fetch -q origin && git log --oneline origin/main -1
