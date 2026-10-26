#!/bin/bash
# Is one GHCR tag already published? (M11-14 item 2 "Immutable tags"; ADR-0060 §1-§2)
#   bash scripts/release/ghcr_tag_state.sh ghcr.io/<owner>/<pkg>-<tool>:<tag>
# Runs `docker buildx imagetools inspect REF` (log in first: a private package answers 401/403):
#   - inspect exits 0                                            -> prints `exists`, exit 0
#   - error matching `manifest unknown|not found|NAME_UNKNOWN`  -> prints `absent`, exit 0
#   - any other error (401/403, network, rate limit, no output)  -> the error on stderr, exit 1
# containers.yml calls it in the plan job and again right before the push, so both checks
# classify alike. An authentication failure is never read as "absent".
set -uo pipefail

if [ "$#" -ne 1 ] || [ -z "$1" ]; then
  echo "usage: ghcr_tag_state.sh REF" >&2
  exit 2
fi
ref="$1"
out="$(docker buildx imagetools inspect "$ref" 2>&1)"
rc=$?
if [ "$rc" -eq 0 ]; then
  echo "exists"
  exit 0
fi
if printf '%s\n' "$out" | grep -qE 'manifest unknown|not found|NAME_UNKNOWN'; then
  echo "absent"
  exit 0
fi
{
  echo "ghcr_tag_state: cannot tell whether $ref exists (docker exit $rc):"
  printf '%s\n' "$out"
} >&2
exit 1
