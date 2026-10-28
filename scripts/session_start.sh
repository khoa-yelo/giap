#!/bin/bash
# Runs at every Claude Code session start or resume on FarmShare.
cd "${CLAUDE_PROJECT_DIR:-.}"
source scripts/farmshare_env.sh
if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "WARNING: not inside a Slurm allocation. Do not run tests or tools on a login node (R14)."
fi
if [ -f pyproject.toml ] && command -v uv >/dev/null; then
  uv sync --frozen --all-extras >"$GIPPA_LOG_DIR/uv-sync.log" 2>&1 || echo "uv sync failed (see $GIPPA_LOG_DIR/uv-sync.log)"
fi
if [ -f scripts/probe_env.py ]; then uv run --no-sync python scripts/probe_env.py --quiet || true; fi
# Remind the lead of GIPPA jobs still queued or running
squeue --me --noheader -o '%i %j %T %M %l' 2>/dev/null | grep -E ' (giap|gippa)-' || true
exit 0
