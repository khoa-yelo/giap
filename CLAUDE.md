@AGENTS.md

# Claude Code notes
- At session start: read PROGRESS.md, then .giap/capabilities.json (run scripts/probe_env.py if missing).
- The main session is the lead; use the subagents in .claude/agents/ for the other roles (spec §14.1).
- You run on Stanford FarmShare inside a Slurm interactive allocation (spec §3). Source scripts/farmshare_env.sh.
- GitHub: `gh pr`, `gh issue`, `gh run` work (spec §3.10). Never push tags, never push to main.
- Commands longer than ~10 min or bigger than the allocation: sbatch from scripts/slurm/, record the job ID in PROGRESS.md, poll with squeue/sacct (spec §3.8).
- Never run work on a login node; cancel only GIAP jobs by ID; respect the R14 caps.
- Unreachable source: log it under PROGRESS.md › Blocked and continue elsewhere; never use unofficial copies (R6).
- End every work unit: update PROGRESS.md, commit, push.
